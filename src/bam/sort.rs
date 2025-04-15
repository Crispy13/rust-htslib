use rayon::{
    iter::{
        IndexedParallelIterator, IntoParallelRefIterator, ParallelDrainRange, ParallelIterator,
    },
    slice::ParallelSliceMut,
};
use regex::{CaptureMatches, Matches, Regex, bytes};
use tempfile::NamedTempFile;

use super::{Header, Read as _, record::Record};
use anyhow::{Error, anyhow};
use std::{
    borrow::Borrow,
    cell::{LazyCell, OnceCell},
    cmp::Ordering,
    collections::{BinaryHeap, HashSet},
    num::NonZero,
    sync::OnceLock,
};

/*
TODO:
1. use memory block and file simultaneously (currently, either one of them is used.)
2. performance? currently samtools vs this: 5.6s vs 9s.
*/
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SortBy {
    QueryNameLexi,
    QueryNameNatural,
    Coordinate,
}

impl SortBy {
    #[allow(non_snake_case)]
    const fn SO_str(&self) -> &'static str {
        match self {
            SortBy::QueryNameLexi | SortBy::QueryNameNatural => "queryname",
            SortBy::Coordinate => "coordinate",
        }
    }
}

impl SortBy {
    #[inline]
    fn compare(&self, a: &Record, b: &Record) -> Ordering {
        match self {
            SortBy::QueryNameLexi => {
                match a.qname().cmp(b.qname()) {
                    Ordering::Equal => {}
                    oth => return oth,
                }
                Self::compare_read_num(a, b)
            }
            SortBy::QueryNameNatural => {
                match natural_sort_cmp(
                    std::str::from_utf8(a.qname()).unwrap(),
                    std::str::from_utf8(b.qname()).unwrap(),
                ) {
                    Ordering::Equal => {}
                    oth => return oth,
                }

                Self::compare_read_num(a, b)
            }
            SortBy::Coordinate => {
                let tid_cmp = a.tid().cmp(&b.tid());
                if tid_cmp == Ordering::Equal {
                    a.pos().cmp(&b.pos())
                } else {
                    tid_cmp
                }
            }
        }
    }

    #[inline]
    fn compare_read_num(a: &Record, b: &Record) -> Ordering {
        {
            if a.is_first_in_template() {
                1
            } else if a.is_last_in_template() {
                2
            } else {
                unreachable!()
            }
        }
        .cmp(&{
            if b.is_first_in_template() {
                1
            } else if b.is_last_in_template() {
                2
            } else {
                unreachable!()
            }
        })
    }

    #[inline]
    fn equal(&self, a: &Record, b: &Record) -> bool {
        match self {
            SortBy::QueryNameLexi | SortBy::QueryNameNatural => a.qname().eq(b.qname()),
            SortBy::Coordinate => a.tid().eq(&b.tid()) && a.pos().eq(&b.pos()),
        }
    }
}

/// ## NOTE:  
/// This is self-referential struct.  
/// So Never move this struct or modify the inner field.
// struct RecordNaturalSortData<'a> {
//     /// never modify after initialized.
//     inner: Box<Record>,
//     sort_data: Vec<NaturalSortElem<'a>>,
// }

// impl<'a> RecordNaturalSortData<'a> {
//     unsafe fn new(inner: Record) -> Result<Self, Error> {
//         let mut s = Self {
//             inner:inner.into(),
//             sort_data: vec![],
//         };

//         // let p = &s.inner as *const _;

//         s.sort_data = NaturalSortElemIter::new(std::str::from_utf8({
//             let p: *const [u8] = s.inner.qname();

//             unsafe { &*p }
//         })?)
//         .collect();

//         Ok(s)
//     }
// }

struct OptionArrIter<'a, T, const N: usize> {
    inner: &'a [Option<T>; N],
    idx: usize,
}

impl<'a, T, const N: usize> OptionArrIter<'a, T, N> {
    fn new(inner: &'a [Option<T>; N]) -> Self {
        Self { inner, idx: 0 }
    }
}

impl<'a, T, const N: usize> Iterator for OptionArrIter<'a, T, N> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx < N {
            let r = match &self.inner[self.idx] {
                Some(v) => Some(v),
                None => None,
            };
            self.idx += 1;
            r
        } else {
            None
        }
    }
}

// struct NaturalSortKey([Option<NaturalSortElem>; 256], usize);
struct NaturalSortKey(Vec<NaturalSortElem>, usize);

impl NaturalSortKey {
    fn from_record(record: &Record) -> Result<Self, Error> {
        let v = NaturalSortElemIter::new(std::str::from_utf8(record.qname())?).collect::<Vec<_>>();
        // let v = NaturalSortElemExtractor::new(record.qname()).collect();

        // let v = {
        //     let mut v = [const { None }; 256];
        //     let iter = NaturalSortElemIter::new(std::str::from_utf8(record.qname())?);
        //     iter.enumerate().for_each(|(i, e)| {
        //         if i < 256 {
        //             v[i] = Some(e);
        //         } else {
        //             panic!("get parsed elems more than 256.")
        //         }
        //     });

        //     v
        // };

        Ok(Self(v, get_read_num(record)))
    }
}

impl PartialOrd for NaturalSortKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for NaturalSortKey {
    fn cmp(&self, other: &Self) -> Ordering {
        match _natural_sort_cmp(self.0.iter(), other.0.iter()) {
            Ordering::Equal => {}
            oth => return oth,
        }

        self.1.cmp(&other.1)
    }
}

impl PartialEq for NaturalSortKey {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0 && self.1 == other.1
    }
}

impl Eq for NaturalSortKey {}

fn get_read_num(record: &Record) -> usize {
    if record.is_first_in_template() {
        1
    } else if record.is_last_in_template() {
        2
    } else {
        3
        // panic!("not first and not last: {}", String::from_utf8_lossy(record.qname()))
    }
}

struct NaturalSortElemExtractor<'a> {
    hay: &'a [u8],
    idx: usize,
}

impl<'a> NaturalSortElemExtractor<'a> {
    fn new(hay: &'a [u8]) -> Self {
        Self { hay, idx: 0 }
    }
}

impl<'a> Iterator for NaturalSortElemExtractor<'a> {
    type Item = NaturalSortElem;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.hay.len() {
            return None;
        }

        let hay = &self.hay[self.idx..];

        // let start = self.idx;
        let mut end = 0;
        let r = if hay[0].is_ascii_digit() {
            end += 1;

            if hay.len() > 1 {
                for b in &hay[1..] {
                    if b.is_ascii_digit() {
                        end += 1;
                    } else {
                        break;
                    }
                }
            }

            NaturalSortElem::Num(
                std::str::from_utf8(&hay[0..end])
                    .unwrap()
                    .parse::<i64>()
                    .unwrap_or_else(|err| {
                        panic!("{err:?} {}", std::str::from_utf8(&hay[0..end]).unwrap())
                    }),
            )
        } else {
            end += 1;

            if hay.len() > 1 {
                for b in &hay[1..] {
                    if !b.is_ascii_digit() {
                        end += 1;
                    } else {
                        break;
                    }
                }
            }

            // NaturalSortElem::Text(String::from_utf8_lossy(&hay[0..end]).into_owned())
            NaturalSortElem::Text(String::from_utf8(hay[0..end].to_vec()).unwrap())
        };

        self.idx += end;

        Some(r)
    }
}

/// ## Examples
/// ### 1. Normal
/// ```ignore
/// let mut rs = RecordSorter::new(100000, SortBy::QueryNameNatural, 4);
///
/// for record in indexed_reader.records() {
///     rs.add_record(record); // if inner vector is full, then sort and dump stored records into a temp file.
/// }
///
/// rs.sort()?;
///
/// let iter = rs.into_iter()?.collect::<Vec<_>>();
///
/// ```
/// ### 2. No chunking (All records in memory)
/// ```ignore
/// let mut rs = RecordSorter::new(0, SortBy::QueryNameNatural, 4); // note that chunk_size is 0.
/// ...
///
/// // we did not use chunking. so all records are in memory(vector).
/// let sorted_records = rs.into_record_vec();
///
/// ```
pub struct RecordSorter {
    chunk_size: usize,              // How many records to process per chunk.
    temp_files: Vec<NamedTempFile>, // Temporary file paths for sorted chunks.
    sort_by: SortBy,                // Your comparator enum.
    record_vec: Vec<Record>,
    hts_tp: Option<crate::tpool::ThreadPool>,
    rayon_tp: Option<rayon::ThreadPool>,
}

impl RecordSorter {
    pub fn new(chunk_size: usize, sort_by: SortBy, threads: usize) -> Result<Self, anyhow::Error> {
        Ok(Self {
            chunk_size,
            temp_files: Vec::new(),
            sort_by,
            record_vec: Vec::with_capacity(chunk_size),
            hts_tp: {
                if threads > 1 {
                    Some(crate::tpool::ThreadPool::new(threads as u32)?)
                } else {
                    None
                }
            },
            rayon_tp: if threads > 1 {
                Some(
                    rayon::ThreadPoolBuilder::new()
                        .num_threads(threads)
                        .build()?,
                )
            } else {
                None
            },
        })
    }

    fn hts_tp(&self) -> Option<&crate::tpool::ThreadPool> {
        self.hts_tp.as_ref()
    }

    #[inline]
    pub fn add_record(&mut self, record: Record) -> Result<(), anyhow::Error> {
        if self.chunk_size != 0 && self.record_vec.len() == self.chunk_size {
            self.sort_and_dump()?;
        }

        self.record_vec.push(record);

        Ok(())
    }

    // fn check_len_and_dump(&mut self) -> Result<(), anyhow::Error> {
    //     if self.record_vec.len() == self.chunk_size {
    //         self.sort_and_dump()?;
    //     }

    //     Ok(())
    // }

    pub fn into_record_vec(self) -> Vec<Record> {
        self.record_vec
    }

    pub fn into_iter(self) -> Result<SortedRecordIterator, Error> {
        SortedRecordIterator::new(self)
    }

    pub fn sort(&mut self) -> Result<(), Error> {
        match &self.rayon_tp {
            Some(tp) => tp.scope(|_s| {
                match self.sort_by {
                    SortBy::QueryNameNatural => {
                        eprintln!("making sort_data");

                        let mut sort_data = self
                            .record_vec
                            .par_drain(..)
                            // .enumerate()
                            .map(|e| Ok((NaturalSortKey::from_record(&e)?, e)))
                            .collect::<Result<Vec<_>, Error>>()?;

                        eprintln!("made sort_data");

                        sort_data.par_sort_by(|a, b| a.0.cmp(&b.0));
                        eprintln!("sorted.");

                        self.record_vec.extend(sort_data.into_iter().map(|v| v.1));
                        // self.record_vec
                        //     .par_sort_by(|a, b| self.sort_by.compare(a, b));

                        // self.record_vec.par_sort_by_cached_key(|r| NaturalSortKey::from_record(r));
                    }
                    SortBy::QueryNameLexi | SortBy::Coordinate => {
                        self.record_vec
                            .par_sort_by(|a, b| self.sort_by.compare(&a, &b));
                    }
                }
                Result::<_, Error>::Ok(())
            })?,
            None => {
                self.record_vec.sort_by(|a, b| self.sort_by.compare(&a, &b));
            }
        }

        Ok(())

        // self.record_vec.sort_by(|a, b| self.sort_by.compare(a, b));
    }

    fn dump(&mut self) -> Result<(), anyhow::Error> {
        let tmp_file = NamedTempFile::new()?;
        let header = super::Header::default();
        // let header = self
        //     .record_vec
        //     .first()
        //     .ok_or_else(|| anyhow!("dump called when having no record."))?
        //     .header
        //     .as_ref()
        //     .map_or(super::Header::default(), |v| {
        //         super::Header::from_template(&v)
        //     });

        let mut writer = super::Writer::from_path(tmp_file.path(), &header, super::Format::Bam)?;

        if let Some(ref tp) = self.hts_tp {
            writer.set_thread_pool(tp)?;
        }

        for record in self.record_vec.drain(..) {
            writer.write(&record)?;
        }

        self.temp_files.push(tmp_file);

        Ok(())
    }

    fn sort_and_dump(&mut self) -> Result<(), anyhow::Error> {
        self.sort()?;

        self.dump()?;

        Ok(())
    }

    // fn finalize(&mut self) -> Result<(), anyhow::Error> {
    //     match (self.record_vec.len() > 0, self.temp_files.len() > 0) {
    //         (true, true) => {
    //             self.sort_and_dump()?;
    //         }
    //         (true, false) => {
    //             // do nothing
    //             if self.chunk_size != 0 {
    //                 assert!(self.record_vec.len() <= self.chunk_size);
    //             }
    //         }
    //         (false, true) => {
    //             // do nothing
    //         }
    //         (false, false) => {
    //             // do nothing
    //         }
    //     }

    //     Ok(())
    // }

    // pub fn sort_and_into_iter(mut self) -> Result<SortedRecordIterator, Error> {
    //     self.finalize()?;

    //     let bam_readers = self
    //         .temp_files
    //         .iter()
    //         .map(|p| super::Reader::from_path(p.path()).map_err(Error::from))
    //         .collect::<Result<Vec<_>, Error>>()?;

    //     let from_mem = bam_readers.is_empty();

    //     assert!(!(self.record_vec.len() > 0 && bam_readers.len() > 0));

    //     let iter = self.record_vec.into_iter();

    //     Ok(SortedRecordIterator {
    //         mem_iter: iter,
    //         chunk_readers: bam_readers,
    //         from_mem,
    //     })
    // }
}

struct HeapItem {
    record: Record,
    source_idx: usize,
    sort_by: SortBy,
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse the ordering (since BinaryHeap is a max-heap) for ascending order.
        self.sort_by.compare(&other.record, &self.record)
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.sort_by.equal(&self.record, &other.record)
    }
}

impl Eq for HeapItem {}

enum RecordIteratorFromChunk {
    Mem(std::vec::IntoIter<Record>),
    File(super::Reader),
}

impl Iterator for RecordIteratorFromChunk {
    type Item = Result<Record, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            RecordIteratorFromChunk::Mem(it) => it.next().map(Result::Ok),
            RecordIteratorFromChunk::File(reader) => {
                let mut record = Record::default();
                reader
                    .read(&mut record)
                    .map(|e| e.map(|_| record).map_err(Error::from))
            }
        }
    }
}

pub struct SortedRecordIterator {
    mem_iter: std::vec::IntoIter<Record>,
    chunk_readers: Vec<RecordIteratorFromChunk>,
    binary_heap: BinaryHeap<HeapItem>,
    in_mem_only: bool,
    sort_by: SortBy,
}

impl SortedRecordIterator {
    pub fn new(record_sorter: RecordSorter) -> Result<Self, Error> {
        // record_sorter.finalize()?;

        let in_mem_only = record_sorter.temp_files.is_empty();

        let mut mem_iter = vec![].into_iter();
        let mut binary_heap = BinaryHeap::new();

        let mut chunk_readers = vec![];

        if !in_mem_only {
            let iter = record_sorter.record_vec.into_iter();

            chunk_readers.push(RecordIteratorFromChunk::Mem(iter));

            for tmpfile in record_sorter.temp_files {
                chunk_readers.push(RecordIteratorFromChunk::File(
                    super::Reader::from_path(tmpfile.path()).map_err(Error::from)?,
                ));
            }

            // let mut chunk_readers = record_sorter
            //     .temp_files
            //     .iter()
            //     .map(|p| super::Reader::from_path(p.path()).map_err(Error::from))
            //     .collect::<Result<Vec<_>, Error>>()?;

            // assert!(!(record_sorter.record_vec.len() > 0 && chunk_readers.len() > 0));

            binary_heap.reserve(chunk_readers.len());

            for (chunk_idx, reader) in chunk_readers.iter_mut().enumerate() {
                if let Some(rr) = reader.next() {
                    let next_record = rr?;

                    binary_heap.push(HeapItem {
                        record: next_record,
                        source_idx: chunk_idx,
                        sort_by: record_sorter.sort_by,
                    });
                }
            }
        } else {
            mem_iter = record_sorter.record_vec.into_iter();
        }

        Ok(Self {
            mem_iter,
            chunk_readers,
            binary_heap,
            in_mem_only,
            sort_by: record_sorter.sort_by,
        })
    }
}

impl Iterator for SortedRecordIterator {
    type Item = Result<Record, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let elem = if self.in_mem_only {
            self.mem_iter.next()
        } else {
            if let Some(HeapItem {
                record,
                source_idx: chunk_idx,
                sort_by,
            }) = self.binary_heap.pop()
            {
                let reader = &mut self.chunk_readers[chunk_idx];

                if let Some(rr) = reader.next() {
                    let next_record = match rr {
                        Ok(r) => r,
                        Err(err) => return Some(Err(err).map_err(Error::from)),
                    };

                    self.binary_heap.push(HeapItem {
                        record: next_record,
                        source_idx: chunk_idx,
                        sort_by,
                    });
                }

                Some(record)
            } else {
                None
            }
        };

        elem.map(Result::Ok)
    }
}

/*
TODO: Make a reader wrapper, read records with thread pool.
it has inner buffer, to give readed record to the caller.
*/

//— MergedSortedRecordIterator ------------------------------------------------

// This iterator merges multiple SortedRecordIterator sources.
// Each SortedRecordIterator yields items of type `Result<Record, Error>`.
// The merged iterator will yield the overall next record in sorted order.
pub struct MergedSortedRecordIterator {
    sorted_iters: Vec<SortedRecordIterator>,
    merge_heap: BinaryHeap<HeapItem>,
}

impl MergedSortedRecordIterator {
    /// Build a merged iterator from a vector of SortedRecordIterator.
    /// This will fetch the first record from each iterator and place it
    /// into a binary heap so that the minimum (per sort_by) is always available.
    pub fn new(mut sorted_iters: Vec<SortedRecordIterator>) -> Result<Self, Error> {
        let mut merge_heap = BinaryHeap::new();

        if sorted_iters
            .iter()
            .map(|e| e.sort_by)
            .collect::<HashSet<_>>()
            .len()
            > 1
        {
            panic!("Multiple sort_by detected")
        }

        // For each SortedRecordIterator, call try_next() and push the record onto the heap.
        for (idx, iter) in sorted_iters.iter_mut().enumerate() {
            if let Some(rr) = iter.next() {
                let record = rr?;

                merge_heap.push(HeapItem {
                    record,
                    source_idx: idx,
                    // Assume that each SortedRecordIterator knows which sort order
                    // it used. In our case, we carry the sort_by along from the chunk sorter.
                    sort_by: iter.sort_by,
                });
            }
        }

        Ok(Self {
            sorted_iters,
            merge_heap,
        })
    }
}

impl Iterator for MergedSortedRecordIterator {
    type Item = Result<Record, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        // If our heap is empty, all iterators are exhausted.
        let top = self.merge_heap.pop()?;
        let result_record = top.record;
        let iter_idx = top.source_idx;
        let sort_by = top.sort_by; // same sort order

        // Now, fetch the next record from the corresponding SortedRecordIterator.
        match self.sorted_iters[iter_idx].next() {
            Some(Err(err)) => return Some(Err(err)),
            Some(Ok(next_record)) => {
                self.merge_heap.push(HeapItem {
                    record: next_record,
                    source_idx: iter_idx,
                    sort_by,
                });
            }
            None => {} // This sub-iterator is exhausted.
        }
        Some(Ok(result_record))
    }
}

struct NaturalSorter<'a> {
    lhs_elems: Vec<&'a [u8]>,
    rhs_elems: Vec<&'a [u8]>,
}

impl<'a> NaturalSorter<'a> {
    /// ```python
    /// import re
    /// def natural_sort(l):
    ///     convert = lambda text: int(text) if text.isdigit() else text.lower()
    ///     alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    ///     return sorted(l, key=alphanum_key)
    /// ```
    fn natural_sort_cmp(&mut self, lhs: &str, rhs: &str) {
        static NUM_REP: OnceLock<Regex> = OnceLock::new();

        let num_rep = NUM_REP.get_or_init(|| Regex::new(r"\d+").unwrap());

        for caps in num_rep.captures_iter(lhs) {
            let g0 = caps.get(0).unwrap();
            let start = g0.start();
            let end = g0.end();

            if start > 0 {}
        }
    }
}

struct NaturalSortElemIter<'r, 'h> {
    cap_matches: Matches<'r, 'h>,
    next_idx: [Option<(usize, usize)>; 2],
    hay: &'h str,
}

impl<'r, 'h> NaturalSortElemIter<'r, 'h> {
    // thread_local! {
    //     static NUM_REP: LazyCell<Regex> = LazyCell::new(|| Regex::new(r"\d+").unwrap());
    // }

    fn num_rep() -> &'static Regex {
        static NUM_REP: OnceLock<Regex> = OnceLock::new();

        NUM_REP.get_or_init(|| Regex::new(r"\d+").unwrap())
    }

    fn new(hay: &'h str) -> Self {
        let cap_matches = Self::num_rep().find_iter(&hay);

        let mut s = Self {
            cap_matches,
            next_idx: [None; 2],
            hay,
        };

        if let Some(g0) = s.cap_matches.next() {
            // let g0 = c1aps.unwrap();
            let n_start = g0.start();
            let n_end = g0.end();

            let _ = s.next_idx[1].insert((n_start, n_end));

            if 0 < n_start {
                // a value before this num
                let _ = s.next_idx[0].insert((0, n_start));
            }
        } else {
            let _ = s.next_idx[1].insert((0, hay.len()));
        }

        s
    }
}

impl<'r, 'h> Iterator for NaturalSortElemIter<'r, 'h> {
    type Item = NaturalSortElem;

    fn next(&mut self) -> Option<Self::Item> {
        let mut elem = None;
        for (i, idx_opt) in self.next_idx.iter_mut().enumerate() {
            if let Some((start, end)) = idx_opt.take() {
                let _ = elem.insert(&self.hay[start..end]);

                if i == 1 {
                    // get next idxs from cap_matches
                    if let Some(g0) = self.cap_matches.next() {
                        // let g0 = caps.get(0).unwrap();
                        let n_start = g0.start();
                        let n_end = g0.end();

                        let _ = self.next_idx[1].insert((n_start, n_end));

                        if end < n_start {
                            // a value before this num
                            let _ = self.next_idx[0].insert((end, n_start));
                        }
                    } else {
                        if end < self.hay.len() {
                            let _ = self.next_idx[1].insert((end, self.hay.len()));
                        }
                        break;
                    }
                }

                break;
            }
        }

        if let Some(elem) = elem {
            if elem.as_bytes().get(0).map_or(false, |b| b.is_ascii_digit()) {
                Some(NaturalSortElem::Num(elem.parse::<i64>().unwrap()))
            } else {
                Some(NaturalSortElem::Text(elem.to_owned()))
            }
        } else {
            None
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone)]
enum NaturalSortElem {
    Num(i64),
    // Text(&'a str),
    Text(String),
}

impl PartialOrd for NaturalSortElem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for NaturalSortElem {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (NaturalSortElem::Num(a), NaturalSortElem::Num(b)) => a.cmp(b),
            (NaturalSortElem::Text(a), NaturalSortElem::Text(b)) => a.cmp(b),
            (NaturalSortElem::Num(_), NaturalSortElem::Text(_)) => Ordering::Less,
            (NaturalSortElem::Text(_), NaturalSortElem::Num(_)) => Ordering::Greater,
        }
    }
}

fn _natural_sort_cmp<'a>(
    mut a_iter: impl Iterator<Item = impl Borrow<NaturalSortElem>>,
    mut b_iter: impl Iterator<Item = impl Borrow<NaturalSortElem>>,
) -> std::cmp::Ordering {
    loop {
        match (a_iter.next(), b_iter.next()) {
            (Some(a), Some(b)) => match a.borrow().cmp(b.borrow()) {
                Ordering::Equal => {}
                oth => return oth,
            },
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
            (None, None) => {
                break;
            }
        }
    }

    Ordering::Equal
}

fn natural_sort_cmp(lhs: &str, rhs: &str) -> std::cmp::Ordering {
    let a_iter = NaturalSortElemIter::new(lhs);
    let b_iter = NaturalSortElemIter::new(rhs);

    _natural_sort_cmp(a_iter, b_iter)
}

/// Modify sort-related tags in bam header.
/// this modifies:
/// 1. SO tag
/// 2. SS tag (if sorted by query name)
/// 3. GO tag (remove if exists)
pub fn modify_sort_tag_in_header(header: &mut Header, sort_by: SortBy) {
    // The @HD line should be present, with either the SO tag or the GO tag (but not both) specified.

    if let Some(hr) = header.records_mut().first_mut() {
        let rep = regex::bytes::Regex::new(r"^@HD((?:\t[^\t\n]+)*)(\n)?([\S\s]*)*").unwrap();

        let mut new_hd_line = vec![b"@HD".to_vec()];
        match rep.captures(&hr) {
            Some(caps) => {
                let mut found_so_tag = false;

                macro_rules! modify_sort_tag {
                    () => {
                        new_hd_line.push(format!("SO:{}", sort_by.SO_str()).into_bytes());

                        match sort_by {
                            SortBy::QueryNameLexi => {
                                new_hd_line.push("SS:queryname:lexicographical".into())
                            }
                            SortBy::QueryNameNatural => {
                                new_hd_line.push("SS:queryname:natural".into())
                            }
                            _ => {}
                        }
                    };
                }

                if let Some(values) = caps.get(1) {
                    for value in values.as_bytes().trim_ascii().split(|&b| b == b'\t') {
                        if value.starts_with(b"SO") {
                            found_so_tag = true;

                            modify_sort_tag!();
                        } else if value.starts_with(b"SS") | value.starts_with(b"GO") {
                            eprintln!(
                                "Dropping this tag from original bam: {}",
                                String::from_utf8_lossy(value)
                            );
                            continue;
                        } else {
                            new_hd_line.push(value.to_vec())
                        }
                    }
                }

                if !found_so_tag {
                    modify_sort_tag!();
                }

                let mut new_hd_line2 = new_hd_line.join(&b'\t');

                if caps.get(2).is_some() {
                    new_hd_line2.push(b'\n');
                };

                if let Some(values) = caps.get(3) {
                    // eprintln!("Added rest header: {} bytes", values.as_bytes().len());

                    new_hd_line2.extend(values.as_bytes());
                }

                let _ = std::mem::replace(hr, new_hd_line2);
            }
            None => panic!("Invalid header: {}", String::from_utf8_lossy(hr)),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::rc::Rc;

    use regex::Regex;

    use super::super::Reader;
    use super::*;
    // Assume TEST_BAM is the path to a small BAM file for testing.
    const TEST_BAM: &str = "/home/eck/workspace/common_resources/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam";
    const THREADS: usize = 4;

    #[test]
    fn test_sorted_record_iterator_from_memory() -> Result<(), Box<dyn std::error::Error>> {
        #[global_allocator]
        static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

        // Obtain header from the test BAM.
        let header = {
            let mut reader = Reader::from_path(TEST_BAM).expect("Could not open test BAM");
            reader.header().to_owned()
        };

        let threads = THREADS;

        // Create a RecordSorter with a small chunk size.
        let mut record_sorter = RecordSorter::new(0, SortBy::QueryNameNatural, threads)?;

        // Read the BAM file and add records to the sorter.
        let mut reader = Reader::from_path(TEST_BAM).expect("Failed to open test BAM");

        let tp = if threads > 1 {
            Some(crate::tpool::ThreadPool::new(threads as u32)?)
        } else {
            None
        };

        if let Some(ref tp) = tp {
            reader.set_thread_pool(tp)?;
        }
        {
            for result in reader.records() {
                let record = result.expect("Error reading record");
                record_sorter
                    .add_record(record)
                    .expect("Failed to add record");
            }
        }

        // Finalize so that any remaining in-memory records are flushed.
        record_sorter.sort()?;
        // record_sorter.finalize().expect("Finalize failed");

        // Create a SortedRecordIterator from the sorter.
        let mut sorted_iter = record_sorter.into_iter()?;

        // Iterate over sorted records and check sort order by query name.
        let mut count = 0;
        let mut prev_record: Option<Record> = None;

        let mut header = super::super::Header::from_template(reader.header());
        // The @HD line should be present, with either the SO tag or the GO tag (but not both) specified.

        modify_sort_tag_in_header(&mut header, sorted_iter.sort_by);
        // let mut header = super::super::Header::from_template(reader.header());

        let header_view = Rc::new(super::super::HeaderView::from_header(&header));

        let mut writer = super::super::Writer::from_path(
            TEST_BAM[TEST_BAM.rfind("/").map_or(0, |v| v + 1)..].replace(".bam", ".sort.cus.bam"),
            &header,
            super::super::Format::Bam,
        )?;

        if let Some(ref tp) = tp {
            writer.set_thread_pool(tp)?;
        }

        // let def_record = Record::default();

        while let Some(result) = sorted_iter.next() {
            let mut record = result.expect("Error iterating sorted record");
            // let cur_name = record.qname().to_vec();
            if let Some(prev) = prev_record {
                // let cmp = sorted_iter.sort_by.compare(&prev, &record);
                // assert!(
                //     cmp != Ordering::Greater,
                //     "i:{}, Records are not sorted by query name: {} > {} {:?}",
                //     count,
                //     String::from_utf8_lossy(prev.qname()),
                //     String::from_utf8_lossy(record.qname()),
                //     cmp,
                // );
            }
            // assert_ne!(def_record, record);

            writer.write(&record)?;

            prev_record = Some(record);
            count += 1;

            // record.set_header(header_view.clone());
        }
        println!("Total sorted records (in memory): {}", count);
        assert!(count > 0, "No records found");
        Ok(())
    }

    #[test]
    fn test_merged_sorted_record_iterator() -> Result<(), Box<dyn std::error::Error>> {
        // Obtain header from the test BAM.
        let header = {
            let mut reader = Reader::from_path(TEST_BAM).expect("Could not open test BAM");
            reader.header().to_owned()
        };

        // Create two RecordSorter instances to simulate separate sorted iterators.
        let mut record_sorter1 = RecordSorter::new(500, SortBy::QueryNameLexi, THREADS)?;
        let mut record_sorter2 = RecordSorter::new(500, SortBy::QueryNameLexi, THREADS)?;

        {
            let mut reader = Reader::from_path(TEST_BAM).expect("Failed to open test BAM");
            let mut toggle = false;
            for result in reader.records() {
                let record = result.expect("Error reading record");
                if toggle {
                    record_sorter1
                        .add_record(record)
                        .expect("Failed to add record in sorter1");
                } else {
                    record_sorter2
                        .add_record(record)
                        .expect("Failed to add record in sorter2");
                }
                toggle = !toggle;
            }
        }

        // record_sorter1.finalize().expect("Finalize sorter1 failed");
        // record_sorter2.finalize().expect("Finalize sorter2 failed");

        let sorted_iter1 =
            SortedRecordIterator::new(record_sorter1).expect("Failed to build sorted_iter1");
        let sorted_iter2 =
            SortedRecordIterator::new(record_sorter2).expect("Failed to build sorted_iter2");

        let mut merged_iter = MergedSortedRecordIterator::new(vec![sorted_iter1, sorted_iter2])
            .expect("Failed to build MergedSortedRecordIterator");

        let mut count = 0;
        let mut prev_name: Option<Vec<u8>> = None;

        while let Some(result) = merged_iter.next() {
            let record = result.expect("Error in merged iterator");
            let cur_name = record.qname().to_vec();
            if let Some(prev) = prev_name {
                assert!(
                    prev <= cur_name,
                    "Merged records are not sorted by query name"
                );
            }
            prev_name = Some(cur_name);
            count += 1;
        }
        println!("Total merged sorted records: {}", count);
        assert!(count > 0, "No records found in merged iterator");
        Ok(())
    }

    #[test]
    fn test_split() -> Result<(), Box<dyn std::error::Error>> {
        let rep = Regex::new(r"(\d+)")?;

        let hay1 = "asdf1234.a111234";

        eprintln!("{:?}", rep.split(&hay1).collect::<Vec<_>>());

        Ok(())
    }

    #[test]
    fn test_natural_elem_iter() {
        let hay = "aaa123aabc.234";

        assert_eq!(
            vec![
                NaturalSortElem::Text("aaa".into()),
                NaturalSortElem::Num(123),
                NaturalSortElem::Text("aabc.".into()),
                NaturalSortElem::Num(234)
            ],
            NaturalSortElemIter::new(hay).collect::<Vec<_>>()
        );

        let hay = "aaa";

        assert_eq!(
            vec![NaturalSortElem::Text("aaa".into()),],
            NaturalSortElemIter::new(hay).collect::<Vec<_>>()
        );

        let hay = "1aaa";

        assert_eq!(
            vec![NaturalSortElem::Num(1), NaturalSortElem::Text("aaa".into()),],
            NaturalSortElemIter::new(hay).collect::<Vec<_>>()
        );

        let hay = "11";

        assert_eq!(
            vec![NaturalSortElem::Num(11),],
            NaturalSortElemIter::new(hay).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_natural_sort_cmp() {
        let a = "aaa1.33";
        let b = "aaa1.1234";

        assert_eq!(natural_sort_cmp(a, b), Ordering::Less);

        let a = "aaa1.3333";
        let b = "aaa1.1234";

        assert_eq!(natural_sort_cmp(a, b), Ordering::Greater);

        let a = "baa1.33";
        let b = "aaa1.1234";

        assert_eq!(natural_sort_cmp(a, b), Ordering::Greater);

        let a = "1baa1.33";
        let b = "aaa1.1234";

        assert_eq!(natural_sort_cmp(a, b), Ordering::Less);
    }

    #[test]
    fn test_equal_strings() {
        assert_eq!(natural_sort_cmp("a12b3", "a12b3"), Ordering::Equal);
    }

    #[test]
    fn test_prefix() {
        // "a12b" is a prefix of "a12b3", so should come before.
        assert_eq!(natural_sort_cmp("a12b", "a12b3"), Ordering::Less);
        assert_eq!(natural_sort_cmp("a12b3", "a12b"), Ordering::Greater);
    }

    #[test]
    fn test_different_tokens() {
        // Compare where first token (text) is equal, then numeric
        // "a7b" vs "a12b" should compare by the numeric token (7 < 12)
        assert_eq!(natural_sort_cmp("a7b", "a12b"), Ordering::Less);
    }

    #[test]
    fn test_complex() {
        // A mixed example comparing different tokens
        assert_eq!(
            natural_sort_cmp("file2version10", "file2version2"),
            Ordering::Greater
        );
        // Here, the common "file2version" is equal,
        // and then numerical tokens 10 > 2, so ordering is Greater.
    }

    #[test]
    fn test_single_number() {
        let input = b"123";
        let mut extractor = NaturalSortElemExtractor::new(input);
        assert_eq!(extractor.next(), Some(NaturalSortElem::Num(123)));
        assert_eq!(extractor.next(), None);
    }

    #[test]
    fn test_single_text() {
        let input = b"abc";
        let mut extractor = NaturalSortElemExtractor::new(input);
        assert_eq!(
            extractor.next(),
            Some(NaturalSortElem::Text("abc".to_string()))
        );
        assert_eq!(extractor.next(), None);
    }

    #[test]
    fn test_mixed_elements() {
        // Input: "abc123def"
        // expected:
        //   NaturalSortElem::Text("abc")
        //   NaturalSortElem::Num(123)
        //   NaturalSortElem::Text("def")
        let input = b"abc123def";
        let elems: Vec<_> = NaturalSortElemExtractor::new(input).collect();
        let expected = vec![
            NaturalSortElem::Text("abc".to_string()),
            NaturalSortElem::Num(123),
            NaturalSortElem::Text("def".to_string()),
        ];
        assert_eq!(elems, expected);
    }

    #[test]
    fn test_leading_zeros_number() {
        // Input with leading zeros: "007abc"
        // expected:
        //   NaturalSortElem::Num(7)    // "007" parses to 7
        //   NaturalSortElem::Text("abc")
        let input = b"007abc";
        let elems: Vec<_> = NaturalSortElemExtractor::new(input).collect();
        let expected = vec![
            NaturalSortElem::Num(7),
            NaturalSortElem::Text("abc".to_string()),
        ];
        assert_eq!(elems, expected);
    }

    #[test]
    fn test_empty_input() {
        let input = b"";
        let elems: Vec<_> = NaturalSortElemExtractor::new(input).collect();
        assert_eq!(elems.len(), 0);
    }
}
