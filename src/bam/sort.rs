use rayon::slice::ParallelSliceMut;
use regex::{CaptureMatches, Regex, bytes};
use tempfile::NamedTempFile;

use super::{Read as _, record::Record};
use anyhow::{Error, anyhow};
use std::{
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
    QueryName,
    Coordinate,
}

impl SortBy {
    #[allow(non_snake_case)]
    const fn SO_str(&self) -> &'static str {
        match self {
            SortBy::QueryName => "queryname",
            SortBy::Coordinate => "coordinate",
        }
    }
}

impl SortBy {
    #[inline]
    fn compare(&self, a: &Record, b: &Record) -> Ordering {
        match self {
            SortBy::QueryName => a.qname().cmp(b.qname()),
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

    fn equal(&self, a: &Record, b: &Record) -> bool {
        match self {
            SortBy::QueryName => a.qname().eq(b.qname()),
            SortBy::Coordinate => a.tid().eq(&b.tid()) && a.pos().eq(&b.pos()),
        }
    }
}

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

    pub fn sort(&mut self) {
        match &self.rayon_tp {
            Some(tp) => tp.scope(|_s| {
                self.record_vec
                    .par_sort_by(|a, b| self.sort_by.compare(a, b));
            }),
            None => {
                self.record_vec.sort_by(|a, b| self.sort_by.compare(a, b));
            }
        }

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
        self.sort();

        self.dump()?;

        Ok(())
    }

    fn finalize(&mut self) -> Result<(), anyhow::Error> {
        match (self.record_vec.len() > 0, self.temp_files.len() > 0) {
            (true, true) => {
                self.sort_and_dump()?;
            }
            (true, false) => {
                // do nothing
                if self.chunk_size != 0 {
                    assert!(self.record_vec.len() <= self.chunk_size);
                }
            }
            (false, true) => {
                // do nothing
            }
            (false, false) => {
                // do nothing
            }
        }

        Ok(())
    }

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
    cap_matches: CaptureMatches<'r, 'h>,
    next_idx: [Option<(usize, usize)>; 2],
    hay: &'h str,
}

impl<'r, 'h> NaturalSortElemIter<'r, 'h> {
    fn num_rep() -> &'static Regex {
        static NUM_REP: OnceLock<Regex> = OnceLock::new();

        NUM_REP.get_or_init(|| Regex::new(r"\d+").unwrap())
    }

    fn new(hay: &'h str) -> Self {
        let cap_matches = Self::num_rep().captures_iter(&hay);

        let mut s = Self {
            cap_matches,
            next_idx: [None; 2],
            hay,
        };

        if let Some(caps) = s.cap_matches.next() {
            let g0 = caps.get(0).unwrap();
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
    type Item = NaturalSortElem<'h>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut elem = None;
        for (i, idx_opt) in self.next_idx.iter_mut().enumerate() {
            if let Some((start, end)) = idx_opt.take() {
                let _ = elem.insert(&self.hay[start..end]);

                if i == 1 {
                    // get next idxs from cap_matches
                    if let Some(caps) = self.cap_matches.next() {
                        let g0 = caps.get(0).unwrap();
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
                Some(NaturalSortElem::Text(elem))
            }
        } else {
            None
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
enum NaturalSortElem<'a> {
    Num(i64),
    Text(&'a str),
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
        let mut record_sorter = RecordSorter::new(0, SortBy::QueryName, threads)?;

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
        record_sorter.sort();
        record_sorter.finalize().expect("Finalize failed");

        // Create a SortedRecordIterator from the sorter.
        let mut sorted_iter =
            SortedRecordIterator::new(record_sorter).expect("Failed to build SortedRecordIterator");

        // Iterate over sorted records and check sort order by query name.
        let mut count = 0;
        let mut prev_name: Option<Vec<u8>> = None;

        let mut header = super::super::Header::from_template(reader.header());
        // The @HD line should be present, with either the SO tag or the GO tag (but not both) specified.

        if let Some(hr) = header.records_mut().first_mut() {
            // let rep = regex::bytes::Regex::new(r"@HD(\t[^\t\n]+)*(\n)?(.*)").unwrap();

            // match rep.captures(&hr) {
            //     Some(caps) => {
            //         let values = caps.get(1).map(|e| e.as_bytes()).unwrap_or(b"");

            //         // replace SO tag
            //         let rep2 = regex::bytes::Regex::new(r"\tSO:([^\t\n]+)").unwrap();
            //         let rep3 = regex::bytes::Regex::new(r"\tSS:([^\t\n]+)").unwrap();
            //         let rep4 = regex::bytes::Regex::new(r"\tGO:([^\t\n]+)").unwrap();

            //         if matches!(sorted_iter.sort_by, SortBy::QueryName) {
            //             // let ss_exists = values.windows(2).position(|b| b == b"SS").is_some();
            //             let ss_exists = rep3.captures(values).is_some();

            //             if !ss_exists {
            //                 let values2 = rep2.replace(
            //                     values,
            //                     format!(
            //                         "\tSO:{}\tSS:{}:natural",
            //                         sorted_iter.sort_by.SO_str(),
            //                         sorted_iter.sort_by.SO_str()
            //                     )
            //                     .as_bytes(),
            //                 );
            //             } else {
            //                 let values2 = rep2.replace(
            //                     values,
            //                     format!(
            //                         "\tSO:{}",
            //                         sorted_iter.sort_by.SO_str(),
            //                     )
            //                     .as_bytes(),
            //                 );

            //                 let values3 = rep3.replace(
            //                     &*values2,
            //                     format!(
            //                         "\tSS:{}:natural",
            //                         sorted_iter.sort_by.SO_str(),
            //                     )
            //                     .as_bytes(),
            //                 );
            //             }
            //         }
            //     }
            //     None => panic!("Invalid header: {}", String::from_utf8_lossy(hr)),
            // }
            let rep = regex::bytes::Regex::new(r"^@HD((?:\t[^\t\n]+)*)(\n)?([\S\s]*)*").unwrap();
            // eprintln!("header: {}", String::from_utf8_lossy(hr));

            let mut new_hd_line = vec![b"@HD".to_vec()];
            match rep.captures(&hr) {
                Some(caps) => {
                    if let Some(values) = caps.get(1) {
                        for value in values.as_bytes().trim_ascii().split(|&b| b == b'\t') {
                            if value.starts_with(b"SO") {
                                new_hd_line.push(
                                    format!("SO:{}", sorted_iter.sort_by.SO_str()).into_bytes(),
                                );

                                if matches!(sorted_iter.sort_by, SortBy::QueryName) {
                                    new_hd_line.push("SS:queryname:natural".into())
                                }
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

                    let mut new_hd_line2 = new_hd_line.join(&b'\t');

                    if caps.get(2).is_some() {
                        new_hd_line2.push(b'\n');
                    };

                    if let Some(values) = caps.get(3) {
                        eprintln!("Added rest header: {} bytes", values.as_bytes().len());

                        new_hd_line2.extend(values.as_bytes());
                    }

                    // eprintln!("Modified header:{}", String::from_utf8_lossy(&new_hd_line2));
                    let _ = std::mem::replace(hr, new_hd_line2);
                }
                None => panic!("Invalid header: {}", String::from_utf8_lossy(hr)),
            }
        }
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

        let def_record = Record::default();

        while let Some(result) = sorted_iter.next() {
            let mut record = result.expect("Error iterating sorted record");
            let cur_name = record.qname().to_vec();
            if let Some(prev) = prev_name {
                assert!(
                    prev <= cur_name,
                    "i:{}, Records are not sorted by query name: {} > {}",
                    count,
                    String::from_utf8_lossy(&prev),
                    String::from_utf8_lossy(&cur_name)
                );
            }
            prev_name = Some(cur_name);
            count += 1;

            assert_ne!(def_record, record);
            // record.set_header(header_view.clone());
            writer.write(&record)?;
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
        let mut record_sorter1 = RecordSorter::new(500, SortBy::QueryName, THREADS)?;
        let mut record_sorter2 = RecordSorter::new(500, SortBy::QueryName, THREADS)?;

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

        record_sorter1.finalize().expect("Finalize sorter1 failed");
        record_sorter2.finalize().expect("Finalize sorter2 failed");

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
                NaturalSortElem::Text("aaa"),
                NaturalSortElem::Num(123),
                NaturalSortElem::Text("aabc."),
                NaturalSortElem::Num(234)
            ],
            NaturalSortElemIter::new(hay).collect::<Vec<_>>()
        );

        let hay = "aaa";

        assert_eq!(
            vec![NaturalSortElem::Text("aaa"),],
            NaturalSortElemIter::new(hay).collect::<Vec<_>>()
        );

        let hay = "1aaa";

        assert_eq!(
            vec![NaturalSortElem::Num(1), NaturalSortElem::Text("aaa"),],
            NaturalSortElemIter::new(hay).collect::<Vec<_>>()
        );

        let hay = "11";

        assert_eq!(
            vec![NaturalSortElem::Num(11),],
            NaturalSortElemIter::new(hay).collect::<Vec<_>>()
        );
    }
}
