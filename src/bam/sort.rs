use tempfile::NamedTempFile;

use super::{Read as _, record::Record};
use anyhow::{Error, anyhow};
use std::{cmp::Ordering, collections::BinaryHeap, marker::PhantomData};

#[derive(Debug, Clone, Copy)]
pub enum SortBy {
    QueryName,
    Coordinate,
}

impl SortBy {
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
}

impl RecordSorter {
    pub fn new(chunk_size: usize, sort_by: SortBy) -> Self {
        Self {
            chunk_size,
            temp_files: Vec::new(),
            sort_by,
            record_vec: Vec::with_capacity(chunk_size),
        }
    }

    pub fn push_record(&mut self, record: Record) -> Result<(), anyhow::Error> {
        if self.record_vec.len() == self.chunk_size {
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
        self.record_vec.sort_by(|a, b| self.sort_by.compare(a, b));
    }

    fn dump(&mut self) -> Result<(), anyhow::Error> {
        let tmp_file = NamedTempFile::new()?;
        let header = super::Header::new();
        let mut writer = super::Writer::from_path(tmp_file.path(), &header, super::Format::Bam)?;

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
                assert!(self.record_vec.len() <= self.chunk_size);
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

pub trait SortKey {
    fn compare(lhs: &Record, rhs: &Record) -> Ordering;

    fn equal(lhs: &Record, rhs: &Record) -> bool;
}

struct SortKeyCoord;

impl SortKey for SortKeyCoord {
    fn compare(lhs: &Record, rhs: &Record) -> Ordering {
        lhs.qname().cmp(rhs.qname())
    }

    fn equal(lhs: &Record, rhs: &Record) -> bool {
        lhs.qname() == rhs.qname()
    }
}

struct SortKeyName;

impl SortKey for SortKeyName {
    fn compare(lhs: &Record, rhs: &Record) -> Ordering {
        match lhs.tid().cmp(&rhs.tid()) {
            Ordering::Equal => {}
            oth => return oth,
        }

        lhs.pos().cmp(&rhs.pos())
    }

    fn equal(lhs: &Record, rhs: &Record) -> bool {
        lhs.tid() == rhs.tid() && lhs.pos() == rhs.pos()
    }
}

struct HeapItem<S>
where
    S: SortKey,
{
    record: Record,
    chunk_idx: usize,
    _marker: PhantomData<S>,
}

impl<S> Ord for HeapItem<S>
where
    S: SortKey,
{
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse the ordering (since BinaryHeap is a max-heap) for ascending order.
        S::compare(&other.record, &self.record)
    }
}

impl<S> PartialOrd for HeapItem<S>
where
    S: SortKey,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(other.cmp(self))
    }
}

impl<S> PartialEq for HeapItem<S>
where
    S: SortKey,
{
    fn eq(&self, other: &Self) -> bool {
        S::equal(&other.record, &self.record)
    }
}

impl<S> Eq for HeapItem<S> where S: SortKey {}

pub struct SortedRecordIterator<S: SortKey> {
    mem_iter: std::vec::IntoIter<Record>,
    chunk_readers: Vec<super::Reader>,
    binary_heap: BinaryHeap<HeapItem<S>>,
    from_mem: bool,
}

impl<S: SortKey> SortedRecordIterator<S> {
    pub fn new(mut record_sorter: RecordSorter) -> Result<Self, Error> {
        record_sorter.finalize()?;

        let mut chunk_readers = record_sorter
            .temp_files
            .iter()
            .map(|p| super::Reader::from_path(p.path()).map_err(Error::from))
            .collect::<Result<Vec<_>, Error>>()?;

        let from_mem = chunk_readers.is_empty();

        assert!(!(record_sorter.record_vec.len() > 0 && chunk_readers.len() > 0));

        let iter = record_sorter.record_vec.into_iter();

        let mut binary_heap = BinaryHeap::with_capacity(chunk_readers.len());

        for (chunk_idx, reader) in chunk_readers.iter_mut().enumerate() {
            let mut next_record = Record::new();

            if let Some(rr) = reader.read(&mut next_record) {
                rr?;

                binary_heap.push(HeapItem {
                    record: next_record,
                    chunk_idx,
                    _marker: Default::default(),
                });
            }
        }

        Ok(Self {
            mem_iter: iter,
            chunk_readers,
            binary_heap,
            from_mem,
        })
    }
}

impl<S: SortKey> Iterator for SortedRecordIterator<S> {
    type Item = Result<Record, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let elem = if self.from_mem {
            self.mem_iter.next()
        } else {
            if let Some(HeapItem {
                record,
                chunk_idx,
                ..
                // _marker,
                // sort_by,
            }) = self.binary_heap.pop()
            {
                let reader = &mut self.chunk_readers[chunk_idx];
                let mut next_record = Record::new();

                if let Some(rr) = reader.read(&mut next_record) {
                    match rr {
                        Ok(_) => {}
                        Err(err) => return Some(Err(err).map_err(Error::from)),
                    }

                    self.binary_heap.push(HeapItem {
                        record: next_record,
                        chunk_idx,
                        _marker: Default::default(),
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
