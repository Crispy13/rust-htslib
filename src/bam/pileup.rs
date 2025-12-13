// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::fmt;
use std::iter;
use std::slice;

use crate::htslib;

use crate::bam;
use crate::bam::record;
use crate::errors::{Error, Result};

/// Iterator over alignments of a pileup.
pub type Alignments<'a> = iter::Map<
    slice::Iter<'a, htslib::bam_pileup1_t>,
    fn(&'a htslib::bam_pileup1_t) -> Alignment<'a>,
>;

/// A pileup over one genomic position.
#[derive(Debug)]
pub struct Pileup {
    inner: *const htslib::bam_pileup1_t,
    depth: u32,
    tid: u32,
    pos: u32,
}

impl Pileup {
    pub fn tid(&self) -> u32 {
        self.tid
    }

    pub fn pos(&self) -> u32 {
        self.pos
    }

    pub fn depth(&self) -> u32 {
        self.depth
    }

    pub fn alignments(&self) -> Alignments<'_> {
        self.inner().iter().map(Alignment::new)
    }

    fn inner(&self) -> &[htslib::bam_pileup1_t] {
        unsafe {
            slice::from_raw_parts(
                self.inner as *mut htslib::bam_pileup1_t,
                self.depth as usize,
            )
        }
    }
}

/// An aligned read in a pileup.
pub struct Alignment<'a> {
    inner: &'a htslib::bam_pileup1_t,
}

impl fmt::Debug for Alignment<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Alignment")
    }
}

impl<'a> Alignment<'a> {
    pub fn new(inner: &'a htslib::bam_pileup1_t) -> Self {
        Alignment { inner }
    }

    /// Position within the read. None if either `is_del` or `is_refskip`.
    pub fn qpos(&self) -> Option<usize> {
        if self.is_del() || self.is_refskip() {
            // there is no alignment position in such a case
            None
        } else {
            Some(self.inner.qpos as usize)
        }
    }

    /// position of the read base at the pileup site, 0-based.
    ///
    /// If the current position is a deletion, returns the next
    /// aligned base.
    pub fn qpos_or_next(&self) -> usize {
        self.inner.qpos as usize
    }

    /// Insertion, deletion (with length) if indel starts at next base or None otherwise.
    pub fn indel(&self) -> Indel {
        match self.inner.indel {
            len if len < 0 => Indel::Del(-len as u32),
            len if len > 0 => Indel::Ins(len as u32),
            _ => Indel::None,
        }
    }

    /// Whether there is a deletion in the alignment at this position.
    pub fn is_del(&self) -> bool {
        self.inner.is_del() != 0
    }

    /// Whether the alignment starts at this position.
    pub fn is_head(&self) -> bool {
        self.inner.is_head() != 0
    }

    /// Whether the alignment ends at this position.
    pub fn is_tail(&self) -> bool {
        self.inner.is_tail() != 0
    }

    /// Whether this position is marked as refskip in the CIGAR string.
    pub fn is_refskip(&self) -> bool {
        self.inner.is_refskip() != 0
    }

    /// The corresponding record.
    pub fn record(&self) -> record::Record {
        record::Record::from_inner(self.inner.b)
    }
}

#[derive(PartialEq, Eq, Debug, Copy, Clone, Hash)]
pub enum Indel {
    Ins(u32),
    Del(u32),
    None,
}

#[derive(Copy, Clone)]
pub struct PileupOption {
    /// Maximum value= `i32::MAX`
    pub max_depth: i32,
    pub ignore_overlaps: bool,
    //TODO: stepper? of pysam.
}

impl Default for PileupOption {
    fn default() -> Self {
        Self {
            max_depth: 8000,
            ignore_overlaps: true,
        }
    }
}

/// Iterator over pileups.
#[derive(Debug)]
pub struct Pileups<'a, R: bam::Read> {
    #[allow(dead_code)]
    reader: &'a mut R,
    itr: htslib::bam_plp_t,
}

impl<'a, R: bam::Read> Pileups<'a, R> {
    pub fn new(reader: &'a mut R, itr: htslib::bam_plp_t) -> Self {
        Pileups { reader, itr }
    }

    pub fn with_option(reader: &'a mut R, itr: htslib::bam_plp_t, option: PileupOption) -> Self {
        let mut s = Pileups { reader, itr };

        let PileupOption {
            max_depth,
            ignore_overlaps,
        } = option;
        s.set_max_depth(max_depth as u32);

        if ignore_overlaps {
            s.ignore_overlaps();
        }

        s
    }

    /// Warning: because htslib internally uses signed integer for depth this method
    /// will panic if `depth` exceeds `i32::MAX`.
    pub fn set_max_depth(&mut self, depth: u32) {
        if depth > i32::MAX as u32 {
            panic!(
                "Maximum value for pileup depth is {} but {} was provided",
                i32::MAX,
                depth
            )
        }
        let intdepth = depth as i32;
        unsafe {
            htslib::bam_plp_set_maxcnt(self.itr, intdepth);
        }
    }

    pub fn ignore_overlaps(&mut self) {
        unsafe {
            if htslib::bam_plp_init_overlaps(self.itr) != 0 {
                panic!("htslib::bam_plp_init_overlaps function failed.")
            }
        }
    }
}

impl<R: bam::Read> Iterator for Pileups<'_, R> {
    type Item = Result<Pileup>;

    #[allow(clippy::match_bool)]
    fn next(&mut self) -> Option<Result<Pileup>> {
        let (mut tid, mut pos, mut depth) = (0i32, 0i32, 0i32);
        let inner = unsafe { htslib::bam_plp_auto(self.itr, &mut tid, &mut pos, &mut depth) };

        match inner.is_null() {
            true if depth == -1 => Some(Err(Error::BamPileup)),
            true => None,
            false => Some(Ok(Pileup {
                inner,
                depth: depth as u32,
                tid: tid as u32,
                pos: pos as u32,
            })),
        }
    }
}

impl<R: bam::Read> Drop for Pileups<'_, R> {
    fn drop(&mut self) {
        unsafe {
            htslib::bam_plp_reset(self.itr);
            htslib::bam_plp_destroy(self.itr);
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::bam;
    use crate::bam::Read;

    #[test]
    fn test_max_pileup() {
        let mut bam = bam::Reader::from_path("test/test.bam").unwrap();
        let mut p = bam.pileup();
        p.set_max_depth(0u32);
        p.set_max_depth(800u32);
    }

    #[test]
    fn test_ignore_overlaps() {
        let mut bam = bam::Reader::from_path("test/test.bam").unwrap();
        let mut p = bam.pileup();
        p.ignore_overlaps();
    }

    #[test]
    #[should_panic]
    fn test_max_pileup_to_high() {
        let mut bam = bam::Reader::from_path("test/test.bam").unwrap();
        let mut p = bam.pileup();
        p.set_max_depth((i32::MAX as u32) + 1);
    }
}
