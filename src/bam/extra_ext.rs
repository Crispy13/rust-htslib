use std::collections::HashMap;

use linear_map::LinearMap;
use crate::errors::Error;

use super::{record::{Aux, Cigar, CigarStringView}, HeaderView, Record};

/**
 * Header information about a read group.
 *
 * This struct is based on SAMReadGroupRecord of htsjdk(Java).
 */
#[derive(Debug, Clone)]
pub struct SAMReadGroupRecord<'m> {
    m_read_group_id: &'m str,
    inner_map: &'m LinearMap<String, String>,
}

/* Platform values for the @RG-PL tag */
enum PlatformValue {
    /// Deprecated. Use `PlatformValue::DNBSEQ` instead.
    BGI,

    /** Capillary */
    CAPILLARY,

    /** MGI/BGI */
    DNBSEQ,

    /** Element Biosciences */
    ELEMENT,

    /** Helicos Biosciences */
    HELICOS,

    /** Illumina */
    ILLUMINA,

    /** Iontorrent */
    IONTORRENT,

    /** 454 Life Sciences */
    LS454,

    /** Oxford Nanopore */
    ONT,

    /** Deprecated. OTHER is not an official value.  It is recommended to omit PL if it is not in this list or is unknown. */
    OTHER,

    /** Pacific Biotechnology */
    PACBIO,

    /** Singular Genomics */
    SINGULAR,

    /** Life Technologies */
    SOLID,

    /** Ultima Genomics */
    ULTIMA,
}
impl<'m> SAMReadGroupRecord<'m> {
    pub(crate) const RG: &'static str = "RG";
    pub(super) const READ_GROUP_ID_TAG: &'static str = "ID";
    const SEQUENCING_CENTER_TAG: &'static str = "CN";
    const DESCRIPTION_TAG: &'static str = "DS";
    const DATE_RUN_PRODUCED_TAG: &'static str = "DT";
    const FLOW_ORDER_TAG: &'static str = "FO";
    const KEY_SEQUENCE_TAG: &'static str = "KS";
    const LIBRARY_TAG: &'static str = "LB";
    const PROGRAM_GROUP_TAG: &'static str = "PG";
    const PREDICTED_MEDIAN_INSERT_SIZE_TAG: &'static str = "PI";
    const PLATFORM_TAG: &'static str = "PL";
    const PLATFORM_MODEL_TAG: &'static str = "PM";
    const PLATFORM_UNIT_TAG: &'static str = "PU";
    const READ_GROUP_SAMPLE_TAG: &'static str = "SM";
    const BARCODE_TAG: &'static str = "BC";

    pub fn from_header_map(
        read_group_map: &'m LinearMap<String, String>,
    ) -> SAMReadGroupRecord<'m> {
        Self {
            m_read_group_id: read_group_map.get(Self::READ_GROUP_ID_TAG).unwrap(),
            inner_map: read_group_map,
        }
    }

    pub fn get_read_group_id(&self) -> &'m str {
        self.m_read_group_id
    }

    pub fn get_sample(&self) -> Option<&'m str> {
        self.inner_map
            .get(Self::READ_GROUP_SAMPLE_TAG)
            .map(|e| e.as_str())
    }

    pub fn get_library(&self) -> Option<&'m str> {
        self.inner_map.get(Self::LIBRARY_TAG).map(|e| e.as_str())
    }

    pub fn get_platform_unit(&self) -> Option<&'m str> {
        self.inner_map
            .get(Self::PLATFORM_UNIT_TAG)
            .map(|e| e.as_str())
    }

    pub fn get_platform(&self) -> Option<&'m str> {
        self.inner_map.get(Self::PLATFORM_TAG).map(|e| e.as_str())
    }

    pub fn get_flow_order(&self) -> Option<&'m str> {
        self.inner_map.get(Self::FLOW_ORDER_TAG).map(|e| e.as_str())
    }

    pub fn get_key_sequence(&self) -> Option<&'m str> {
        self.inner_map
            .get(Self::KEY_SEQUENCE_TAG)
            .map(|e| e.as_str())
    }

    pub fn id(&self) -> Option<&'m str> {
        self.inner_map
            .get(Self::READ_GROUP_ID_TAG)
            .map(|e| e.as_str())
    }
}

pub trait RecordExt {
    fn get_read_group(&self) -> Result<SAMReadGroupRecord, Error>;

    fn get_str_aux(&self, tag: &[u8]) -> Result<&str, Error>;

    fn get_unclipped_start(&self) -> u32;

    fn get_unclipped_end(&self) -> u32;

    fn get_mate_unclipped_start(&self) -> Result<u32, Error>;

    fn get_mate_unclipped_end(&self) -> Result<u32, Error>;

    /**
     * the read is either a PCR duplicate or an optical duplicate.
     */
    fn set_duplicate_read_flag(&mut self, flag:bool);

    fn is_secondary_or_supplementary(&self) -> bool;
}

fn get_unclipped_start(alignment_start: i64, cigar: CigarStringView) -> u32 {
    let mut unclipped_start = alignment_start as u32;

    for cigar_element in cigar.iter() {
        match cigar_element {
            Cigar::SoftClip(len) | Cigar::HardClip(len) => {
                unclipped_start -= *len;
            }
            _ => break,
        }
    }

    unclipped_start
}

fn get_unclipped_end(alignment_end: i64, cigar: CigarStringView) -> u32 {
    let mut unclipped_end = alignment_end as u32;

    for cigar_element in cigar.iter().rev() {
        match cigar_element {
            Cigar::SoftClip(len) | Cigar::HardClip(len) => {
                unclipped_end += *len;
            }
            _ => break,
        }
    }

    unclipped_end
}

impl RecordExt for Record {
    /// # Errors
    /// Failed when:
    /// - `self.mate_cigar()` fails.
    fn get_mate_unclipped_start(&self) -> Result<u32, Error> {
        let mate_cigar = self.mate_cigar()?;

        Ok(get_unclipped_start(mate_cigar.pos(), mate_cigar))
    }

    /// # Errors
    /// Failed when:
    /// - `self.mate_cigar()` fails.
    fn get_mate_unclipped_end(&self) -> Result<u32, Error> {
        let mate_cigar = self.mate_cigar()?;

        Ok(get_unclipped_end(mate_cigar.end_pos(), mate_cigar))
    }

    fn get_unclipped_start(&self) -> u32 {
        let cigar = self.cigar();
        get_unclipped_start(cigar.pos(), cigar)
    }

    fn get_unclipped_end(&self) -> u32 {
        let cigar = self.cigar();

        get_unclipped_end(cigar.end_pos(), cigar)
    }

    /// .
    ///
    /// # Errors
    ///
    /// This function will return an error if:
    /// 1. "RG" tag's type of record is other than string type.
    /// 2. inner `Record::read_aux_field()` failed.
    /// 3. No "RG" in header.
    /// 4. multiple read group data for a single read group id in the header.
    fn get_read_group(&self) -> Result<SAMReadGroupRecord, Error> {
        let rg = match self.aux(SAMReadGroupRecord::RG.as_bytes()) {
            Ok(aux) => match aux {
                Aux::String(v) => v,
                _ => Err(Error::BamAuxParsingError)?,
            },
            Err(err) => Err(err)?,
        };

        let rg_info_map = self
            .header
            .as_ref()
            .unwrap()
            .header_map()
            .inner
            .get(SAMReadGroupRecord::RG)
            .ok_or_else(|| Error::BamUndefinedTag {
                tag: SAMReadGroupRecord::RG.to_string(),
            })?
            .iter()
            .filter_map(|e| {
                e.get(SAMReadGroupRecord::READ_GROUP_ID_TAG).and_then(|v| {
                    if v == rg {
                        Some(e)
                    } else {
                        None
                    }
                })
            })
            .collect::<Vec<_>>();

        // check exactly one rg info exists.
        match rg_info_map.len() {
            1 => (),
            _ => Err(Error::BamUndefinedTag {
                tag: SAMReadGroupRecord::RG.to_string(),
            })?,
        }

        let rg_info = rg_info_map.first().unwrap().clone();

        Ok(SAMReadGroupRecord {
            m_read_group_id: rg,
            inner_map: rg_info,
        })
    }

    fn get_str_aux(&self, tag: &[u8]) -> Result<&str, Error> {
        self.aux(tag).and_then(|aux| match aux {
            Aux::String(s) => Ok(s),
            _ => Err(Error::BamAuxParsingError)?,
        })
    }
    
    fn set_duplicate_read_flag(&mut self, flag:bool) {
        const DUP_FLAG:u16=0x400;
        
        if flag {
            self.set_flags(self.flags() | DUP_FLAG)
        } else {
            self.set_flags(self.flags() & !DUP_FLAG)
        }
    }

    fn is_secondary_or_supplementary(&self) -> bool {
        self.is_secondary() || self.is_supplementary()
    }
}

#[derive(Debug, Clone)]
pub struct HeaderMap {
    pub(super) inner: HashMap<String, Vec<LinearMap<String, String>>>,
}

impl HeaderMap {
    const HD: &'static str = "HD";
    const SO: &'static str = "SO";

    pub fn get_read_groups(&self) -> Option<Vec<SAMReadGroupRecord>> {
        self.inner.get(SAMReadGroupRecord::RG).and_then(|v| {
            Some(
                v.into_iter()
                    .map(|e| SAMReadGroupRecord::from_header_map(e))
                    .collect::<Vec<_>>(),
            )
        })
    }

    pub fn get_sort_order(&self) -> Option<&str> {
        self.inner
            .get(Self::HD)
            .and_then(|e| e.first())
            .and_then(|e| e.get(Self::SO))
            .and_then(|e| Some(e.as_str()))
    }
}

pub trait HeaderViewExt {
    /// Get `SAMReadGroupRecord` per read group ID.
    ///
    /// # Errors
    /// Return Error when the header does not have "RG" tag.
    ///
    fn get_read_groups(&self) -> Result<Vec<SAMReadGroupRecord>, Error>;
}
impl HeaderViewExt for HeaderView {
    fn get_read_groups(&self) -> Result<Vec<SAMReadGroupRecord>, Error> {
        let rg_info_map = self
            .header_map
            .inner
            .get(SAMReadGroupRecord::RG)
            .ok_or_else(|| Error::BamUndefinedTag {
                tag: SAMReadGroupRecord::RG.to_string(),
            })?
            .iter()
            .map(|lm| SAMReadGroupRecord::from_header_map(lm))
            .collect::<Vec<_>>();

        Ok(rg_info_map)
    }
}
    