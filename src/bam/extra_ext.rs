use std::{collections::HashMap, convert::TryFrom};

use crate::{bam::record::AuxArray, errors::Error};
use linear_map::LinearMap;

use super::{
    record::{Aux, Cigar, CigarString, CigarStringView},
    HeaderView, Record,
};

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

const TAGS_TO_REVERSE_COMPLEMENT: [&'static str; 2] = ["E2", "SQ"];
const TAGS_TO_REVERSE: [&'static str; 2] = ["OQ", "U2"];

/**
 * If a read has this reference index, it is unaligned, but not all unaligned reads have
 * this reference index (see above).
 */
const NO_ALIGNMENT_REFERENCE_INDEX: i32 = -1;

const NO_ALIGNMENT_START: i64 = 0;

const NO_ALIGNMENT_CIGAR: &str = "*";

const NO_MAPPING_QUALITY: i32 = 0;

pub trait RecordExt {
    fn setter(&mut self) -> RecordModifier;
    fn reverse_complement(
        &mut self,
        tags_to_rev_comp: &[&str],
        tags_to_reverse: &[&str],
        inplace: bool,
    );
    fn is_fr_pair(&self) -> bool;
    fn get_read_group(&self) -> Result<SAMReadGroupRecord, Error>;

    fn get_str_aux(&self, tag: &[u8]) -> Result<&str, Error>;

    fn get_unclipped_start(&self) -> u32;

    fn get_unclipped_end(&self) -> u32;

    fn get_mate_unclipped_start(&self) -> Result<u32, Error>;

    fn get_mate_unclipped_end(&self) -> Result<u32, Error>;

    /**
     * the read is either a PCR duplicate or an optical duplicate.
     */
    fn set_duplicate_read_flag(&mut self, flag: bool);

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

    fn set_duplicate_read_flag(&mut self, flag: bool) {
        const DUP_FLAG: u16 = 0x400;

        if flag {
            self.set_flags(self.flags() | DUP_FLAG)
        } else {
            self.set_flags(self.flags() & !DUP_FLAG)
        }
    }

    fn is_secondary_or_supplementary(&self) -> bool {
        self.is_secondary() || self.is_supplementary()
    }

    fn is_fr_pair(&self) -> bool {
        match self.read_pair_orientation() {
            bio_types::sequence::SequenceReadPairOrientation::F1R2
            | bio_types::sequence::SequenceReadPairOrientation::F2R1 => true,
            _ => false,
        }
    }

    fn reverse_complement(
        &mut self,
        tags_to_rev_comp: &[&str],
        tags_to_reverse: &[&str],
        inplace: bool,
    ) {
        let mut read_bases = self.seq().encoded.to_vec();
        reverse_complement(&mut read_bases);

        let mut qualities = self.qual().to_vec();
        qualities.reverse();

        self.set(
            &self.qname().to_vec(),
            Some(&self.cigar()),
            &read_bases,
            &qualities,
        );

        // Deal with tags that need to be reverse complemented
        for &tag in tags_to_rev_comp {
            if let Ok(mut value) = self.aux(tag.as_bytes()) {
                if let Aux::ArrayU8(arr) = &value {
                    let mut arr = arr.iter().collect::<Vec<_>>();
                    reverse_complement(&mut arr);

                    value = Aux::ArrayU8(AuxArray::from_bytes(&arr));

                    remove_and_push!(self, tag, value);
                } else if let Aux::String(s) = &value {
                    //SequenceUtil.reverseComplement is in-place for bytes but copies Strings since they are immutable.
                    let mut s = (*s).to_owned().into_bytes();
                    reverse_complement(&mut s);

                    value = Aux::String(std::str::from_utf8(&s).unwrap());
                    remove_and_push!(self, tag, value);
                } else {
                    panic!("Don't know how to reverse complement: {:?}", value);
                }
            }
        }

        // Deal with tags that needed to just be reversed
        for &tag in tags_to_reverse {
            if let Ok(mut value) = self.aux(tag.as_bytes()) {
                match value {
                    Aux::String(s) => {
                        let vs = String::from_utf8(s.bytes().rev().collect::<Vec<_>>()).unwrap();
                        value = Aux::String(&vs);

                        remove_and_push!(self, tag, value);
                    }
                    Aux::ArrayU8(arr) => {
                        reverse_aux_array_and_set!(Aux::ArrayU8, arr, value, self, tag)
                    }
                    Aux::ArrayI16(arr) => {
                        reverse_aux_array_and_set!(Aux::ArrayI16, arr, value, self, tag)
                    }
                    Aux::ArrayI32(arr) => {
                        reverse_aux_array_and_set!(Aux::ArrayI32, arr, value, self, tag)
                    }
                    Aux::ArrayFloat(arr) => {
                        reverse_aux_array_and_set!(Aux::ArrayFloat, arr, value, self, tag)
                    }
                    _ => panic!("Don't know how to reverse: {:?}", value),
                }
            }
        }
    }

    fn setter(&mut self) -> RecordModifier {
        RecordModifier::new(self)
    }
}

pub struct RecordModifier<'r> {
    rec: &'r mut Record,
    qname: Option<Vec<u8>>,
    cigar: Option<CigarString>,
    seq: Option<Vec<u8>>,
    quals: Option<Vec<u8>>,
}

macro_rules! record_modifier_setter {
    ($func_name:ident, $self:ident, $field:ident, $ty:ty) => {
        pub fn $func_name(mut $self, value:$ty) -> Self {
            $self.$field = Some(value);
            $self
        }
    };
}

impl<'r> RecordModifier<'r> {
    fn new(rec: &'r mut Record) -> Self {
        Self {
            rec,
            qname: None,
            cigar: None,
            seq: None,
            quals: None,
        }
    }

    record_modifier_setter!(set_qname, self, qname, Vec<u8>);
    record_modifier_setter!(set_cigar, self, cigar, CigarString);
    record_modifier_setter!(set_seq, self, seq, Vec<u8>);
    record_modifier_setter!(set_quals, self, quals, Vec<u8>);

    pub fn modify_record(self) {
        let rec = self.rec;

        let new_qname = self.qname.unwrap_or_else(|| rec.qname().to_vec());
        let new_cigar = self.cigar;
        let new_seq = self.seq.unwrap_or_else(|| rec.seq().encoded.to_vec());
        let new_quals = self.quals.unwrap_or_else(|| rec.qual().to_vec());

        rec.set(
            &new_qname,
            new_cigar.as_ref(),
            &new_seq,
            &new_quals
        );
    }
}

/// remove tag and its original value and push a new value to the tag.
macro_rules! remove_and_push {
    ($rec:ident, $tag:ident, $value:ident) => {
        $rec.remove_aux($tag.as_bytes()).unwrap();
        $rec.push_aux($tag.as_bytes(), $value).unwrap();
    };
}

pub(crate) use remove_and_push;

macro_rules! reverse_aux_array_and_set {
    ($array_variant:path, $arr:ident, $value:ident, $self:ident, $tag:ident) => {{
        let mut arr = $arr.iter().collect::<Vec<_>>();
        arr.reverse();

        $value = $array_variant(AuxArray::from(&arr));

        remove_and_push!($self, $tag, $value);
    }};
}

pub(crate) use reverse_aux_array_and_set;

/// convert a sequence into its reverse complement. (in-place)
pub(crate) fn reverse_complement(seq: &mut [u8]) {
    // let mut s = seq.to_string().into_bytes(); //seq.chars().collect::<Vec<char>>();
    let mut s_iter = seq.iter_mut();

    loop {
        match (s_iter.next(), s_iter.next_back()) {
            (Some(f), Some(b)) => {
                (*f, *b) = (get_complement_base(b), get_complement_base(f));
            }
            (Some(f), None) => {
                *f = get_complement_base(f);
                break;
            }
            (None, None) => {
                break;
            }
            (None, Some(b)) => {
                // this is unreachable.
                unreachable!();
                break;
            }
        }
    }

    // String::from_utf8(seq).unwrap()
}

#[inline]
fn get_complement_base(base: &u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }
}

/**
 * Strip mapping information from a SAMRecord.
 * <p>
 * WARNING: by clearing the secondary and supplementary flags,
 * this may have the affect of producing multiple distinct records with the
 * same read name and flags, which may lead to invalid SAM/BAM output.
 * Callers of this method should make sure to deal with this issue.
 */
pub fn make_read_unmapped(rec: &mut Record) {
    if rec.is_reverse() {
        rec.reverse_complement(&TAGS_TO_REVERSE_COMPLEMENT, &TAGS_TO_REVERSE, true);
        rec.unset_reverse();
    }

    rec.unset_duplicate();
    rec.set_tid(NO_ALIGNMENT_REFERENCE_INDEX);
    rec.set_pos(NO_ALIGNMENT_START);
    rec.set(
        &rec.qname().to_vec(),
        Some(&CigarString::try_from(NO_ALIGNMENT_CIGAR.as_bytes()).unwrap()),
        &rec.seq().encoded.to_vec(),
        &rec.qual().to_vec(),
    );

    rec.set_insert_size(0);
    rec.unset_secondary();
    rec.unset_supplementary();
    rec.unset_proper_pair();
    rec.set_unmapped();
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

pub trait CigarExt {
    fn is_clipping(&self) -> bool;
    fn modify_len(&mut self, len: u32);
    fn consumes_reference_bases(&self) -> bool;
    fn consumes_read_bases(&self) -> bool;
    /** Returns true if the operator is a M, a X or a EQ */
    fn is_alignment(&self) -> bool;
}

impl CigarExt for Cigar {
    fn is_alignment(&self) -> bool {
        match self {
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => true,
            _ => false,
        }
    }

    fn consumes_read_bases(&self) -> bool {
        match self {
            Cigar::Match(_)
            | Cigar::Ins(_)
            | Cigar::SoftClip(_)
            | Cigar::Equal(_)
            | Cigar::Diff(_) => true,
            _ => false,
        }
    }

    fn consumes_reference_bases(&self) -> bool {
        match self {
            Cigar::Match(_)
            | Cigar::Del(_)
            | Cigar::RefSkip(_)
            | Cigar::Equal(_)
            | Cigar::Diff(_) => true,
            _ => false,
        }
    }

    fn modify_len(&mut self, len: u32) {
        match self {
            Cigar::Match(ref mut l) => *l = len,
            Cigar::Ins(ref mut l) => *l = len,
            Cigar::Del(ref mut l) => *l = len,
            Cigar::RefSkip(ref mut l) => *l = len,
            Cigar::SoftClip(ref mut l) => *l = len,
            Cigar::HardClip(ref mut l) => *l = len,
            Cigar::Pad(ref mut l) => *l = len,
            Cigar::Equal(ref mut l) => *l = len,
            Cigar::Diff(ref mut l) => *l = len,
        }
    }

    fn is_clipping(&self) -> bool {
        match self {
            Cigar::SoftClip(_) | Cigar::HardClip(_) => true,
            _ => false,
        }
    }
}

pub trait CigarStringExt {
    /** Coalesces adjacent operators of the same type into single operators. */
    fn coalesce(&self) -> CigarString;
}

impl CigarStringExt for CigarString {
    fn coalesce(&self) -> CigarString {
        let mut iter = self.iter();
        let mut builder = Vec::with_capacity(self.len());

        while let Some(elem) = iter.next() {
            let same_len = iter
                .by_ref()
                .take_while(|c| std::mem::discriminant(elem) == std::mem::discriminant(c))
                .fold(elem.len(), |a, x| a + x.len());

            let mut same_elem = elem.clone();
            same_elem.modify_len(same_len);

            builder.push(same_elem);
        }

        CigarString::from(builder)
    }
}
