use std::cmp::Ordering;

use crate::bam::{self, Read, Record, Writer, Header};
use crate::errors::Result;

/// Sort order supported by the in-memory sorter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SortOrder {
    /// Coordinate sort: reference id, position, then strand. Unmapped reads are placed last
    /// if the reference count is known.
    Coordinate,
    /// Query-name sort using lexicographic (byte) comparison.
    QueryNameLexicographic,
    /// Query-name sort using natural comparison ("a2" < "a10").
    QueryNameNatural,
}

/// Options that control how records are ordered.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SortOptions {
    /// Sort order to apply.
    pub order: SortOrder,
    /// Number of reference sequences. When set, unmapped reads (tid == -1) sort after all
    /// mapped reads, mirroring samtools behaviour.
    pub reference_count: Option<u32>,
}

impl SortOptions {
    /// Create options with the given order and no reference-count hint.
    pub const fn new(order: SortOrder) -> Self {
        Self {
            order,
            reference_count: None,
        }
    }

    /// Create options using the reference-count hint from an existing header view.
    pub fn with_reference_count(order: SortOrder, reference_count: u32) -> Self {
        Self {
            order,
            reference_count: Some(reference_count),
        }
    }

    /// Get the SO: tag value for the @HD line.
    fn so_tag(&self) -> &'static str {
        match self.order {
            SortOrder::Coordinate => "coordinate",
            SortOrder::QueryNameLexicographic | SortOrder::QueryNameNatural => "queryname",
        }
    }

    /// Get the SS: tag value for the @HD line, if applicable.
    fn ss_tag(&self) -> Option<&'static str> {
        match self.order {
            SortOrder::QueryNameNatural => Some("queryname:natural"),
            SortOrder::QueryNameLexicographic => Some("queryname:lexicographical"),
            SortOrder::Coordinate => None,
        }
    }
}

/// Sort all records from `reader` into `writer`. This currently loads all records into memory;
/// for large inputs prefer chunked/external sorting built on top of this API.
///
/// **Note:** This function does not update the header. Use `sort_and_write` if you need
/// automatic header updates.
pub fn sort_reader_to_writer<R: Read>(reader: &mut R, writer: &mut Writer, opts: SortOptions) -> Result<()> {
    let mut records = Vec::new();
    let mut record = Record::new();

    while let Some(res) = reader.read(&mut record) {
        res?;
        records.push(record.clone());
    }

    sort_records_in_place(&mut records, opts);

    for rec in records {
        writer.write(&rec)?;
    }

    Ok(())
}

/// Sort records in memory and write to output with proper header updates.
///
/// This function:
/// 1. Sorts records in place according to `opts`
/// 2. Updates the output header's @HD SO: and SS: tags
/// 3. Writes sorted records to `output_path`
///
/// # Arguments
///
/// * `records` - Mutable slice of records to sort and write
/// * `header` - Template header to use (will be updated with SO:/SS: tags)
/// * `output_path` - Path for the output BAM/SAM/CRAM file
/// * `format` - Output format
/// * `opts` - Sort options
///
/// # Example
///
/// ```no_run
/// use rust_htslib::bam::{self, Read, Format, Header, Record};
/// use rust_htslib::bam::sort::{SortOrder, SortOptions, sort_and_write};
///
/// let mut reader = bam::Reader::from_path("input.bam").unwrap();
/// let header = Header::from_template(reader.header());
/// 
/// // Read and modify records
/// let mut records = Vec::new();
/// let mut record = Record::new();
/// while let Some(r) = reader.read(&mut record) {
///     r.unwrap();
///     // Modify record...
///     record.set_pos(record.pos() + 100);
///     records.push(record.clone());
/// }
///
/// let opts = SortOptions::with_reference_count(
///     SortOrder::Coordinate,
///     reader.header().target_count()
/// );
///
/// sort_and_write(&mut records, &header, "sorted.bam", Format::Bam, opts).unwrap();
/// ```
pub fn sort_and_write<P: AsRef<std::path::Path>>(
    records: &mut [Record],
    header: &Header,
    output_path: P,
    format: bam::Format,
    opts: SortOptions,
) -> Result<()> {
    // Sort in place
    sort_records_in_place(records, opts);

    // Create output header with updated sort order
    let mut out_header = header.clone();
    update_header_sort_order(&mut out_header, opts);

    // Write sorted records
    let mut writer = Writer::from_path(output_path, &out_header, format)?;
    for rec in records.iter() {
        writer.write(rec)?;
    }

    Ok(())
}

/// Update header with appropriate SO: and SS: tags for the sort order.
fn update_header_sort_order(header: &mut Header, opts: SortOptions) {
    // Build new @HD line with SO: and optionally SS:
    let mut hd_record = crate::bam::header::HeaderRecord::new(b"HD");
    hd_record.push_tag(b"VN", "1.6");
    hd_record.push_tag(b"SO", opts.so_tag());
    
    if let Some(ss) = opts.ss_tag() {
        hd_record.push_tag(b"SS", ss);
    }

    // Remove any existing @HD line(s) and add the new one while preserving
    // all other header lines. We split on newlines so headers stored as a
    // single blob (the `from_template` case) are handled correctly.
    let mut lines: Vec<Vec<u8>> = header
        .to_bytes()
        .split(|b| *b == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| line.to_vec())
        .collect();

    lines.retain(|line| !line.starts_with(b"@HD"));
    lines.insert(0, hd_record.to_bytes());

    let records = header.records_mut();
    records.clear();
    records.extend(lines);
}

/// Sort a vector of records in memory, returning the sorted vector.
///
/// This is a convenience wrapper around `sort_records_in_place` that
/// consumes and returns the vector.
pub fn sort_records_in_memory(mut records: Vec<Record>, opts: SortOptions) -> Vec<Record> {
    sort_records_in_place(&mut records, opts);
    records
}

/// Sort the provided records in-place according to the requested options.
pub fn sort_records_in_place(records: &mut [Record], opts: SortOptions) {
    match opts.order {
        SortOrder::Coordinate => {
            let n_ref = opts.reference_count;
            records.sort_by(|a, b| coordinate_cmp(a, b, n_ref));
        }
        SortOrder::QueryNameLexicographic => {
            records.sort_by(|a, b| query_name_cmp(a, b, false));
        }
        SortOrder::QueryNameNatural => {
            records.sort_by(|a, b| query_name_cmp(a, b, true));
        }
    }
}

fn coordinate_cmp(a: &Record, b: &Record, n_ref: Option<u32>) -> Ordering {
    let ta = map_tid(a.tid(), n_ref);
    let tb = map_tid(b.tid(), n_ref);

    match ta.cmp(&tb) {
        Ordering::Equal => {}
        ord => return ord,
    }

    match a.pos().cmp(&b.pos()) {
        Ordering::Equal => {}
        ord => return ord,
    }

    a.is_reverse().cmp(&b.is_reverse())
}

fn map_tid(tid: i32, n_ref: Option<u32>) -> i64 {
    if tid >= 0 {
        tid as i64
    } else if let Some(n) = n_ref {
        n as i64
    } else {
        i64::MAX
    }
}

fn query_name_cmp(a: &Record, b: &Record, natural: bool) -> Ordering {
    let cmp = if natural {
        natural_cmp_bytes(a.qname(), b.qname())
    } else {
        a.qname().cmp(b.qname())
    };

    if cmp != Ordering::Equal {
        return cmp;
    }

    flag_rank(a.flags()).cmp(&flag_rank(b.flags()))
}

fn flag_rank(flag: u16) -> u16 {
    // Matches the ranking used in samtools: READ1, READ2, PRIMARY, SUPPLEMENTARY, SECONDARY.
    ((flag & 0x00c0) << 8) | ((flag & 0x0100) << 3) | ((flag & 0x0800) >> 3)
}

fn natural_cmp_bytes(a: &[u8], b: &[u8]) -> Ordering {
    let mut ia = 0;
    let mut ib = 0;

    while ia < a.len() && ib < b.len() {
        let ca = a[ia];
        let cb = b[ib];

        let da = ca.is_ascii_digit();
        let db = cb.is_ascii_digit();

        if da && db {
            let (num_a, next_a) = read_digits(a, ia);
            let (num_b, next_b) = read_digits(b, ib);

            let cmp = compare_digit_chunks(num_a, num_b);
            if cmp != Ordering::Equal {
                return cmp;
            }

            ia = next_a;
            ib = next_b;
        } else {
            match ca.cmp(&cb) {
                Ordering::Equal => {
                    ia += 1;
                    ib += 1;
                }
                ord => return ord,
            }
        }
    }

    match (ia == a.len(), ib == b.len()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        (false, false) => a.len().cmp(&b.len()),
    }
}

fn read_digits(bytes: &[u8], start: usize) -> (&[u8], usize) {
    let mut end = start;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
    }
    (&bytes[start..end], end)
}

fn compare_digit_chunks(a: &[u8], b: &[u8]) -> Ordering {
    let a_trimmed = trim_leading_zeros(a);
    let b_trimmed = trim_leading_zeros(b);

    match a_trimmed.len().cmp(&b_trimmed.len()) {
        Ordering::Equal => {}
        ord => return ord,
    }

    match a_trimmed.cmp(b_trimmed) {
        Ordering::Equal => Ordering::Equal,
        ord => ord,
    }
}

fn trim_leading_zeros(digits: &[u8]) -> &[u8] {
    let mut idx = 0;
    while idx < digits.len() && digits[idx] == b'0' {
        idx += 1;
    }

    if idx == digits.len() {
        // All zeros -> treat as single zero to keep ordering stable.
        &digits[digits.len().saturating_sub(1)..]
    } else {
        &digits[idx..]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::htslib;

    #[test]
    fn natural_cmp_matches_numeric_intuition() {
        assert_eq!(natural_cmp_bytes(b"a2", b"a10"), Ordering::Less);
        assert_eq!(natural_cmp_bytes(b"a10", b"a2"), Ordering::Greater);
        assert_eq!(natural_cmp_bytes(b"a01", b"a1"), Ordering::Equal);
        assert_eq!(natural_cmp_bytes(b"a001", b"a1"), Ordering::Equal);
        assert_eq!(natural_cmp_bytes(b"x9y", b"x9z"), Ordering::Less);
    }

    #[test]
    fn flag_rank_orders_read1_before_read2() {
        let r1 = flag_rank(htslib::BAM_FREAD1 as u16);
        let r2 = flag_rank(htslib::BAM_FREAD2 as u16);
        assert!(r1 < r2);
    }

    #[test]
    fn coordinate_sort_sends_unmapped_to_end_when_reference_count_known() {
        let opts = SortOptions::with_reference_count(SortOrder::Coordinate, 1);

        let mut mapped = Record::new();
        mapped.set_tid(0);
        mapped.set_pos(5);

        let mut unmapped = Record::new();
        unmapped.set_tid(-1);

        let mut records = vec![unmapped, mapped.clone()];
        sort_records_in_place(&mut records, opts);

        assert_eq!(records[0].tid(), mapped.tid());
        assert_eq!(records[1].tid(), -1);
    }

    #[test]
    fn sort_records_in_memory_convenience() {
        let mut r1 = Record::new();
        r1.set_tid(0);
        r1.set_pos(100);
        r1.set_qname(b"read1");

        let mut r2 = Record::new();
        r2.set_tid(0);
        r2.set_pos(50);
        r2.set_qname(b"read2");

        let records = vec![r1.clone(), r2.clone()];
        let sorted = sort_records_in_memory(records, SortOptions::new(SortOrder::Coordinate));

        assert_eq!(sorted[0].pos(), 50);
        assert_eq!(sorted[1].pos(), 100);
    }
}
