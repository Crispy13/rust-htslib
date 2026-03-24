#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use rust_htslib::bam::{self, Read};

#[test]
fn test_bam_indexed_reader_with_mimalloc() {
    let path = format!("{}/test/test.bam", env!("CARGO_MANIFEST_DIR"));
    let mut reader = bam::IndexedReader::from_path(&path).expect("Failed to open indexed BAM");

    reader.fetch(("CHROMOSOME_I", 0u32, 2u32)).unwrap();
    let mut count = 0u32;
    for result in reader.records() {
        let record = result.expect("Failed to read record");
        assert!(!record.qname().is_empty());
        count += 1;
    }
    assert!(count > 0, "Expected to find records in the region");
}

#[test]
fn test_bam_reader_with_mimalloc() {
    let path = format!("{}/test/test.bam", env!("CARGO_MANIFEST_DIR"));
    let mut reader = bam::Reader::from_path(&path).expect("Failed to open BAM");

    let mut count = 0u32;
    for result in reader.records() {
        let record = result.expect("Failed to read record");
        assert!(!record.qname().is_empty());
        count += 1;
    }
    assert!(count > 0, "Expected to find records");
}

#[test]
fn test_bam_record_clone_with_mimalloc() {
    let path = format!("{}/test/test.bam", env!("CARGO_MANIFEST_DIR"));
    let mut reader = bam::Reader::from_path(&path).expect("Failed to open BAM");

    let mut record = bam::Record::new();
    reader.read(&mut record).expect("Failed to read").expect("No records");

    let cloned = record.clone();
    assert_eq!(record.qname(), cloned.qname());
    assert_eq!(record.seq().as_bytes(), cloned.seq().as_bytes());
}