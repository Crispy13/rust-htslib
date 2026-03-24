#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use rust_htslib::faidx;

#[test]
fn test_faidx_fetch_seq_with_mimalloc() {
    let path = format!("{}/test/test_cram.fa", env!("CARGO_MANIFEST_DIR"));
    let reader = faidx::Reader::from_path(&path).expect("Failed to open faidx");
    let seq = reader.fetch_seq("chr1", 0, 9).unwrap();
    assert_eq!(seq.len(), 10);
    assert_eq!(seq, b"GGGCACAGCC");
}

#[test]
fn test_faidx_fetch_seq_string_with_mimalloc() {
    let path = format!("{}/test/test_cram.fa", env!("CARGO_MANIFEST_DIR"));
    let reader = faidx::Reader::from_path(&path).expect("Failed to open faidx");
    let seq = reader.fetch_seq_string("chr1", 0, 9).unwrap();
    assert_eq!(seq.len(), 10);
    assert_eq!(seq, "GGGCACAGCC");
}

#[test]
fn test_faidx_multiple_fetches_with_mimalloc() {
    let path = format!("{}/test/test_cram.fa", env!("CARGO_MANIFEST_DIR"));
    let reader = faidx::Reader::from_path(&path).expect("Failed to open faidx");

    let seq1 = reader.fetch_seq("chr1", 0, 9).unwrap();
    assert_eq!(seq1, b"GGGCACAGCC");

    let seq2 = reader.fetch_seq("chr1", 110, 119).unwrap();
    assert_eq!(seq2.len(), 10);
    assert_eq!(seq2, b"CCCCTCCGTG");
}

/// Test with a large chromosome sequence (~249MB) to reproduce the SIGSEGV
/// reported when mimalloc tries to free large mmap-backed allocations.
/// Run with: cargo test --test mimalloc_faidx -- --ignored
#[test]
#[ignore]
fn test_faidx_large_seq_with_mimalloc() {
    let path = "/home/eck/workspace/vardict_rs/testdata/hs37d5.fa";
    if !std::path::Path::new(path).exists() {
        eprintln!("Skipping: large FASTA not found at {}", path);
        return;
    }

    let reader = faidx::Reader::from_path(path).expect("Failed to open faidx");

    // Fetch entire chr1 (~249MB) - this triggers mimalloc's large allocation path.
    let seq = reader.fetch_seq("1", 0, 249_250_620).unwrap();
    assert_eq!(seq.len(), 249_250_621);

    // Verify first few bases.
    assert_eq!(&seq[..10], b"NNNNNNNNNN");

    // seq is dropped here - with the old Vec::from_raw_parts bug, this would SIGSEGV
    // because mimalloc tries to free a pointer allocated by htslib's C malloc.
}