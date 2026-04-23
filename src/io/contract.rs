//! Shared streaming fragment-reader contract for FASTQ, BAM, and CRAM.
//!
//! The pipeline only depends on this contract:
//! - one `open_reader()` call per screening pass,
//! - pull-based iteration via `Iterator<Item = Result<FragmentRecord>>`,
//! - bounded buffering at reader level,
//! - pair-preserving `FragmentRecord` output across formats,
//! - natural backpressure propagation because callers control `next()` pacing.
//!
//! Internal buffering granularity still differs by format: FASTQ and BAM stream record-by-record,
//! while CRAM remains container-streaming because that is format-inherent in `noodles`.
