pub mod bam;
pub mod cram;
pub mod fastq;

pub use bam::BamFragmentReader;
pub use cram::CramFragmentReader;
pub use fastq::{normalize_qname, FastqPairReader, FragmentRecord};
