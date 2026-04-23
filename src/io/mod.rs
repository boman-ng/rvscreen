pub mod bam;
pub mod cram;
pub mod fastq;
mod contract;
mod fragment_reader;
mod input;

pub use bam::BamFragmentReader;
pub use cram::CramFragmentReader;
pub use fastq::{normalize_qname, FastqPairReader};
pub use fragment_reader::{FragmentReaderFactory, FragmentRecord, FragmentStream};
pub use input::ScreenInput;
