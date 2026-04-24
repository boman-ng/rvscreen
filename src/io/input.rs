use super::{BamFragmentReader, CramFragmentReader, FastqPairReader};
use super::{FragmentReaderFactory, FragmentStream};
use crate::error::{Result, RvScreenError};
use std::path::{Path, PathBuf};

const COMPOSITE_FASTA: &str = "composite.fa";

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ScreenInput {
    FastqPair {
        r1: PathBuf,
        r2: PathBuf,
    },
    Bam {
        path: PathBuf,
    },
    Cram {
        path: PathBuf,
        reference_fasta: PathBuf,
    },
}

impl ScreenInput {
    pub fn from_cli_inputs(inputs: &[PathBuf], reference_bundle_dir: &Path) -> Result<Self> {
        match inputs {
            [r1, r2] => {
                if !is_fastq_like(r1) || !is_fastq_like(r2) {
                    return Err(RvScreenError::validation(
                        "input",
                        format!(
                            "two-input mode currently supports only paired FASTQ files; got `{}` and `{}`",
                            r1.display(),
                            r2.display()
                        ),
                    ));
                }

                Ok(Self::FastqPair {
                    r1: r1.clone(),
                    r2: r2.clone(),
                })
            }
            [path] if is_bam_like(path) => Ok(Self::Bam { path: path.clone() }),
            [path] if is_cram_like(path) => Ok(Self::Cram {
                path: path.clone(),
                reference_fasta: reference_bundle_dir.join(COMPOSITE_FASTA),
            }),
            [path] if is_fastq_like(path) => Err(RvScreenError::validation(
                "input",
                format!(
                    "single-file interleaved FASTQ input `{}` is not yet wired in Task 18; provide two paired FASTQ files, one BAM/uBAM, or one CRAM",
                    path.display()
                ),
            )),
            [path] => Err(RvScreenError::validation(
                "input",
                format!(
                    "unsupported single input `{}`; supported inputs are paired FASTQ, BAM/uBAM, or CRAM",
                    path.display()
                ),
            )),
            [] => Err(RvScreenError::validation(
                "input",
                "at least one input file is required",
            )),
            _ => Err(RvScreenError::validation(
                "input",
                "screen currently accepts exactly two paired FASTQs or exactly one BAM/uBAM/CRAM input",
            )),
        }
    }
}

impl FragmentReaderFactory for ScreenInput {
    fn open_reader(&self) -> Result<FragmentStream> {
        match self {
            Self::FastqPair { r1, r2 } => Ok(Box::new(FastqPairReader::open(r1, r2)?)),
            Self::Bam { path } => Ok(Box::new(BamFragmentReader::open(path)?)),
            Self::Cram {
                path,
                reference_fasta,
            } => Ok(Box::new(CramFragmentReader::open(
                path,
                Some(reference_fasta),
            )?)),
        }
    }

    fn open_fastq_pair_reader(&self) -> Option<Result<FastqPairReader>> {
        match self {
            Self::FastqPair { r1, r2 } => Some(FastqPairReader::open(r1, r2)),
            _ => None,
        }
    }

    fn kind_label(&self) -> &'static str {
        match self {
            Self::FastqPair { r1, .. } if is_gzip_path(r1) => "fastq.gz",
            Self::FastqPair { .. } => "fastq",
            Self::Bam { path } if is_ubam_like(path) => "ubam",
            Self::Bam { .. } => "bam",
            Self::Cram { .. } => "cram",
        }
    }
}

fn is_fastq_like(path: &Path) -> bool {
    let path = path.to_string_lossy().to_ascii_lowercase();
    path.ends_with(".fastq")
        || path.ends_with(".fq")
        || path.ends_with(".fastq.gz")
        || path.ends_with(".fq.gz")
}

fn is_bam_like(path: &Path) -> bool {
    let path = path.to_string_lossy().to_ascii_lowercase();
    path.ends_with(".bam") || path.ends_with(".ubam")
}

fn is_ubam_like(path: &Path) -> bool {
    path.to_string_lossy()
        .to_ascii_lowercase()
        .ends_with(".ubam")
}

fn is_cram_like(path: &Path) -> bool {
    path.to_string_lossy()
        .to_ascii_lowercase()
        .ends_with(".cram")
}

fn is_gzip_path(path: &Path) -> bool {
    path.to_string_lossy().to_ascii_lowercase().ends_with(".gz")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paired_fastq_inputs_stay_owned_by_io_boundary() {
        let input = ScreenInput::from_cli_inputs(
            &[
                PathBuf::from("sample_R1.fastq.gz"),
                PathBuf::from("sample_R2.fastq.gz"),
            ],
            Path::new("bundle"),
        )
        .expect("paired FASTQ inputs should be accepted");

        assert!(matches!(input, ScreenInput::FastqPair { .. }));
        assert_eq!(input.kind_label(), "fastq.gz");
    }

    #[test]
    fn bam_and_cram_inputs_keep_kind_specific_metadata() {
        let bam = ScreenInput::from_cli_inputs(&[PathBuf::from("reads.ubam")], Path::new("bundle"))
            .expect("uBAM input should be accepted");
        assert!(matches!(bam, ScreenInput::Bam { .. }));
        assert_eq!(bam.kind_label(), "ubam");

        let cram =
            ScreenInput::from_cli_inputs(&[PathBuf::from("reads.cram")], Path::new("bundle"))
                .expect("CRAM input should be accepted");
        match cram {
            ScreenInput::Cram {
                reference_fasta, ..
            } => {
                assert_eq!(reference_fasta, PathBuf::from("bundle/composite.fa"));
            }
            other => panic!("expected CRAM input, got {other:?}"),
        }
    }

    #[test]
    fn single_fastq_is_rejected_with_boundary_guidance() {
        let error =
            ScreenInput::from_cli_inputs(&[PathBuf::from("reads.fastq")], Path::new("bundle"))
                .expect_err("single FASTQ input should be rejected");

        assert!(error
            .to_string()
            .contains("single-file interleaved FASTQ input"));
    }
}
