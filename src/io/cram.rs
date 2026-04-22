use super::fastq::{normalize_qname, FragmentRecord};
use crate::error::{Result, RvScreenError};
use noodles::cram;
use std::collections::{BTreeSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};

#[derive(Debug)]
pub struct CramFragmentReader {
    path: PathBuf,
    reference_path: PathBuf,
    buffered_pairs: VecDeque<FragmentRecord>,
    decode_attempted: bool,
}

impl CramFragmentReader {
    pub fn open<P, R>(path: P, reference_path: Option<R>) -> Result<Self>
    where
        P: AsRef<Path>,
        R: AsRef<Path>,
    {
        let path = path.as_ref().to_path_buf();
        File::open(&path).map_err(|source| RvScreenError::io(&path, source))?;

        let reference_path = reference_path
            .map(|path| path.as_ref().to_path_buf())
            .ok_or_else(|| {
                RvScreenError::validation(
                    "reference",
                    format!(
                        "CRAM input `{}` requires an explicit reference FASTA path before decoding",
                        path.display()
                    ),
                )
            })?;

        File::open(&reference_path).map_err(|source| RvScreenError::io(&reference_path, source))?;
        ensure_reference_fasta_has_sequence_definitions(&reference_path)?;

        Ok(Self {
            path,
            reference_path,
            buffered_pairs: VecDeque::new(),
            decode_attempted: false,
        })
    }

    pub fn next_pair(&mut self) -> Option<Result<FragmentRecord>> {
        if let Some(pair) = self.buffered_pairs.pop_front() {
            return Some(Ok(pair));
        }

        if self.decode_attempted {
            return None;
        }

        if let Err(err) = self.decode_all_pairs() {
            self.decode_attempted = true;
            return Some(Err(err));
        }

        self.buffered_pairs.pop_front().map(Ok)
    }

    fn decode_all_pairs(&mut self) -> Result<()> {
        let reference_names = read_reference_sequence_names(&self.reference_path)?;
        let file =
            File::open(&self.path).map_err(|source| RvScreenError::io(&self.path, source))?;
        let mut reader = cram::io::Reader::new(file);
        let header = reader
            .read_header()
            .map_err(|source| RvScreenError::io(&self.path, source))?;

        let required_reference_names = header
            .reference_sequences()
            .keys()
            .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
            .collect::<Vec<_>>();

        validate_reference_coverage(
            &self.path,
            &self.reference_path,
            &reference_names,
            required_reference_names.clone(),
        )?;

        let records = catch_unwind(AssertUnwindSafe(|| {
            reader.records(&header).collect::<std::io::Result<Vec<_>>>()
        }));

        let records = match records {
            Ok(Ok(records)) => records,
            Ok(Err(source)) => return Err(RvScreenError::io(&self.path, source)),
            Err(_) => {
                return Err(reference_decode_limitation_error(
                    &self.path,
                    &self.reference_path,
                    required_reference_names,
                ))
            }
        };

        let mut pending_mate = None;

        for (index, record) in records.iter().enumerate() {
            let mate = extract_mate_record(
                &self.path,
                index as u64 + 1,
                record
                    .name()
                    .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned()),
                record.sequence().as_ref().to_vec(),
                record.quality_scores().as_ref().to_vec(),
            )?;

            if let Some(first) = pending_mate.take() {
                self.buffered_pairs
                    .push_back(build_fragment_record(&self.path, &first, &mate)?);
            } else {
                pending_mate = Some(mate);
            }
        }

        if let Some(orphan) = pending_mate {
            return Err(RvScreenError::parse(
                &self.path,
                orphan.record_number,
                format!(
                    "orphaned CRAM record {} `{}`: adjacent mate is missing; provide name-sorted input",
                    orphan.record_number, orphan.qname,
                ),
            ));
        }

        self.decode_attempted = true;

        Ok(())
    }
}

impl Iterator for CramFragmentReader {
    type Item = Result<FragmentRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_pair()
    }
}

#[derive(Debug, Clone)]
struct MateRecord {
    record_number: u64,
    qname: String,
    fragment_key: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

fn extract_mate_record(
    path: &Path,
    record_number: u64,
    qname: Option<String>,
    seq: Vec<u8>,
    qual: Vec<u8>,
) -> Result<MateRecord> {
    let qname = qname.ok_or_else(|| {
        RvScreenError::parse(
            path,
            record_number,
            format!("CRAM record {record_number} is missing QNAME"),
        )
    })?;

    Ok(MateRecord {
        record_number,
        fragment_key: normalize_qname(&qname),
        qname,
        seq,
        qual,
    })
}

fn build_fragment_record(
    path: &Path,
    first: &MateRecord,
    second: &MateRecord,
) -> Result<FragmentRecord> {
    if first.fragment_key != second.fragment_key {
        return Err(RvScreenError::parse(
            path,
            first.record_number,
            format!(
                "adjacent CRAM records are not a name-sorted pair: record {} `{}` -> `{}`, record {} `{}` -> `{}`",
                first.record_number,
                first.qname,
                first.fragment_key,
                second.record_number,
                second.qname,
                second.fragment_key,
            ),
        ));
    }

    Ok(FragmentRecord {
        fragment_key: first.fragment_key.clone(),
        r1_seq: first.seq.clone(),
        r1_qual: first.qual.clone(),
        r2_seq: second.seq.clone(),
        r2_qual: second.qual.clone(),
    })
}

fn ensure_reference_fasta_has_sequence_definitions(path: &Path) -> Result<()> {
    read_reference_sequence_names(path).map(|_| ())
}

fn read_reference_sequence_names(path: &Path) -> Result<BTreeSet<String>> {
    let file = File::open(path).map_err(|source| RvScreenError::io(path, source))?;
    let reader = BufReader::new(file);
    let mut names = BTreeSet::new();

    for (index, line) in reader.lines().enumerate() {
        let line_number = index as u64 + 1;
        let line = line.map_err(|source| RvScreenError::io(path, source))?;

        if let Some(definition) = line.strip_prefix('>') {
            let name = definition.split_ascii_whitespace().next().unwrap_or("");
            if name.is_empty() {
                return Err(RvScreenError::parse(
                    path,
                    line_number,
                    "FASTA reference header is missing a sequence name",
                ));
            }

            names.insert(name.to_owned());
        }
    }

    if names.is_empty() {
        return Err(RvScreenError::validation(
            "reference",
            format!(
                "reference FASTA `{}` does not contain any sequence definitions",
                path.display()
            ),
        ));
    }

    Ok(names)
}

fn validate_reference_coverage(
    cram_path: &Path,
    reference_path: &Path,
    reference_names: &BTreeSet<String>,
    required_reference_names: Vec<String>,
) -> Result<()> {
    if required_reference_names.is_empty() {
        return Ok(());
    }

    let missing: Vec<_> = required_reference_names
        .into_iter()
        .filter(|name| !reference_names.contains(name.as_str()))
        .collect();

    if missing.is_empty() {
        Ok(())
    } else {
        Err(RvScreenError::validation(
            "reference",
            format!(
                "reference FASTA `{}` is missing CRAM header reference sequences required by `{}`: {}",
                reference_path.display(),
                cram_path.display(),
                missing.join(", "),
            ),
        ))
    }
}

fn reference_decode_limitation_error(
    cram_path: &Path,
    reference_path: &Path,
    required_reference_names: Vec<String>,
) -> RvScreenError {
    let required = if required_reference_names.is_empty() {
        "<none declared in CRAM header>".to_owned()
    } else {
        required_reference_names.join(", ")
    };

    RvScreenError::validation(
        "reference",
        format!(
            "CRAM input `{}` could not be fully decoded with provided FASTA `{}`. This build validates the FASTA path and header reference names, but cannot wire the FASTA into noodles CRAM decoding because the direct `noodles` dependency omits the `fasta` feature; only CRAMs that decode without an external repository are supported. Header references: [{}]",
            cram_path.display(),
            reference_path.display(),
            required,
        ),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn test_missing_reference_returns_validation_error() {
        let tempdir = tempdir().expect("tempdir should be created");
        let cram_path = tempdir.path().join("sample.cram");
        fs::write(&cram_path, b"CRAM-placeholder").expect("placeholder CRAM path should exist");

        let err = CramFragmentReader::open(&cram_path, Option::<&Path>::None)
            .expect_err("missing CRAM reference must fail");

        match err {
            RvScreenError::ValidationError { field, reason } => {
                assert_eq!(field, "reference");
                assert!(reason.contains("reference FASTA"));
                assert!(reason.contains("CRAM"));
            }
            other => panic!("expected ValidationError, got {other:?}"),
        }
    }

    #[test]
    fn test_open_rejects_reference_fasta_without_headers() {
        let tempdir = tempdir().expect("tempdir should be created");
        let cram_path = tempdir.path().join("sample.cram");
        let reference_path = tempdir.path().join("reference.fa");
        fs::write(&cram_path, b"CRAM-placeholder").expect("placeholder CRAM path should exist");
        fs::write(&reference_path, b"ACGT\nACGT\n").expect("reference FASTA should be written");

        let err = CramFragmentReader::open(&cram_path, Some(&reference_path))
            .expect_err("reference FASTA without headers must fail");

        match err {
            RvScreenError::ValidationError { field, reason } => {
                assert_eq!(field, "reference");
                assert!(reason.contains("does not contain any sequence definitions"));
            }
            other => panic!("expected ValidationError, got {other:?}"),
        }
    }

    #[test]
    fn test_validate_reference_coverage_reports_missing_header_sequence() {
        let tempdir = tempdir().expect("tempdir should be created");
        let cram_path = tempdir.path().join("sample.cram");
        let reference_path = tempdir.path().join("reference.fa");
        let reference_names = BTreeSet::from([String::from("chr1")]);

        let err = validate_reference_coverage(
            &cram_path,
            &reference_path,
            &reference_names,
            vec![String::from("chr1"), String::from("chr2")],
        )
        .expect_err("missing CRAM header reference should fail validation");

        match err {
            RvScreenError::ValidationError { field, reason } => {
                assert_eq!(field, "reference");
                assert!(reason.contains("chr2"));
                assert!(reason.contains("sample.cram"));
            }
            other => panic!("expected ValidationError, got {other:?}"),
        }
    }
}
