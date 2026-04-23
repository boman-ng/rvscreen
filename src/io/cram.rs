use super::fastq::normalize_qname;
use super::fragment_reader::{
    build_name_sorted_fragment_record, mate_slot_from_flags, orphaned_name_sorted_record_error,
    FragmentMateRecord,
};
use super::FragmentRecord;
use crate::error::{Result, RvScreenError};
use noodles::{
    cram,
    fasta::{self, repository::adapters::IndexedReader},
    sam,
};
use std::collections::BTreeSet;
use std::fs::File;
use std::io::ErrorKind;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};

pub struct CramFragmentReader {
    path: PathBuf,
    reference_path: PathBuf,
    reference_repository: fasta::Repository,
    reader: cram::io::Reader<File>,
    header: sam::Header,
    container: cram::io::reader::Container,
    container_records: std::vec::IntoIter<sam::alignment::RecordBuf>,
    pending_mate: Option<FragmentMateRecord>,
    record_number: u64,
    finished: bool,
}

impl std::fmt::Debug for CramFragmentReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CramFragmentReader")
            .field("path", &self.path)
            .field("reference_path", &self.reference_path)
            .field("has_pending_mate", &self.pending_mate.is_some())
            .field("record_number", &self.record_number)
            .field("finished", &self.finished)
            .finish()
    }
}

impl CramFragmentReader {
    pub fn open<P, R>(path: P, reference_path: Option<R>) -> Result<Self>
    where
        P: AsRef<Path>,
        R: AsRef<Path>,
    {
        let path = path.as_ref().to_path_buf();
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

        let (reference_repository, _) = build_reference_sequence_repository(&reference_path)?;
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(reference_repository.clone())
            .build_from_path(&path)
            .map_err(|source| RvScreenError::io(&path, source))?;
        let header = reader
            .read_header()
            .map_err(|source| RvScreenError::io(&path, source))?;

        Ok(Self {
            path,
            reference_path,
            reference_repository,
            reader,
            header,
            container: cram::io::reader::Container::default(),
            container_records: Vec::new().into_iter(),
            pending_mate: None,
            record_number: 0,
            finished: false,
        })
    }

    pub fn next_pair(&mut self) -> Option<Result<FragmentRecord>> {
        if self.finished {
            return None;
        }

        loop {
            let record = match self.next_record_buf() {
                Ok(Some(record)) => record,
                Ok(None) => {
                    self.finished = true;
                    return self.pending_mate.take().map(|orphan| Err(orphan_cram_record_error(
                        &self.path,
                        &orphan,
                    )));
                }
                Err(err) => {
                    self.finished = true;
                    return Some(Err(err));
                }
            };

            self.record_number += 1;

            let mate = match extract_mate_record(&self.path, self.record_number, &record) {
                Ok(mate) => mate,
                Err(err) => {
                    self.finished = true;
                    return Some(Err(err));
                }
            };

            if let Some(first) = self.pending_mate.take() {
                let pair = build_fragment_record(&self.path, &first, &mate);
                if pair.is_err() {
                    self.finished = true;
                }
                return Some(pair);
            }

            self.pending_mate = Some(mate);
        }
    }

    fn next_record_buf(&mut self) -> Result<Option<sam::alignment::RecordBuf>> {
        loop {
            if let Some(record) = self.container_records.next() {
                return Ok(Some(record));
            }

            if !self.read_next_container_records()? {
                return Ok(None);
            }
        }
    }

    fn read_next_container_records(&mut self) -> Result<bool> {
        let container_len = self
            .reader
            .read_container(&mut self.container)
            .map_err(|source| RvScreenError::io(&self.path, source))?;

        if container_len == 0 {
            return Ok(false);
        }

        let container = &self.container;
        let header = &self.header;
        let reference_repository = self.reference_repository.clone();

        let records = catch_unwind(AssertUnwindSafe(|| {
            decode_container_records(container, header, reference_repository)
        }));

        let records = match records {
            Ok(Ok(records)) => records,
            Ok(Err(source)) => return Err(RvScreenError::io(&self.path, source)),
            Err(_) => {
                return Err(reference_decode_failure_error(
                    &self.path,
                    &self.reference_path,
                    &self.header,
                ))
            }
        };

        self.container_records = records.into_iter();

        Ok(true)
    }
}

impl Iterator for CramFragmentReader {
    type Item = Result<FragmentRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_pair()
    }
}

fn extract_mate_record(
    path: &Path,
    record_number: u64,
    record: &sam::alignment::RecordBuf,
) -> Result<FragmentMateRecord> {
    let qname = record.name().map(|name| name.to_string()).ok_or_else(|| {
        RvScreenError::parse(
            path,
            record_number,
            format!("CRAM record {record_number} is missing QNAME"),
        )
    })?;

    let fragment_key = normalize_qname(&qname);

    Ok(FragmentMateRecord::new(
        record_number,
        qname,
        fragment_key,
        record.sequence().as_ref().to_vec(),
        record.quality_scores().as_ref().to_vec(),
        mate_slot_from_flags(record.flags()),
    ))
}

fn build_fragment_record(
    path: &Path,
    first: &FragmentMateRecord,
    second: &FragmentMateRecord,
) -> Result<FragmentRecord> {
    build_name_sorted_fragment_record(path, "CRAM", first, second)
}

fn decode_container_records(
    container: &cram::io::reader::Container,
    header: &sam::Header,
    reference_repository: fasta::Repository,
) -> std::io::Result<Vec<sam::alignment::RecordBuf>> {
    let compression_header = container.compression_header()?;

    container
        .slices()
        .map(|result| {
            let slice = result?;
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            slice
                .records(
                    reference_repository.clone(),
                    header,
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                )
                .and_then(|records| {
                    records
                        .into_iter()
                        .map(|record| {
                            sam::alignment::RecordBuf::try_from_alignment_record(header, &record)
                        })
                        .collect::<std::io::Result<Vec<_>>>()
                })
        })
        .collect::<std::io::Result<Vec<_>>>()
        .map(|record_batches| record_batches.into_iter().flatten().collect())
}

fn build_reference_sequence_repository(path: &Path) -> Result<(fasta::Repository, BTreeSet<String>)> {
    let fai_path = reference_fai_path(path);
    let index = fasta::fai::fs::read(&fai_path)
        .map_err(|source| map_reference_index_error(path, &fai_path, source))?;

    let reference_names: BTreeSet<_> = index
        .as_ref()
        .iter()
        .map(|record| String::from_utf8_lossy(record.name().as_ref()).into_owned())
        .collect();

    if reference_names.is_empty() {
        return Err(RvScreenError::validation(
            "reference",
            format!(
                "reference FASTA index `{}` for `{}` does not contain any sequence definitions",
                fai_path.display(),
                path.display()
            ),
        ));
    }

    let indexed_reader = fasta::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(path)
        .map_err(|source| RvScreenError::io(path, source))?;
    let repository = fasta::Repository::new(IndexedReader::new(indexed_reader));

    Ok((repository, reference_names))
}

fn reference_fai_path(reference_path: &Path) -> PathBuf {
    let mut fai = reference_path.as_os_str().to_owned();
    fai.push(".fai");
    PathBuf::from(fai)
}

fn map_reference_index_error(
    reference_path: &Path,
    fai_path: &Path,
    source: std::io::Error,
) -> RvScreenError {
    if source.kind() == ErrorKind::NotFound {
        return RvScreenError::validation(
            "reference",
            format!(
                "CRAM input requires FASTA index `{}` alongside `{}`; run `samtools faidx {}` before decoding",
                fai_path.display(),
                reference_path.display(),
                reference_path.display()
            ),
        );
    }

    RvScreenError::io(fai_path, source)
}

fn reference_decode_failure_error(
    cram_path: &Path,
    reference_path: &Path,
    header: &sam::Header,
) -> RvScreenError {
    let required_reference_names = header
        .reference_sequences()
        .keys()
        .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
        .collect::<Vec<_>>();
    let required = if required_reference_names.is_empty() {
        "<none declared in CRAM header>".to_owned()
    } else {
        required_reference_names.join(", ")
    };

    RvScreenError::validation(
        "reference",
        format!(
            "CRAM input `{}` could not be decoded with FASTA `{}` despite validated FASTA/FAI reference wiring. The decoder panicked while materializing record data, which usually indicates a malformed CRAM or a FASTA/FAI that does not match the CRAM payload. Header references: [{}]",
            cram_path.display(),
            reference_path.display(),
            required,
        ),
    )
}

fn orphan_cram_record_error(path: &Path, orphan: &FragmentMateRecord) -> RvScreenError {
    orphaned_name_sorted_record_error(
        path,
        "CRAM",
        orphan,
        "provide name-sorted input",
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::{
        cram,
        fasta::{self, repository::adapters::IndexedReader},
        sam::{
            self,
            alignment::{RecordBuf, io::Write as _},
            header::record::value::{Map, map::ReferenceSequence},
        },
    };
    use std::fs;
    use std::io::{self, Write};
    use std::num::NonZeroUsize;
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
    fn test_open_rejects_reference_fasta_without_fai_index() {
        let tempdir = tempdir().expect("tempdir should be created");
        let cram_path = tempdir.path().join("sample.cram");
        let reference_path = tempdir.path().join("reference.fa");
        fs::write(&cram_path, b"CRAM-placeholder").expect("placeholder CRAM path should exist");
        fs::write(&reference_path, b">chr1\nACGT\n").expect("reference FASTA should be written");

        let err = CramFragmentReader::open(&cram_path, Some(&reference_path))
            .expect_err("reference FASTA without FAI must fail");

        match err {
            RvScreenError::ValidationError { field, reason } => {
                assert_eq!(field, "reference");
                assert!(reason.contains("samtools faidx"));
                assert!(reason.contains("reference.fa.fai"));
            }
            other => panic!("expected ValidationError, got {other:?}"),
        }
    }

    #[test]
    fn test_open_accepts_unused_header_reference_absent_from_fasta() {
        let tempdir = tempdir().expect("tempdir should be created");
        let reader_reference_path = tempdir.path().join("reference.fa");
        let writer_reference_path = tempdir.path().join("writer-reference.fa");
        let cram_path = tempdir.path().join("sample.cram");
        let reference_sequence = "ACGT".repeat(3000);

        write_test_reference(&reader_reference_path, "chr1", &reference_sequence)
            .expect("reference FASTA and FAI should be written");
        write_multi_reference(&writer_reference_path, &[("chr1", &reference_sequence), ("chr2", &reference_sequence)])
            .expect("writer reference FASTA and FAI should be written");
        write_test_cram_with_header(
            &cram_path,
            &writer_reference_path,
            &[("chr1", reference_sequence.len()), ("chr2", reference_sequence.len())],
            &[SyntheticCramRecord::first(
                "READ:1:FCX123:1:1101:1000:2000",
                b"ACGT",
                &[30, 31, 32, 33],
            )],
        )
        .expect("CRAM fixture should be written with unused extra header reference");

        let mut reader = CramFragmentReader::open(&cram_path, Some(&reader_reference_path))
            .expect("unused header-only contigs should not fail open-time validation");
        let err = reader
            .next_pair()
            .expect("single mapped record should yield one terminal result")
            .expect_err("single mapped record should remain orphaned without mate");

        match err {
            RvScreenError::ParseError { reason, .. } => {
                assert!(reason.contains("orphaned CRAM record"));
            }
            other => panic!("expected ParseError, got {other:?}"),
        }
    }

    #[test]
    fn test_reads_name_sorted_cram_pairs_across_multiple_containers() {
        let tempdir = tempdir().expect("tempdir should be created");
        let reference_path = tempdir.path().join("reference.fa");
        let cram_path = tempdir.path().join("sample.cram");
        let reference_sequence = "ACGT".repeat(3000);

        write_test_reference(&reference_path, "chr1", &reference_sequence)
            .expect("reference FASTA and FAI should be written");

        let pair_count = 5_121usize;
        let mut records = Vec::with_capacity(pair_count * 2);

        for pair_index in 0..pair_count {
            let qname = format!("READ:1:FCX123:1:1101:{}:{}", 1000 + pair_index, 2000 + pair_index);
            records.push(SyntheticCramRecord::first(&qname, b"ACGT", &[30, 31, 32, 33]));
            records.push(SyntheticCramRecord::last(&qname, b"TGCA", &[34, 35, 36, 37]));
        }

        write_test_cram(&cram_path, &reference_path, "chr1", reference_sequence.len(), &records)
            .expect("CRAM fixture should be written");

        let mut reader = CramFragmentReader::open(&cram_path, Some(&reference_path))
            .expect("reader should open CRAM with FASTA/FAI");
        let mut pair_index = 0usize;

        while let Some(record) = reader.next_pair() {
            let record = record.expect("CRAM pair should decode cleanly");
            if pair_index == 0 {
                assert_eq!(record.fragment_key, "READ:1:FCX123:1:1101:1000:2000");
            }

            pair_index += 1;
        }

        assert_eq!(pair_index, pair_count);
    }

    #[test]
    fn test_reader_orders_flagged_cram_mates_into_r1_r2_slots() {
        let tempdir = tempdir().expect("tempdir should be created");
        let reference_path = tempdir.path().join("reference.fa");
        let cram_path = tempdir.path().join("slot-ordered.cram");
        let reference_sequence = "ACGT".repeat(3000);

        write_test_reference(&reference_path, "chr1", &reference_sequence)
            .expect("reference FASTA and FAI should be written");
        write_test_cram(
            &cram_path,
            &reference_path,
            "chr1",
            reference_sequence.len(),
            &[
                SyntheticCramRecord::last(
                    "READ:1:FCX123:1:1101:1000:2000",
                    b"TGCA",
                    &[34, 35, 36, 37],
                ),
                SyntheticCramRecord::first(
                    "READ:1:FCX123:1:1101:1000:2000",
                    b"ACGT",
                    &[30, 31, 32, 33],
                ),
            ],
        )
        .expect("CRAM fixture should be written");

        let record = CramFragmentReader::open(&cram_path, Some(&reference_path))
            .expect("reader should open CRAM")
            .next_pair()
            .expect("pair should be available")
            .expect("pair should decode");

        assert_eq!(record.r1_seq, b"ACGT");
        assert_eq!(record.r1_qual, vec![30, 31, 32, 33]);
        assert_eq!(record.r2_seq, b"TGCA");
        assert_eq!(record.r2_qual, vec![34, 35, 36, 37]);
    }

    #[derive(Debug, Clone)]
    struct SyntheticCramRecord {
        qname: String,
        flags: sam::alignment::record::Flags,
        sequence: Vec<u8>,
        quality_scores: Vec<u8>,
    }

    impl SyntheticCramRecord {
        fn first(qname: &str, sequence: &[u8], quality_scores: &[u8]) -> Self {
            Self::new(
                qname,
                sam::alignment::record::Flags::SEGMENTED
                    | sam::alignment::record::Flags::UNMAPPED
                    | sam::alignment::record::Flags::MATE_UNMAPPED
                    | sam::alignment::record::Flags::FIRST_SEGMENT,
                sequence,
                quality_scores,
            )
        }

        fn last(qname: &str, sequence: &[u8], quality_scores: &[u8]) -> Self {
            Self::new(
                qname,
                sam::alignment::record::Flags::SEGMENTED
                    | sam::alignment::record::Flags::UNMAPPED
                    | sam::alignment::record::Flags::MATE_UNMAPPED
                    | sam::alignment::record::Flags::LAST_SEGMENT,
                sequence,
                quality_scores,
            )
        }

        fn new(
            qname: &str,
            flags: sam::alignment::record::Flags,
            sequence: &[u8],
            quality_scores: &[u8],
        ) -> Self {
            Self {
                qname: qname.to_owned(),
                flags,
                sequence: sequence.to_vec(),
                quality_scores: quality_scores.to_vec(),
            }
        }

        fn to_record_buf(&self) -> RecordBuf {
            RecordBuf::builder()
                .set_name(self.qname.clone())
                .set_flags(self.flags)
                .set_sequence(self.sequence.clone().into())
                .set_quality_scores(self.quality_scores.clone().into())
                .build()
        }
    }

    fn write_test_reference(path: &Path, name: &str, sequence: &str) -> io::Result<()> {
        write_multi_reference(path, &[(name, sequence)])
    }

    fn write_multi_reference(path: &Path, references: &[(&str, &str)]) -> io::Result<()> {
        let mut fasta = File::create(path)?;
        let mut offset = 0u64;
        let mut fai = String::new();

        for (name, sequence) in references {
            writeln!(fasta, ">{name}")?;
            writeln!(fasta, "{sequence}")?;

            let sequence_offset = offset + (name.len() + 2) as u64;
            fai.push_str(&format!(
                "{name}\t{}\t{sequence_offset}\t{}\t{}\n",
                sequence.len(),
                sequence.len(),
                sequence.len() + 1,
            ));
            offset = sequence_offset + sequence.len() as u64 + 1;
        }

        let fai_path = reference_fai_path(path);
        fs::write(fai_path, fai)
    }

    fn write_test_cram(
        path: &Path,
        reference_path: &Path,
        reference_name: &str,
        reference_len: usize,
        records: &[SyntheticCramRecord],
    ) -> io::Result<()> {
        write_test_cram_with_header(
            path,
            reference_path,
            &[(reference_name, reference_len)],
            records,
        )
    }

    fn write_test_cram_with_header(
        path: &Path,
        reference_path: &Path,
        header_references: &[(&str, usize)],
        records: &[SyntheticCramRecord],
    ) -> io::Result<()> {
        let repository = fasta::io::indexed_reader::Builder::default()
            .build_from_path(reference_path)
            .map(IndexedReader::new)
            .map(fasta::Repository::new)?;
        let header = header_references
            .iter()
            .fold(sam::Header::builder(), |builder, (reference_name, reference_len)| {
                builder.add_reference_sequence(
                    *reference_name,
                    Map::<ReferenceSequence>::new(NonZeroUsize::new(*reference_len).unwrap()),
                )
            })
            .build();
        let file = File::create(path)?;
        let mut writer = cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_writer(file);

        writer.write_header(&header)?;

        for record in records {
            writer.write_alignment_record(&header, &record.to_record_buf())?;
        }

        writer.try_finish(&header)
    }
}
