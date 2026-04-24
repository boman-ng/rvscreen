pub use super::FragmentRecord;
use crate::error::{Result, RvScreenError};
use flate2::read::MultiGzDecoder;
use seq_io::fastq::{Error as FastqError, Reader, Record};
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::time::Duration;

const FASTQ_RECORD_LINES: u64 = 4;

#[derive(Debug, Clone)]
pub(crate) struct BorrowedFragmentRecord<'a> {
    pub fragment_key: Cow<'a, str>,
    pub r1_seq: &'a [u8],
    pub r1_qual: &'a [u8],
    pub r2_seq: &'a [u8],
    pub r2_qual: &'a [u8],
}

impl<'a> BorrowedFragmentRecord<'a> {
    pub(crate) fn into_owned_with_qualities(self) -> Result<FragmentRecord> {
        Ok(FragmentRecord::new(
            self.fragment_key.into_owned(),
            self.r1_seq.to_vec(),
            decode_fastq_quality_scores(self.r1_qual)?,
            self.r2_seq.to_vec(),
            decode_fastq_quality_scores(self.r2_qual)?,
        ))
    }

    pub(crate) fn into_owned_without_qualities(self) -> FragmentRecord {
        FragmentRecord::new(
            self.fragment_key.into_owned(),
            self.r1_seq.to_vec(),
            Vec::new(),
            self.r2_seq.to_vec(),
            Vec::new(),
        )
    }
}

pub struct FastqPairReader {
    r1_path: PathBuf,
    r2_path: PathBuf,
    r1_reader: Reader<Box<dyn Read>>,
    r2_reader: Reader<Box<dyn Read>>,
    pair_index: u64,
}

impl FastqPairReader {
    pub fn open<P1, P2>(r1_path: P1, r2_path: P2) -> Result<Self>
    where
        P1: AsRef<Path>,
        P2: AsRef<Path>,
    {
        let r1_path = r1_path.as_ref().to_path_buf();
        let r2_path = r2_path.as_ref().to_path_buf();

        Ok(Self {
            r1_reader: open_fastq_reader(&r1_path)?,
            r2_reader: open_fastq_reader(&r2_path)?,
            r1_path,
            r2_path,
            pair_index: 0,
        })
    }

    pub fn next_pair(&mut self) -> Option<Result<FragmentRecord>> {
        self.map_next_pair(|fragment| fragment.into_owned_with_qualities())
    }

    pub(crate) fn for_each_borrowed_pair_with_timing<F>(&mut self, mut visit: F) -> Result<Duration>
    where
        F: for<'a> FnMut(BorrowedFragmentRecord<'a>) -> Result<()>,
    {
        let mut io_wall_duration = Duration::default();

        loop {
            let started = std::time::Instant::now();
            let item = self.map_next_pair(|fragment| visit(fragment));
            io_wall_duration = io_wall_duration.saturating_add(started.elapsed());

            match item {
                Some(result) => result?,
                None => break,
            }
        }

        Ok(io_wall_duration)
    }

    fn map_next_pair<T, F>(&mut self, on_pair: F) -> Option<Result<T>>
    where
        F: for<'a> FnOnce(BorrowedFragmentRecord<'a>) -> Result<T>,
    {
        let line = fastq_record_line(self.pair_index);
        let r1_path = self.r1_path.clone();
        let r2_path = self.r2_path.clone();
        let r1_item = self.r1_reader.next();
        let r2_item = self.r2_reader.next();

        let item = match (r1_item, r2_item) {
            (None, None) => None,
            (Some(Err(err)), _) => Some(Err(map_fastq_error(&r1_path, err))),
            (_, Some(Err(err))) => Some(Err(map_fastq_error(&r2_path, err))),
            (Some(Ok(r1)), Some(Ok(r2))) => Some(
                build_borrowed_fragment_record(&r1_path, &r2_path, line, &r1, &r2)
                    .and_then(on_pair),
            ),
            (Some(Ok(r1)), None) => Some(Err(orphan_error(
                &r1_path,
                &r2_path,
                line,
                &header_text(&r1),
            ))),
            (None, Some(Ok(r2))) => Some(Err(orphan_error(
                &r2_path,
                &r1_path,
                line,
                &header_text(&r2),
            ))),
        };

        if matches!(item, Some(Ok(_))) {
            self.pair_index += 1;
        }

        item
    }
}

impl Iterator for FastqPairReader {
    type Item = Result<FragmentRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_pair()
    }
}

pub fn normalize_qname(qname: &str) -> String {
    normalize_qname_view(qname).to_owned()
}

fn build_borrowed_fragment_record<'a>(
    r1_path: &Path,
    r2_path: &Path,
    line: u64,
    r1: &'a impl Record,
    r2: &'a impl Record,
) -> Result<BorrowedFragmentRecord<'a>> {
    let r1_header = header_text(r1);
    let r2_header = header_text(r2);
    let r1_key = normalize_record_key(r1.head());
    let r2_key = normalize_record_key(r2.head());

    if r1_key != r2_key {
        return Err(RvScreenError::parse(
            r1_path,
            line,
            format!(
                "mismatched paired-end fragment keys between `{}` and `{}`: R1 header `{}` -> `{}`, R2 header `{}` -> `{}`",
                r1_path.display(),
                r2_path.display(),
                r1_header,
                r1_key,
                r2_header,
                r2_key,
            ),
        ));
    }

    Ok(BorrowedFragmentRecord {
        fragment_key: r1_key,
        r1_seq: r1.seq(),
        r1_qual: r1.qual(),
        r2_seq: r2.seq(),
        r2_qual: r2.qual(),
    })
}

fn orphan_error(
    orphan_path: &Path,
    mate_path: &Path,
    line: u64,
    orphan_header: &str,
) -> RvScreenError {
    RvScreenError::parse(
        orphan_path,
        line,
        format!(
            "orphaned read `{orphan_header}` in `{}`: mate FASTQ `{}` ended early",
            orphan_path.display(),
            mate_path.display(),
        ),
    )
}

fn open_fastq_reader(path: &Path) -> Result<Reader<Box<dyn Read>>> {
    let file = File::open(path).map_err(|source| RvScreenError::io(path, source))?;
    let stream: Box<dyn Read> = if is_gzip_path(path) {
        Box::new(MultiGzDecoder::new(BufReader::new(file)))
    } else {
        Box::new(file)
    };

    Ok(Reader::new(stream))
}

fn is_gzip_path(path: &Path) -> bool {
    path.extension().is_some_and(|extension| extension == "gz")
}

fn normalize_qname_view(qname: &str) -> &str {
    let id = qname.split_ascii_whitespace().next().unwrap_or(qname);
    id.strip_suffix("/1")
        .or_else(|| id.strip_suffix("/2"))
        .unwrap_or(id)
}

fn normalize_record_key(head: &[u8]) -> Cow<'_, str> {
    match std::str::from_utf8(head) {
        Ok(header) => Cow::Borrowed(normalize_qname_view(header)),
        Err(_) => Cow::Owned(normalize_qname(String::from_utf8_lossy(head).as_ref())),
    }
}

fn header_text(record: &impl Record) -> Cow<'_, str> {
    String::from_utf8_lossy(record.head())
}

fn decode_fastq_quality_scores(qualities: &[u8]) -> Result<Vec<u8>> {
    qualities
        .iter()
        .copied()
        .map(decode_fastq_quality_score)
        .collect()
}

fn decode_fastq_quality_score(symbol: u8) -> Result<u8> {
    symbol.checked_sub(33).ok_or_else(|| {
        RvScreenError::validation(
            "input",
            format!("FASTQ quality byte {symbol} is below Phred+33 encoding range"),
        )
    })
}

fn fastq_record_line(pair_index: u64) -> u64 {
    pair_index * FASTQ_RECORD_LINES + 1
}

fn map_fastq_error(path: &Path, error: FastqError) -> RvScreenError {
    let reason = error.to_string();
    match error {
        FastqError::Io(source) => RvScreenError::io(path, source),
        FastqError::UnequalLengths { pos, .. }
        | FastqError::InvalidStart { pos, .. }
        | FastqError::InvalidSep { pos, .. }
        | FastqError::UnexpectedEnd { pos } => RvScreenError::parse(path, pos.line, reason),
        FastqError::BufferLimit => RvScreenError::parse(path, 0, reason),
    }
}

#[cfg(test)]
#[path = "../../tests/testutil/mod.rs"]
mod testutil;

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{write::GzEncoder, Compression};
    use std::fs;
    use std::io::{self, Write};
    use std::path::{Path, PathBuf};
    use tempfile::tempdir;
    use testutil::{generate_fastq_pair, FastqPairConfig, ReadComponent, SyntheticSource};

    #[test]
    fn test_reads_normal_paired_fastq() {
        let tempdir = tempdir().expect("tempdir should be created");
        let (r1_path, r2_path) = generate_fastq_pair(
            &FastqPairConfig::new(100, 75, 42)
                .with_output_dir(tempdir.path().join("generated"))
                .with_components(vec![ReadComponent::new(SyntheticSource::Human, 1.0)]),
        )
        .expect("FASTQ generation should succeed");

        let mut reader = FastqPairReader::open(&r1_path, &r2_path)
            .expect("paired FASTQ reader should open generated files");
        let mut pair_count = 0usize;

        while let Some(record) = reader.next_pair() {
            let record = record.expect("generated FASTQ should read cleanly");
            assert!(!record.fragment_key.is_empty());
            assert_eq!(record.r1_seq.len(), 75);
            assert_eq!(record.r2_seq.len(), 75);
            assert_eq!(record.r1_qual.len(), 75);
            assert_eq!(record.r2_qual.len(), 75);
            assert!(record.r1_qual.iter().all(|score| *score <= 93));
            assert!(record.r2_qual.iter().all(|score| *score <= 93));
            pair_count += 1;
        }

        assert_eq!(pair_count, 100);
    }

    #[test]
    fn test_normalize_qname_handles_suffix_and_illumina_markers() {
        assert_eq!(
            normalize_qname("INST:1:FC:1:1:100:200/1"),
            "INST:1:FC:1:1:100:200"
        );
        assert_eq!(
            normalize_qname("INST:1:FC:1:1:100:200/2"),
            "INST:1:FC:1:1:100:200"
        );
        assert_eq!(
            normalize_qname("INST:1:FC:1:1:100:200 1:N:0:ACGT"),
            "INST:1:FC:1:1:100:200"
        );
        assert_eq!(
            normalize_qname("INST:1:FC:1:1:100:200 2:N:0:ACGT"),
            "INST:1:FC:1:1:100:200"
        );
    }

    #[test]
    fn test_reader_accepts_normalized_suffix_and_illumina_pairs() {
        let tempdir = tempdir().expect("tempdir should be created");
        let generated_dir = tempdir.path().join("generated");
        let (r1_path, r2_path) = generate_fastq_pair(
            &FastqPairConfig::new(2, 75, 7)
                .with_output_dir(&generated_dir)
                .with_components(vec![ReadComponent::new(SyntheticSource::Human, 1.0)]),
        )
        .expect("FASTQ generation should succeed");

        let suffix_dir = tempdir.path().join("suffix");
        let (suffix_r1, suffix_r2) = rewrite_pair_headers(
            &r1_path,
            &r2_path,
            &suffix_dir,
            |_, header| format!("{header}/1"),
            |_, header| format!("{header}/2"),
        )
        .expect("suffix FASTQ rewrite should succeed");
        let suffix_records = read_all_pairs(&suffix_r1, &suffix_r2)
            .expect("reader should accept /1 and /2 normalized pairs");
        assert_eq!(suffix_records.len(), 2);
        assert!(suffix_records
            .iter()
            .all(|record| record.fragment_key.split(':').count() == 7));

        let illumina_dir = tempdir.path().join("illumina");
        let (illumina_r1, illumina_r2) = rewrite_pair_headers(
            &r1_path,
            &r2_path,
            &illumina_dir,
            |_, header| format!("{header} 1:N:0:ACGT"),
            |_, header| format!("{header} 2:N:0:ACGT"),
        )
        .expect("Illumina FASTQ rewrite should succeed");
        let illumina_records = read_all_pairs(&illumina_r1, &illumina_r2)
            .expect("reader should accept Illumina mate markers");
        assert_eq!(illumina_records.len(), 2);
        assert!(illumina_records
            .iter()
            .all(|record| record.fragment_key.split(':').count() == 7));
    }

    #[test]
    fn test_pair_mismatch_returns_parse_error() {
        let tempdir = tempdir().expect("tempdir should be created");
        let generated_dir = tempdir.path().join("generated");
        let (r1_path, r2_path) = generate_fastq_pair(
            &FastqPairConfig::new(1, 75, 99)
                .with_output_dir(&generated_dir)
                .with_components(vec![ReadComponent::new(SyntheticSource::Human, 1.0)]),
        )
        .expect("FASTQ generation should succeed");

        let mismatch_dir = tempdir.path().join("mismatch");
        let (_rewritten_r1, rewritten_r2) = rewrite_pair_headers(
            &r1_path,
            &r2_path,
            &mismatch_dir,
            |_, header| header.to_owned(),
            |_, _| "DIFFERENT:1:FCX123:1:1101:9999:9999".to_owned(),
        )
        .expect("mismatch FASTQ rewrite should succeed");

        let mut reader = FastqPairReader::open(&r1_path, &rewritten_r2)
            .expect("reader should open mismatch pair");
        let err = reader
            .next_pair()
            .expect("mismatched pair should return an error item")
            .expect_err("mismatched pair should fail");

        match err {
            RvScreenError::ParseError { line, reason, .. } => {
                assert_eq!(line, 1);
                assert!(reason.contains("mismatched paired-end fragment keys"));
                assert!(reason.contains("DIFFERENT:1:FCX123:1:1101:9999:9999"));
            }
            other => panic!("expected ParseError for mismatch, got {other:?}"),
        }
    }

    #[test]
    fn test_orphaned_reads_are_reported() {
        let tempdir = tempdir().expect("tempdir should be created");
        let generated_dir = tempdir.path().join("generated");
        let (r1_path, r2_path) = generate_fastq_pair(
            &FastqPairConfig::new(2, 75, 123)
                .with_output_dir(&generated_dir)
                .with_components(vec![ReadComponent::new(SyntheticSource::Human, 1.0)]),
        )
        .expect("FASTQ generation should succeed");

        let orphan_dir = tempdir.path().join("orphan");
        fs::create_dir_all(&orphan_dir).expect("orphan dir should be created");
        let orphan_r2 = orphan_dir.join("synthetic_R2.fastq");
        truncate_to_records(&r2_path, 1, &orphan_r2).expect("R2 truncation should succeed");

        let mut reader = FastqPairReader::open(&r1_path, &orphan_r2)
            .expect("reader should open orphan test pair");
        reader
            .next_pair()
            .expect("first pair should exist")
            .expect("first pair should read cleanly");
        let err = reader
            .next_pair()
            .expect("second pair should surface orphan error")
            .expect_err("second pair should fail as orphan");

        match err {
            RvScreenError::ParseError { line, reason, .. } => {
                assert_eq!(line, 5);
                assert!(reason.contains("orphaned read"));
                assert!(reason.contains("ended early"));
            }
            other => panic!("expected ParseError for orphan, got {other:?}"),
        }
    }

    #[test]
    fn test_reads_gzip_fastq_pairs() {
        let tempdir = tempdir().expect("tempdir should be created");
        let generated_dir = tempdir.path().join("generated");
        let (r1_path, r2_path) = generate_fastq_pair(
            &FastqPairConfig::new(8, 75, 2026)
                .with_output_dir(&generated_dir)
                .with_components(vec![ReadComponent::new(SyntheticSource::Human, 1.0)]),
        )
        .expect("FASTQ generation should succeed");

        let gzip_dir = tempdir.path().join("gzip");
        fs::create_dir_all(&gzip_dir).expect("gzip dir should be created");
        let r1_gz = gzip_dir.join("synthetic_R1.fastq.gz");
        let r2_gz = gzip_dir.join("synthetic_R2.fastq.gz");
        gzip_copy(&r1_path, &r1_gz).expect("R1 gzip conversion should succeed");
        gzip_copy(&r2_path, &r2_gz).expect("R2 gzip conversion should succeed");

        let records = read_all_pairs(&r1_gz, &r2_gz).expect("gzip FASTQ pair should read cleanly");
        assert_eq!(records.len(), 8);
    }

    #[test]
    fn test_decodes_fastq_qualities_to_phred_scores() {
        let tempdir = tempdir().expect("tempdir should be created");
        let r1_path = tempdir.path().join("r1.fastq");
        let r2_path = tempdir.path().join("r2.fastq");

        fs::write(&r1_path, "@read/1\nACGT\n+\n?@AB\n").expect("R1 FASTQ should be written");
        fs::write(&r2_path, "@read/2\nTGCA\n+\nBCDE\n").expect("R2 FASTQ should be written");

        let record = FastqPairReader::open(&r1_path, &r2_path)
            .expect("reader should open")
            .next_pair()
            .expect("pair should exist")
            .expect("pair should decode");

        assert_eq!(record.r1_qual, vec![30, 31, 32, 33]);
        assert_eq!(record.r2_qual, vec![33, 34, 35, 36]);
    }

    #[test]
    fn test_borrowed_pair_iteration_skips_quality_decode_until_requested() {
        let tempdir = tempdir().expect("tempdir should be created");
        let r1_path = tempdir.path().join("r1.fastq");
        let r2_path = tempdir.path().join("r2.fastq");

        fs::write(&r1_path, b"@read/1\nACGT\n+\n\x1f\x1f\x1f\x1f\n")
            .expect("R1 FASTQ should be written");
        fs::write(&r2_path, b"@read/2\nTGCA\n+\n\x1f\x1f\x1f\x1f\n")
            .expect("R2 FASTQ should be written");

        let mut borrowed_reader = FastqPairReader::open(&r1_path, &r2_path)
            .expect("reader should open for borrowed iteration");
        let mut visited = 0usize;
        let io_wall_duration = borrowed_reader
            .for_each_borrowed_pair_with_timing(|fragment| {
                assert_eq!(fragment.fragment_key.as_ref(), "read");
                assert_eq!(fragment.r1_seq, b"ACGT");
                assert_eq!(fragment.r2_seq, b"TGCA");
                visited += 1;
                Ok(())
            })
            .expect("borrowed iteration should not decode FASTQ qualities");
        assert_eq!(visited, 1);
        assert!(io_wall_duration > Duration::default());

        let err = FastqPairReader::open(&r1_path, &r2_path)
            .expect("reader should reopen")
            .next_pair()
            .expect("pair should exist")
            .expect_err("owned pair construction should still decode qualities");
        assert!(err
            .to_string()
            .contains("FASTQ quality byte 31 is below Phred+33 encoding range"));
    }

    fn read_all_pairs(r1_path: &Path, r2_path: &Path) -> Result<Vec<FragmentRecord>> {
        let reader = FastqPairReader::open(r1_path, r2_path)?;
        reader.collect()
    }

    fn rewrite_pair_headers<F1, F2>(
        r1_source: &Path,
        r2_source: &Path,
        output_dir: &Path,
        r1_transform: F1,
        r2_transform: F2,
    ) -> io::Result<(PathBuf, PathBuf)>
    where
        F1: Fn(usize, &str) -> String,
        F2: Fn(usize, &str) -> String,
    {
        fs::create_dir_all(output_dir)?;
        let r1_dest = output_dir.join("synthetic_R1.fastq");
        let r2_dest = output_dir.join("synthetic_R2.fastq");
        rewrite_fastq_headers(r1_source, &r1_dest, r1_transform)?;
        rewrite_fastq_headers(r2_source, &r2_dest, r2_transform)?;
        Ok((r1_dest, r2_dest))
    }

    fn rewrite_fastq_headers<F>(source: &Path, dest: &Path, transform: F) -> io::Result<()>
    where
        F: Fn(usize, &str) -> String,
    {
        let text = fs::read_to_string(source)?;
        let mut rewritten = String::with_capacity(text.len() + 64);

        for (line_index, line) in text.lines().enumerate() {
            if line_index % 4 == 0 {
                let header = line.strip_prefix('@').ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "FASTQ header missing '@'")
                })?;
                rewritten.push('@');
                rewritten.push_str(&transform(line_index / 4, header));
            } else {
                rewritten.push_str(line);
            }
            rewritten.push('\n');
        }

        fs::write(dest, rewritten)
    }

    fn truncate_to_records(source: &Path, record_count: usize, dest: &Path) -> io::Result<()> {
        let text = fs::read_to_string(source)?;
        let lines: Vec<_> = text.lines().collect();
        let keep = record_count.checked_mul(4).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "record count overflowed")
        })?;
        if keep > lines.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "cannot keep more FASTQ records than exist",
            ));
        }

        let mut truncated = lines[..keep].join("\n");
        truncated.push('\n');
        fs::write(dest, truncated)
    }

    fn gzip_copy(source: &Path, dest: &Path) -> io::Result<()> {
        let mut encoder = GzEncoder::new(File::create(dest)?, Compression::default());
        encoder.write_all(&fs::read(source)?)?;
        encoder.finish()?;
        Ok(())
    }
}
