use super::fastq::normalize_qname;
use super::fragment_reader::{
    build_name_sorted_fragment_record, mate_slot_from_flags, orphaned_name_sorted_record_error,
    FragmentMateRecord,
};
use super::FragmentRecord;
use crate::error::{Result, RvScreenError};
use noodles::bam;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

const BAM_MAGIC_NUMBER: [u8; 4] = *b"BAM\x01";

pub struct BamFragmentReader {
    path: PathBuf,
    reader: bam::io::Reader<Box<dyn Read>>,
    record: bam::Record,
    pending_mate: Option<FragmentMateRecord>,
    record_number: u64,
    exhausted: bool,
}

impl std::fmt::Debug for BamFragmentReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BamFragmentReader")
            .field("path", &self.path)
            .field("record", &self.record)
            .field("pending_mate", &self.pending_mate)
            .field("record_number", &self.record_number)
            .field("exhausted", &self.exhausted)
            .finish_non_exhaustive()
    }
}

impl BamFragmentReader {
    pub fn open<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let path = path.as_ref().to_path_buf();
        let mut reader = open_bam_reader(&path)?;
        let raw_header = read_raw_bam_header(&mut reader, &path)?;

        reject_coordinate_sorted_input(&path, &raw_header)?;

        Ok(Self {
            path,
            reader,
            record: bam::Record::default(),
            pending_mate: None,
            record_number: 0,
            exhausted: false,
        })
    }

    pub fn next_pair(&mut self) -> Option<Result<FragmentRecord>> {
        if self.exhausted {
            return None;
        }

        loop {
            match self.read_next_mate() {
                Ok(Some(mate)) => {
                    if let Some(first) = self.pending_mate.take() {
                        let fragment = match build_fragment_record(&self.path, &first, &mate) {
                            Ok(fragment) => fragment,
                            Err(err) => {
                                self.exhausted = true;
                                return Some(Err(err));
                            }
                        };

                        return Some(Ok(fragment));
                    }

                    self.pending_mate = Some(mate);
                }
                Ok(None) => {
                    self.exhausted = true;

                    return self
                        .pending_mate
                        .take()
                        .map(|orphan| Err(orphan_bam_record_error(&self.path, &orphan)));
                }
                Err(err) => {
                    self.exhausted = true;
                    return Some(Err(err));
                }
            }
        }
    }

    fn read_next_mate(&mut self) -> Result<Option<FragmentMateRecord>> {
        match self.reader.read_record(&mut self.record) {
            Ok(0) => Ok(None),
            Ok(_) => {
                self.record_number = self.record_number.saturating_add(1);
                extract_mate_record(&self.path, self.record_number, &self.record).map(Some)
            }
            Err(source) => Err(RvScreenError::io(&self.path, source)),
        }
    }
}

impl Iterator for BamFragmentReader {
    type Item = Result<FragmentRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_pair()
    }
}

fn read_raw_bam_header<R>(reader: &mut bam::io::Reader<R>, path: &Path) -> Result<String>
where
    R: Read,
{
    let mut header_reader = reader.header_reader();
    let magic = header_reader
        .read_magic_number()
        .map_err(|source| RvScreenError::io(path, source))?;

    if magic != BAM_MAGIC_NUMBER {
        return Err(RvScreenError::parse(
            path,
            0,
            format!("invalid BAM magic number: expected {BAM_MAGIC_NUMBER:?}, got {magic:?}"),
        ));
    }

    let mut raw_sam_header_reader = header_reader
        .raw_sam_header_reader()
        .map_err(|source| RvScreenError::io(path, source))?;
    let mut raw_header = String::new();
    raw_sam_header_reader
        .read_to_string(&mut raw_header)
        .map_err(|source| RvScreenError::io(path, source))?;
    raw_sam_header_reader
        .discard_to_end()
        .map_err(|source| RvScreenError::io(path, source))?;
    header_reader
        .read_reference_sequences()
        .map_err(|source| RvScreenError::io(path, source))?;

    Ok(raw_header)
}

fn open_bam_reader(path: &Path) -> Result<bam::io::Reader<Box<dyn Read>>> {
    let file = File::open(path).map_err(|source| RvScreenError::io(path, source))?;
    let decoder: Box<dyn Read> = Box::new(bam::io::Reader::new(file).into_inner());
    Ok(bam::io::Reader::from(decoder))
}

fn reject_coordinate_sorted_input(path: &Path, raw_header: &str) -> Result<()> {
    if matches!(parse_sort_order(raw_header), Some("coordinate")) {
        return Err(RvScreenError::validation(
            "input",
            format!(
                "BAM input `{}` is coordinate-sorted; provide name-sorted BAM/uBAM or FASTQ instead",
                path.display()
            ),
        ));
    }

    Ok(())
}

fn parse_sort_order(raw_header: &str) -> Option<&str> {
    raw_header
        .lines()
        .find(|line| line.starts_with("@HD\t"))
        .and_then(|line| line.split('\t').find_map(|field| field.strip_prefix("SO:")))
}

fn extract_mate_record(
    path: &Path,
    record_number: u64,
    record: &bam::Record,
) -> Result<FragmentMateRecord> {
    let qname = record
        .name()
        .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
        .ok_or_else(|| {
            RvScreenError::parse(
                path,
                record_number,
                format!("BAM record {record_number} is missing QNAME"),
            )
        })?;

    let fragment_key = normalize_qname(&qname);

    Ok(FragmentMateRecord::new(
        record_number,
        qname,
        fragment_key,
        record.sequence().iter().collect(),
        record.quality_scores().iter().collect(),
        mate_slot_from_flags(record.flags()),
    ))
}

fn build_fragment_record(
    path: &Path,
    first: &FragmentMateRecord,
    second: &FragmentMateRecord,
) -> Result<FragmentRecord> {
    build_name_sorted_fragment_record(path, "BAM", first, second)
}

fn orphan_bam_record_error(path: &Path, orphan: &FragmentMateRecord) -> RvScreenError {
    orphaned_name_sorted_record_error(
        path,
        "BAM",
        orphan,
        "provide name-sorted BAM/uBAM input",
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam;
    use std::fs;
    use std::io::{self, Write};
    use tempfile::tempdir;

    #[test]
    fn test_reads_name_sorted_adjacent_bam_pairs() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bam_path = tempdir.path().join("queryname.bam");

        write_test_bam(
            &bam_path,
            "queryname",
            &[
                SyntheticBamRecord::new(
                    "READ:1:FCX123:1:1101:1000:2000",
                    b"ACGT",
                    &[30, 31, 32, 33],
                ),
                SyntheticBamRecord::new(
                    "READ:1:FCX123:1:1101:1000:2000",
                    b"TGCA",
                    &[34, 35, 36, 37],
                ),
                SyntheticBamRecord::new(
                    "READ:1:FCX123:1:1101:1001:2001",
                    b"GGTT",
                    &[20, 21, 22, 23],
                ),
                SyntheticBamRecord::new(
                    "READ:1:FCX123:1:1101:1001:2001",
                    b"AACC",
                    &[24, 25, 26, 27],
                ),
            ],
        )
        .expect("BAM fixture should be written");

        let records: Vec<_> = BamFragmentReader::open(&bam_path)
            .expect("reader should open name-sorted BAM")
            .collect::<Result<_>>()
            .expect("name-sorted BAM should read cleanly");

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].fragment_key, "READ:1:FCX123:1:1101:1000:2000");
        assert_eq!(records[0].r1_seq, b"ACGT");
        assert_eq!(records[0].r1_qual, vec![30, 31, 32, 33]);
        assert_eq!(records[0].r2_seq, b"TGCA");
        assert_eq!(records[0].r2_qual, vec![34, 35, 36, 37]);
        assert_eq!(records[1].fragment_key, "READ:1:FCX123:1:1101:1001:2001");
    }

    #[test]
    fn test_rejects_coordinate_sorted_bam() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bam_path = tempdir.path().join("coordinate.bam");

        write_test_bam(
            &bam_path,
            "coordinate",
            &[
                SyntheticBamRecord::new(
                    "READ:1:FCX123:1:1101:1000:2000",
                    b"ACGT",
                    &[30, 31, 32, 33],
                ),
                SyntheticBamRecord::new(
                    "READ:1:FCX123:1:1101:1000:2000",
                    b"TGCA",
                    &[34, 35, 36, 37],
                ),
            ],
        )
        .expect("BAM fixture should be written");

        let err =
            BamFragmentReader::open(&bam_path).expect_err("coordinate-sorted BAM must be rejected");

        match err {
            RvScreenError::ValidationError { field, reason } => {
                assert_eq!(field, "input");
                assert!(reason.contains("coordinate-sorted"));
                assert!(reason.contains("name-sorted"));
            }
            other => panic!("expected ValidationError, got {other:?}"),
        }
    }

    #[test]
    fn test_reads_unmapped_bam_through_bam_path() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bam_path = tempdir.path().join("sample.ubam");

        write_test_bam(
            &bam_path,
            "queryname",
            &[
                SyntheticBamRecord::new(
                    "UBAM:1:FCX123:1:1101:1000:2000",
                    b"NNNN",
                    &[10, 11, 12, 13],
                ),
                SyntheticBamRecord::new(
                    "UBAM:1:FCX123:1:1101:1000:2000",
                    b"TTTT",
                    &[14, 15, 16, 17],
                ),
            ],
        )
        .expect("uBAM fixture should be written");

        let record = BamFragmentReader::open(&bam_path)
            .expect("reader should open uBAM through BAM path")
            .next_pair()
            .expect("uBAM should yield one pair")
            .expect("uBAM pair should decode");

        assert_eq!(record.fragment_key, "UBAM:1:FCX123:1:1101:1000:2000");
        assert_eq!(record.r1_seq, b"NNNN");
        assert_eq!(record.r2_seq, b"TTTT");
    }

    #[test]
    fn test_streaming_reader_defers_orphan_detection_until_iteration() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bam_path = tempdir.path().join("orphan-late.bam");

        write_test_bam(
            &bam_path,
            "queryname",
            &[
                SyntheticBamRecord::new("READ:1:FCX123:1:1101:1000:2000", b"ACGT", &[30, 31, 32, 33]),
                SyntheticBamRecord::new("READ:1:FCX123:1:1101:1000:2000", b"TGCA", &[34, 35, 36, 37]),
                SyntheticBamRecord::new("READ:1:FCX123:1:1101:1001:2001", b"GGTT", &[20, 21, 22, 23]),
            ],
        )
        .expect("BAM fixture should be written");

        let mut reader = BamFragmentReader::open(&bam_path)
            .expect("reader should open without preloading the whole BAM");

        let first = reader
            .next_pair()
            .expect("first pair should be available")
            .expect("first pair should decode");
        assert_eq!(first.fragment_key, "READ:1:FCX123:1:1101:1000:2000");

        let err = reader
            .next_pair()
            .expect("orphan should surface when iteration reaches EOF")
            .expect_err("orphaned BAM record should fail");

        match err {
            RvScreenError::ParseError { line, reason, .. } => {
                assert_eq!(line, 3);
                assert!(reason.contains("orphaned BAM record 3"));
            }
            other => panic!("expected ParseError, got {other:?}"),
        }

        assert!(reader.next_pair().is_none());
    }

    #[test]
    fn test_streaming_reader_reports_non_adjacent_pair_mismatch() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bam_path = tempdir.path().join("mismatch-late.bam");

        write_test_bam(
            &bam_path,
            "queryname",
            &[
                SyntheticBamRecord::new("READ:1:FCX123:1:1101:1000:2000", b"ACGT", &[30, 31, 32, 33]),
                SyntheticBamRecord::new("READ:1:FCX123:1:1101:1000:2000", b"TGCA", &[34, 35, 36, 37]),
                SyntheticBamRecord::new("READ:1:FCX123:1:1101:1001:2001", b"GGTT", &[20, 21, 22, 23]),
                SyntheticBamRecord::new("READ:1:FCX123:1:1101:1002:2002", b"AACC", &[24, 25, 26, 27]),
            ],
        )
        .expect("BAM fixture should be written");

        let mut reader = BamFragmentReader::open(&bam_path)
            .expect("reader should open without pre-validating later records");

        let first = reader
            .next_pair()
            .expect("first pair should be available")
            .expect("first pair should decode");
        assert_eq!(first.fragment_key, "READ:1:FCX123:1:1101:1000:2000");

        let err = reader
            .next_pair()
            .expect("mismatched adjacent pair should surface during iteration")
            .expect_err("mismatched fragment keys must fail");

        match err {
            RvScreenError::ParseError { line, reason, .. } => {
                assert_eq!(line, 3);
                assert!(reason.contains("adjacent BAM records are not a name-sorted pair"));
                assert!(reason.contains("1001:2001"));
                assert!(reason.contains("1002:2002"));
            }
            other => panic!("expected ParseError, got {other:?}"),
        }

        assert!(reader.next_pair().is_none());
    }

    #[test]
    fn test_streaming_reader_orders_flagged_mates_into_r1_r2_slots() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bam_path = tempdir.path().join("slot-ordered.bam");

        write_test_bam(
            &bam_path,
            "queryname",
            &[
                SyntheticBamRecord::last(
                    "READ:1:FCX123:1:1101:1000:2000",
                    b"TGCA",
                    &[34, 35, 36, 37],
                ),
                SyntheticBamRecord::first(
                    "READ:1:FCX123:1:1101:1000:2000",
                    b"ACGT",
                    &[30, 31, 32, 33],
                ),
            ],
        )
        .expect("BAM fixture should be written");

        let record = BamFragmentReader::open(&bam_path)
            .expect("reader should open BAM")
            .next_pair()
            .expect("pair should be available")
            .expect("pair should decode");

        assert_eq!(record.r1_seq, b"ACGT");
        assert_eq!(record.r1_qual, vec![30, 31, 32, 33]);
        assert_eq!(record.r2_seq, b"TGCA");
        assert_eq!(record.r2_qual, vec![34, 35, 36, 37]);
    }

    #[derive(Debug, Clone)]
    struct SyntheticBamRecord {
        qname: String,
        flags: sam::alignment::record::Flags,
        sequence: Vec<u8>,
        quality_scores: Vec<u8>,
    }

    impl SyntheticBamRecord {
        fn new(qname: &str, sequence: &[u8], quality_scores: &[u8]) -> Self {
            Self::with_flags(
                qname,
                sam::alignment::record::Flags::SEGMENTED
                    | sam::alignment::record::Flags::UNMAPPED
                    | sam::alignment::record::Flags::MATE_UNMAPPED,
                sequence,
                quality_scores,
            )
        }

        fn first(qname: &str, sequence: &[u8], quality_scores: &[u8]) -> Self {
            Self::with_flags(
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
            Self::with_flags(
                qname,
                sam::alignment::record::Flags::SEGMENTED
                    | sam::alignment::record::Flags::UNMAPPED
                    | sam::alignment::record::Flags::MATE_UNMAPPED
                    | sam::alignment::record::Flags::LAST_SEGMENT,
                sequence,
                quality_scores,
            )
        }

        fn with_flags(
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
    }

    fn write_test_bam(
        path: &Path,
        sort_order: &str,
        records: &[SyntheticBamRecord],
    ) -> io::Result<()> {
        let file = File::create(path)?;
        let mut writer = bam::io::Writer::new(file);

        writer.get_mut().write_all(&encode_bam_header(sort_order))?;
        for record in records {
            writer.get_mut().write_all(&encode_bam_record(record)?)?;
        }

        writer.try_finish()?;

        Ok(())
    }

    fn encode_bam_header(sort_order: &str) -> Vec<u8> {
        let raw_header = format!("@HD\tVN:1.6\tSO:{sort_order}\n");
        let mut encoded = Vec::new();
        encoded.extend(BAM_MAGIC_NUMBER);
        encoded.extend((raw_header.len() as u32).to_le_bytes());
        encoded.extend(raw_header.as_bytes());
        encoded.extend(0u32.to_le_bytes());
        encoded
    }

    fn encode_bam_record(record: &SyntheticBamRecord) -> io::Result<Vec<u8>> {
        if record.sequence.len() != record.quality_scores.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "sequence length must match quality-score length",
            ));
        }

        let read_name = format!("{}\0", record.qname).into_bytes();
        let seq_len = i32::try_from(record.sequence.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let read_name_len = u8::try_from(read_name.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let bin_mq_nl =
            u32::from(read_name_len) | (u32::from(255u8) << 8) | (u32::from(4680u16) << 16);
        let flags = u16::from(record.flags);
        let flag_nc = u32::from(flags) << 16;
        let packed_sequence = pack_sequence(&record.sequence)?;

        let mut body = Vec::new();
        body.extend((-1i32).to_le_bytes());
        body.extend((-1i32).to_le_bytes());
        body.extend(bin_mq_nl.to_le_bytes());
        body.extend(flag_nc.to_le_bytes());
        body.extend(seq_len.to_le_bytes());
        body.extend((-1i32).to_le_bytes());
        body.extend((-1i32).to_le_bytes());
        body.extend(0i32.to_le_bytes());
        body.extend(read_name);
        body.extend(packed_sequence);
        body.extend(&record.quality_scores);

        let block_size = u32::try_from(body.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let mut encoded = Vec::with_capacity(body.len() + 4);
        encoded.extend(block_size.to_le_bytes());
        encoded.extend(body);
        Ok(encoded)
    }

    fn pack_sequence(sequence: &[u8]) -> io::Result<Vec<u8>> {
        let mut packed = Vec::with_capacity(sequence.len().div_ceil(2));

        for chunk in sequence.chunks(2) {
            let left = encode_base(chunk[0])? << 4;
            let right = chunk
                .get(1)
                .copied()
                .map(encode_base)
                .transpose()?
                .unwrap_or(0);
            packed.push(left | right);
        }

        Ok(packed)
    }

    fn encode_base(base: u8) -> io::Result<u8> {
        match base.to_ascii_uppercase() {
            b'=' => Ok(0),
            b'A' => Ok(1),
            b'C' => Ok(2),
            b'G' => Ok(4),
            b'T' => Ok(8),
            b'N' => Ok(15),
            other => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("unsupported BAM test base `{}`", char::from(other)),
            )),
        }
    }

    #[test]
    fn test_pack_sequence_round_trip_layout() {
        assert_eq!(
            pack_sequence(b"ACGT").expect("packing should succeed"),
            vec![0x12, 0x48]
        );
    }

    #[test]
    fn test_test_bam_writer_creates_nonempty_file() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bam_path = tempdir.path().join("fixture.bam");

        write_test_bam(
            &bam_path,
            "queryname",
            &[SyntheticBamRecord::new(
                "READ:1",
                b"ACGT",
                &[30, 31, 32, 33],
            )],
        )
        .expect("fixture BAM should be written");

        assert!(
            fs::metadata(&bam_path)
                .expect("metadata should be readable")
                .len()
                > 0
        );
    }
}
