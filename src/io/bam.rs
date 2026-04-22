use super::fastq::{normalize_qname, FragmentRecord};
use crate::error::{Result, RvScreenError};
use noodles::bam;
use std::collections::VecDeque;
use std::fs::File;
use std::io::Read;
use std::path::Path;

const BAM_MAGIC_NUMBER: [u8; 4] = *b"BAM\x01";

#[derive(Debug)]
pub struct BamFragmentReader {
    buffered_pairs: VecDeque<FragmentRecord>,
}

impl BamFragmentReader {
    pub fn open<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path).map_err(|source| RvScreenError::io(&path, source))?;
        let mut reader = bam::io::Reader::new(file);
        let raw_header = read_raw_bam_header(&mut reader, &path)?;

        reject_coordinate_sorted_input(&path, &raw_header)?;

        let buffered_pairs = decode_all_pairs(&path, &mut reader)?;

        Ok(Self { buffered_pairs })
    }

    pub fn next_pair(&mut self) -> Option<Result<FragmentRecord>> {
        self.buffered_pairs.pop_front().map(Ok)
    }
}

impl Iterator for BamFragmentReader {
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

fn decode_all_pairs<R>(
    path: &Path,
    reader: &mut bam::io::Reader<R>,
) -> Result<VecDeque<FragmentRecord>>
where
    R: Read,
{
    let mut buffered_pairs = VecDeque::new();
    let mut pending_mate = None;

    for (index, result) in reader.records().enumerate() {
        let record = result.map_err(|source| RvScreenError::io(path, source))?;
        let mate = extract_mate_record(path, index as u64 + 1, &record)?;

        if let Some(first) = pending_mate.take() {
            buffered_pairs.push_back(build_fragment_record(path, &first, &mate)?);
        } else {
            pending_mate = Some(mate);
        }
    }

    if let Some(orphan) = pending_mate {
        return Err(orphan_bam_record_error(path, &orphan));
    }

    Ok(buffered_pairs)
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
) -> Result<MateRecord> {
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

    Ok(MateRecord {
        record_number,
        fragment_key: normalize_qname(&qname),
        qname,
        seq: record.sequence().iter().collect(),
        qual: record.quality_scores().iter().collect(),
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
                "adjacent BAM records are not a name-sorted pair: record {} `{}` -> `{}`, record {} `{}` -> `{}`",
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

fn orphan_bam_record_error(path: &Path, orphan: &MateRecord) -> RvScreenError {
    RvScreenError::parse(
        path,
        orphan.record_number,
        format!(
            "orphaned BAM record {} `{}`: adjacent mate is missing; provide name-sorted BAM/uBAM input",
            orphan.record_number, orphan.qname,
        ),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
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

    #[derive(Debug, Clone)]
    struct SyntheticBamRecord {
        qname: String,
        sequence: Vec<u8>,
        quality_scores: Vec<u8>,
    }

    impl SyntheticBamRecord {
        fn new(qname: &str, sequence: &[u8], quality_scores: &[u8]) -> Self {
            Self {
                qname: qname.to_owned(),
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
        let flags = 0x01u16 | 0x04u16 | 0x08u16;
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
