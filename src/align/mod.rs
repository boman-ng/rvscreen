pub mod competitive;

pub use competitive::{
    AlignHit, CompetitiveAligner, CompetitivePreset, FragmentAlignResult, VirusTarget,
};

use crate::error::{Result, RvScreenError};
use minimap2::{Aligner, Built, Mapping, Strand};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProbeStrand {
    Forward,
    Reverse,
}

impl From<Strand> for ProbeStrand {
    fn from(value: Strand) -> Self {
        match value {
            Strand::Forward => Self::Forward,
            Strand::Reverse => Self::Reverse,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProbeMapping {
    pub target_name: Option<String>,
    pub mapq: u32,
    pub strand: ProbeStrand,
    pub target_start: i32,
    pub target_end: i32,
    pub nm: Option<i32>,
    pub alignment_score: Option<i32>,
}

pub fn probe_minimap2_alignment<P>(fasta_path: P, read_seq: &[u8]) -> Result<Vec<ProbeMapping>>
where
    P: AsRef<Path>,
{
    let aligner = build_sr_probe_aligner(fasta_path.as_ref())?;
    let mappings = aligner
        .map(
            read_seq,
            false,
            false,
            None,
            None,
            Some(b"probe-single-read"),
        )
        .map_err(|reason| RvScreenError::alignment("probe-single-read", reason))?;

    Ok(extract_probe_mappings(mappings))
}

pub fn probe_paired_alignment<P>(
    fasta_path: P,
    read1_seq: &[u8],
    read2_seq: &[u8],
) -> Result<(Vec<ProbeMapping>, Vec<ProbeMapping>)>
where
    P: AsRef<Path>,
{
    let aligner = build_sr_probe_aligner(fasta_path.as_ref())?;
    let (read1_mappings, read2_mappings) = aligner
        .map_pair(
            read1_seq,
            read2_seq,
            false,
            false,
            None,
            None,
            Some(b"probe-read-pair"),
        )
        .map_err(|reason| RvScreenError::alignment("probe-read-pair", reason))?;

    Ok((
        extract_probe_mappings(read1_mappings),
        extract_probe_mappings(read2_mappings),
    ))
}

fn build_sr_probe_aligner(fasta_path: &Path) -> Result<Aligner<Built>> {
    Aligner::builder()
        .sr()
        .with_cigar()
        .with_index(fasta_path, None)
        .map_err(|reason| {
            RvScreenError::validation(
                "reference_fasta",
                format!(
                    "failed to index `{}` with minimap2 sr preset: {reason}",
                    fasta_path.display()
                ),
            )
        })
}

fn extract_probe_mappings(mappings: Vec<Mapping>) -> Vec<ProbeMapping> {
    mappings.into_iter().map(extract_probe_mapping).collect()
}

fn extract_probe_mapping(mapping: Mapping) -> ProbeMapping {
    let target_name = mapping
        .target_name
        .as_ref()
        .map(|target_name| target_name.as_str().to_owned());
    let nm = mapping.alignment.as_ref().map(|alignment| alignment.nm);
    let alignment_score = mapping
        .alignment
        .as_ref()
        .and_then(|alignment| alignment.alignment_score);

    ProbeMapping {
        target_name,
        mapq: mapping.mapq,
        strand: mapping.strand.into(),
        target_start: mapping.target_start,
        target_end: mapping.target_end,
        nm,
        alignment_score,
    }
}

#[cfg(test)]
#[path = "../../tests/testutil/mod.rs"]
mod testutil;

#[cfg(test)]
mod tests {
    use super::*;
    use seq_io::fastq::{Reader, Record};
    use std::io;
    use std::path::Path;
    use testutil::{
        generate_fastq_pair, generate_mini_reference, mini_reference_fasta_path, FastqPairConfig,
        ReadComponent, SyntheticSource, VirusSelector,
    };

    #[test]
    fn test_single_read_alignment() {
        let reference_dir =
            generate_mini_reference().expect("mini reference generation should succeed");
        let fasta_path = mini_reference_fasta_path(&reference_dir);
        let (r1_path, _) = generate_fastq_pair(
            &FastqPairConfig::new(1, 75, 11)
                .with_components(vec![ReadComponent::new(SyntheticSource::Human, 1.0)]),
        )
        .expect("FASTQ generation should succeed");
        let read_seq = read_first_fastq_sequence(&r1_path).expect("R1 FASTQ should be readable");

        let mappings = probe_minimap2_alignment(&fasta_path, &read_seq)
            .expect("single-read minimap2 probe should succeed");
        let best_mapping =
            best_mapping(&mappings).expect("single read should produce at least one mapping");

        assert_eq!(best_mapping.target_name.as_deref(), Some("mini_human_chr1"));
        assert!(
            best_mapping.mapq > 0,
            "expected positive MAPQ, got {best_mapping:?}"
        );
    }

    #[test]
    fn test_paired_alignment() {
        let reference_dir =
            generate_mini_reference().expect("mini reference generation should succeed");
        let fasta_path = mini_reference_fasta_path(&reference_dir);
        let (r1_path, r2_path) =
            generate_fastq_pair(&FastqPairConfig::new(1, 75, 22).with_components(vec![
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    1.0,
                ),
            ]))
            .expect("paired FASTQ generation should succeed");
        let read1_seq = read_first_fastq_sequence(&r1_path).expect("R1 FASTQ should be readable");
        let read2_seq = read_first_fastq_sequence(&r2_path).expect("R2 FASTQ should be readable");

        let (read1_mappings, read2_mappings) =
            probe_paired_alignment(&fasta_path, &read1_seq, &read2_seq)
                .expect("paired minimap2 probe should succeed");

        assert!(
            !read1_mappings.is_empty(),
            "expected non-empty R1 mappings from minimap2"
        );
        assert!(
            !read2_mappings.is_empty(),
            "expected non-empty R2 mappings from minimap2"
        );

        let read1_best = best_mapping(&read1_mappings).expect("R1 should contain a best mapping");
        let read2_best = best_mapping(&read2_mappings).expect("R2 should contain a best mapping");

        assert_eq!(read1_best.target_name.as_deref(), Some("mini_virus_alpha"));
        assert_eq!(read2_best.target_name.as_deref(), Some("mini_virus_alpha"));
    }

    #[test]
    fn test_probe_metadata_includes_mapq_and_alignment_metrics() {
        let reference_dir =
            generate_mini_reference().expect("mini reference generation should succeed");
        let fasta_path = mini_reference_fasta_path(&reference_dir);
        let (r1_path, _) =
            generate_fastq_pair(&FastqPairConfig::new(1, 75, 33).with_components(vec![
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV2.1".to_string())),
                    1.0,
                ),
            ]))
            .expect("FASTQ generation should succeed");
        let read_seq = read_first_fastq_sequence(&r1_path).expect("R1 FASTQ should be readable");

        let mappings = probe_minimap2_alignment(&fasta_path, &read_seq)
            .expect("single-read minimap2 probe should succeed");
        let best_mapping = best_mapping(&mappings).expect("expected at least one mapping");

        assert_eq!(best_mapping.target_name.as_deref(), Some("mini_virus_beta"));
        assert!(
            best_mapping.mapq > 0,
            "expected positive MAPQ, got {best_mapping:?}"
        );
        assert!(
            best_mapping.target_start >= 0,
            "expected non-negative start, got {best_mapping:?}"
        );
        assert!(
            best_mapping.target_end > best_mapping.target_start,
            "expected end after start, got {best_mapping:?}"
        );
        assert!(
            best_mapping.nm.is_some(),
            "expected NM to be extracted, got {best_mapping:?}"
        );
        assert!(
            best_mapping.alignment_score.is_some(),
            "expected alignment score to be extracted, got {best_mapping:?}"
        );
    }

    fn read_first_fastq_sequence(path: &Path) -> io::Result<Vec<u8>> {
        let mut reader = Reader::from_path(path)?;
        let record = reader
            .next()
            .transpose()
            .map_err(|err| io::Error::new(io::ErrorKind::InvalidData, err.to_string()))?
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "FASTQ did not contain a record",
                )
            })?;
        Ok(record.seq().to_vec())
    }

    fn best_mapping(mappings: &[ProbeMapping]) -> Option<&ProbeMapping> {
        mappings
            .iter()
            .max_by_key(|mapping| (mapping.mapq, mapping.alignment_score.unwrap_or(i32::MIN)))
    }
}
