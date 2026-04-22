use rvscreen::types::{BundleManifest, ContigEntry};
use seq_io::fastq::Record;
use sha2::{Digest, Sha256};
use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::{self, BufWriter, ErrorKind, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
use std::time::{SystemTime, UNIX_EPOCH};

const MINI_REFERENCE_PREFIX: &str = "rvscreen-mini-reference";
const FASTQ_PREFIX: &str = "rvscreen-fastq-pair";
const MINI_REFERENCE_FASTA: &str = "mini_reference.fa";
const MINI_REFERENCE_MANIFEST: &str = "manifest.json";
const ADAPTER_SEQUENCE: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACGATCGATCGTACGTCAGTCAGTCAGT";
const LOW_COMPLEXITY_MOTIFS: [&str; 4] = ["AAAAAAAA", "ATATATAT", "CCCCGGGG", "GCGCGCGC"];
static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

#[derive(Debug, Clone, PartialEq)]
pub struct FastqPairConfig {
    pub read_pairs: usize,
    pub read_length: usize,
    pub phred_score: u8,
    pub seed: u64,
    pub output_dir: Option<PathBuf>,
    pub components: Vec<ReadComponent>,
}

impl Default for FastqPairConfig {
    fn default() -> Self {
        Self {
            read_pairs: 32,
            read_length: 75,
            phred_score: 35,
            seed: 42,
            output_dir: None,
            components: vec![ReadComponent::new(SyntheticSource::Human, 1.0)],
        }
    }
}

impl FastqPairConfig {
    pub fn new(read_pairs: usize, read_length: usize, seed: u64) -> Self {
        Self {
            read_pairs,
            read_length,
            seed,
            ..Self::default()
        }
    }

    pub fn with_output_dir(mut self, output_dir: impl Into<PathBuf>) -> Self {
        self.output_dir = Some(output_dir.into());
        self
    }

    pub fn with_phred_score(mut self, phred_score: u8) -> Self {
        self.phred_score = phred_score;
        self
    }

    pub fn with_components(mut self, components: Vec<ReadComponent>) -> Self {
        self.components = components;
        self
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReadComponent {
    pub source: SyntheticSource,
    pub fraction: f64,
}

impl ReadComponent {
    pub fn new(source: SyntheticSource, fraction: f64) -> Self {
        Self { source, fraction }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum SyntheticSource {
    Human,
    Virus(VirusSelector),
    LowComplexity,
    Adapter,
}

#[derive(Debug, Clone, PartialEq)]
pub enum VirusSelector {
    Accession(String),
    SourceType(String),
}

#[derive(Debug, Clone)]
struct MiniReferenceRecord {
    entry: ContigEntry,
    sequence: String,
}

pub fn generate_fastq_pair(config: &FastqPairConfig) -> io::Result<(PathBuf, PathBuf)> {
    validate_fastq_config(config)?;

    let reference = mini_reference_records();
    let counts = allocate_component_counts(config.read_pairs, &config.components)?;
    let quality = quality_string(config.read_length, config.phred_score)?;
    let output_dir = prepare_output_dir(config.output_dir.as_deref(), FASTQ_PREFIX)?;
    let r1_path = output_dir.join("synthetic_R1.fastq");
    let r2_path = output_dir.join("synthetic_R2.fastq");

    let mut rng = SplitMix64::new(config.seed);
    let mut r1_writer = BufWriter::new(File::create(&r1_path)?);
    let mut r2_writer = BufWriter::new(File::create(&r2_path)?);
    let mut pair_index = 0usize;

    for (component, count) in config.components.iter().zip(counts.into_iter()) {
        for _ in 0..count {
            let fragment =
                synthesize_fragment(component, config.read_length, &reference, &mut rng)?;
            let qname = illumina_qname(&component.source, pair_index);
            let (r1_seq, r2_seq) = split_paired_fragment(&fragment, config.read_length)?;

            write_fastq_record(&mut r1_writer, &qname, &r1_seq, &quality)?;
            write_fastq_record(&mut r2_writer, &qname, &r2_seq, &quality)?;
            pair_index += 1;
        }
    }

    r1_writer.flush()?;
    r2_writer.flush()?;

    Ok((r1_path, r2_path))
}

pub fn generate_mini_reference() -> io::Result<PathBuf> {
    let reference_dir = prepare_output_dir(None, MINI_REFERENCE_PREFIX)?;
    let fasta_path = mini_reference_fasta_path(&reference_dir);
    let manifest_path = mini_reference_manifest_path(&reference_dir);
    let records = mini_reference_records();

    let mut fasta_writer = BufWriter::new(File::create(&fasta_path)?);
    for record in &records {
        writeln!(fasta_writer, ">{}", record.entry.contig)?;
        for chunk in record.sequence.as_bytes().chunks(80) {
            fasta_writer.write_all(chunk)?;
            fasta_writer.write_all(b"\n")?;
        }
    }
    fasta_writer.flush()?;

    let manifest = BundleManifest(records.into_iter().map(|record| record.entry).collect());
    let manifest_json = serde_json::to_string_pretty(&manifest).map_err(io::Error::other)?;
    fs::write(&manifest_path, manifest_json)?;

    Ok(reference_dir)
}

pub fn mini_reference_fasta_path(reference_dir: &Path) -> PathBuf {
    reference_dir.join(MINI_REFERENCE_FASTA)
}

pub fn mini_reference_manifest_path(reference_dir: &Path) -> PathBuf {
    reference_dir.join(MINI_REFERENCE_MANIFEST)
}

fn validate_fastq_config(config: &FastqPairConfig) -> io::Result<()> {
    if config.read_pairs == 0 {
        return Err(invalid_input("read_pairs must be greater than zero"));
    }

    if config.read_length == 0 {
        return Err(invalid_input("read_length must be greater than zero"));
    }

    if config.components.is_empty() {
        return Err(invalid_input("at least one read component is required"));
    }

    let total_fraction: f64 = config
        .components
        .iter()
        .map(|component| component.fraction)
        .sum();
    if config
        .components
        .iter()
        .any(|component| !component.fraction.is_finite() || component.fraction <= 0.0)
    {
        return Err(invalid_input(
            "all component fractions must be finite and > 0",
        ));
    }

    if (total_fraction - 1.0).abs() > 1e-9 {
        return Err(invalid_input("component fractions must sum to 1.0"));
    }

    if config.phred_score > 93 {
        return Err(invalid_input(
            "phred_score must be in 0..=93 for Phred+33 FASTQ",
        ));
    }

    let reference = mini_reference_records();
    for component in &config.components {
        match &component.source {
            SyntheticSource::Human => {
                let human = reference
                    .iter()
                    .find(|record| record.entry.group == "human")
                    .expect("mini reference must include a human contig");
                ensure_supports_read_length(human.sequence.len(), config.read_length, "human")?;
            }
            SyntheticSource::Virus(selector) => {
                let virus = select_virus_record(&reference, selector)?;
                ensure_supports_read_length(
                    virus.sequence.len(),
                    config.read_length,
                    &format!("virus selector {selector:?}"),
                )?;
            }
            SyntheticSource::LowComplexity | SyntheticSource::Adapter => {}
        }
    }

    Ok(())
}

fn ensure_supports_read_length(
    sequence_len: usize,
    read_length: usize,
    label: &str,
) -> io::Result<()> {
    if read_length * 2 > sequence_len {
        return Err(invalid_input(&format!(
            "read_length {read_length} is too long for {label} sequence length {sequence_len}"
        )));
    }

    Ok(())
}

fn allocate_component_counts(
    read_pairs: usize,
    components: &[ReadComponent],
) -> io::Result<Vec<usize>> {
    let mut counts = vec![0usize; components.len()];
    let mut remainders = Vec::with_capacity(components.len());
    let mut assigned = 0usize;

    for (index, component) in components.iter().enumerate() {
        let exact = component.fraction * read_pairs as f64;
        let whole = exact.floor() as usize;
        counts[index] = whole;
        assigned += whole;
        remainders.push((index, exact - whole as f64));
    }

    let remaining = read_pairs
        .checked_sub(assigned)
        .ok_or_else(|| invalid_input("component allocation overflowed total read pairs"))?;

    remainders.sort_by(|left, right| match right.1.total_cmp(&left.1) {
        Ordering::Equal => left.0.cmp(&right.0),
        order => order,
    });

    for (index, _) in remainders.into_iter().take(remaining) {
        counts[index] += 1;
    }

    Ok(counts)
}

fn synthesize_fragment(
    component: &ReadComponent,
    read_length: usize,
    reference: &[MiniReferenceRecord],
    rng: &mut SplitMix64,
) -> io::Result<String> {
    let fragment_length = read_length * 2;

    match &component.source {
        SyntheticSource::Human => {
            let human = reference
                .iter()
                .find(|record| record.entry.group == "human")
                .expect("mini reference must include a human contig");
            sample_reference_fragment(&human.sequence, fragment_length, rng)
        }
        SyntheticSource::Virus(selector) => {
            let virus = select_virus_record(reference, selector)?;
            sample_reference_fragment(&virus.sequence, fragment_length, rng)
        }
        SyntheticSource::LowComplexity => Ok(repeat_pattern(
            LOW_COMPLEXITY_MOTIFS[rng.next_usize(LOW_COMPLEXITY_MOTIFS.len())],
            fragment_length,
        )),
        SyntheticSource::Adapter => Ok(repeat_pattern(ADAPTER_SEQUENCE, fragment_length)),
    }
}

fn sample_reference_fragment(
    sequence: &str,
    fragment_length: usize,
    rng: &mut SplitMix64,
) -> io::Result<String> {
    if fragment_length > sequence.len() {
        return Err(invalid_input(
            "requested fragment exceeds source sequence length",
        ));
    }

    let max_start = sequence.len() - fragment_length;
    let start = if max_start == 0 {
        0
    } else {
        rng.next_usize(max_start + 1)
    };

    Ok(sequence[start..start + fragment_length].to_string())
}

fn select_virus_record<'a>(
    reference: &'a [MiniReferenceRecord],
    selector: &VirusSelector,
) -> io::Result<&'a MiniReferenceRecord> {
    reference
        .iter()
        .find(|record| {
            record.entry.group == "virus"
                && match selector {
                    VirusSelector::Accession(accession) => record.entry.accession == *accession,
                    VirusSelector::SourceType(source_type) => {
                        record.entry.source_type == *source_type
                    }
                }
        })
        .ok_or_else(|| invalid_input(&format!("no virus contig matched selector {selector:?}")))
}

fn split_paired_fragment(fragment: &str, read_length: usize) -> io::Result<(String, String)> {
    if fragment.len() < read_length * 2 {
        return Err(invalid_input(
            "fragment must be at least 2x the read length",
        ));
    }

    let r1 = fragment[..read_length].to_string();
    let r2 = reverse_complement(&fragment[fragment.len() - read_length..]);
    Ok((r1, r2))
}

fn illumina_qname(source: &SyntheticSource, pair_index: usize) -> String {
    let instrument = match source {
        SyntheticSource::Human => "HUMAN",
        SyntheticSource::Virus(_) => "VIRUS",
        SyntheticSource::LowComplexity => "LOWCX",
        SyntheticSource::Adapter => "ADAPT",
    };
    let tile = 1101 + (pair_index % 50);
    let x = 1_000 + pair_index;
    let y = 2_000 + pair_index;

    format!("{instrument}:1:FCX123:1:{tile}:{x}:{y}")
}

fn quality_string(read_length: usize, phred_score: u8) -> io::Result<String> {
    let ascii = phred_score
        .checked_add(33)
        .ok_or_else(|| invalid_input("phred_score overflowed Phred+33 encoding"))?;
    if ascii > 126 {
        return Err(invalid_input("Phred+33 quality must map to ASCII <= 126"));
    }

    Ok(std::iter::repeat_n(ascii as char, read_length).collect())
}

fn write_fastq_record(
    writer: &mut impl Write,
    qname: &str,
    sequence: &str,
    quality: &str,
) -> io::Result<()> {
    writeln!(writer, "@{qname}")?;
    writeln!(writer, "{sequence}")?;
    writeln!(writer, "+")?;
    writeln!(writer, "{quality}")?;
    Ok(())
}

fn prepare_output_dir(output_dir: Option<&Path>, prefix: &str) -> io::Result<PathBuf> {
    let path = match output_dir {
        Some(path) => path.to_path_buf(),
        None => unique_temp_dir(prefix),
    };
    fs::create_dir_all(&path)?;
    Ok(path)
}

fn unique_temp_dir(prefix: &str) -> PathBuf {
    let suffix = TEMP_COUNTER.fetch_add(1, AtomicOrdering::Relaxed);
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_nanos())
        .unwrap_or_default();
    std::env::temp_dir().join(format!("{prefix}-{}-{nanos}-{suffix}", std::process::id()))
}

fn mini_reference_records() -> Vec<MiniReferenceRecord> {
    vec![
        MiniReferenceRecord {
            entry: ContigEntry {
                contig: "mini_human_chr1".to_string(),
                accession: "CHM13-MINI-001".to_string(),
                taxid: 9606,
                virus_name: "human_background".to_string(),
                segment: None,
                group: "human".to_string(),
                genome_length: 1024,
                source_release: "mini-2026-04".to_string(),
                source_type: "host-genome".to_string(),
                masked_regions: vec![],
            },
            sequence: random_dna(1024, 0xA1B2_C3D4_E5F6_0101),
        },
        MiniReferenceRecord {
            entry: ContigEntry {
                contig: "mini_virus_alpha".to_string(),
                accession: "NC_SYNTHV1.1".to_string(),
                taxid: 2_697_049,
                virus_name: "Synthetic virus alpha".to_string(),
                segment: Some("segment-A".to_string()),
                group: "virus".to_string(),
                genome_length: 512,
                source_release: "mini-2026-04".to_string(),
                source_type: "refseq-virus".to_string(),
                masked_regions: vec![],
            },
            sequence: random_dna(512, 0x1100_2200_3300_4400),
        },
        MiniReferenceRecord {
            entry: ContigEntry {
                contig: "mini_virus_beta".to_string(),
                accession: "NC_SYNTHV2.1".to_string(),
                taxid: 11_103,
                virus_name: "Synthetic virus beta".to_string(),
                segment: Some("segment-B".to_string()),
                group: "virus".to_string(),
                genome_length: 544,
                source_release: "mini-2026-04".to_string(),
                source_type: "spike-in-virus".to_string(),
                masked_regions: vec![],
            },
            sequence: random_dna(544, 0x5500_6600_7700_8800),
        },
        MiniReferenceRecord {
            entry: ContigEntry {
                contig: "mini_decoy_adapter".to_string(),
                accession: "DECOY-MINI-001".to_string(),
                taxid: 32_630,
                virus_name: "adapter_decoy".to_string(),
                segment: None,
                group: "decoy".to_string(),
                genome_length: 224,
                source_release: "mini-2026-04".to_string(),
                source_type: "adapter-decoy".to_string(),
                masked_regions: vec![],
            },
            sequence: repeat_pattern("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTAACCGGTTCCAAGGTT", 224),
        },
    ]
}

fn random_dna(length: usize, seed: u64) -> String {
    let mut rng = SplitMix64::new(seed);
    let mut sequence = String::with_capacity(length);
    for _ in 0..length {
        let base = match rng.next_usize(4) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        };
        sequence.push(base);
    }
    sequence
}

fn repeat_pattern(pattern: &str, length: usize) -> String {
    pattern.chars().cycle().take(length).collect()
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

fn invalid_input(message: &str) -> io::Error {
    io::Error::new(ErrorKind::InvalidInput, message.to_string())
}

fn sha256_hex(path: &Path) -> io::Result<String> {
    let data = fs::read(path)?;
    let digest = Sha256::digest(data);
    Ok(format!("{digest:x}"))
}

fn normalized_fragment_id(qname: &str) -> &str {
    qname
        .strip_suffix("/1")
        .or_else(|| qname.strip_suffix("/2"))
        .unwrap_or(qname)
}

#[derive(Debug, Clone)]
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn next_usize(&mut self, upper_bound: usize) -> usize {
        if upper_bound <= 1 {
            0
        } else {
            (self.next_u64() % upper_bound as u64) as usize
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use seq_io::fastq::Reader;

    #[test]
    fn test_fastq_pair_reproducible() {
        let output_a = unique_temp_dir("fastq-repro-a");
        let output_b = unique_temp_dir("fastq-repro-b");
        let config = FastqPairConfig::new(48, 75, 42)
            .with_phred_score(37)
            .with_components(vec![
                ReadComponent::new(SyntheticSource::Human, 0.75),
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    0.25,
                ),
            ]);

        let (r1_first, r2_first) = generate_fastq_pair(&config.clone().with_output_dir(&output_a))
            .expect("first FASTQ generation should succeed");
        let (r1_second, r2_second) = generate_fastq_pair(&config.with_output_dir(&output_b))
            .expect("second FASTQ generation should succeed");

        assert_eq!(
            sha256_hex(&r1_first).expect("R1 digest should compute"),
            sha256_hex(&r1_second).expect("R1 digest should compute")
        );
        assert_eq!(
            sha256_hex(&r2_first).expect("R2 digest should compute"),
            sha256_hex(&r2_second).expect("R2 digest should compute")
        );
    }

    #[test]
    fn test_fastq_pair_seq_io_parseable() {
        let output_dir = unique_temp_dir("fastq-parseable");
        let config = FastqPairConfig::new(24, 80, 20260420)
            .with_output_dir(&output_dir)
            .with_components(vec![
                ReadComponent::new(SyntheticSource::Human, 0.50),
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::SourceType("spike-in-virus".to_string())),
                    0.25,
                ),
                ReadComponent::new(SyntheticSource::LowComplexity, 0.125),
                ReadComponent::new(SyntheticSource::Adapter, 0.125),
            ]);

        let (r1_path, r2_path) =
            generate_fastq_pair(&config).expect("FASTQ generation should succeed");

        let mut r1_reader = Reader::from_path(&r1_path).expect("R1 FASTQ should open");
        let mut r2_reader = Reader::from_path(&r2_path).expect("R2 FASTQ should open");
        let mut seen_pairs = 0usize;

        loop {
            match (r1_reader.next(), r2_reader.next()) {
                (Some(Ok(r1)), Some(Ok(r2))) => {
                    let qname_r1 = r1.id().expect("R1 id should be UTF-8");
                    let qname_r2 = r2.id().expect("R2 id should be UTF-8");

                    assert_eq!(
                        normalized_fragment_id(qname_r1),
                        normalized_fragment_id(qname_r2)
                    );
                    assert_eq!(qname_r1, qname_r2);
                    assert_eq!(qname_r1.split(':').count(), 7);
                    assert_eq!(r1.seq().len(), config.read_length);
                    assert_eq!(r2.seq().len(), config.read_length);
                    assert_eq!(r1.qual().len(), config.read_length);
                    assert_eq!(r2.qual().len(), config.read_length);
                    seen_pairs += 1;
                }
                (None, None) => break,
                (Some(Err(err)), _) | (_, Some(Err(err))) => {
                    panic!("seq_io should parse synthetic FASTQ cleanly: {err}")
                }
                _ => panic!("R1 and R2 must contain the same number of records"),
            }
        }

        assert_eq!(seen_pairs, config.read_pairs);
    }

    #[test]
    fn test_mixed_sample_ratio() {
        let output_dir = unique_temp_dir("fastq-ratio");
        let config = FastqPairConfig::new(1_000, 75, 7)
            .with_output_dir(&output_dir)
            .with_components(vec![
                ReadComponent::new(SyntheticSource::Human, 0.99),
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    0.01,
                ),
            ]);

        let (r1_path, _) = generate_fastq_pair(&config).expect("FASTQ generation should succeed");
        let mut reader = Reader::from_path(&r1_path).expect("generated FASTQ should open");
        let mut human_reads = 0usize;
        let mut virus_reads = 0usize;

        while let Some(record) = reader.next() {
            let record = record.expect("generated FASTQ should parse");
            let qname = record.id().expect("QNAME should be UTF-8");
            if qname.starts_with("HUMAN:") {
                human_reads += 1;
            } else if qname.starts_with("VIRUS:") {
                virus_reads += 1;
            }
        }

        assert_eq!(human_reads + virus_reads, 1_000);
        assert_eq!(virus_reads, 10);
        assert_eq!(human_reads, 990);
        assert!(((virus_reads as f64 / 1_000.0) - 0.01).abs() < f64::EPSILON);
    }

    #[test]
    fn test_generate_mini_reference_contains_required_layers() {
        let reference_dir =
            generate_mini_reference().expect("mini reference generation should succeed");
        let fasta = fs::read_to_string(mini_reference_fasta_path(&reference_dir))
            .expect("mini reference FASTA should be readable");
        let manifest: BundleManifest = serde_json::from_str(
            &fs::read_to_string(mini_reference_manifest_path(&reference_dir))
                .expect("mini reference manifest should be readable"),
        )
        .expect("manifest JSON should deserialize");

        assert!(fasta.contains(">mini_human_chr1"));
        assert!(fasta.contains(">mini_virus_alpha"));
        assert!(fasta.contains(">mini_virus_beta"));
        assert!(fasta.contains(">mini_decoy_adapter"));

        let groups: Vec<_> = manifest
            .0
            .iter()
            .map(|entry| entry.group.as_str())
            .collect();
        assert!(groups.contains(&"human"));
        assert!(groups.contains(&"virus"));
        assert!(groups.contains(&"decoy"));

        let human = manifest
            .0
            .iter()
            .find(|entry| entry.group == "human")
            .expect("human contig should exist");
        let decoy = manifest
            .0
            .iter()
            .find(|entry| entry.group == "decoy")
            .expect("decoy contig should exist");
        let virus_count = manifest
            .0
            .iter()
            .filter(|entry| entry.group == "virus")
            .count();

        assert_eq!(human.genome_length, 1024);
        assert_eq!(decoy.genome_length, 224);
        assert_eq!(virus_count, 2);
    }
}
