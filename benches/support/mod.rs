use crate::testutil::generate_mini_reference;
use rvscreen::reference::{build_reference_bundle, BuildReferenceBundleRequest};
use rvscreen::types::{BundleManifest, ContigEntry};
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

pub struct ReferenceBundleFixture {
    pub bundle_dir: PathBuf,
    pub version: String,
}

pub struct CalibrationProfile<'a> {
    pub profile_id: &'a str,
    pub seed: u64,
    pub rounds: &'a [usize],
    pub max_rounds: usize,
    pub max_background_ratio: f64,
}

pub fn prepare_reference_bundle(base_dir: &Path) -> io::Result<ReferenceBundleFixture> {
    let mini_reference_dir = generate_mini_reference()?;
    let manifest_path = mini_reference_dir.join("manifest.json");
    let fasta_path = mini_reference_dir.join("mini_reference.fa");
    let manifest: BundleManifest =
        serde_json::from_str(&fs::read_to_string(&manifest_path)?).map_err(io_other)?;
    let sequences = read_fasta_records(&fasta_path)?;

    let inputs_dir = base_dir.join("reference-inputs");
    fs::create_dir_all(&inputs_dir)?;
    let host_fasta = inputs_dir.join("host.fa");
    let virus_fasta = inputs_dir.join("virus.fa");
    let decoy_fasta = inputs_dir.join("decoy.fa");
    let manifest_out = inputs_dir.join("manifest.json");
    let taxonomy = inputs_dir.join("taxonomy.tsv");

    write_group_fasta(&host_fasta, &manifest.0, &sequences, "human")?;
    write_group_fasta(&virus_fasta, &manifest.0, &sequences, "virus")?;
    write_group_fasta(&decoy_fasta, &manifest.0, &sequences, "decoy")?;
    fs::copy(&manifest_path, &manifest_out)?;
    write_taxonomy(&taxonomy, &manifest.0)?;

    let bundle_dir = base_dir.join("reference-bundle");
    let outcome = build_reference_bundle(&BuildReferenceBundleRequest {
        host_fasta,
        virus_fasta,
        decoy_fasta: Some(decoy_fasta),
        manifest: manifest_out,
        taxonomy,
        out_dir: bundle_dir.clone(),
    })
    .map_err(io_other)?;

    Ok(ReferenceBundleFixture {
        bundle_dir,
        version: outcome.bundle.version,
    })
}

pub fn write_calibration_profile(
    dir: PathBuf,
    reference_bundle: &str,
    config: &CalibrationProfile<'_>,
) -> io::Result<PathBuf> {
    let rounds = config
        .rounds
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    fs::create_dir_all(&dir)?;
    fs::write(
        dir.join("profile.toml"),
        format!(
            "profile_id = \"{}\"\nstatus = \"release_candidate\"\nreference_bundle = \"{reference_bundle}\"\nbackend = \"minimap2\"\npreset = \"sr-conservative\"\nseed = {}\nsupported_input = [\"fastq\", \"fastq.gz\", \"bam\", \"ubam\", \"cram\"]\nsupported_read_type = [\"illumina_pe_shortread\"]\nnegative_control_required = false\n\n[sampling]\nmode = \"representative\"\nrounds = [{rounds}]\nmax_rounds = {}\n\n[fragment_rules]\nmin_mapq = 0\nmin_as_diff = 0\nmax_nm = 100\nrequire_pair_consistency = true\n\n[candidate_rules]\nmin_nonoverlap_fragments = 1\nmin_breadth = 0.0\nmax_background_ratio = {}\n\n[decision_rules]\ntheta_pos = 0.01\ntheta_neg = 0.0001\nallow_indeterminate = true\n",
            config.profile_id, config.seed, config.max_rounds, config.max_background_ratio
        ),
    )?;
    Ok(dir)
}

fn read_fasta_records(path: &Path) -> io::Result<Vec<(String, String)>> {
    let mut records = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_sequence = String::new();

    for line in fs::read_to_string(path)?.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = current_header.replace(header.to_string()) {
                records.push((previous, std::mem::take(&mut current_sequence)));
            }
        } else if !line.trim().is_empty() {
            current_sequence.push_str(line.trim());
        }
    }

    if let Some(header) = current_header {
        records.push((header, current_sequence));
    }

    Ok(records)
}

fn write_group_fasta(
    path: &Path,
    entries: &[ContigEntry],
    sequences: &[(String, String)],
    group: &str,
) -> io::Result<()> {
    let mut body = String::new();

    for entry in entries.iter().filter(|entry| entry.group == group) {
        let sequence = sequences
            .iter()
            .find_map(|(header, sequence)| (header == &entry.contig).then_some(sequence))
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("missing FASTA sequence for {}", entry.contig),
                )
            })?;
        body.push('>');
        body.push_str(&entry.contig);
        body.push('\n');
        body.push_str(sequence);
        body.push('\n');
    }

    fs::write(path, body)
}

fn write_taxonomy(path: &Path, entries: &[ContigEntry]) -> io::Result<()> {
    let mut rows = String::from("taxid\tname\n");
    let mut seen = std::collections::BTreeSet::new();

    for entry in entries {
        if seen.insert((entry.taxid, entry.virus_name.clone())) {
            rows.push_str(&format!("{}\t{}\n", entry.taxid, entry.virus_name));
        }
    }

    fs::write(path, rows)
}

fn io_other(error: impl ToString) -> io::Error {
    io::Error::other(error.to_string())
}
