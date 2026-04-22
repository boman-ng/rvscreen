use rvscreen::pipeline::PreparedScreenRunner;
use rvscreen::reference::{build_reference_bundle, BuildReferenceBundleRequest};
use rvscreen::types::{BundleManifest, ContigEntry};
use rvscreen::{cli::ScreenArgs, cli::ScreenMode};
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

#[allow(dead_code, unused_imports)]
#[path = "../../tests/testutil/mod.rs"]
mod testutil;

use testutil::{
    generate_fastq_pair, generate_mini_reference, FastqPairConfig, ReadComponent, SyntheticSource,
    VirusSelector,
};

const RSS_FRAGMENT_COUNT: usize = 32_000;

fn main() -> io::Result<()> {
    let threads = std::env::args()
        .nth(1)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing thread count arg"))?
        .parse::<usize>()
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error.to_string()))?;

    let fixture = ProbeFixture::new()?;
    let runner = fixture.prepared_runner(threads)?;
    let output_dir = fixture.base_dir.join(format!("rss-out-{threads}"));
    let outcome = runner
        .run_without_report(&fixture.screen_args(threads, output_dir.clone()))
        .map_err(io_other)?;
    let rss_kb = current_rss_kb()?;
    eprintln!(
        "threads={threads}\tsampled_fragments={}\trss_kb={rss_kb}",
        outcome.summary.sampled_fragments
    );
    fs::remove_dir_all(output_dir).ok();
    Ok(())
}

struct ProbeFixture {
    base_dir: PathBuf,
    bundle_dir: PathBuf,
    calibration_dir: PathBuf,
    r1: PathBuf,
    r2: PathBuf,
}

impl ProbeFixture {
    fn new() -> io::Result<Self> {
        let base_dir = unique_temp_dir("task23-rss-probe")?;
        let bundle = prepare_reference_bundle(&base_dir)?;
        let calibration_dir =
            write_calibration_profile(base_dir.join("calibration"), &bundle.version)?;
        let (r1, r2) = generate_fastq_pair(
            &FastqPairConfig::new(RSS_FRAGMENT_COUNT, 100, 2501)
                .with_output_dir(base_dir.join("fastq"))
                .with_components(vec![
                    ReadComponent::new(SyntheticSource::Human, 0.95),
                    ReadComponent::new(
                        SyntheticSource::Virus(VirusSelector::Accession(
                            "NC_SYNTHV1.1".to_string(),
                        )),
                        0.05,
                    ),
                ]),
        )?;

        Ok(Self {
            base_dir,
            bundle_dir: bundle.bundle_dir,
            calibration_dir,
            r1,
            r2,
        })
    }

    fn screen_args(&self, threads: usize, out: PathBuf) -> ScreenArgs {
        ScreenArgs {
            input: vec![self.r1.clone(), self.r2.clone()],
            reference_bundle: self.bundle_dir.clone(),
            calibration_profile: self.calibration_dir.clone(),
            negative_control: None,
            out,
            mode: ScreenMode::Representative,
            threads,
        }
    }

    fn prepared_runner(&self, threads: usize) -> io::Result<PreparedScreenRunner> {
        PreparedScreenRunner::prepare(&self.screen_args(
            threads,
            self.base_dir.join(format!("prepare-out-{threads}")),
        ))
        .map_err(io_other)
    }
}

struct ReferenceBundleFixture {
    bundle_dir: PathBuf,
    version: String,
}

fn prepare_reference_bundle(base_dir: &Path) -> io::Result<ReferenceBundleFixture> {
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

fn write_calibration_profile(dir: PathBuf, reference_bundle: &str) -> io::Result<PathBuf> {
    fs::create_dir_all(&dir)?;
    fs::write(
        dir.join("profile.toml"),
        format!(
            "profile_id = \"rvscreen_calib_task23_rss\"\nstatus = \"release_candidate\"\nreference_bundle = \"{reference_bundle}\"\nbackend = \"minimap2\"\npreset = \"sr-conservative\"\nseed = 20260420\nsupported_input = [\"fastq\", \"fastq.gz\", \"bam\", \"ubam\", \"cram\"]\nsupported_read_type = [\"illumina_pe_shortread\"]\nnegative_control_required = false\n\n[sampling]\nmode = \"representative\"\nrounds = [{RSS_FRAGMENT_COUNT}]\nmax_rounds = 1\n\n[fragment_rules]\nmin_mapq = 0\nmin_as_diff = 0\nmax_nm = 100\nrequire_pair_consistency = true\n\n[candidate_rules]\nmin_nonoverlap_fragments = 1\nmin_breadth = 0.0\nmax_background_ratio = 0.0\n\n[decision_rules]\ntheta_pos = 0.01\ntheta_neg = 0.0001\nallow_indeterminate = true\n"
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

fn current_rss_kb() -> io::Result<u64> {
    let status = fs::read_to_string("/proc/self/status")?;
    status
        .lines()
        .find_map(|line| {
            line.strip_prefix("VmHWM:")
                .or_else(|| line.strip_prefix("VmRSS:"))
        })
        .and_then(|value| value.split_whitespace().next())
        .and_then(|value| value.parse::<u64>().ok())
        .ok_or_else(|| io::Error::other("VmRSS/VmHWM not found in /proc/self/status"))
}

fn io_other(error: impl ToString) -> io::Error {
    io::Error::other(error.to_string())
}

fn unique_temp_dir(prefix: &str) -> io::Result<PathBuf> {
    let dir = std::env::temp_dir().join(format!("{prefix}-{}", unique_run_id()));
    fs::create_dir_all(&dir)?;
    Ok(dir)
}

fn unique_run_id() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system time should be after epoch")
        .as_nanos()
}
