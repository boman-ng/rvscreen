use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rvscreen::pipeline::PreparedScreenRunner;
use rvscreen::reference::{build_reference_bundle, BuildReferenceBundleRequest};
use rvscreen::types::{BundleManifest, ContigEntry};
use rvscreen::{cli::ScreenArgs, cli::ScreenMode};
use std::env;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};
use tempfile::tempdir;

#[allow(dead_code, unused_imports)]
#[path = "../tests/testutil/mod.rs"]
mod testutil;

use testutil::{
    generate_fastq_pair, generate_mini_reference, FastqPairConfig, ReadComponent, SyntheticSource,
    VirusSelector,
};

const BENCH_FRAGMENT_COUNT: usize = 32_000;
const DEFAULT_THREADS: [usize; 2] = [1, 4];
const DEFAULT_SUMMARY_REPEATS: usize = 5;
const CI_SUMMARY_REPEATS: usize = 3;
const DEFAULT_SUMMARY_PATH: &str = "target/criterion/task23_parallel_screen_summary.tsv";
const MAX_THREAD_BUDGET: usize = 16;

fn configured_criterion() -> Criterion {
    let mut criterion = Criterion::default();
    if ci_mode() {
        criterion = criterion
            .sample_size(10)
            .warm_up_time(Duration::from_secs(1))
            .measurement_time(Duration::from_secs(1));
    }
    criterion.configure_from_args()
}

fn parallel_screen_benchmark(criterion: &mut Criterion) {
    let fixture = BenchmarkFixture::new().expect("benchmark fixture should build");
    let threads = benchmark_threads();
    let mut group = criterion.benchmark_group("task23_parallel_screen");
    group.throughput(Throughput::Elements(BENCH_FRAGMENT_COUNT as u64));

    for threads in &threads {
        let runner = fixture
            .prepared_runner(*threads)
            .expect("prepared runner should build");
        group.bench_with_input(
            BenchmarkId::from_parameter(*threads),
            threads,
            |bencher, &threads| {
                bencher.iter(|| {
                    let output_dir = fixture
                        .base_dir
                        .join(format!("bench-out-{threads}-{}", unique_run_id()));
                    let outcome = runner
                        .run_without_report(&fixture.screen_args(threads, output_dir.clone()))
                        .expect("screen benchmark run should succeed");
                    assert_eq!(
                        outcome.summary.sampled_fragments,
                        BENCH_FRAGMENT_COUNT as u64
                    );
                    fs::remove_dir_all(output_dir).ok();
                });
            },
        );
    }

    group.finish();

    let summary_rows = collect_summary(&fixture, &threads).expect("task23 summary should build");
    write_summary(&summary_rows).expect("task23 summary should be written");
}

criterion_group! {
    name = benches;
    config = configured_criterion();
    targets = parallel_screen_benchmark
}
criterion_main!(benches);

struct SummaryRow {
    metric: &'static str,
    scenario: &'static str,
    threads: usize,
    repeats: usize,
    mean: f64,
    unit: &'static str,
    cv_pct: f64,
}

struct BenchmarkFixture {
    base_dir: PathBuf,
    bundle_dir: PathBuf,
    calibration_dir: PathBuf,
    r1: PathBuf,
    r2: PathBuf,
}

impl BenchmarkFixture {
    fn new() -> io::Result<Self> {
        let base_dir = tempdir()?.keep();
        let bundle = prepare_reference_bundle(&base_dir)?;
        let calibration_dir =
            write_calibration_profile(base_dir.join("calibration"), &bundle.version)?;
        let (r1, r2) = generate_fastq_pair(
            &FastqPairConfig::new(BENCH_FRAGMENT_COUNT, 100, 2401)
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
            "profile_id = \"rvscreen_calib_task23_bench\"\nstatus = \"release_candidate\"\nreference_bundle = \"{reference_bundle}\"\nbackend = \"minimap2\"\npreset = \"sr-conservative\"\nseed = 20260420\nsupported_input = [\"fastq\", \"fastq.gz\", \"bam\", \"ubam\", \"cram\"]\nsupported_read_type = [\"illumina_pe_shortread\"]\nnegative_control_required = false\n\n[sampling]\nmode = \"representative\"\nrounds = [{BENCH_FRAGMENT_COUNT}]\nmax_rounds = 1\n\n[fragment_rules]\nmin_mapq = 0\nmin_as_diff = 0\nmax_nm = 100\nrequire_pair_consistency = true\n\n[candidate_rules]\nmin_nonoverlap_fragments = 1\nmin_breadth = 0.0\nmax_background_ratio = 0.0\n\n[decision_rules]\ntheta_pos = 0.01\ntheta_neg = 0.0001\nallow_indeterminate = true\n"
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

fn collect_summary(fixture: &BenchmarkFixture, threads: &[usize]) -> io::Result<Vec<SummaryRow>> {
    let mut rows = Vec::with_capacity(threads.len());
    let repeats = summary_repeats();

    for &threads in threads {
        let runner = fixture.prepared_runner(threads)?;
        rows.push(measure_fragments_per_second(threads, repeats, || {
            let output_dir = fixture
                .base_dir
                .join(format!("summary-out-{threads}-{}", unique_run_id()));
            let started = Instant::now();
            let outcome = runner
                .run_without_report(&fixture.screen_args(threads, output_dir.clone()))
                .map_err(io_other)?;
            assert_eq!(
                outcome.summary.sampled_fragments,
                BENCH_FRAGMENT_COUNT as u64,
                "task23 summary should process all benchmark fragments"
            );
            fs::remove_dir_all(output_dir).ok();
            Ok(started.elapsed())
        })?);
    }

    Ok(rows)
}

fn measure_fragments_per_second<F>(threads: usize, repeats: usize, mut run: F) -> io::Result<SummaryRow>
where
    F: FnMut() -> io::Result<Duration>,
{
    let mut samples = Vec::with_capacity(repeats);
    for _ in 0..repeats {
        let elapsed = run()?.as_secs_f64();
        samples.push(BENCH_FRAGMENT_COUNT as f64 / elapsed.max(f64::MIN_POSITIVE));
    }

    Ok(SummaryRow {
        metric: "throughput",
        scenario: "screen_run_without_report",
        threads,
        repeats,
        mean: mean(&samples),
        unit: "fragments/s",
        cv_pct: coefficient_of_variation_pct(&samples),
    })
}

fn write_summary(rows: &[SummaryRow]) -> io::Result<()> {
    let path = summary_path();
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }

    let mut body = String::from("metric\tscenario\tthreads\trepeats\tmean\tunit\tcv_pct\n");
    for row in rows {
        body.push_str(&format!(
            "{}\t{}\t{}\t{}\t{:.3}\t{}\t{:.3}\n",
            row.metric, row.scenario, row.threads, row.repeats, row.mean, row.unit, row.cv_pct
        ));
    }
    fs::write(path, body)
}

fn summary_path() -> PathBuf {
    env::var_os("RVSCREEN_TASK23_SUMMARY_PATH")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(DEFAULT_SUMMARY_PATH))
}

fn summary_repeats() -> usize {
    env::var("RVSCREEN_TASK23_SUMMARY_REPEATS")
        .ok()
        .map(|value| parse_positive_usize(&value, "RVSCREEN_TASK23_SUMMARY_REPEATS"))
        .unwrap_or_else(|| {
            if ci_mode() {
                CI_SUMMARY_REPEATS
            } else {
                DEFAULT_SUMMARY_REPEATS
            }
        })
}

fn benchmark_threads() -> Vec<usize> {
    parse_thread_list("RVSCREEN_TASK23_THREADS", &DEFAULT_THREADS[..])
}

fn parse_thread_list(name: &str, default: &[usize]) -> Vec<usize> {
    let threads = env::var(name)
        .ok()
        .map(|value| {
            value
                .split(',')
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| parse_positive_usize(value, name))
                .collect::<Vec<_>>()
        })
        .unwrap_or_else(|| default.to_vec());

    if threads.is_empty() {
        panic!("{name} must contain at least one thread count");
    }

    let max_threads = max_thread_budget();
    for thread in &threads {
        assert!(
            *thread <= max_threads,
            "{name} cannot exceed the {max_threads}-thread repo budget"
        );
    }
    threads
}

fn max_thread_budget() -> usize {
    env::var("RVSCREEN_BENCH_MAX_THREADS")
        .ok()
        .map(|value| parse_positive_usize(&value, "RVSCREEN_BENCH_MAX_THREADS"))
        .map(|value| {
            assert!(
                value <= MAX_THREAD_BUDGET,
                "RVSCREEN_BENCH_MAX_THREADS cannot exceed {MAX_THREAD_BUDGET}"
            );
            value
        })
        .unwrap_or(MAX_THREAD_BUDGET)
}

fn parse_positive_usize(value: &str, name: &str) -> usize {
    let parsed = value
        .parse::<usize>()
        .unwrap_or_else(|error| panic!("{name} must be an integer: {error}"));
    assert!(parsed > 0, "{name} must be greater than zero");
    parsed
}

fn ci_mode() -> bool {
    env::var("RVSCREEN_BENCH_CI")
        .map(|value| {
            let normalized = value.trim();
            !normalized.is_empty() && normalized != "0" && !normalized.eq_ignore_ascii_case("false")
        })
        .unwrap_or(false)
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn coefficient_of_variation_pct(values: &[f64]) -> f64 {
    if values.len() <= 1 {
        return 0.0;
    }

    let mean = mean(values);
    if mean == 0.0 {
        return 0.0;
    }

    let variance = values
        .iter()
        .map(|value| {
            let delta = *value - mean;
            delta * delta
        })
        .sum::<f64>()
        / (values.len() as f64 - 1.0);
    variance.sqrt() / mean * 100.0
}

fn unique_run_id() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system time should be after epoch")
        .as_nanos()
}
