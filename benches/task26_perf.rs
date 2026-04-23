use criterion::{black_box, BenchmarkId, Criterion, Throughput};
use rayon::prelude::*;
use rvscreen::align::{CompetitiveAligner, CompetitivePreset, FragmentAlignResult};
use rvscreen::io::{FastqPairReader, FragmentRecord};
use rvscreen::pipeline::{run_screen, PreparedScreenRunner};
use rvscreen::reference::{build_reference_bundle, BuildReferenceBundleRequest};
use rvscreen::types::{BundleManifest, ContigEntry};
use rvscreen::{cli::ScreenArgs, cli::ScreenMode};
use std::env;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};
use tempfile::TempDir;

#[allow(dead_code, unused_imports)]
#[path = "../tests/testutil/mod.rs"]
mod testutil;

use testutil::{
    generate_fastq_pair, generate_mini_reference, FastqPairConfig, ReadComponent, SyntheticSource,
    VirusSelector,
};

const IO_FRAGMENT_COUNT: usize = 32_000;
const ALIGN_FRAGMENT_COUNT: usize = 4_096;
const E2E_FRAGMENT_COUNT: usize = 10_000;
const DEFAULT_SUMMARY_REPEATS: usize = 5;
const DEFAULT_E2E_SUMMARY_REPEATS: usize = 3;
const CI_SUMMARY_REPEATS: usize = 3;
const CI_E2E_SUMMARY_REPEATS: usize = 2;
const DEFAULT_THREAD_CURVE: [usize; 4] = [1, 2, 4, 8];
const CI_THREAD_CURVE: [usize; 2] = [1, 4];
const DEFAULT_E2E_THREADS: [usize; 2] = [1, 4];
const SUMMARY_PATH: &str = "target/criterion/task26_summary.tsv";
const MAX_THREAD_BUDGET: usize = 16;

fn main() {
    real_main().expect("Task 26 benchmark suite should succeed");
}

fn real_main() -> io::Result<()> {
    let settings = BenchSettings::from_env();
    let mut criterion = Criterion::default()
        .sample_size(settings.sample_size)
        .warm_up_time(settings.warm_up_time)
        .measurement_time(settings.measurement_time)
        .configure_from_args();

    let io_fixture = BenchmarkFixture::new(IO_FRAGMENT_COUNT, 2601)?;
    let alignment_fixture = BenchmarkFixture::new(ALIGN_FRAGMENT_COUNT, 2602)?;
    let e2e_fixture = BenchmarkFixture::new(E2E_FRAGMENT_COUNT, 2603)?;
    let alignment_fragments = alignment_fixture.load_fragments()?;

    benchmark_fastq_input(&mut criterion, &io_fixture);
    benchmark_alignment(
        &mut criterion,
        &alignment_fixture,
        &alignment_fragments,
        &settings,
    );
    benchmark_end_to_end(&mut criterion, &e2e_fixture, &settings);

    let summary_rows = collect_summary(
        &io_fixture,
        &alignment_fixture,
        &alignment_fragments,
        &e2e_fixture,
        &settings,
    )?;
    write_summary(&summary_rows, &settings.summary_path)?;

    criterion.final_summary();
    Ok(())
}

fn benchmark_fastq_input(criterion: &mut Criterion, fixture: &BenchmarkFixture) {
    let mut group = criterion.benchmark_group("task26_input_throughput");
    group.throughput(Throughput::Elements(fixture.fragment_count as u64));
    group.bench_function("fastq_pair_read", |bencher| {
        bencher.iter(|| {
            let pair_count = count_fastq_pairs(&fixture.r1, &fixture.r2)
                .expect("FASTQ throughput benchmark should read all pairs");
            black_box(pair_count)
        });
    });
    group.finish();
}

fn benchmark_alignment(
    criterion: &mut Criterion,
    fixture: &BenchmarkFixture,
    fragments: &[FragmentRecord],
    settings: &BenchSettings,
) {
    let mut group = criterion.benchmark_group("task26_alignment_throughput");
    group.throughput(Throughput::Elements(fragments.len() as u64));

    for &threads in &settings.alignment_threads {
        let aligner = fixture
            .competitive_aligner(threads)
            .expect("alignment benchmark aligner should build");
        let thread_pool = thread_pool_for(threads).expect("alignment thread pool should build");

        group.bench_with_input(
            BenchmarkId::from_parameter(threads),
            &threads,
            |bencher, _| {
                bencher.iter(|| {
                    let checksum = align_all_fragments(fragments, &aligner, thread_pool.as_ref())
                        .expect("alignment throughput benchmark should succeed");
                    black_box(checksum)
                });
            },
        );
    }

    group.finish();
}

fn benchmark_end_to_end(criterion: &mut Criterion, fixture: &BenchmarkFixture, settings: &BenchSettings) {
    let mut group = criterion.benchmark_group("task26_end_to_end_latency");
    group.sample_size(settings.sample_size);
    group.measurement_time(settings.end_to_end_measurement_time);
    group.throughput(Throughput::Elements(fixture.fragment_count as u64));

    for &threads in &settings.end_to_end_threads {
        group.bench_with_input(
            BenchmarkId::from_parameter(threads),
            &threads,
            |bencher, &threads| {
                bencher.iter(|| {
                    let output_dir = fixture.output_dir("criterion-e2e", threads);
                    let outcome = run_screen(&fixture.screen_args(threads, output_dir.clone()))
                        .expect("end-to-end latency benchmark should succeed");
                    assert_eq!(
                        outcome.summary.sampled_fragments,
                        fixture.fragment_count as u64
                    );
                    fs::remove_dir_all(output_dir).ok();
                    black_box(outcome.summary.rounds_run)
                });
            },
        );
    }

    group.finish();

    let mut cold_group = criterion.benchmark_group("task26_prepare_cold_start");
    cold_group.sample_size(settings.sample_size);
    cold_group.measurement_time(settings.end_to_end_measurement_time);

    for &threads in &settings.end_to_end_threads {
        cold_group.bench_with_input(
            BenchmarkId::from_parameter(threads),
            &threads,
            |bencher, &threads| {
                bencher.iter(|| {
                    let output_dir = fixture.output_dir("criterion-prepare", threads);
                    let args = fixture.screen_args(threads, output_dir.clone());
                    let prepared = PreparedScreenRunner::prepare(&args)
                        .expect("prepare cold-start benchmark should succeed");
                    black_box(prepared.aligner_reference_input().to_path_buf());
                    fs::remove_dir_all(output_dir).ok();
                });
            },
        );
    }

    cold_group.finish();
}

fn collect_summary(
    io_fixture: &BenchmarkFixture,
    alignment_fixture: &BenchmarkFixture,
    alignment_fragments: &[FragmentRecord],
    e2e_fixture: &BenchmarkFixture,
    settings: &BenchSettings,
) -> io::Result<Vec<SummaryRow>> {
    let mut rows = Vec::new();
    rows.push(measure_pairs_per_second(
        "input_throughput",
        "fastq_pair_read",
        1,
        settings.summary_repeats,
        || count_fastq_pairs(&io_fixture.r1, &io_fixture.r2).map(|value| value as u64),
    )?);

    for &threads in &settings.alignment_threads {
        let aligner = alignment_fixture.competitive_aligner(threads)?;
        let thread_pool = thread_pool_for(threads)?;
        rows.push(measure_pairs_per_second(
            "alignment_throughput",
            "competitive_align",
            threads,
            settings.summary_repeats,
            || align_all_fragments(alignment_fragments, &aligner, thread_pool.as_ref()),
        )?);

        rows.push(measure_latency_ms(
            "startup_latency",
            "competitive_aligner_prepare",
            threads,
            settings.summary_repeats,
            || {
                let started = Instant::now();
                let aligner = alignment_fixture.competitive_aligner(threads)?;
                black_box(aligner.uses_prebuilt_index());
                Ok(started.elapsed())
            },
        )?);
    }

    for &threads in &settings.end_to_end_threads {
        rows.push(measure_latency_ms(
            "end_to_end_latency",
            "screen_to_report_bundle",
            threads,
            settings.end_to_end_summary_repeats,
            || {
                let output_dir = e2e_fixture.output_dir("summary-e2e", threads);
                let started = Instant::now();
                let outcome = run_screen(&e2e_fixture.screen_args(threads, output_dir.clone()))
                    .map_err(io_other)?;
                assert_eq!(
                    outcome.summary.sampled_fragments,
                    e2e_fixture.fragment_count as u64
                );
                fs::remove_dir_all(output_dir).ok();
                Ok(started.elapsed())
            },
        )?);

        rows.push(measure_latency_ms(
            "startup_latency",
            "screen_prepare_only",
            threads,
            settings.end_to_end_summary_repeats,
            || {
                let output_dir = e2e_fixture.output_dir("summary-prepare", threads);
                let args = e2e_fixture.screen_args(threads, output_dir.clone());
                let started = Instant::now();
                let prepared = PreparedScreenRunner::prepare(&args).map_err(io_other)?;
                black_box(prepared.aligner_uses_prebuilt_index());
                fs::remove_dir_all(output_dir).ok();
                Ok(started.elapsed())
            },
        )?);

        rows.push(measure_latency_ms(
            "end_to_end_latency",
            "prepared_runner_execute_only",
            threads,
            settings.end_to_end_summary_repeats,
            || {
                let prepare_output_dir = e2e_fixture.output_dir("summary-prepared", threads);
                let execute_output_dir = e2e_fixture.output_dir("summary-prepared-run", threads);
                let args = e2e_fixture.screen_args(threads, prepare_output_dir.clone());
                let prepared = PreparedScreenRunner::prepare(&args).map_err(io_other)?;
                let mut execute_args = args.clone();
                execute_args.out = execute_output_dir.clone();
                let started = Instant::now();
                let outcome = prepared.run_without_report(&execute_args).map_err(io_other)?;
                assert_eq!(
                    outcome.summary.sampled_fragments,
                    e2e_fixture.fragment_count as u64
                );
                fs::remove_dir_all(prepare_output_dir).ok();
                fs::remove_dir_all(execute_output_dir).ok();
                Ok(started.elapsed())
            },
        )?);
    }

    let peak_rss_kb = current_peak_rss_kb()?;
    rows.push(SummaryRow {
        metric: "peak_rss",
        scenario: "cargo_bench_process_hwm".to_string(),
        threads: "-".to_string(),
        repeats: 1,
        mean: peak_rss_kb as f64,
        unit: "KiB",
        cv_pct: 0.0,
    });

    Ok(rows)
}

fn measure_pairs_per_second<F>(
    metric: &'static str,
    scenario: &'static str,
    threads: usize,
    repeats: usize,
    mut run: F,
) -> io::Result<SummaryRow>
where
    F: FnMut() -> io::Result<u64>,
{
    let mut samples = Vec::with_capacity(repeats);
    for _ in 0..repeats {
        let started = Instant::now();
        let processed_pairs = run()?;
        let elapsed = started.elapsed().as_secs_f64();
        samples.push(processed_pairs as f64 / elapsed.max(f64::MIN_POSITIVE));
    }

    Ok(SummaryRow {
        metric,
        scenario: scenario.to_string(),
        threads: threads.to_string(),
        repeats,
        mean: mean(&samples),
        unit: "pairs/s",
        cv_pct: coefficient_of_variation_pct(&samples),
    })
}

fn measure_latency_ms<F>(
    metric: &'static str,
    scenario: &'static str,
    threads: usize,
    repeats: usize,
    mut run: F,
) -> io::Result<SummaryRow>
where
    F: FnMut() -> io::Result<Duration>,
{
    let mut samples = Vec::with_capacity(repeats);
    for _ in 0..repeats {
        samples.push(run()?.as_secs_f64() * 1_000.0);
    }

    Ok(SummaryRow {
        metric,
        scenario: scenario.to_string(),
        threads: threads.to_string(),
        repeats,
        mean: mean(&samples),
        unit: "ms",
        cv_pct: coefficient_of_variation_pct(&samples),
    })
}

fn count_fastq_pairs(r1: &Path, r2: &Path) -> io::Result<usize> {
    let mut reader = FastqPairReader::open(r1, r2).map_err(io_other)?;
    let mut count = 0usize;
    while let Some(record) = reader.next() {
        record.map_err(io_other)?;
        count = count.saturating_add(1);
    }
    Ok(count)
}

fn align_all_fragments(
    fragments: &[FragmentRecord],
    aligner: &CompetitiveAligner,
    thread_pool: Option<&rayon::ThreadPool>,
) -> io::Result<u64> {
    let mapped = if let Some(thread_pool) = thread_pool {
        thread_pool.install(|| {
            fragments
                .par_iter()
                .map(|fragment| aligner.align_fragment(fragment).map(alignment_checksum))
                .collect::<Vec<_>>()
        })
    } else {
        fragments
            .iter()
            .map(|fragment| aligner.align_fragment(fragment).map(alignment_checksum))
            .collect::<Vec<_>>()
    };

    mapped
        .into_iter()
        .try_fold(0u64, |accumulator, value| {
            value.map(|value| accumulator.saturating_add(value))
        })
        .map_err(io_other)
}

fn alignment_checksum(result: FragmentAlignResult) -> u64 {
    let primary_mapq = u64::from(
        result
            .best_virus_hit
            .as_ref()
            .or(result.best_host_hit.as_ref())
            .map(|hit| hit.mapq)
            .unwrap_or(0),
    );
    primary_mapq + result.secondary_hits.len() as u64
}

fn thread_pool_for(threads: usize) -> io::Result<Option<rayon::ThreadPool>> {
    if threads <= 1 {
        return Ok(None);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map(Some)
        .map_err(io_other)
}

fn write_summary(rows: &[SummaryRow], path: &Path) -> io::Result<()> {
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

fn current_peak_rss_kb() -> io::Result<u64> {
    let status = fs::read_to_string("/proc/self/status")?;
    status
        .lines()
        .find_map(|line| {
            line.strip_prefix("VmHWM:")
                .or_else(|| line.strip_prefix("VmRSS:"))
        })
        .and_then(|value| value.split_whitespace().next())
        .and_then(|value| value.parse::<u64>().ok())
        .ok_or_else(|| io::Error::other("VmHWM/VmRSS not found in /proc/self/status"))
}

struct SummaryRow {
    metric: &'static str,
    scenario: String,
    threads: String,
    repeats: usize,
    mean: f64,
    unit: &'static str,
    cv_pct: f64,
}

struct BenchSettings {
    alignment_threads: Vec<usize>,
    end_to_end_threads: Vec<usize>,
    summary_repeats: usize,
    end_to_end_summary_repeats: usize,
    sample_size: usize,
    warm_up_time: Duration,
    measurement_time: Duration,
    end_to_end_measurement_time: Duration,
    summary_path: PathBuf,
}

impl BenchSettings {
    fn from_env() -> Self {
        let ci = ci_mode();
        Self {
            alignment_threads: parse_thread_list(
                "RVSCREEN_TASK26_ALIGNMENT_THREADS",
                if ci {
                    &CI_THREAD_CURVE[..]
                } else {
                    &DEFAULT_THREAD_CURVE[..]
                },
            ),
            end_to_end_threads: parse_thread_list(
                "RVSCREEN_TASK26_E2E_THREADS",
                &DEFAULT_E2E_THREADS[..],
            ),
            summary_repeats: env::var("RVSCREEN_TASK26_SUMMARY_REPEATS")
                .ok()
                .map(|value| parse_positive_usize(&value, "RVSCREEN_TASK26_SUMMARY_REPEATS"))
                .unwrap_or(if ci {
                    CI_SUMMARY_REPEATS
                } else {
                    DEFAULT_SUMMARY_REPEATS
                }),
            end_to_end_summary_repeats: env::var("RVSCREEN_TASK26_E2E_SUMMARY_REPEATS")
                .ok()
                .map(|value| {
                    parse_positive_usize(&value, "RVSCREEN_TASK26_E2E_SUMMARY_REPEATS")
                })
                .unwrap_or(if ci {
                    CI_E2E_SUMMARY_REPEATS
                } else {
                    DEFAULT_E2E_SUMMARY_REPEATS
                }),
            sample_size: env::var("RVSCREEN_TASK26_SAMPLE_SIZE")
                .ok()
                .map(|value| parse_positive_usize(&value, "RVSCREEN_TASK26_SAMPLE_SIZE"))
                .unwrap_or(10),
            warm_up_time: Duration::from_secs(1),
            measurement_time: Duration::from_secs(if ci { 1 } else { 2 }),
            end_to_end_measurement_time: Duration::from_secs(if ci { 1 } else { 3 }),
            summary_path: env::var_os("RVSCREEN_TASK26_SUMMARY_PATH")
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from(SUMMARY_PATH)),
        }
    }
}

struct BenchmarkFixture {
    _tempdir: TempDir,
    base_dir: PathBuf,
    bundle_dir: PathBuf,
    calibration_dir: PathBuf,
    r1: PathBuf,
    r2: PathBuf,
    fragment_count: usize,
}

impl BenchmarkFixture {
    fn new(fragment_count: usize, seed: u64) -> io::Result<Self> {
        let tempdir = tempfile::tempdir()?;
        let base_dir = tempdir.path().to_path_buf();
        let bundle = prepare_reference_bundle(&base_dir)?;
        let calibration_dir = write_calibration_profile(
            base_dir.join("calibration"),
            &bundle.version,
            fragment_count,
        )?;
        let (r1, r2) = generate_fastq_pair(
            &FastqPairConfig::new(fragment_count, 100, seed)
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
            _tempdir: tempdir,
            base_dir,
            bundle_dir: bundle.bundle_dir,
            calibration_dir,
            r1,
            r2,
            fragment_count,
        })
    }

    fn load_fragments(&self) -> io::Result<Vec<FragmentRecord>> {
        FastqPairReader::open(&self.r1, &self.r2)
            .map_err(io_other)?
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(io_other)
    }

    fn competitive_aligner(&self, threads: usize) -> io::Result<CompetitiveAligner> {
        CompetitiveAligner::new_with_threads(
            &self.bundle_dir,
            CompetitivePreset::SrConservative,
            threads,
        )
        .map_err(io_other)
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
    fn output_dir(&self, prefix: &str, threads: usize) -> PathBuf {
        self.base_dir
            .join(format!("{prefix}-{threads}-{}", unique_run_id()))
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

fn write_calibration_profile(
    dir: PathBuf,
    reference_bundle: &str,
    fragment_count: usize,
) -> io::Result<PathBuf> {
    fs::create_dir_all(&dir)?;
    fs::write(
        dir.join("profile.toml"),
        format!(
            "profile_id = \"rvscreen_calib_task26_perf\"\nstatus = \"release_candidate\"\nreference_bundle = \"{reference_bundle}\"\nbackend = \"minimap2\"\npreset = \"sr-conservative\"\nseed = 20260421\nsupported_input = [\"fastq\", \"fastq.gz\", \"bam\", \"ubam\", \"cram\"]\nsupported_read_type = [\"illumina_pe_shortread\"]\nnegative_control_required = false\n\n[sampling]\nmode = \"representative\"\nrounds = [{fragment_count}]\nmax_rounds = 1\n\n[fragment_rules]\nmin_mapq = 0\nmin_as_diff = 0\nmax_nm = 100\nrequire_pair_consistency = true\n\n[candidate_rules]\nmin_nonoverlap_fragments = 1\nmin_breadth = 0.0\nmax_background_ratio = 0.0\n\n[decision_rules]\ntheta_pos = 0.01\ntheta_neg = 0.0001\nallow_indeterminate = true\n"
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

fn unique_run_id() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system time should be after epoch")
        .as_nanos()
}
