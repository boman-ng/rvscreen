use criterion::{black_box, BenchmarkId, Criterion, Throughput};
use rayon::prelude::*;
use rvscreen::align::{CompetitiveAligner, CompetitivePreset, FragmentAlignResult};
use rvscreen::io::{FastqPairReader, FragmentRecord};
use rvscreen::pipeline::{
    run_screen, PreparedScreenRunner, ScreenPerfMetrics, REPRESENTATIVE_STARTUP_STAGE_LABELS,
};
use rvscreen::{cli::ScreenArgs, cli::ScreenMode};
use serde_json::json;
use std::env;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};
use tempfile::TempDir;

#[allow(dead_code, unused_imports)]
#[path = "../tests/testutil/mod.rs"]
mod testutil;

mod support;

use support::{prepare_reference_bundle, write_calibration_profile, CalibrationProfile};
use testutil::{
    generate_fastq_pair, FastqPairConfig, ReadComponent, SyntheticSource, VirusSelector,
};

const IO_FRAGMENT_COUNT: usize = 32_000;
const ALIGN_FRAGMENT_COUNT: usize = 4_096;
const E2E_FRAGMENT_COUNT: usize = 10_000;
const DEFAULT_SUMMARY_REPEATS: usize = 5;
const DEFAULT_E2E_SUMMARY_REPEATS: usize = 3;
const CI_SUMMARY_REPEATS: usize = 3;
const CI_E2E_SUMMARY_REPEATS: usize = 2;
const DEFAULT_THREAD_CURVE: [usize; 4] = [1, 2, 4, 8];
const CI_THREAD_CURVE: [usize; 2] = [1, 16];
const DEFAULT_E2E_THREADS: [usize; 2] = [1, 16];
const SUMMARY_PATH: &str = "target/criterion/representative_perf_summary.tsv";
const BASELINE_EVIDENCE_PATH: &str = ".sisyphus/evidence/task-1-startup-baseline.json";
const TASK7_EVIDENCE_PATH: &str = ".sisyphus/evidence/task-7-startup-before-after.json";
const MAX_THREAD_BUDGET: usize = 16;

fn main() {
    real_main().expect("representative perf benchmark suite should succeed");
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
    write_baseline_evidence(&e2e_fixture, &settings)?;
    write_task7_startup_evidence(&e2e_fixture, &settings)?;

    criterion.final_summary();
    Ok(())
}

fn benchmark_fastq_input(criterion: &mut Criterion, fixture: &BenchmarkFixture) {
    let mut group = criterion.benchmark_group("representative_perf_input_throughput");
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
    let mut group = criterion.benchmark_group("representative_perf_alignment_throughput");
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

fn benchmark_end_to_end(
    criterion: &mut Criterion,
    fixture: &BenchmarkFixture,
    settings: &BenchSettings,
) {
    let mut group = criterion.benchmark_group("representative_perf_end_to_end_latency");
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

    let mut cold_group = criterion.benchmark_group("representative_perf_prepare_cold_start");
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
            "representative_retained_sample_align",
            threads,
            settings.summary_repeats,
            || align_all_fragments(alignment_fragments, &aligner, thread_pool.as_ref()),
        )?);

        rows.push(measure_latency_ms(
            "startup_latency",
            "representative_aligner_prepare_only",
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
            "representative_scan_align_report_bundle",
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
            "representative_prepare_only",
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
            "representative_execute_only_no_scan_report",
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
                let outcome = prepared
                    .run_without_report(&execute_args)
                    .map_err(io_other)?;
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

    let perf_rows = collect_representative_stage_summary(e2e_fixture, settings)?;
    rows.extend(perf_rows);

    Ok(rows)
}

fn collect_representative_stage_summary(
    e2e_fixture: &BenchmarkFixture,
    settings: &BenchSettings,
) -> io::Result<Vec<SummaryRow>> {
    let mut rows = Vec::new();

    for &threads in &settings.end_to_end_threads {
        rows.push(measure_perf_metric(
            "input_fragments",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.input_fragments as f64,
            "fragments",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-scan"),
        )?);
        rows.push(measure_perf_metric(
            "qc_passing_fragments",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.qc_passing_fragments as f64,
            "fragments",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-align"),
        )?);
        rows.push(measure_perf_metric(
            "representative_sampled_final",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.representative_sampled_final as f64,
            "fragments",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-aggregate"),
        )?);
        rows.push(measure_perf_metric(
            "representative_round_sampled",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.representative_round_sampled as f64,
            "fragments",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-sample"),
        )?);
        rows.push(measure_perf_metric(
            "minimap2_calls",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.minimap2_calls as f64,
            "calls",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-scan-latency"),
        )?);
        rows.push(measure_perf_metric(
            "io_wall_ms",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.io_wall_duration.as_secs_f64() * 1_000.0,
            "ms",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-align-latency"),
        )?);
        rows.push(measure_perf_metric(
            "qc_hash_wall_ms",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.qc_hash_wall_duration.as_secs_f64() * 1_000.0,
            "ms",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-aggregate-latency"),
        )?);
        rows.push(measure_perf_metric(
            "align_wall_ms",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.align_wall_duration.as_secs_f64() * 1_000.0,
            "ms",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-align-wall"),
        )?);
        rows.push(measure_perf_metric(
            "aggregate_wall_ms",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.aggregate_wall_duration.as_secs_f64() * 1_000.0,
            "ms",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-aggregate-wall"),
        )?);
        rows.push(measure_perf_metric(
            "executor_queue_wait_ms",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.executor_queue_wait_duration.as_secs_f64() * 1_000.0,
            "ms",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-queue-wait"),
        )?);
        rows.push(measure_perf_metric(
            "peak_heap_bytes",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.peak_heap_bytes as f64,
            "bytes",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-heap"),
        )?);
        rows.push(measure_perf_metric(
            "aggregate_candidate_count",
            "representative_scan_qc_align_aggregate",
            threads,
            settings.end_to_end_summary_repeats,
            |metrics| metrics.aggregate_candidate_count as f64,
            "candidates",
            || run_perf_screen(e2e_fixture, threads, "summary-stage-candidate-count"),
        )?);

        for &label in &REPRESENTATIVE_STARTUP_STAGE_LABELS {
            rows.push(measure_perf_metric(
                label,
                "representative_startup_stage_timing",
                threads,
                settings.end_to_end_summary_repeats,
                |metrics| stage_timing_ms(metrics, label),
                "ms",
                || run_perf_screen(e2e_fixture, threads, "summary-stage-startup"),
            )?);
        }
    }

    Ok(rows)
}

fn stage_timing_ms(metrics: &ScreenPerfMetrics, label: &str) -> f64 {
    metrics
        .stage_timings
        .iter()
        .find(|timing| timing.label == label)
        .map(|timing| timing.wall_millis())
        .unwrap_or(0.0)
}

fn write_baseline_evidence(fixture: &BenchmarkFixture, settings: &BenchSettings) -> io::Result<()> {
    let threads = settings.baseline_threads;
    let output_dir = fixture.output_dir("task-1-startup-baseline", threads);
    let args = fixture.screen_args(threads, output_dir.clone());
    let representative_command = fixture.representative_command(threads, &output_dir);
    let started = Instant::now();
    let outcome = run_screen(&args).map_err(io_other)?;
    let walltime_seconds = started.elapsed().as_secs_f64();
    let stage_timings = REPRESENTATIVE_STARTUP_STAGE_LABELS
        .iter()
        .map(|label| {
            let wall_ms = stage_timing_ms(&outcome.perf_metrics, label);
            json!({
                "label": label,
                "wall_ms": wall_ms,
            })
        })
        .collect::<Vec<_>>();
    let peak_rss_kb = current_peak_rss_kb().ok();
    let peak_rss_bytes = peak_rss_kb.map(|value| value.saturating_mul(1024));

    let payload = json!({
        "command": &settings.baseline_command,
        "representative_command": representative_command,
        "cache_state": "cold-ish first run or warm rerun; OS page cache was not flushed",
        "input_fixture": {
            "name": "representative_perf_e2e_fixture",
            "fragment_count": fixture.fragment_count,
            "generator_seed": fixture.seed,
            "read_length": 100,
            "components": [
                {"source": "human", "fraction": 0.95},
                {"source": "NC_SYNTHV1.1", "fraction": 0.05}
            ],
            "r1": fixture.r1.display().to_string(),
            "r2": fixture.r2.display().to_string(),
            "reference_bundle": fixture.bundle_dir.display().to_string(),
            "calibration_profile": fixture.calibration_dir.display().to_string(),
        },
        "threads": threads,
        "walltime_seconds": walltime_seconds,
        "cpu_percent": null,
        "cpu_percent_note": "unavailable inside benchmark process; capture command is intended to be wrapped with /usr/bin/time -v and may be enriched from its output",
        "max_rss_bytes": peak_rss_bytes,
        "max_rss_note": if peak_rss_bytes.is_some() { "captured from /proc/self/status VmHWM" } else { "unavailable: /proc/self/status VmHWM/VmRSS not found" },
        "stage_timings": stage_timings,
    });

    if let Some(parent) = settings.baseline_evidence_path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(
        &settings.baseline_evidence_path,
        serde_json::to_vec_pretty(&payload).map_err(io_other)?,
    )?;
    fs::remove_dir_all(output_dir).ok();
    Ok(())
}

fn write_task7_startup_evidence(
    fixture: &BenchmarkFixture,
    settings: &BenchSettings,
) -> io::Result<()> {
    let before = fs::read_to_string(&settings.task7_before_evidence_path)
        .ok()
        .and_then(|body| serde_json::from_str::<serde_json::Value>(&body).ok());
    let after_runs = settings
        .end_to_end_threads
        .iter()
        .copied()
        .map(|threads| task7_after_run(fixture, threads))
        .collect::<io::Result<Vec<_>>>()?;

    let payload = json!({
        "command": &settings.baseline_command,
        "cache_state": "cold-ish first run or warm rerun; OS page cache was not flushed",
        "optimization": {
            "implemented": "FASTQ representative retention uses ordered batch workers for rvscreen-owned cheap QC/hash scoring when --threads > 1",
            "determinism_guard": "batch results are merged by original read order before existing retention and round-prefix derivation",
            "minimap2_index_load": "unchanged serial upstream .mmi load; no unsafe Send/Sync wrappers or multithreaded .mmi claim",
        },
        "thread_metadata": {
            "requested_after_threads": &settings.end_to_end_threads,
            "max_worker_threads": MAX_THREAD_BUDGET,
            "retention_workers": "at most normalized --threads; one serial FASTQ reader/decompressor feeds owned batches",
        },
        "before": before,
        "after": after_runs,
        "external_time": {
            "cpu_percent": null,
            "user_time_seconds": null,
            "system_time_seconds": null,
            "walltime_seconds": null,
            "max_rss_kbytes": null,
            "note": "populate from /usr/bin/time -v wrapper for the benchmark command when available",
        },
    });

    if let Some(parent) = settings.task7_evidence_path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(
        &settings.task7_evidence_path,
        serde_json::to_vec_pretty(&payload).map_err(io_other)?,
    )
}

fn task7_after_run(fixture: &BenchmarkFixture, threads: usize) -> io::Result<serde_json::Value> {
    let output_dir = fixture.output_dir("task-7-startup-after", threads);
    let args = fixture.screen_args(threads, output_dir.clone());
    let representative_command = fixture.representative_command(threads, &output_dir);
    let started = Instant::now();
    let outcome = run_screen(&args).map_err(io_other)?;
    let walltime_seconds = started.elapsed().as_secs_f64();
    let stage_timings = REPRESENTATIVE_STARTUP_STAGE_LABELS
        .iter()
        .map(|label| {
            json!({
                "label": label,
                "wall_ms": stage_timing_ms(&outcome.perf_metrics, label),
            })
        })
        .collect::<Vec<_>>();
    let peak_rss_kb = current_peak_rss_kb().ok();
    fs::remove_dir_all(output_dir).ok();

    Ok(json!({
        "threads": threads,
        "representative_command": representative_command,
        "walltime_seconds": walltime_seconds,
        "cpu_percent": null,
        "user_time_seconds": null,
        "system_time_seconds": null,
        "max_rss_bytes": peak_rss_kb.map(|value| value.saturating_mul(1024)),
        "input_fragments": outcome.perf_metrics.input_fragments,
        "qc_passing_fragments": outcome.perf_metrics.qc_passing_fragments,
        "representative_sampled_final": outcome.perf_metrics.representative_sampled_final,
        "representative_round_sampled": outcome.perf_metrics.representative_round_sampled,
        "round_sampled_fragments": outcome
            .rounds
            .iter()
            .map(|round| round.sampled_fragments)
            .collect::<Vec<_>>(),
        "stage_timings": stage_timings,
    }))
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

fn measure_perf_metric<FRun, FExtract>(
    metric: &'static str,
    scenario: &'static str,
    threads: usize,
    repeats: usize,
    mut extract: FExtract,
    unit: &'static str,
    mut run: FRun,
) -> io::Result<SummaryRow>
where
    FRun: FnMut() -> io::Result<ScreenPerfMetrics>,
    FExtract: FnMut(&ScreenPerfMetrics) -> f64,
{
    let mut samples = Vec::with_capacity(repeats);
    for _ in 0..repeats {
        let metrics = run()?;
        samples.push(extract(&metrics));
    }

    Ok(SummaryRow {
        metric,
        scenario: scenario.to_string(),
        threads: threads.to_string(),
        repeats,
        mean: mean(&samples),
        unit,
        cv_pct: coefficient_of_variation_pct(&samples),
    })
}

fn run_perf_screen(
    fixture: &BenchmarkFixture,
    threads: usize,
    prefix: &str,
) -> io::Result<ScreenPerfMetrics> {
    let output_dir = fixture.output_dir(prefix, threads);
    let outcome =
        run_screen(&fixture.screen_args(threads, output_dir.clone())).map_err(io_other)?;
    fs::remove_dir_all(output_dir).ok();
    Ok(outcome.perf_metrics)
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
            "{}\t{}\t{}\t{}\t{:.6}\t{}\t{:.3}\n",
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
    baseline_evidence_path: PathBuf,
    task7_before_evidence_path: PathBuf,
    task7_evidence_path: PathBuf,
    baseline_threads: usize,
    baseline_command: String,
}

impl BenchSettings {
    fn from_env() -> Self {
        let ci = ci_mode();
        Self {
            alignment_threads: parse_thread_list(
                "RVSCREEN_REPRESENTATIVE_PERF_ALIGNMENT_THREADS",
                Some("RVSCREEN_TASK26_ALIGNMENT_THREADS"),
                if ci {
                    &CI_THREAD_CURVE[..]
                } else {
                    &DEFAULT_THREAD_CURVE[..]
                },
            ),
            end_to_end_threads: parse_thread_list(
                "RVSCREEN_REPRESENTATIVE_PERF_E2E_THREADS",
                Some("RVSCREEN_TASK26_E2E_THREADS"),
                &DEFAULT_E2E_THREADS[..],
            ),
            summary_repeats: env::var("RVSCREEN_REPRESENTATIVE_PERF_SUMMARY_REPEATS")
                .or_else(|_| env::var("RVSCREEN_TASK26_SUMMARY_REPEATS"))
                .ok()
                .map(|value| {
                    parse_positive_usize(&value, "RVSCREEN_REPRESENTATIVE_PERF_SUMMARY_REPEATS")
                })
                .unwrap_or(if ci {
                    CI_SUMMARY_REPEATS
                } else {
                    DEFAULT_SUMMARY_REPEATS
                }),
            end_to_end_summary_repeats: env::var(
                "RVSCREEN_REPRESENTATIVE_PERF_E2E_SUMMARY_REPEATS",
            )
            .or_else(|_| env::var("RVSCREEN_TASK26_E2E_SUMMARY_REPEATS"))
            .ok()
            .map(|value| {
                parse_positive_usize(&value, "RVSCREEN_REPRESENTATIVE_PERF_E2E_SUMMARY_REPEATS")
            })
            .unwrap_or(if ci {
                CI_E2E_SUMMARY_REPEATS
            } else {
                DEFAULT_E2E_SUMMARY_REPEATS
            }),
            sample_size: env::var("RVSCREEN_REPRESENTATIVE_PERF_SAMPLE_SIZE")
                .or_else(|_| env::var("RVSCREEN_TASK26_SAMPLE_SIZE"))
                .ok()
                .map(|value| {
                    parse_positive_usize(&value, "RVSCREEN_REPRESENTATIVE_PERF_SAMPLE_SIZE")
                })
                .unwrap_or(10),
            warm_up_time: Duration::from_secs(1),
            measurement_time: Duration::from_secs(if ci { 1 } else { 2 }),
            end_to_end_measurement_time: Duration::from_secs(if ci { 1 } else { 3 }),
            summary_path: env::var_os("RVSCREEN_REPRESENTATIVE_PERF_SUMMARY_PATH")
                .or_else(|| env::var_os("RVSCREEN_TASK26_SUMMARY_PATH"))
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from(SUMMARY_PATH)),
            baseline_evidence_path: env::var_os("RVSCREEN_REPRESENTATIVE_PERF_BASELINE_EVIDENCE_PATH")
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from(BASELINE_EVIDENCE_PATH)),
            task7_before_evidence_path: env::var_os("RVSCREEN_REPRESENTATIVE_PERF_TASK7_BEFORE_EVIDENCE_PATH")
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from(BASELINE_EVIDENCE_PATH)),
            task7_evidence_path: env::var_os("RVSCREEN_REPRESENTATIVE_PERF_TASK7_EVIDENCE_PATH")
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from(TASK7_EVIDENCE_PATH)),
            baseline_threads: env::var("RVSCREEN_REPRESENTATIVE_PERF_BASELINE_THREADS")
                .ok()
                .map(|value| parse_positive_usize(&value, "RVSCREEN_REPRESENTATIVE_PERF_BASELINE_THREADS"))
                .unwrap_or(16)
                .min(max_thread_budget()),
            baseline_command: env::var("RVSCREEN_REPRESENTATIVE_PERF_BASELINE_COMMAND")
                .unwrap_or_else(|_| {
                    "env RVSCREEN_BENCH_CI=1 RVSCREEN_REPRESENTATIVE_PERF_ALIGNMENT_THREADS=1,16 RVSCREEN_REPRESENTATIVE_PERF_E2E_THREADS=1,16 RVSCREEN_REPRESENTATIVE_PERF_SUMMARY_REPEATS=1 RVSCREEN_REPRESENTATIVE_PERF_E2E_SUMMARY_REPEATS=1 cargo bench --bench representative_perf -- --quiet".to_string()
                }),
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
    seed: u64,
}

impl BenchmarkFixture {
    fn new(fragment_count: usize, seed: u64) -> io::Result<Self> {
        let tempdir = tempfile::tempdir()?;
        let base_dir = tempdir.path().to_path_buf();
        let bundle = prepare_reference_bundle(&base_dir)?;
        let first_round = (fragment_count / 2).max(1);
        let rounds = [first_round, fragment_count];
        let calibration_dir = write_calibration_profile(
            base_dir.join("calibration"),
            &bundle.version,
            &CalibrationProfile::absolute_override(
                "rvscreen_calib_representative_perf",
                20260421,
                &rounds,
                1.0,
            ),
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
            seed,
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
            rounds: None,
            round_proportions: None,
            allow_sampling_threshold_override: false,
            threads,
        }
    }

    fn representative_command(&self, threads: usize, out: &Path) -> String {
        format!(
            "rvscreen screen --input {} {} --reference-bundle {} --calibration-profile {} --out {} --mode representative --threads {}",
            self.r1.display(),
            self.r2.display(),
            self.bundle_dir.display(),
            self.calibration_dir.display(),
            out.display(),
            threads,
        )
    }

    fn output_dir(&self, prefix: &str, threads: usize) -> PathBuf {
        self.base_dir
            .join(format!("{prefix}-{threads}-{}", unique_run_id()))
    }
}

fn io_other(error: impl ToString) -> io::Error {
    io::Error::other(error.to_string())
}

fn parse_thread_list(name: &str, legacy_name: Option<&str>, default: &[usize]) -> Vec<usize> {
    let threads = env::var(name)
        .or_else(|_| {
            legacy_name
                .map(env::var)
                .unwrap_or_else(|| Err(env::VarError::NotPresent))
        })
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
