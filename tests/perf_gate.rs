use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::tempdir;

const TASK23_THROUGHPUT_DROP_PCT: f64 = 25.0;
const TASK23_RSS_INCREASE_PCT: f64 = 35.0;
const TASK26_REPRESENTATIVE_ALIGNMENT_DROP_PCT: f64 = 30.0;
const TASK26_E2E_LATENCY_INCREASE_PCT: f64 = 35.0;
const TASK26_RSS_INCREASE_PCT: f64 = 35.0;
const TASK23_MIN_THROUGHPUT_1: f64 = 20_000.0;
const TASK23_MIN_THROUGHPUT_4: f64 = 80_000.0;
const TASK23_MAX_RSS_KIB: f64 = 120_000.0;
const TASK26_MIN_INPUT_THROUGHPUT: f64 = 500_000.0;
const TASK26_MIN_REPRESENTATIVE_ALIGNMENT_THROUGHPUT_4: f64 = 8_000_000.0;
const TASK26_MAX_STARTUP_LATENCY_4_MS: f64 = 20.0;
const TASK26_MAX_E2E_LATENCY_4_MS: f64 = 200.0;
const TASK26_MAX_RSS_KIB: f64 = 120_000.0;

#[test]
fn perf_gate_accepts_current_artifacts_and_small_metric_deltas() {
    let tempdir = tempdir().expect("tempdir should build");
    let current_dir = tempdir.path().join("current");
    let base_dir = tempdir.path().join("base");

    write_task23_summary(
        &current_dir,
        &[(1, 58_000.0), (4, 210_000.0)][..],
        48_000.0,
        3_200_000.0,
        18_000_000.0,
        55.0,
        0.95,
        50_000.0,
    );
    write_task23_summary(
        &base_dir,
        &[(1, 60_000.0), (4, 220_000.0)][..],
        45_000.0,
        3_100_000.0,
        19_000_000.0,
        50.0,
        0.90,
        48_000.0,
    );

    let current = load_evidence(&current_dir).expect("current evidence should parse");
    let base = load_evidence(&base_dir).expect("base evidence should parse");
    let sanity_errors = validate_current_absolute_guard(&current);
    assert!(
        sanity_errors.is_empty(),
        "expected current artifact sanity checks to pass, got {sanity_errors:?}"
    );
    let errors = compare_against_base(&base, &current);
    assert!(
        errors.is_empty(),
        "expected no gate failures, got {errors:?}"
    );
}

#[test]
fn perf_gate_rejects_clear_regressions() {
    let tempdir = tempdir().expect("tempdir should build");
    let current_dir = tempdir.path().join("current");
    let base_dir = tempdir.path().join("base");

    write_task23_summary(
        &current_dir,
        &[(1, 15_000.0), (4, 70_000.0)][..],
        80_000.0,
        900_000.0,
        7_000_000.0,
        250.0,
        25.0,
        80_000.0,
    );
    write_task23_summary(
        &base_dir,
        &[(1, 60_000.0), (4, 220_000.0)][..],
        45_000.0,
        3_100_000.0,
        19_000_000.0,
        50.0,
        0.90,
        48_000.0,
    );

    let current = load_evidence(&current_dir).expect("current evidence should parse");
    let base = load_evidence(&base_dir).expect("base evidence should parse");
    let sanity_errors = validate_current_absolute_guard(&current);
    let errors = compare_against_base(&base, &current);

    assert!(
        sanity_errors
            .iter()
            .any(|error| error.contains("task23 throughput @4")),
        "expected task23 current-only guard to fail, got {sanity_errors:?}"
    );
    assert!(
        errors
            .iter()
            .any(|error| error.contains("task23 throughput @4")),
        "expected task23 throughput regression, got {errors:?}"
    );
    assert!(
        errors
            .iter()
            .any(|error| error.contains("task26 representative scan+align end-to-end latency @4")),
        "expected task26 latency regression, got {errors:?}"
    );
    assert!(
        errors
            .iter()
            .any(|error| error.contains("task23 peak RSS @4")),
        "expected task23 RSS regression, got {errors:?}"
    );
}

#[test]
fn perf_gate_ci_artifacts_are_valid_when_configured() {
    let Some(current_dir) = env::var_os("RVSCREEN_PERF_GATE_CURRENT_DIR") else {
        eprintln!("skipping perf gate CI artifact test: RVSCREEN_PERF_GATE_CURRENT_DIR is unset");
        return;
    };

    let current_dir = PathBuf::from(current_dir);
    let current = load_evidence(&current_dir)
        .unwrap_or_else(|error| panic!("failed to validate current perf evidence: {error}"));
    let sanity_errors = validate_current_absolute_guard(&current);
    assert!(
        sanity_errors.is_empty(),
        "current perf evidence failed sanity guard:\n{}",
        sanity_errors.join("\n")
    );

    let Some(base_dir) = env::var_os("RVSCREEN_PERF_GATE_BASE_DIR") else {
        eprintln!("no base perf directory configured; validated current evidence shape only");
        return;
    };

    let base_dir = PathBuf::from(base_dir);
    if !base_dir.exists() {
        eprintln!(
            "base perf directory {:?} does not exist; validated current evidence shape only",
            base_dir
        );
        return;
    }

    let base = match load_evidence(&base_dir) {
        Ok(base) => base,
        Err(error) => {
            eprintln!(
                "base perf evidence in {:?} is unavailable or not yet on the new artifact format ({error}); validated current evidence shape only",
                base_dir
            );
            return;
        }
    };
    let errors = compare_against_base(&base, &current);
    assert!(
        errors.is_empty(),
        "obvious performance regression guard failed:\n{}",
        errors.join("\n")
    );
}

#[derive(Clone, Debug)]
struct SummaryRow {
    metric: String,
    scenario: String,
    threads: String,
    mean: f64,
    unit: String,
}

#[derive(Clone, Debug)]
struct RssRow {
    metric: String,
    scenario: String,
    threads: String,
    sampled_fragments: u64,
    value: f64,
    unit: String,
}

#[derive(Clone, Debug)]
struct EvidenceSet {
    task23_summary: BTreeMap<(String, String, String), SummaryRow>,
    task23_rss: BTreeMap<(String, String, String), RssRow>,
    task26_summary: BTreeMap<(String, String, String), SummaryRow>,
}

fn load_evidence(dir: &Path) -> Result<EvidenceSet, String> {
    let task23_summary = load_summary_map(&dir.join("task23_summary.tsv"))?;
    let task23_rss = load_rss_map(&dir.join("task23_rss.tsv"))?;
    let task26_summary = load_summary_map(&dir.join("task26_summary.tsv"))?;

    require_summary(
        &task23_summary,
        "throughput",
        "screen_run_without_report",
        "1",
    )?;
    require_summary(
        &task23_summary,
        "throughput",
        "screen_run_without_report",
        "4",
    )?;
    require_rss(&task23_rss, "peak_rss", "parallel_screen_run", "4")?;

    require_summary(&task26_summary, "input_throughput", "fastq_pair_read", "1")?;
    require_summary(
        &task26_summary,
        "alignment_throughput",
        "representative_retained_sample_align",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "startup_latency",
        "representative_aligner_prepare_only",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "startup_latency",
        "representative_prepare_only",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "end_to_end_latency",
        "representative_scan_align_report_bundle",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "end_to_end_latency",
        "representative_execute_only_no_scan_report",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "input_fragments",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "qc_passing_fragments",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "representative_sampled_final",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "representative_round_sampled",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "minimap2_calls",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "io_wall_ms",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "qc_hash_wall_ms",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "align_wall_ms",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "aggregate_wall_ms",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "executor_queue_wait_ms",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "peak_heap_bytes",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(
        &task26_summary,
        "aggregate_candidate_count",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_summary(&task26_summary, "peak_rss", "cargo_bench_process_hwm", "-")?;
    require_positive_metric(
        &task26_summary,
        "io_wall_ms",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;
    require_positive_metric(
        &task26_summary,
        "peak_heap_bytes",
        "representative_scan_qc_align_aggregate",
        "4",
    )?;

    Ok(EvidenceSet {
        task23_summary,
        task23_rss,
        task26_summary,
    })
}

fn compare_against_base(base: &EvidenceSet, current: &EvidenceSet) -> Vec<String> {
    let mut errors = Vec::new();

    compare_throughput(
        &mut errors,
        "task23 throughput @1",
        &base.task23_summary,
        &current.task23_summary,
        "throughput",
        "screen_run_without_report",
        "1",
        TASK23_THROUGHPUT_DROP_PCT,
    );
    compare_throughput(
        &mut errors,
        "task23 throughput @4",
        &base.task23_summary,
        &current.task23_summary,
        "throughput",
        "screen_run_without_report",
        "4",
        TASK23_THROUGHPUT_DROP_PCT,
    );
    compare_rss(
        &mut errors,
        "task23 peak RSS @4",
        &base.task23_rss,
        &current.task23_rss,
        "peak_rss",
        "parallel_screen_run",
        "4",
        TASK23_RSS_INCREASE_PCT,
    );
    compare_throughput(
        &mut errors,
        "task26 representative retained-sample alignment throughput @4",
        &base.task26_summary,
        &current.task26_summary,
        "alignment_throughput",
        "representative_retained_sample_align",
        "4",
        TASK26_REPRESENTATIVE_ALIGNMENT_DROP_PCT,
    );
    compare_latency(
        &mut errors,
        "task26 representative scan+align end-to-end latency @4",
        &base.task26_summary,
        &current.task26_summary,
        "end_to_end_latency",
        "representative_scan_align_report_bundle",
        "4",
        TASK26_E2E_LATENCY_INCREASE_PCT,
    );
    compare_summary_rss(
        &mut errors,
        "task26 peak RSS",
        &base.task26_summary,
        &current.task26_summary,
        "peak_rss",
        "cargo_bench_process_hwm",
        "-",
        TASK26_RSS_INCREASE_PCT,
    );

    errors
}

fn validate_current_absolute_guard(current: &EvidenceSet) -> Vec<String> {
    let mut errors = Vec::new();

    enforce_minimum_summary(
        &mut errors,
        "task23 throughput @1",
        &current.task23_summary,
        "throughput",
        "screen_run_without_report",
        "1",
        TASK23_MIN_THROUGHPUT_1,
    );
    enforce_minimum_summary(
        &mut errors,
        "task23 throughput @4",
        &current.task23_summary,
        "throughput",
        "screen_run_without_report",
        "4",
        TASK23_MIN_THROUGHPUT_4,
    );
    enforce_maximum_rss(
        &mut errors,
        "task23 peak RSS @4",
        &current.task23_rss,
        "peak_rss",
        "parallel_screen_run",
        "4",
        TASK23_MAX_RSS_KIB,
    );
    enforce_minimum_summary(
        &mut errors,
        "task26 input throughput @1",
        &current.task26_summary,
        "input_throughput",
        "fastq_pair_read",
        "1",
        TASK26_MIN_INPUT_THROUGHPUT,
    );
    enforce_minimum_summary(
        &mut errors,
        "task26 representative retained-sample alignment throughput @4",
        &current.task26_summary,
        "alignment_throughput",
        "representative_retained_sample_align",
        "4",
        TASK26_MIN_REPRESENTATIVE_ALIGNMENT_THROUGHPUT_4,
    );
    enforce_maximum_summary(
        &mut errors,
        "task26 representative prepare-only startup latency @4",
        &current.task26_summary,
        "startup_latency",
        "representative_prepare_only",
        "4",
        TASK26_MAX_STARTUP_LATENCY_4_MS,
    );
    enforce_maximum_summary(
        &mut errors,
        "task26 representative scan+align end-to-end latency @4",
        &current.task26_summary,
        "end_to_end_latency",
        "representative_scan_align_report_bundle",
        "4",
        TASK26_MAX_E2E_LATENCY_4_MS,
    );
    enforce_maximum_summary(
        &mut errors,
        "task26 peak RSS",
        &current.task26_summary,
        "peak_rss",
        "cargo_bench_process_hwm",
        "-",
        TASK26_MAX_RSS_KIB,
    );

    errors
}

fn load_summary_map(path: &Path) -> Result<BTreeMap<(String, String, String), SummaryRow>, String> {
    let body = fs::read_to_string(path)
        .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
    let mut lines = body.lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("{} is empty", path.display()))?;
    if header != "metric\tscenario\tthreads\trepeats\tmean\tunit\tcv_pct" {
        return Err(format!(
            "unexpected summary header in {}: {header}",
            path.display()
        ));
    }

    let mut rows = BTreeMap::new();
    for (index, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<_> = line.split('\t').collect();
        if parts.len() != 7 {
            return Err(format!(
                "{} line {} expected 7 tab-separated columns, got {}",
                path.display(),
                index + 2,
                parts.len()
            ));
        }

        let row = SummaryRow {
            metric: parts[0].to_string(),
            scenario: parts[1].to_string(),
            threads: parts[2].to_string(),
            mean: parts[4].parse::<f64>().map_err(|error| {
                format!(
                    "{} line {} has invalid mean {:?}: {error}",
                    path.display(),
                    index + 2,
                    parts[4]
                )
            })?,
            unit: parts[5].to_string(),
        };
        rows.insert(
            (
                row.metric.clone(),
                row.scenario.clone(),
                row.threads.clone(),
            ),
            row,
        );
    }
    Ok(rows)
}

fn load_rss_map(path: &Path) -> Result<BTreeMap<(String, String, String), RssRow>, String> {
    let body = fs::read_to_string(path)
        .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
    let mut lines = body.lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("{} is empty", path.display()))?;
    if header != "metric\tscenario\tthreads\tsampled_fragments\tvalue\tunit" {
        return Err(format!(
            "unexpected RSS header in {}: {header}",
            path.display()
        ));
    }

    let mut rows = BTreeMap::new();
    for (index, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<_> = line.split('\t').collect();
        if parts.len() != 6 {
            return Err(format!(
                "{} line {} expected 6 tab-separated columns, got {}",
                path.display(),
                index + 2,
                parts.len()
            ));
        }
        let row = RssRow {
            metric: parts[0].to_string(),
            scenario: parts[1].to_string(),
            threads: parts[2].to_string(),
            sampled_fragments: parts[3].parse::<u64>().map_err(|error| {
                format!(
                    "{} line {} has invalid sampled_fragments {:?}: {error}",
                    path.display(),
                    index + 2,
                    parts[3]
                )
            })?,
            value: parts[4].parse::<f64>().map_err(|error| {
                format!(
                    "{} line {} has invalid value {:?}: {error}",
                    path.display(),
                    index + 2,
                    parts[4]
                )
            })?,
            unit: parts[5].to_string(),
        };
        rows.insert(
            (
                row.metric.clone(),
                row.scenario.clone(),
                row.threads.clone(),
            ),
            row,
        );
    }
    Ok(rows)
}

fn require_summary(
    rows: &BTreeMap<(String, String, String), SummaryRow>,
    metric: &str,
    scenario: &str,
    threads: &str,
) -> Result<(), String> {
    let key = (
        metric.to_string(),
        scenario.to_string(),
        threads.to_string(),
    );
    let row = rows
        .get(&key)
        .ok_or_else(|| format!("missing summary row {:?}", key))?;
    if !(row.mean.is_finite() && row.mean > 0.0) {
        return Err(format!(
            "summary row {:?} must have a positive finite mean",
            key
        ));
    }
    Ok(())
}

fn require_positive_metric(
    rows: &BTreeMap<(String, String, String), SummaryRow>,
    metric: &str,
    scenario: &str,
    threads: &str,
) -> Result<(), String> {
    let key = (
        metric.to_string(),
        scenario.to_string(),
        threads.to_string(),
    );
    let row = rows
        .get(&key)
        .ok_or_else(|| format!("missing summary row {:?}", key))?;
    if !(row.mean.is_finite() && row.mean > 0.0) {
        return Err(format!(
            "summary row {:?} must be strictly positive to prove the metric is wired meaningfully",
            key
        ));
    }
    Ok(())
}

fn require_rss(
    rows: &BTreeMap<(String, String, String), RssRow>,
    metric: &str,
    scenario: &str,
    threads: &str,
) -> Result<(), String> {
    let key = (
        metric.to_string(),
        scenario.to_string(),
        threads.to_string(),
    );
    let row = rows
        .get(&key)
        .ok_or_else(|| format!("missing RSS row {:?}", key))?;
    if row.sampled_fragments == 0 {
        return Err(format!("RSS row {:?} must record sampled fragments", key));
    }
    if !(row.value.is_finite() && row.value > 0.0) {
        return Err(format!(
            "RSS row {:?} must have a positive finite value",
            key
        ));
    }
    Ok(())
}

fn compare_throughput(
    errors: &mut Vec<String>,
    label: &str,
    base: &BTreeMap<(String, String, String), SummaryRow>,
    current: &BTreeMap<(String, String, String), SummaryRow>,
    metric: &str,
    scenario: &str,
    threads: &str,
    max_drop_pct: f64,
) {
    let key = &(
        metric.to_string(),
        scenario.to_string(),
        threads.to_string(),
    );
    let base_row = &base[key];
    let current_row = &current[key];
    let min_allowed = base_row.mean * (1.0 - max_drop_pct / 100.0);
    if current_row.mean < min_allowed {
        errors.push(format!(
            "{label} regressed from {:.3} {} to {:.3} {} (allowed drop {:.1}%)",
            base_row.mean, base_row.unit, current_row.mean, current_row.unit, max_drop_pct
        ));
    }
}

fn compare_latency(
    errors: &mut Vec<String>,
    label: &str,
    base: &BTreeMap<(String, String, String), SummaryRow>,
    current: &BTreeMap<(String, String, String), SummaryRow>,
    metric: &str,
    scenario: &str,
    threads: &str,
    max_increase_pct: f64,
) {
    let key = &(
        metric.to_string(),
        scenario.to_string(),
        threads.to_string(),
    );
    let base_row = &base[key];
    let current_row = &current[key];
    let max_allowed = base_row.mean * (1.0 + max_increase_pct / 100.0);
    if current_row.mean > max_allowed {
        errors.push(format!(
            "{label} regressed from {:.3} {} to {:.3} {} (allowed increase {:.1}%)",
            base_row.mean, base_row.unit, current_row.mean, current_row.unit, max_increase_pct
        ));
    }
}

fn compare_summary_rss(
    errors: &mut Vec<String>,
    label: &str,
    base: &BTreeMap<(String, String, String), SummaryRow>,
    current: &BTreeMap<(String, String, String), SummaryRow>,
    metric: &str,
    scenario: &str,
    threads: &str,
    max_increase_pct: f64,
) {
    compare_latency(
        errors,
        label,
        base,
        current,
        metric,
        scenario,
        threads,
        max_increase_pct,
    );
}

fn compare_rss(
    errors: &mut Vec<String>,
    label: &str,
    base: &BTreeMap<(String, String, String), RssRow>,
    current: &BTreeMap<(String, String, String), RssRow>,
    metric: &str,
    scenario: &str,
    threads: &str,
    max_increase_pct: f64,
) {
    let key = &(
        metric.to_string(),
        scenario.to_string(),
        threads.to_string(),
    );
    let base_row = &base[key];
    let current_row = &current[key];
    let max_allowed = base_row.value * (1.0 + max_increase_pct / 100.0);
    if current_row.value > max_allowed {
        errors.push(format!(
            "{label} regressed from {:.3} {} to {:.3} {} (allowed increase {:.1}%)",
            base_row.value, base_row.unit, current_row.value, current_row.unit, max_increase_pct
        ));
    }
}

fn enforce_minimum_summary(
    errors: &mut Vec<String>,
    label: &str,
    rows: &BTreeMap<(String, String, String), SummaryRow>,
    metric: &str,
    scenario: &str,
    threads: &str,
    minimum: f64,
) {
    let key = (
        metric.to_string(),
        scenario.to_string(),
        threads.to_string(),
    );
    let row = &rows[&key];
    if row.mean < minimum {
        errors.push(format!(
            "{label} fell below absolute floor: {:.3} {} < {:.3} {}",
            row.mean, row.unit, minimum, row.unit
        ));
    }
}

fn enforce_maximum_summary(
    errors: &mut Vec<String>,
    label: &str,
    rows: &BTreeMap<(String, String, String), SummaryRow>,
    metric: &str,
    scenario: &str,
    threads: &str,
    maximum: f64,
) {
    let key = (
        metric.to_string(),
        scenario.to_string(),
        threads.to_string(),
    );
    let row = &rows[&key];
    if row.mean > maximum {
        errors.push(format!(
            "{label} exceeded absolute ceiling: {:.3} {} > {:.3} {}",
            row.mean, row.unit, maximum, row.unit
        ));
    }
}

fn enforce_maximum_rss(
    errors: &mut Vec<String>,
    label: &str,
    rows: &BTreeMap<(String, String, String), RssRow>,
    metric: &str,
    scenario: &str,
    threads: &str,
    maximum: f64,
) {
    let key = (
        metric.to_string(),
        scenario.to_string(),
        threads.to_string(),
    );
    let row = &rows[&key];
    if row.value > maximum {
        errors.push(format!(
            "{label} exceeded absolute ceiling: {:.3} {} > {:.3} {}",
            row.value, row.unit, maximum, row.unit
        ));
    }
}

#[allow(clippy::too_many_arguments)]
fn write_task23_summary(
    dir: &Path,
    task23_throughput: &[(usize, f64)],
    task23_rss: f64,
    input_throughput: f64,
    alignment_throughput: f64,
    e2e_latency: f64,
    startup_latency: f64,
    task26_rss: f64,
) {
    fs::create_dir_all(dir).expect("fixture directory should build");
    let execute_only_latency = (e2e_latency - startup_latency).max(1.0);

    let mut task23_body = String::from("metric\tscenario\tthreads\trepeats\tmean\tunit\tcv_pct\n");
    for (threads, mean) in task23_throughput {
        task23_body.push_str(&format!(
            "throughput\tscreen_run_without_report\t{threads}\t3\t{mean:.3}\tfragments/s\t1.000\n"
        ));
    }
    fs::write(dir.join("task23_summary.tsv"), task23_body)
        .expect("task23 summary fixture should write");

    fs::write(
        dir.join("task23_rss.tsv"),
        format!(
            "metric\tscenario\tthreads\tsampled_fragments\tvalue\tunit\npeak_rss\tparallel_screen_run\t4\t32000\t{task23_rss:.3}\tKiB\n"
        ),
    )
    .expect("task23 rss fixture should write");

    fs::write(
        dir.join("task26_summary.tsv"),
        format!(
            "metric\tscenario\tthreads\trepeats\tmean\tunit\tcv_pct\n\
input_throughput\tfastq_pair_read\t1\t3\t{input_throughput:.3}\tpairs/s\t1.000\n\
alignment_throughput\trepresentative_retained_sample_align\t4\t3\t{alignment_throughput:.3}\tpairs/s\t1.000\n\
startup_latency\trepresentative_aligner_prepare_only\t4\t2\t{startup_latency:.3}\tms\t1.000\n\
startup_latency\trepresentative_prepare_only\t4\t2\t{startup_latency:.3}\tms\t1.000\n\
end_to_end_latency\trepresentative_scan_align_report_bundle\t4\t2\t{e2e_latency:.3}\tms\t1.000\n\
end_to_end_latency\trepresentative_execute_only_no_scan_report\t4\t2\t{execute_only_latency:.3}\tms\t1.000\n\
peak_rss\tcargo_bench_process_hwm\t-\t1\t{task26_rss:.3}\tKiB\t0.000\n\
input_fragments\trepresentative_scan_qc_align_aggregate\t4\t2\t10000.000\tfragments\t0.000\n\
qc_passing_fragments\trepresentative_scan_qc_align_aggregate\t4\t2\t10000.000\tfragments\t0.000\n\
representative_sampled_final\trepresentative_scan_qc_align_aggregate\t4\t2\t10000.000\tfragments\t0.000\n\
representative_round_sampled\trepresentative_scan_qc_align_aggregate\t4\t2\t10000.000\tfragments\t0.000\n\
minimap2_calls\trepresentative_scan_qc_align_aggregate\t4\t2\t10000.000\tcalls\t0.000\n\
io_wall_ms\trepresentative_scan_qc_align_aggregate\t4\t2\t10.000\tms\t1.000\n\
qc_hash_wall_ms\trepresentative_scan_qc_align_aggregate\t4\t2\t6.000\tms\t1.000\n\
align_wall_ms\trepresentative_scan_qc_align_aggregate\t4\t2\t20.000\tms\t1.000\n\
aggregate_wall_ms\trepresentative_scan_qc_align_aggregate\t4\t2\t5.000\tms\t1.000\n\
executor_queue_wait_ms\trepresentative_scan_qc_align_aggregate\t4\t2\t1.000\tms\t1.000\n\
peak_heap_bytes\trepresentative_scan_qc_align_aggregate\t4\t2\t1024.000\tbytes\t0.000\n\
aggregate_candidate_count\trepresentative_scan_qc_align_aggregate\t4\t2\t2.000\tcandidates\t0.000\n"
        ),
    )
    .expect("task26 summary fixture should write");
}
