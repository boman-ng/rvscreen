#!/usr/bin/env bash
set -euo pipefail

repo_root="${RVSCREEN_CAPTURE_REPO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}"
output_dir="${1:?usage: capture_ci_perf.sh <output-dir> [rss-threads]}"
rss_threads="${2:-4}"

thread_budget="${RVSCREEN_THREAD_BUDGET:-16}"
build_jobs="${CARGO_BUILD_JOBS:-$thread_budget}"
rayon_threads="${RAYON_NUM_THREADS:-$thread_budget}"

if (( thread_budget < 1 || thread_budget > 16 )); then
    printf 'RVSCREEN_THREAD_BUDGET must be between 1 and 16\n' >&2
    exit 1
fi
if (( build_jobs < 1 || build_jobs > 16 )); then
    printf 'CARGO_BUILD_JOBS must be between 1 and 16\n' >&2
    exit 1
fi
if (( rayon_threads < 1 || rayon_threads > 16 )); then
    printf 'RAYON_NUM_THREADS must be between 1 and 16\n' >&2
    exit 1
fi
if (( rss_threads < 1 || rss_threads > 16 )); then
    printf 'rss_threads must be between 1 and 16\n' >&2
    exit 1
fi

mkdir -p "${output_dir}"

export CARGO_BUILD_JOBS="${build_jobs}"
export RAYON_NUM_THREADS="${rayon_threads}"
export RVSCREEN_BENCH_CI="${RVSCREEN_BENCH_CI:-1}"
export RVSCREEN_BENCH_MAX_THREADS="${RVSCREEN_BENCH_MAX_THREADS:-16}"
export RVSCREEN_REPRESENTATIVE_SCREEN_SUMMARY_PATH="${output_dir}/task23_summary.tsv"
export RVSCREEN_TASK23_RSS_PATH="${output_dir}/task23_rss.tsv"
export RVSCREEN_REPRESENTATIVE_PERF_SUMMARY_PATH="${output_dir}/task26_summary.tsv"
export CARGO_TERM_COLOR="never"

cd "${repo_root}"

cargo bench --bench representative_screen 2>&1 | tee "${output_dir}/representative_screen_bench.log"
cargo run --profile bench --bin task23_rss_probe -- "${rss_threads}" 2>&1 | tee "${output_dir}/task23_rss.log"
cargo bench --bench representative_perf 2>&1 | tee "${output_dir}/representative_perf_bench.log"
