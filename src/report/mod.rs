pub(crate) mod contract;
mod writer;

use crate::error::Result;
use crate::types::{CandidateCall, RoundRecord, RunManifest, SampleSummary};
use std::path::Path;

pub use writer::ReportWriter;

pub fn write_report_bundle(
    output_dir: impl AsRef<Path>,
    summary: &SampleSummary,
    candidates: &[CandidateCall],
    rounds: &[RoundRecord],
    manifest: &RunManifest,
) -> Result<()> {
    ReportWriter::write(output_dir, summary, candidates, rounds, manifest)
}
