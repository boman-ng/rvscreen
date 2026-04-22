use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

/// RVScreen — rapid, conservative, auditable viral signal detection in human sequencing data.
#[derive(Parser)]
#[command(name = "rvscreen", version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

/// Top-level subcommands.
#[derive(Subcommand)]
pub enum Commands {
    /// Build and manage reference bundles.
    #[command(name = "ref")]
    Ref {
        #[command(subcommand)]
        action: RefAction,
    },
    /// Generate a calibration profile from benchmark data.
    Calibrate(CalibrateArgs),
    /// Screen a sample for viral signal.
    Screen(ScreenArgs),
    /// Audit pipeline artifacts for integrity and release readiness.
    Audit {
        #[command(subcommand)]
        action: AuditAction,
    },
}

/// Arguments for `rvscreen calibrate`.
#[derive(Args, Debug, Clone)]
pub struct CalibrateArgs {
    /// Reference bundle directory bound to the generated profile.
    #[arg(long, value_name = "DIR")]
    pub reference_bundle: PathBuf,

    /// Benchmark manifest YAML/JSON describing datasets and expected outcomes.
    #[arg(long, value_name = "FILE")]
    pub benchmark_manifest: PathBuf,

    /// Output calibration profile directory.
    #[arg(long, value_name = "DIR")]
    pub out: PathBuf,
}

/// Sampling modes for `rvscreen screen`.
#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScreenMode {
    Representative,
    Streaming,
}

/// Arguments for `rvscreen screen`.
#[derive(Args, Debug, Clone)]
pub struct ScreenArgs {
    /// Input files: two paired FASTQs, one BAM/uBAM, or one CRAM.
    #[arg(long, value_name = "FILE", required = true, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Reference bundle directory.
    #[arg(long, value_name = "DIR")]
    pub reference_bundle: PathBuf,

    /// Calibration profile directory.
    #[arg(long, value_name = "DIR")]
    pub calibration_profile: PathBuf,

    /// Optional negative-control result JSON.
    #[arg(long, value_name = "FILE")]
    pub negative_control: Option<PathBuf>,

    /// Output report bundle directory.
    #[arg(long, value_name = "DIR")]
    pub out: PathBuf,

    /// Sampling mode.
    #[arg(long, value_enum, default_value_t = ScreenMode::Representative)]
    pub mode: ScreenMode,

    /// Requested thread count.
    #[arg(long, value_name = "N", default_value_t = default_threads())]
    pub threads: usize,
}

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
}

/// Subcommands under `ref`.
#[derive(Subcommand)]
pub enum RefAction {
    /// Build a reference bundle from host and virus FASTA inputs.
    Build(RefBuildArgs),
}

/// Arguments for `rvscreen ref build`.
#[derive(Args, Debug, Clone)]
pub struct RefBuildArgs {
    /// Host FASTA layer (e.g. GRCh38 or T2T-CHM13).
    #[arg(long, value_name = "PATH")]
    pub host_fasta: PathBuf,

    /// Virus FASTA panel.
    #[arg(long, value_name = "PATH")]
    pub virus_fasta: PathBuf,

    /// Optional decoy/background FASTA layer.
    #[arg(long, value_name = "PATH")]
    pub decoy_fasta: Option<PathBuf>,

    /// Pre-authored manifest JSON for all bundle contigs.
    #[arg(long, value_name = "PATH")]
    pub manifest: PathBuf,

    /// Taxonomy table to copy into the bundle.
    #[arg(long, value_name = "PATH")]
    pub taxonomy: PathBuf,

    /// Output bundle directory.
    #[arg(long, value_name = "DIR")]
    pub out: PathBuf,
}

/// Subcommands under `audit`.
#[derive(Subcommand)]
pub enum AuditAction {
    /// Verify a report bundle's integrity and release status.
    Verify(AuditVerifyArgs),
}

/// Arguments for `rvscreen audit verify`.
#[derive(Args, Debug, Clone)]
pub struct AuditVerifyArgs {
    /// Report bundle directory to validate.
    #[arg(long, value_name = "DIR")]
    pub report_bundle: PathBuf,

    /// Optional reference bundle directory for version cross-checks.
    #[arg(long, value_name = "DIR")]
    pub reference_bundle: Option<PathBuf>,

    /// Optional calibration profile directory for version and release-status cross-checks.
    #[arg(long, value_name = "DIR")]
    pub calibration_profile: Option<PathBuf>,
}
