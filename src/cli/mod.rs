use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

const DEFAULT_THREAD_BUDGET: usize = 16;

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
        .min(DEFAULT_THREAD_BUDGET)
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

#[cfg(test)]
mod tests {
    use super::*;
    use clap::{CommandFactory, Parser};
    use std::ffi::OsString;
    use std::path::PathBuf;

    #[test]
    fn clap_contract_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn parses_ref_build_contract() {
        let cli = parse([
            "rvscreen",
            "ref",
            "build",
            "--host-fasta",
            "host.fa",
            "--virus-fasta",
            "virus.fa",
            "--decoy-fasta",
            "decoy.fa",
            "--manifest",
            "manifest.json",
            "--taxonomy",
            "taxonomy.tsv",
            "--out",
            "reference-bundle",
        ]);

        let Commands::Ref { action } = cli.command else {
            panic!("expected ref command");
        };
        let RefAction::Build(args) = action;

        assert_eq!(args.host_fasta, PathBuf::from("host.fa"));
        assert_eq!(args.virus_fasta, PathBuf::from("virus.fa"));
        assert_eq!(args.decoy_fasta, Some(PathBuf::from("decoy.fa")));
        assert_eq!(args.manifest, PathBuf::from("manifest.json"));
        assert_eq!(args.taxonomy, PathBuf::from("taxonomy.tsv"));
        assert_eq!(args.out, PathBuf::from("reference-bundle"));
    }

    #[test]
    fn parses_calibrate_contract() {
        let cli = parse([
            "rvscreen",
            "calibrate",
            "--reference-bundle",
            "reference-bundle",
            "--benchmark-manifest",
            "bench.yaml",
            "--out",
            "calibration-profile",
        ]);

        let Commands::Calibrate(args) = cli.command else {
            panic!("expected calibrate command");
        };

        assert_eq!(args.reference_bundle, PathBuf::from("reference-bundle"));
        assert_eq!(args.benchmark_manifest, PathBuf::from("bench.yaml"));
        assert_eq!(args.out, PathBuf::from("calibration-profile"));
    }

    #[test]
    fn parses_screen_contract_with_default_representative_mode() {
        let cli = parse([
            "rvscreen",
            "screen",
            "--input",
            "reads_R1.fastq.gz",
            "reads_R2.fastq.gz",
            "--reference-bundle",
            "reference-bundle",
            "--calibration-profile",
            "calibration-profile",
            "--negative-control",
            "negative-control.json",
            "--out",
            "report-bundle",
        ]);

        let Commands::Screen(args) = cli.command else {
            panic!("expected screen command");
        };

        assert_eq!(
            args.input,
            vec![
                PathBuf::from("reads_R1.fastq.gz"),
                PathBuf::from("reads_R2.fastq.gz"),
            ]
        );
        assert_eq!(args.reference_bundle, PathBuf::from("reference-bundle"));
        assert_eq!(
            args.calibration_profile,
            PathBuf::from("calibration-profile")
        );
        assert_eq!(
            args.negative_control,
            Some(PathBuf::from("negative-control.json"))
        );
        assert_eq!(args.out, PathBuf::from("report-bundle"));
        assert_eq!(args.mode, ScreenMode::Representative);
        assert_eq!(args.threads, default_threads());
    }

    #[test]
    fn screen_default_threads_respects_task9_thread_budget() {
        assert!(default_threads() <= DEFAULT_THREAD_BUDGET);

        let mut command = <ScreenArgs as clap::Args>::augment_args(clap::Command::new("screen"));
        let mut help = Vec::new();
        command
            .write_long_help(&mut help)
            .expect("screen help should render");
        let help = String::from_utf8(help).expect("screen help should be UTF-8");

        assert!(
            help.contains(&format!("[default: {}]", default_threads())),
            "screen help should advertise the clamped default thread count: {help}"
        );
    }

    #[test]
    fn parses_screen_contract_with_streaming_mode_and_threads() {
        let cli = parse([
            "rvscreen",
            "screen",
            "--input",
            "reads.cram",
            "--reference-bundle",
            "reference-bundle",
            "--calibration-profile",
            "calibration-profile",
            "--out",
            "report-bundle",
            "--mode",
            "streaming",
            "--threads",
            "4",
        ]);

        let Commands::Screen(args) = cli.command else {
            panic!("expected screen command");
        };

        assert_eq!(args.input, vec![PathBuf::from("reads.cram")]);
        assert_eq!(args.mode, ScreenMode::Streaming);
        assert_eq!(args.threads, 4);
    }

    #[test]
    fn parses_audit_verify_contract() {
        let cli = parse([
            "rvscreen",
            "audit",
            "verify",
            "--report-bundle",
            "report-bundle",
            "--reference-bundle",
            "reference-bundle",
            "--calibration-profile",
            "calibration-profile",
        ]);

        let Commands::Audit { action } = cli.command else {
            panic!("expected audit command");
        };
        let AuditAction::Verify(args) = action;

        assert_eq!(args.report_bundle, PathBuf::from("report-bundle"));
        assert_eq!(
            args.reference_bundle,
            Some(PathBuf::from("reference-bundle"))
        );
        assert_eq!(
            args.calibration_profile,
            Some(PathBuf::from("calibration-profile"))
        );
    }

    #[test]
    fn rejects_unknown_top_level_command_family() {
        let error = match Cli::try_parse_from(["rvscreen", "report"]) {
            Ok(_) => panic!("unsupported command family should be rejected"),
            Err(error) => error,
        };
        let rendered = error.to_string();
        let mut help = Vec::new();
        Cli::command()
            .write_long_help(&mut help)
            .expect("CLI help should render");
        let help = String::from_utf8(help).expect("CLI help should be UTF-8");

        assert!(rendered.contains("unrecognized subcommand 'report'"), "{rendered}");
        assert!(help.contains("ref"), "{help}");
        assert!(help.contains("calibrate"), "{help}");
        assert!(help.contains("screen"), "{help}");
        assert!(help.contains("audit"), "{help}");
    }

    fn parse(args: impl IntoIterator<Item = impl Into<OsString> + Clone>) -> Cli {
        Cli::try_parse_from(args).expect("CLI arguments should parse")
    }
}
