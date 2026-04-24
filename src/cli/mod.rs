use crate::types::validate_round_proportion_values;
use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

const DEFAULT_THREAD_BUDGET: usize = 16;

/// RVScreen: rapid, conservative, auditable viral signal detection in human sequencing data.
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
    /// Recommended one-command representative screening entrypoint.
    #[command(
        about = "Run the recommended representative screening workflow.",
        long_about = "Run the recommended representative screening workflow. By default, `rvscreen run` uses representative mode, which samples across the input before reporting production-oriented results.\n\nRecommended one-command usage:\n  rvscreen run --input reads_R1.fastq.gz reads_R2.fastq.gz --reference-bundle reference-bundle --calibration-profile calibration-profile --out report-bundle\n\n`run` is a thin alias for `screen` with the same flags."
    )]
    Run(ScreenArgs),
    /// Screen a sample for viral signal.
    #[command(
        about = "Screen a sample for viral signal.",
        long_about = "Screen a sample for viral signal. Use `rvscreen run` for the recommended one-command representative workflow; `screen` remains available for explicit mode and sampling control."
    )]
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

    /// Sampling mode. The default is representative, which samples across the input.
    #[arg(long, value_enum, default_value_t = ScreenMode::Representative)]
    pub mode: ScreenMode,

    /// Absolute fragment-count sampling rounds, comma-separated or repeated.
    #[arg(long, value_name = "N[,N...]", value_delimiter = ',', num_args = 1..)]
    pub rounds: Option<Vec<u64>>,

    /// Proportional sampling rounds as fractions or percentages, with default 10k-10M clamping.
    #[arg(
        long,
        value_name = "P[,P...]",
        value_delimiter = ',',
        num_args = 1..,
        value_parser = parse_round_proportion
    )]
    pub round_proportions: Option<Vec<f64>>,

    /// Disable the default 10k-10M proportional clamp; warnings are still emitted.
    #[arg(long = "allow-sampling-threshold-override", default_value_t = false)]
    pub allow_sampling_threshold_override: bool,

    /// Requested thread count.
    #[arg(long, value_name = "N", default_value_t = default_threads())]
    pub threads: usize,
}

impl ScreenArgs {
    pub fn validate_sampling_overrides(&self) -> std::result::Result<(), String> {
        if self.rounds.is_some() && self.round_proportions.is_some() {
            return Err(
                "--rounds and --round-proportions are mutually exclusive sampling overrides"
                    .to_string(),
            );
        }

        if let Some(proportions) = self.round_proportions.as_deref() {
            validate_round_proportion_values(proportions)?;
        }

        Ok(())
    }
}

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
        .min(DEFAULT_THREAD_BUDGET)
}

fn parse_round_proportion(value: &str) -> std::result::Result<f64, String> {
    let value = value.trim();
    if value.is_empty() {
        return Err("round proportion must not be empty".to_string());
    }

    if let Some(percent) = value.strip_suffix('%') {
        let percent = percent.trim();
        if percent.is_empty() {
            return Err("round percentage must include a numeric value before `%`".to_string());
        }
        return percent
            .parse::<f64>()
            .map(|parsed| parsed / 100.0)
            .map_err(|error| format!("failed to parse round percentage `{value}`: {error}"));
    }

    value
        .parse::<f64>()
        .map_err(|error| format!("failed to parse round proportion `{value}`: {error}"))
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
        assert_eq!(args.rounds, None);
        assert_eq!(args.round_proportions, None);
        assert!(!args.allow_sampling_threshold_override);
        assert_eq!(args.threads, default_threads());
    }

    #[test]
    fn parses_run_contract_with_default_representative_mode() {
        let cli = parse([
            "rvscreen",
            "run",
            "--input",
            "reads_R1.fastq.gz",
            "reads_R2.fastq.gz",
            "--reference-bundle",
            "reference-bundle",
            "--calibration-profile",
            "calibration-profile",
            "--out",
            "report-bundle",
        ]);

        let Commands::Run(args) = cli.command else {
            panic!("expected run command");
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
        assert_eq!(args.out, PathBuf::from("report-bundle"));
        assert_eq!(args.mode, ScreenMode::Representative);
        assert!(!args.allow_sampling_threshold_override);
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
    fn task9_help_copy_explains_run_and_sampling_rounds() {
        let top_level_help = help_text(["rvscreen", "--help"]);
        assert!(top_level_help.contains("run"), "{top_level_help}");

        let run_help = help_text(["rvscreen", "run", "--help"]);
        assert!(
            run_help.contains("recommended representative"),
            "{run_help}"
        );
        assert!(
            run_help.contains("Recommended one-command usage"),
            "{run_help}"
        );

        let screen_help = help_text(["rvscreen", "screen", "--help"]);
        assert!(
            screen_help.contains("Absolute fragment-count sampling rounds"),
            "{screen_help}"
        );
        assert!(
            screen_help.contains("fractions or percentages"),
            "{screen_help}"
        );
        assert!(
            screen_help.contains("default 10k-10M clamping"),
            "{screen_help}"
        );
        assert!(
            screen_help.contains("warnings are still emitted"),
            "{screen_help}"
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
        assert_eq!(args.rounds, None);
        assert_eq!(args.round_proportions, None);
        assert!(!args.allow_sampling_threshold_override);
        assert_eq!(args.threads, 4);
    }

    #[test]
    fn parses_screen_absolute_round_overrides() {
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
            "--rounds",
            "10,20,40",
        ]);

        let Commands::Screen(args) = cli.command else {
            panic!("expected screen command");
        };

        assert_eq!(args.rounds, Some(vec![10, 20, 40]));
        assert_eq!(args.round_proportions, None);
        args.validate_sampling_overrides()
            .expect("absolute rounds should validate");
    }

    #[test]
    fn parses_screen_proportional_round_overrides() {
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
            "--round-proportions",
            "0.25,0.5,1.0",
        ]);

        let Commands::Screen(args) = cli.command else {
            panic!("expected screen command");
        };

        assert_eq!(args.rounds, None);
        assert_eq!(args.round_proportions, Some(vec![0.25, 0.5, 1.0]));
        assert!(!args.allow_sampling_threshold_override);
        args.validate_sampling_overrides()
            .expect("round proportions should validate");
    }

    #[test]
    fn parses_sampling_threshold_override_flag() {
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
            "--round-proportions",
            "0.25,0.5",
            "--allow-sampling-threshold-override",
        ]);

        let Commands::Screen(args) = cli.command else {
            panic!("expected screen command");
        };

        assert_eq!(args.round_proportions, Some(vec![0.25, 0.5]));
        assert!(args.allow_sampling_threshold_override);
        args.validate_sampling_overrides()
            .expect("allow flag should not invalidate proportional rounds");
    }

    #[test]
    fn parses_screen_percent_round_proportion_overrides_as_fractions() {
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
            "--round-proportions",
            "0.2%,10%,100%",
        ]);

        let Commands::Screen(args) = cli.command else {
            panic!("expected screen command");
        };
        let proportions = args
            .round_proportions
            .as_deref()
            .expect("percent proportions should parse");

        assert_approx_eq(proportions[0], 0.002);
        assert_approx_eq(proportions[1], 0.10);
        assert_approx_eq(proportions[2], 1.0);
        args.validate_sampling_overrides()
            .expect("percent round proportions should validate after conversion");
    }

    #[test]
    fn task8_run_alias_accepts_same_sampling_overrides_as_screen() {
        let cli = parse([
            "rvscreen",
            "run",
            "--input",
            "reads.cram",
            "--reference-bundle",
            "reference-bundle",
            "--calibration-profile",
            "calibration-profile",
            "--out",
            "report-bundle",
            "--round-proportions",
            "25%,50%,100%",
            "--allow-sampling-threshold-override",
        ]);

        let Commands::Run(args) = cli.command else {
            panic!("expected run command");
        };

        assert_eq!(args.rounds, None);
        assert_eq!(args.round_proportions, Some(vec![0.25, 0.5, 1.0]));
        assert!(args.allow_sampling_threshold_override);
        args.validate_sampling_overrides()
            .expect("run alias should validate proportional overrides like screen");
    }

    #[test]
    fn rejects_mixed_screen_sampling_overrides() {
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
            "--rounds",
            "10,20",
            "--round-proportions",
            "0.5,1.0",
        ]);

        let Commands::Screen(args) = cli.command else {
            panic!("expected screen command");
        };
        let error = args
            .validate_sampling_overrides()
            .expect_err("mixed sampling overrides should be rejected");

        assert!(error.contains("--rounds"), "{error}");
        assert!(error.contains("--round-proportions"), "{error}");
        assert!(error.contains("mutually exclusive"), "{error}");
    }

    #[test]
    fn rejects_invalid_screen_round_proportions() {
        for (value, expected) in [
            ("0.0", "greater than 0 and less than or equal to 1.0"),
            ("1.2", "greater than 0 and less than or equal to 1.0"),
            ("0.5,0.5", "strictly increasing"),
            ("0.5,0.4", "strictly increasing"),
            ("NaN", "finite"),
            ("inf", "finite"),
            ("120%", "greater than 0 and less than or equal to 1.0"),
        ] {
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
                "--round-proportions",
                value,
            ]);

            let Commands::Screen(args) = cli.command else {
                panic!("expected screen command");
            };
            let error = args
                .validate_sampling_overrides()
                .expect_err("invalid round proportions should be rejected");

            assert!(error.contains(expected), "{value}: {error}");
        }
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

        assert!(
            rendered.contains("unrecognized subcommand 'report'"),
            "{rendered}"
        );
        assert!(help.contains("ref"), "{help}");
        assert!(help.contains("calibrate"), "{help}");
        assert!(help.contains("run"), "{help}");
        assert!(help.contains("screen"), "{help}");
        assert!(help.contains("audit"), "{help}");
    }

    fn parse(args: impl IntoIterator<Item = impl Into<OsString> + Clone>) -> Cli {
        Cli::try_parse_from(args).expect("CLI arguments should parse")
    }

    fn help_text<const N: usize>(args: [&str; N]) -> String {
        match Cli::try_parse_from(args) {
            Ok(_) => panic!("help should render as a clap display error"),
            Err(error) => error.to_string(),
        }
    }

    fn assert_approx_eq(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1e-12,
            "expected {actual} to approximately equal {expected}"
        );
    }
}
