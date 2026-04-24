use anyhow::Result;
use clap::Parser;
use rvscreen::audit::run_audit_verify;
use rvscreen::calibration::run_calibration;
use rvscreen::cli::{AuditAction, Cli, Commands, RefAction};
use rvscreen::pipeline::run_screen;
use rvscreen::reference::{build_reference_bundle, BuildReferenceBundleRequest};
use std::io::{self, Write};
use std::process::ExitCode;
use tracing::info;
use tracing_subscriber::{fmt, EnvFilter};

fn init_tracing() {
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("warn"));
    let _ = fmt().with_env_filter(filter).try_init();
}

fn main() -> Result<ExitCode> {
    init_tracing();

    let cli = Cli::parse();
    match cli.command {
        Commands::Ref { action } => match action {
            RefAction::Build(args) => {
                let outcome = build_reference_bundle(&BuildReferenceBundleRequest {
                    host_fasta: args.host_fasta,
                    virus_fasta: args.virus_fasta,
                    decoy_fasta: args.decoy_fasta,
                    manifest: args.manifest,
                    taxonomy: args.taxonomy,
                    out_dir: args.out,
                })?;
                info!(
                    bundle_version = %outcome.bundle.version,
                    bundle_dir = %outcome.bundle_dir.display(),
                    "built reference bundle"
                );
                Ok(ExitCode::SUCCESS)
            }
        },
        Commands::Calibrate(args) => {
            let outcome = run_calibration(&args)?;
            info!(
                profile_id = %outcome.profile.profile_id,
                profile_status = %outcome.profile.status,
                output_dir = %outcome.output_dir.display(),
                benchmark_runs = outcome.benchmark_runs.len(),
                "calibration profile generated"
            );
            Ok(ExitCode::SUCCESS)
        }
        Commands::Run(args) | Commands::Screen(args) => {
            let outcome = run_screen(&args)?;
            for warning in &outcome.sampling_warnings {
                let round = warning
                    .round
                    .map(|round| round.to_string())
                    .unwrap_or_else(|| "all".to_string());
                writeln!(
                    io::stderr(),
                    "rvscreen_warning kind={} round={} message={:?}",
                    warning.kind.as_str(),
                    round,
                    warning.message
                )?;
            }
            info!(
                sample_id = %outcome.summary.sample_id,
                decision_status = ?outcome.summary.decision_status,
                release_status = ?outcome.summary.release_status,
                rounds_run = outcome.summary.rounds_run,
                "screening pipeline completed"
            );
            Ok(ExitCode::SUCCESS)
        }
        Commands::Audit { action } => match action {
            AuditAction::Verify(args) => {
                let report = run_audit_verify(&args)?;
                if report.passed() {
                    writeln!(io::stdout(), "{}", report.render())?;
                    Ok(ExitCode::SUCCESS)
                } else {
                    writeln!(io::stderr(), "{}", report.render())?;
                    Ok(ExitCode::from(1))
                }
            }
        },
    }
}
