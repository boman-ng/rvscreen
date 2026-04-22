#[path = "testutil/mod.rs"]
mod testutil;

use csv::ReaderBuilder;
use rvscreen::reference::{build_reference_bundle, BuildReferenceBundleRequest};
use rvscreen::report::ReportWriter;
use rvscreen::types::{
    BundleManifest, BundleToml, CandidateCall, ContigEntry, DecisionStatus, EvidenceStrength,
    ReleaseStatus, RoundRecord, RunManifest, SampleSummary, SamplingConfig, StopReason,
};
use rvscreen::{
    calibration::{
        BackendGate, GateStatus, ReferenceGate, ReleaseGate as CalibrationReleaseGate,
        SensitivityGate, SpecificityGate,
    },
    decision::SAMPLING_ONLY_CI_LABEL,
};
use serde::Deserialize;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::tempdir;
use testutil::{
    generate_fastq_pair, generate_mini_reference, FastqPairConfig, ReadComponent, SyntheticSource,
    VirusSelector,
};

#[test]
fn scaffold_compiles() {}

#[test]
fn integration_pure_negative_sample_is_final_negative() {
    let tempdir = tempdir().expect("tempdir should be created");
    let bundle =
        prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
    let calibration_dir = write_calibration_profile(
        tempdir.path().join("calibration-negative"),
        &bundle.version,
        CalibrationProfileOptions {
            write_passing_release_gate: true,
            ..CalibrationProfileOptions::default()
        },
    )
    .expect("calibration profile should be written");
    let negative_control =
        write_negative_control(tempdir.path().join("negative-control.json"), "pass", &[])
            .expect("negative control should be written");
    let (r1, r2) = generate_fastq_pair(
        &FastqPairConfig::new(200, 100, 11)
            .with_output_dir(tempdir.path().join("fastq-negative"))
            .with_components(vec![ReadComponent::new(SyntheticSource::Human, 1.0)]),
    )
    .expect("negative FASTQ should be generated");
    let out_dir = tempdir.path().join("report-negative");

    run_screen_cli([
        os("screen"),
        os("--input"),
        r1.into_os_string(),
        r2.into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.clone().into_os_string(),
        os("--calibration-profile"),
        calibration_dir.into_os_string(),
        os("--negative-control"),
        negative_control.into_os_string(),
        os("--out"),
        out_dir.clone().into_os_string(),
        os("--mode"),
        os("representative"),
    ]);

    let summary = read_summary(&out_dir);
    assert_eq!(
        summary.decision_status,
        rvscreen::types::DecisionStatus::Negative
    );
    assert_eq!(summary.release_status, ReleaseStatus::Final);
    assert_report_bundle_artifacts(&out_dir);
    assert_candidate_accessions(&out_dir, &[]);
    assert!(coverage_file_map(&out_dir).is_empty());
}

#[test]
fn integration_high_titer_positive_sample_is_final_positive() {
    let tempdir = tempdir().expect("tempdir should be created");
    let bundle =
        prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
    let calibration_dir = write_calibration_profile(
        tempdir.path().join("calibration-positive"),
        &bundle.version,
        CalibrationProfileOptions {
            write_passing_release_gate: true,
            ..CalibrationProfileOptions::default()
        },
    )
    .expect("calibration profile should be written");
    let negative_control =
        write_negative_control(tempdir.path().join("negative-control.json"), "pass", &[])
            .expect("negative control should be written");
    let (r1, r2) = generate_fastq_pair(
        &FastqPairConfig::new(200, 100, 22)
            .with_output_dir(tempdir.path().join("fastq-positive"))
            .with_components(vec![
                ReadComponent::new(SyntheticSource::Human, 0.95),
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    0.05,
                ),
            ]),
    )
    .expect("positive FASTQ should be generated");
    let out_dir = tempdir.path().join("report-positive");

    run_screen_cli([
        os("screen"),
        os("--input"),
        r1.into_os_string(),
        r2.into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.clone().into_os_string(),
        os("--calibration-profile"),
        calibration_dir.into_os_string(),
        os("--negative-control"),
        negative_control.into_os_string(),
        os("--out"),
        out_dir.clone().into_os_string(),
        os("--mode"),
        os("representative"),
    ]);

    let summary = read_summary(&out_dir);
    assert_eq!(
        summary.decision_status,
        rvscreen::types::DecisionStatus::Positive
    );
    assert_eq!(summary.release_status, ReleaseStatus::Final);
    assert_report_bundle_artifacts(&out_dir);
    assert_candidate_accessions(&out_dir, &["NC_SYNTHV1.1"]);
    assert_eq!(
        coverage_file_map(&out_dir)
            .keys()
            .cloned()
            .collect::<BTreeSet<_>>(),
        BTreeSet::from(["NC_SYNTHV1.1.coverage.tsv".to_string()])
    );
}

#[test]
fn integration_low_titer_positive_sample_stays_detectable() {
    let tempdir = tempdir().expect("tempdir should be created");
    let bundle =
        prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
    let calibration_dir = write_calibration_profile(
        tempdir.path().join("calibration-low-titer"),
        &bundle.version,
        CalibrationProfileOptions {
            rounds: vec![5_000, 10_000],
            max_rounds: 2,
            theta_pos: 0.00005,
            theta_neg: 0.0,
            write_passing_release_gate: true,
            ..CalibrationProfileOptions::default()
        },
    )
    .expect("low-titer calibration profile should be written");
    let negative_control =
        write_negative_control(tempdir.path().join("negative-control.json"), "pass", &[])
            .expect("negative control should be written");
    let (r1, r2) = generate_fastq_pair(
        &FastqPairConfig::new(10_000, 100, 33)
            .with_output_dir(tempdir.path().join("fastq-low-titer"))
            .with_components(vec![
                ReadComponent::new(SyntheticSource::Human, 0.9999),
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    0.0001,
                ),
            ]),
    )
    .expect("FASTQ should be generated");

    let out_dir = tempdir.path().join("report-low-titer");
    run_screen_cli([
        os("screen"),
        os("--input"),
        r1.into_os_string(),
        r2.into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.clone().into_os_string(),
        os("--calibration-profile"),
        calibration_dir.into_os_string(),
        os("--negative-control"),
        negative_control.into_os_string(),
        os("--out"),
        out_dir.clone().into_os_string(),
        os("--mode"),
        os("representative"),
    ]);

    let summary = read_summary(&out_dir);
    assert!(
        matches!(
            summary.decision_status,
            DecisionStatus::Positive | DecisionStatus::Indeterminate
        ),
        "low-titer sample should remain positive or indeterminate, got {:?}",
        summary.decision_status
    );
    assert_eq!(summary.release_status, ReleaseStatus::Final);
    assert_report_bundle_artifacts(&out_dir);
    let candidate_rows = read_candidate_rows(&out_dir);
    assert!(
        candidate_rows.len() <= 1,
        "expected at most one low-titer candidate row"
    );
    if let Some(candidate) = candidate_rows.first() {
        assert_eq!(candidate.accession_or_group, "NC_SYNTHV1.1");
    }
}

#[test]
fn integration_multi_virus_sample_reports_two_candidates() {
    let tempdir = tempdir().expect("tempdir should be created");
    let bundle =
        prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
    let calibration_dir = write_calibration_profile(
        tempdir.path().join("calibration-multi-virus"),
        &bundle.version,
        CalibrationProfileOptions {
            write_passing_release_gate: true,
            ..CalibrationProfileOptions::default()
        },
    )
    .expect("calibration profile should be written");
    let negative_control =
        write_negative_control(tempdir.path().join("negative-control.json"), "pass", &[])
            .expect("negative control should be written");
    let (r1, r2) = generate_fastq_pair(
        &FastqPairConfig::new(200, 100, 44)
            .with_output_dir(tempdir.path().join("fastq-multi-virus"))
            .with_components(vec![
                ReadComponent::new(SyntheticSource::Human, 0.90),
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    0.05,
                ),
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV2.1".to_string())),
                    0.05,
                ),
            ]),
    )
    .expect("multi-virus FASTQ should be generated");
    let out_dir = tempdir.path().join("report-multi-virus");

    run_screen_cli([
        os("screen"),
        os("--input"),
        r1.into_os_string(),
        r2.into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.clone().into_os_string(),
        os("--calibration-profile"),
        calibration_dir.into_os_string(),
        os("--negative-control"),
        negative_control.into_os_string(),
        os("--out"),
        out_dir.clone().into_os_string(),
        os("--mode"),
        os("representative"),
    ]);

    let summary = read_summary(&out_dir);
    assert_eq!(summary.decision_status, DecisionStatus::Positive);
    assert_eq!(summary.release_status, ReleaseStatus::Final);
    assert_report_bundle_artifacts(&out_dir);
    assert_candidate_accessions(&out_dir, &["NC_SYNTHV1.1", "NC_SYNTHV2.1"]);
    assert_eq!(
        coverage_file_map(&out_dir)
            .keys()
            .cloned()
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "NC_SYNTHV1.1.coverage.tsv".to_string(),
            "NC_SYNTHV2.1.coverage.tsv".to_string(),
        ])
    );
}

#[test]
fn integration_low_complexity_contamination_does_not_trigger_false_positive() {
    let tempdir = tempdir().expect("tempdir should be created");
    let bundle =
        prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
    let calibration_dir = write_calibration_profile(
        tempdir.path().join("calibration-low-complexity"),
        &bundle.version,
        CalibrationProfileOptions {
            write_passing_release_gate: true,
            ..CalibrationProfileOptions::default()
        },
    )
    .expect("calibration profile should be written");
    let negative_control =
        write_negative_control(tempdir.path().join("negative-control.json"), "pass", &[])
            .expect("negative control should be written");
    let (r1, r2) = generate_fastq_pair(
        &FastqPairConfig::new(200, 100, 2301)
            .with_output_dir(tempdir.path().join("fastq-low-complexity"))
            .with_components(vec![
                ReadComponent::new(SyntheticSource::Human, 0.95),
                ReadComponent::new(SyntheticSource::LowComplexity, 0.05),
            ]),
    )
    .expect("FASTQ should be generated");

    let out_dir = tempdir.path().join("report-low-complexity");
    run_screen_cli([
        os("screen"),
        os("--input"),
        r1.into_os_string(),
        r2.into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.clone().into_os_string(),
        os("--calibration-profile"),
        calibration_dir.into_os_string(),
        os("--negative-control"),
        negative_control.into_os_string(),
        os("--out"),
        out_dir.clone().into_os_string(),
        os("--mode"),
        os("representative"),
    ]);

    let summary = read_summary(&out_dir);
    assert_eq!(summary.decision_status, DecisionStatus::Negative);
    assert_eq!(summary.release_status, ReleaseStatus::Final);
    assert!(
        summary.qc_passing_fragments < summary.input_fragments,
        "low-complexity contamination should be filtered at QC"
    );
    assert_report_bundle_artifacts(&out_dir);
    assert_candidate_accessions(&out_dir, &[]);
}

#[test]
fn integration_negative_control_clean_vs_contaminated_degrades_result() {
    let tempdir = tempdir().expect("tempdir should be created");
    let bundle =
        prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
    let calibration_dir = write_calibration_profile(
        tempdir.path().join("calibration-negative-control"),
        &bundle.version,
        CalibrationProfileOptions {
            max_background_ratio: 1.5,
            write_passing_release_gate: true,
            ..CalibrationProfileOptions::default()
        },
    )
    .expect("calibration profile should be written");
    let clean_control = write_negative_control(
        tempdir.path().join("negative-control-clean.json"),
        "pass",
        &[],
    )
    .expect("clean negative control should be written");
    let contaminated_control = write_negative_control(
        tempdir.path().join("negative-control-contaminated.json"),
        "pass",
        &[("NC_SYNTHV1.1", 0.10)],
    )
    .expect("contaminated negative control should be written");
    let (r1, r2) = generate_fastq_pair(
        &FastqPairConfig::new(200, 100, 62)
            .with_output_dir(tempdir.path().join("fastq-negative-control"))
            .with_components(vec![
                ReadComponent::new(SyntheticSource::Human, 0.95),
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    0.05,
                ),
            ]),
    )
    .expect("FASTQ should be generated");

    let clean_out = tempdir.path().join("report-negative-control-clean");
    run_screen_cli([
        os("screen"),
        os("--input"),
        r1.clone().into_os_string(),
        r2.clone().into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.clone().into_os_string(),
        os("--calibration-profile"),
        calibration_dir.clone().into_os_string(),
        os("--negative-control"),
        clean_control.into_os_string(),
        os("--out"),
        clean_out.clone().into_os_string(),
        os("--mode"),
        os("representative"),
    ]);

    let degraded_out = tempdir.path().join("report-negative-control-degraded");
    run_screen_cli([
        os("screen"),
        os("--input"),
        r1.into_os_string(),
        r2.into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.into_os_string(),
        os("--calibration-profile"),
        calibration_dir.into_os_string(),
        os("--negative-control"),
        contaminated_control.into_os_string(),
        os("--out"),
        degraded_out.clone().into_os_string(),
        os("--mode"),
        os("representative"),
    ]);

    assert_eq!(
        read_summary(&clean_out).decision_status,
        DecisionStatus::Positive
    );
    assert_eq!(
        read_summary(&clean_out).release_status,
        ReleaseStatus::Final
    );

    let degraded_summary = read_summary(&degraded_out);
    assert_eq!(
        degraded_summary.decision_status,
        DecisionStatus::Indeterminate
    );
    assert_eq!(degraded_summary.release_status, ReleaseStatus::Final);
    let degraded_candidate = read_candidate_rows(&degraded_out)
        .into_iter()
        .find(|candidate| candidate.accession_or_group == "NC_SYNTHV1.1")
        .expect("degraded run should still emit candidate row for target virus");
    assert_eq!(degraded_candidate.decision, "indeterminate");
    assert!(
        degraded_candidate.background_ratio < 1.5,
        "background ratio should fall below positive threshold when control exceeds sample signal"
    );
}

#[test]
fn integration_streaming_mode_is_provisional() {
    let tempdir = tempdir().expect("tempdir should be created");
    let bundle =
        prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
    let calibration_dir = write_calibration_profile(
        tempdir.path().join("calibration-streaming"),
        &bundle.version,
        CalibrationProfileOptions {
            sampling_mode: "streaming",
            rounds: vec![50, 100],
            max_rounds: 2,
            write_passing_release_gate: true,
            ..CalibrationProfileOptions::default()
        },
    )
    .expect("streaming calibration profile should be written");
    let negative_control =
        write_negative_control(tempdir.path().join("negative-control.json"), "pass", &[])
            .expect("negative control should be written");
    let (r1, r2) = generate_fastq_pair(
        &FastqPairConfig::new(200, 100, 73)
            .with_output_dir(tempdir.path().join("fastq-streaming"))
            .with_components(vec![
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    0.05,
                ),
                ReadComponent::new(SyntheticSource::Human, 0.95),
            ]),
    )
    .expect("FASTQ should be generated");
    let out_dir = tempdir.path().join("report-streaming");

    run_screen_cli([
        os("screen"),
        os("--input"),
        r1.into_os_string(),
        r2.into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.clone().into_os_string(),
        os("--calibration-profile"),
        calibration_dir.into_os_string(),
        os("--negative-control"),
        negative_control.into_os_string(),
        os("--out"),
        out_dir.clone().into_os_string(),
        os("--mode"),
        os("streaming"),
    ]);

    let summary = read_summary(&out_dir);
    assert_eq!(summary.decision_status, DecisionStatus::Positive);
    assert_eq!(summary.release_status, ReleaseStatus::Provisional);
    assert_report_bundle_artifacts(&out_dir);
    assert_candidate_accessions(&out_dir, &["NC_SYNTHV1.1"]);
}

#[test]
fn integration_reproducibility_same_seed_same_outputs() {
    let tempdir = tempdir().expect("tempdir should be created");
    let bundle =
        prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
    let calibration_dir = write_calibration_profile(
        tempdir.path().join("calibration-reproducibility"),
        &bundle.version,
        CalibrationProfileOptions {
            write_passing_release_gate: true,
            ..CalibrationProfileOptions::default()
        },
    )
    .expect("calibration profile should be written");
    let negative_control =
        write_negative_control(tempdir.path().join("negative-control.json"), "pass", &[])
            .expect("negative control should be written");
    let (r1, r2) = generate_fastq_pair(
        &FastqPairConfig::new(200, 100, 84)
            .with_output_dir(tempdir.path().join("fastq-reproducibility"))
            .with_components(vec![
                ReadComponent::new(SyntheticSource::Human, 0.95),
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    0.05,
                ),
            ]),
    )
    .expect("FASTQ should be generated");

    let first_out = tempdir.path().join("report-reproducibility-a");
    run_screen_cli([
        os("screen"),
        os("--input"),
        r1.clone().into_os_string(),
        r2.clone().into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.clone().into_os_string(),
        os("--calibration-profile"),
        calibration_dir.clone().into_os_string(),
        os("--negative-control"),
        negative_control.clone().into_os_string(),
        os("--out"),
        first_out.clone().into_os_string(),
        os("--mode"),
        os("representative"),
    ]);

    let second_out = tempdir.path().join("report-reproducibility-b");
    run_screen_cli([
        os("screen"),
        os("--input"),
        r1.into_os_string(),
        r2.into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.into_os_string(),
        os("--calibration-profile"),
        calibration_dir.into_os_string(),
        os("--negative-control"),
        negative_control.into_os_string(),
        os("--out"),
        second_out.clone().into_os_string(),
        os("--mode"),
        os("representative"),
    ]);

    assert_eq!(read_summary(&first_out), read_summary(&second_out));
    assert_eq!(
        fs::read_to_string(first_out.join("candidate_calls.tsv"))
            .expect("first candidate_calls.tsv should be readable"),
        fs::read_to_string(second_out.join("candidate_calls.tsv"))
            .expect("second candidate_calls.tsv should be readable")
    );
    assert_eq!(
        fs::read_to_string(first_out.join("rounds.tsv"))
            .expect("first rounds.tsv should be readable"),
        fs::read_to_string(second_out.join("rounds.tsv"))
            .expect("second rounds.tsv should be readable")
    );
    assert_eq!(
        fs::read_to_string(first_out.join("run_manifest.json"))
            .expect("first run_manifest.json should be readable"),
        fs::read_to_string(second_out.join("run_manifest.json"))
            .expect("second run_manifest.json should be readable")
    );
    assert_eq!(
        coverage_file_map(&first_out),
        coverage_file_map(&second_out)
    );
}

#[test]
fn integration_full_chain_ref_build_calibrate_screen_audit_verify() {
    let tempdir = tempdir().expect("tempdir should be created");
    let bundle = prepare_reference_bundle_via_cli(tempdir.path())
        .expect("reference bundle should be prepared through CLI");
    let benchmark_manifest =
        write_calibration_benchmark_manifest(tempdir.path(), ManifestFormat::Yaml)
            .expect("benchmark manifest should be written");
    let calibration_dir = tempdir.path().join("calibration-full-chain");

    run_cli([
        os("calibrate"),
        os("--reference-bundle"),
        bundle.bundle_dir.clone().into_os_string(),
        os("--benchmark-manifest"),
        benchmark_manifest.into_os_string(),
        os("--out"),
        calibration_dir.clone().into_os_string(),
    ]);

    for artifact in [
        "profile.toml",
        "thresholds.toml",
        "benchmark_manifest.json",
        "benchmark_summary.tsv",
        "release_gate.json",
        "checksum.sha256",
    ] {
        assert!(
            calibration_dir.join(artifact).exists(),
            "missing calibration artifact {artifact}"
        );
    }

    let profile = fs::read_to_string(calibration_dir.join("profile.toml"))
        .expect("profile.toml should be readable");
    assert!(
        profile.contains(&format!("reference_bundle = \"{}\"", bundle.version)),
        "profile should bind the reference bundle version: {profile}"
    );

    let negative_control =
        write_negative_control(tempdir.path().join("negative-control.json"), "pass", &[])
            .expect("negative control should be written");
    let (r1, r2) = generate_fastq_pair(
        &FastqPairConfig::new(200, 100, 95)
            .with_output_dir(tempdir.path().join("fastq-full-chain"))
            .with_components(vec![
                ReadComponent::new(SyntheticSource::Human, 0.95),
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    0.05,
                ),
            ]),
    )
    .expect("FASTQ should be generated");
    let report_dir = tempdir.path().join("report-full-chain");

    run_screen_cli([
        os("screen"),
        os("--input"),
        r1.into_os_string(),
        r2.into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.clone().into_os_string(),
        os("--calibration-profile"),
        calibration_dir.clone().into_os_string(),
        os("--negative-control"),
        negative_control.into_os_string(),
        os("--out"),
        report_dir.clone().into_os_string(),
        os("--mode"),
        os("representative"),
    ]);

    let summary = read_summary(&report_dir);
    assert_eq!(summary.decision_status, DecisionStatus::Positive);
    assert_eq!(summary.release_status, ReleaseStatus::Final);
    assert_report_bundle_artifacts(&report_dir);
    assert_candidate_accessions(&report_dir, &["NC_SYNTHV1.1"]);

    let audit_output = run_cli_capture([
        os("audit"),
        os("verify"),
        os("--report-bundle"),
        report_dir.into_os_string(),
        os("--reference-bundle"),
        bundle.bundle_dir.into_os_string(),
        os("--calibration-profile"),
        calibration_dir.into_os_string(),
    ]);
    let stdout = String::from_utf8_lossy(&audit_output.stdout);
    assert!(stdout.contains("PASS rvscreen audit verify"), "{stdout}");
}

#[test]
fn audit_verify_passes_valid_report_bundle() {
    let tempdir = tempdir().expect("tempdir should be created");
    let fixture = write_audit_verify_fixture(tempdir.path(), "rvscreen_ref_2026.04.21-r1")
        .expect("audit fixture should be written");

    let output = Command::new(binary_path())
        .args([
            os("audit"),
            os("verify"),
            os("--report-bundle"),
            fixture.report_dir.clone().into_os_string(),
            os("--reference-bundle"),
            fixture.reference_dir.clone().into_os_string(),
            os("--calibration-profile"),
            fixture.profile_dir.clone().into_os_string(),
        ])
        .output()
        .expect("rvscreen audit verify should execute");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(
        output.status.success(),
        "stdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(stdout.contains("PASS rvscreen audit verify"), "{stdout}");
}

#[test]
fn audit_verify_missing_file_returns_exit_code_1() {
    let tempdir = tempdir().expect("tempdir should be created");
    let fixture = write_audit_verify_fixture(tempdir.path(), "rvscreen_ref_2026.04.21-r1")
        .expect("audit fixture should be written");
    fs::remove_file(fixture.report_dir.join("rounds.tsv")).expect("rounds.tsv should be removed");

    let output = Command::new(binary_path())
        .args([
            os("audit"),
            os("verify"),
            os("--report-bundle"),
            fixture.report_dir.into_os_string(),
        ])
        .output()
        .expect("rvscreen audit verify should execute");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert_eq!(
        output.status.code(),
        Some(1),
        "stdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(stderr.contains("rounds.tsv"), "{stderr}");
    assert!(
        stderr.contains("missing required file(s): rounds.tsv"),
        "{stderr}"
    );
}

struct ReferenceBundleFixture {
    bundle_dir: PathBuf,
    version: String,
}

struct ReferenceBundleInputsFixture {
    host_fasta: PathBuf,
    virus_fasta: PathBuf,
    decoy_fasta: PathBuf,
    manifest: PathBuf,
    taxonomy: PathBuf,
}

struct AuditVerifyFixture {
    report_dir: PathBuf,
    reference_dir: PathBuf,
    profile_dir: PathBuf,
}

fn prepare_reference_bundle(base_dir: &Path) -> io::Result<ReferenceBundleFixture> {
    let inputs = prepare_reference_bundle_inputs(base_dir)?;
    let bundle_dir = base_dir.join("reference-bundle");
    let outcome = build_reference_bundle(&BuildReferenceBundleRequest {
        host_fasta: inputs.host_fasta,
        virus_fasta: inputs.virus_fasta,
        decoy_fasta: Some(inputs.decoy_fasta),
        manifest: inputs.manifest,
        taxonomy: inputs.taxonomy,
        out_dir: bundle_dir.clone(),
    })
    .map_err(io_other)?;

    Ok(ReferenceBundleFixture {
        bundle_dir,
        version: outcome.bundle.version,
    })
}

fn prepare_reference_bundle_via_cli(base_dir: &Path) -> io::Result<ReferenceBundleFixture> {
    let inputs = prepare_reference_bundle_inputs(base_dir)?;
    let bundle_dir = base_dir.join("reference-bundle-cli");
    run_cli([
        os("ref"),
        os("build"),
        os("--host-fasta"),
        inputs.host_fasta.into_os_string(),
        os("--virus-fasta"),
        inputs.virus_fasta.into_os_string(),
        os("--decoy-fasta"),
        inputs.decoy_fasta.into_os_string(),
        os("--manifest"),
        inputs.manifest.into_os_string(),
        os("--taxonomy"),
        inputs.taxonomy.into_os_string(),
        os("--out"),
        bundle_dir.clone().into_os_string(),
    ]);

    let bundle: BundleToml = toml::from_str(
        &fs::read_to_string(bundle_dir.join("bundle.toml"))
            .map_err(|error| io::Error::other(error.to_string()))?,
    )
    .map_err(io_other)?;

    Ok(ReferenceBundleFixture {
        bundle_dir,
        version: bundle.version,
    })
}

fn prepare_reference_bundle_inputs(base_dir: &Path) -> io::Result<ReferenceBundleInputsFixture> {
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

    Ok(ReferenceBundleInputsFixture {
        host_fasta,
        virus_fasta,
        decoy_fasta,
        manifest: manifest_out,
        taxonomy,
    })
}

fn write_audit_verify_fixture(
    base_dir: &Path,
    reference_version: &str,
) -> io::Result<AuditVerifyFixture> {
    let report_dir = base_dir.join("audit-report");
    let reference_dir = base_dir.join("audit-reference");
    let profile_dir = base_dir.join("audit-profile");

    let summary = audit_sample_summary(reference_version);
    let manifest = audit_run_manifest(reference_version);
    let candidates = audit_candidate_calls();
    let rounds = audit_rounds();

    ReportWriter::write(&report_dir, &summary, &candidates, &rounds, &manifest)
        .map_err(io_other)?;
    write_reference_bundle_stub(&reference_dir, reference_version)?;
    write_calibration_profile_stub(&profile_dir, reference_version)?;

    Ok(AuditVerifyFixture {
        report_dir,
        reference_dir,
        profile_dir,
    })
}

fn write_reference_bundle_stub(dir: &Path, version: &str) -> io::Result<()> {
    fs::create_dir_all(dir)?;
    fs::write(
        dir.join("bundle.toml"),
        toml::to_string_pretty(&rvscreen::types::BundleToml {
            version: version.to_string(),
            created_at: "2026-04-21T00:00:00Z".to_string(),
            included_layers: vec!["host_backbone".to_string(), "viral_panel".to_string()],
        })
        .map_err(io_other)?,
    )
}

fn write_calibration_profile_stub(dir: &Path, reference_version: &str) -> io::Result<()> {
    fs::create_dir_all(dir)?;
    fs::write(
        dir.join("profile.toml"),
        toml::to_string_pretty(&rvscreen::types::ProfileToml {
            profile_id: "rvscreen_calib_2026.04.21-r1".to_string(),
            status: "release_candidate".to_string(),
            reference_bundle: reference_version.to_string(),
            backend: "minimap2".to_string(),
            preset: "sr-conservative".to_string(),
            seed: 20_260_421,
            supported_input: vec!["fastq".to_string(), "fastq.gz".to_string()],
            supported_read_type: vec!["illumina_pe_shortread".to_string()],
            negative_control_required: false,
            sampling: SamplingConfig {
                mode: "representative".to_string(),
                rounds: vec![50_000, 100_000],
                max_rounds: 2,
            },
            fragment_rules: rvscreen::types::FragmentRules {
                min_mapq: 20,
                min_as_diff: 12,
                max_nm: 8,
                require_pair_consistency: true,
            },
            candidate_rules: rvscreen::types::CandidateRules {
                min_nonoverlap_fragments: 1,
                min_breadth: 0.0,
                max_background_ratio: 0.0,
            },
            decision_rules: rvscreen::types::DecisionRules {
                theta_pos: 0.001,
                theta_neg: 0.0,
                allow_indeterminate: true,
            },
        })
        .map_err(io_other)?,
    )
}

fn audit_sample_summary(reference_version: &str) -> SampleSummary {
    SampleSummary {
        sample_id: "S-AUDIT-001".to_string(),
        reference_bundle_version: reference_version.to_string(),
        calibration_profile_version: "rvscreen_calib_2026.04.21-r1".to_string(),
        backend: "minimap2".to_string(),
        seed: 20_260_421,
        input_fragments: 100_000,
        qc_passing_fragments: 95_000,
        sampled_fragments: 50_000,
        rounds_run: 2,
        stop_reason: StopReason::PositiveBoundaryCrossed,
        decision_status: DecisionStatus::Positive,
        release_status: ReleaseStatus::Provisional,
    }
}

fn audit_run_manifest(reference_version: &str) -> RunManifest {
    RunManifest {
        reference_bundle_version: reference_version.to_string(),
        calibration_profile_version: "rvscreen_calib_2026.04.21-r1".to_string(),
        backend: "minimap2".to_string(),
        seed: 20_260_421,
        input_files: vec![
            "reads_R1.fastq.gz".to_string(),
            "reads_R2.fastq.gz".to_string(),
        ],
    }
}

fn audit_rounds() -> Vec<RoundRecord> {
    vec![
        RoundRecord {
            sampled_fragments: 50_000,
            accepted_virus: 4,
            decision_status: DecisionStatus::Indeterminate,
        },
        RoundRecord {
            sampled_fragments: 100_000,
            accepted_virus: 12,
            decision_status: DecisionStatus::Positive,
        },
    ]
}

fn audit_candidate_calls() -> Vec<CandidateCall> {
    vec![CandidateCall {
        virus_name: "Synthetic audit virus".to_string(),
        taxid: 9_999,
        accession_or_group: "NC_SYNTHV1.1".to_string(),
        accepted_fragments: 12,
        nonoverlap_fragments: 4,
        raw_fraction: 0.00012,
        unique_fraction: 0.00013,
        fraction_ci_95: [0.00008, 0.00018],
        clopper_pearson_upper: 0.0002,
        breadth: 0.01,
        ambiguous_fragments: 1,
        background_ratio: 0.0,
        decision: DecisionStatus::Positive,
        decision_reasons: vec![
            format!(
                "fraction_ci_95_label={}",
                rvscreen::decision::SAMPLING_ONLY_CI_LABEL
            ),
            "unique_fraction_above_theta_pos".to_string(),
        ],
        evidence_strength: EvidenceStrength::High,
    }]
}

#[derive(Debug, Clone, PartialEq)]
struct CalibrationProfileOptions {
    negative_control_required: bool,
    max_background_ratio: f64,
    write_passing_release_gate: bool,
    sampling_mode: &'static str,
    rounds: Vec<u64>,
    max_rounds: u64,
    theta_pos: f64,
    theta_neg: f64,
}

impl Default for CalibrationProfileOptions {
    fn default() -> Self {
        Self {
            negative_control_required: false,
            max_background_ratio: 0.0,
            write_passing_release_gate: false,
            sampling_mode: "representative",
            rounds: vec![50, 100],
            max_rounds: 2,
            theta_pos: 0.01,
            theta_neg: 0.0001,
        }
    }
}

fn write_calibration_profile(
    dir: PathBuf,
    reference_bundle: &str,
    options: CalibrationProfileOptions,
) -> std::io::Result<PathBuf> {
    fs::create_dir_all(&dir)?;
    let rounds = options
        .rounds
        .iter()
        .map(u64::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    fs::write(
        dir.join("profile.toml"),
        format!(
            "profile_id = \"rvscreen_calib_test\"\nstatus = \"release_candidate\"\nreference_bundle = \"{reference_bundle}\"\nbackend = \"minimap2\"\npreset = \"sr-conservative\"\nseed = 20260420\nsupported_input = [\"fastq\", \"fastq.gz\", \"bam\", \"ubam\", \"cram\"]\nsupported_read_type = [\"illumina_pe_shortread\"]\nnegative_control_required = {}\n\n[sampling]\nmode = \"{}\"\nrounds = [{}]\nmax_rounds = {}\n\n[fragment_rules]\nmin_mapq = 0\nmin_as_diff = 0\nmax_nm = 100\nrequire_pair_consistency = true\n\n[candidate_rules]\nmin_nonoverlap_fragments = 1\nmin_breadth = 0.0\nmax_background_ratio = {}\n\n[decision_rules]\ntheta_pos = {}\ntheta_neg = {}\nallow_indeterminate = true\n",
            options.negative_control_required,
            options.sampling_mode,
            rounds,
            options.max_rounds,
            options.max_background_ratio,
            options.theta_pos,
            options.theta_neg,
        ),
    )?;

    if options.write_passing_release_gate {
        fs::write(
            dir.join("release_gate.json"),
            serde_json::to_vec_pretty(&passing_release_gate()).map_err(io_other)?,
        )?;
    }

    Ok(dir)
}

fn passing_release_gate() -> CalibrationReleaseGate {
    CalibrationReleaseGate {
        backend_gate: BackendGate {
            status: GateStatus::Pass,
            details: "integration-test backend gate".to_string(),
        },
        reference_gate: ReferenceGate {
            status: GateStatus::Pass,
            details: "integration-test reference gate".to_string(),
        },
        specificity_gate: SpecificityGate {
            status: GateStatus::Pass,
            details: "integration-test specificity gate".to_string(),
            negative_samples: 1,
            false_positives: 0,
        },
        sensitivity_gate: SensitivityGate {
            status: GateStatus::Pass,
            details: "integration-test sensitivity gate".to_string(),
            spike_in_detected: 1,
            spike_in_total: 1,
        },
    }
}

fn write_negative_control(
    path: PathBuf,
    control_status: &str,
    candidates: &[(&str, f64)],
) -> io::Result<PathBuf> {
    let payload = json!({
        "control_id": "neg-001",
        "control_status": control_status,
        "candidates": candidates
            .iter()
            .map(|(accession_or_group, unique_fraction)| json!({
                "accession_or_group": accession_or_group,
                "unique_fraction": unique_fraction,
            }))
            .collect::<Vec<_>>(),
    });
    fs::write(
        &path,
        serde_json::to_vec_pretty(&payload).map_err(io_other)?,
    )?;
    Ok(path)
}

#[derive(Clone, Copy)]
enum ManifestFormat {
    Json,
    Yaml,
}

#[test]
fn benchmark_manifest_supports_json_helper_path() {
    let tempdir = tempdir().expect("tempdir should be created");
    let manifest_path = write_calibration_benchmark_manifest(tempdir.path(), ManifestFormat::Json)
        .expect("JSON benchmark manifest should be written");

    assert_eq!(
        manifest_path.extension().and_then(|ext| ext.to_str()),
        Some("json")
    );
}

fn write_calibration_benchmark_manifest(
    base_dir: &Path,
    format: ManifestFormat,
) -> io::Result<PathBuf> {
    let benchmark_dir = base_dir.join("benchmarks");
    fs::create_dir_all(&benchmark_dir)?;
    let negative_dir = benchmark_dir.join("negative");
    let spike_dir = benchmark_dir.join("spike");

    generate_fastq_pair(
        &FastqPairConfig::new(200, 100, 44)
            .with_output_dir(&negative_dir)
            .with_components(vec![ReadComponent::new(SyntheticSource::Human, 1.0)]),
    )?;
    generate_fastq_pair(
        &FastqPairConfig::new(200, 100, 55)
            .with_output_dir(&spike_dir)
            .with_components(vec![
                ReadComponent::new(SyntheticSource::Human, 0.95),
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    0.05,
                ),
            ]),
    )?;

    let (manifest_path, manifest_body) = match format {
        ManifestFormat::Json => (
            benchmark_dir.join("manifest.json"),
            r#"{
  "profile": {
    "supported_input": ["fastq"],
    "supported_read_type": ["illumina_pe_shortread"],
    "negative_control_required": false,
    "sampling": {
      "mode": "representative",
      "rounds": [50, 100],
      "max_rounds": 2
    },
    "fragment_rules": {
      "min_mapq": 0,
      "min_as_diff": 0,
      "max_nm": 100,
      "require_pair_consistency": true
    },
    "candidate_rules": {
      "min_nonoverlap_fragments": 1,
      "min_breadth": 0.0,
      "max_background_ratio": 0.0
    },
    "decision_rules": {
      "theta_pos": 0.01,
      "theta_neg": 0.0001,
      "allow_indeterminate": true
    }
  },
  "datasets": [
    {
      "dataset_id": "negative-human",
      "input": ["negative/synthetic_R1.fastq", "negative/synthetic_R2.fastq"],
      "expected_outcome": "negative"
    },
    {
      "dataset_id": "spike-5pct",
      "input": ["spike/synthetic_R1.fastq", "spike/synthetic_R2.fastq"],
      "expected_outcome": "spike_in",
      "expected_target": "NC_SYNTHV1.1"
    }
  ]
}"#,
        ),
        ManifestFormat::Yaml => (
            benchmark_dir.join("manifest.yaml"),
            r#"profile:
  supported_input:
    - fastq
  supported_read_type:
    - illumina_pe_shortread
  negative_control_required: false
  sampling:
    mode: representative
    rounds:
      - 50
      - 100
    max_rounds: 2
  fragment_rules:
    min_mapq: 0
    min_as_diff: 0
    max_nm: 100
    require_pair_consistency: true
  candidate_rules:
    min_nonoverlap_fragments: 1
    min_breadth: 0.0
    max_background_ratio: 0.0
  decision_rules:
    theta_pos: 0.01
    theta_neg: 0.0001
    allow_indeterminate: true
datasets:
  - dataset_id: negative-human
    input:
      - negative/synthetic_R1.fastq
      - negative/synthetic_R2.fastq
    expected_outcome: negative
  - dataset_id: spike-5pct
    input:
      - spike/synthetic_R1.fastq
      - spike/synthetic_R2.fastq
    expected_outcome: spike_in
    expected_target: NC_SYNTHV1.1
"#,
        ),
    };
    fs::write(&manifest_path, manifest_body)?;

    Ok(manifest_path)
}

fn run_screen_cli(args: impl IntoIterator<Item = OsString>) {
    run_cli(args)
}

fn run_cli_capture(args: impl IntoIterator<Item = OsString>) -> std::process::Output {
    let output = Command::new(binary_path())
        .args(args)
        .output()
        .expect("rvscreen binary should execute");

    if !output.status.success() {
        panic!(
            "rvscreen CLI failed\nstatus: {:?}\nstdout:\n{}\nstderr:\n{}",
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
    }

    output
}

fn run_cli(args: impl IntoIterator<Item = OsString>) {
    let _ = run_cli_capture(args);
}

fn binary_path() -> PathBuf {
    std::env::var_os("CARGO_BIN_EXE_rvscreen")
        .map(PathBuf::from)
        .expect("CARGO_BIN_EXE_rvscreen should be set for integration tests")
}

fn read_summary(report_dir: &Path) -> SampleSummary {
    serde_json::from_str(
        &fs::read_to_string(report_dir.join("sample_summary.json"))
            .expect("sample_summary.json should be readable"),
    )
    .expect("sample summary JSON should parse")
}

#[derive(Debug, Deserialize)]
struct CandidateRow {
    accession_or_group: String,
    decision: String,
    background_ratio: f64,
    decision_reasons: String,
}

fn read_candidate_rows(report_dir: &Path) -> Vec<CandidateRow> {
    ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(report_dir.join("candidate_calls.tsv"))
        .expect("candidate_calls.tsv should be readable")
        .deserialize()
        .collect::<Result<Vec<_>, _>>()
        .expect("candidate_calls.tsv rows should deserialize")
}

fn assert_candidate_accessions(report_dir: &Path, expected: &[&str]) {
    let rows = read_candidate_rows(report_dir);
    let actual = rows
        .iter()
        .map(|row| row.accession_or_group.clone())
        .collect::<BTreeSet<_>>();
    let expected = expected.iter().map(|value| (*value).to_string()).collect();
    assert_eq!(actual, expected);

    for row in rows {
        let reasons: Vec<String> =
            serde_json::from_str(&row.decision_reasons).expect("decision reasons should parse");
        assert!(
            reasons.iter().any(|reason| {
                reason == &format!("fraction_ci_95_label={SAMPLING_ONLY_CI_LABEL}")
            }),
            "candidate `{}` should preserve sampling-only CI marker",
            row.accession_or_group
        );
    }
}

fn assert_report_bundle_artifacts(report_dir: &Path) {
    for artifact in [
        "sample_summary.json",
        "candidate_calls.tsv",
        "run_manifest.json",
        "rounds.tsv",
        "logs/run.log",
    ] {
        assert!(
            report_dir.join(artifact).exists(),
            "missing report artifact {artifact}"
        );
    }
    assert!(
        report_dir.join("coverage").is_dir(),
        "coverage dir should exist"
    );
}

fn coverage_file_map(report_dir: &Path) -> BTreeMap<String, String> {
    let coverage_dir = report_dir.join("coverage");
    if !coverage_dir.is_dir() {
        return BTreeMap::new();
    }

    fs::read_dir(&coverage_dir)
        .expect("coverage dir should be readable")
        .map(|entry| {
            let entry = entry.expect("coverage dir entry should be readable");
            let path = entry.path();
            (
                entry.file_name().to_string_lossy().into_owned(),
                fs::read_to_string(path).expect("coverage file should be readable"),
            )
        })
        .collect()
}

fn os(value: &str) -> OsString {
    OsString::from(value)
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
