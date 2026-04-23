#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryTransition {
    PortedAsIs,
    Adapted,
    Replaced,
    PendingRemoval,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryMapping {
    pub boundary: &'static str,
    pub owner_module: &'static str,
    pub primary_api: &'static [&'static str],
    pub transition: BoundaryTransition,
    pub intent: &'static str,
}

pub const V2_BOUNDARY_MAP: &[BoundaryMapping] = &[
    BoundaryMapping {
        boundary: "cli",
        owner_module: "crate::cli",
        primary_api: &["Cli", "Commands", "RefAction", "AuditAction"],
        transition: BoundaryTransition::Adapted,
        intent: "Owns the fixed V2 command family and argument contracts only.",
    },
    BoundaryMapping {
        boundary: "reference bundle",
        owner_module: "crate::reference",
        primary_api: &["build_reference_bundle", "BuildReferenceBundleRequest"],
        transition: BoundaryTransition::Adapted,
        intent: "Owns reference bundle construction and on-disk bundle artifact shape.",
    },
    BoundaryMapping {
        boundary: "calibration profile",
        owner_module: "crate::calibration",
        primary_api: &["run_calibration", "load_profile", "load_reference_bundle"],
        transition: BoundaryTransition::Adapted,
        intent: "Owns calibration profile generation plus release-gate artifacts bound to a reference bundle.",
    },
    BoundaryMapping {
        boundary: "screen pipeline",
        owner_module: "crate::pipeline",
        primary_api: &["run_screen", "PreparedScreenRunner"],
        transition: BoundaryTransition::Adapted,
        intent: "Owns orchestration across I/O, sampling, align, adjudicate, aggregate, decision, and report without redefining domain rules.",
    },
    BoundaryMapping {
        boundary: "i/o readers",
        owner_module: "crate::io",
        primary_api: &["ScreenInput", "FragmentReaderFactory"],
        transition: BoundaryTransition::Adapted,
        intent: "Owns input-kind classification plus fragment-reader construction for FASTQ, BAM/uBAM, and CRAM.",
    },
    BoundaryMapping {
        boundary: "sampling",
        owner_module: "crate::sampling",
        primary_api: &["RepresentativeSampler", "StreamingSampler", "RoundManager"],
        transition: BoundaryTransition::PortedAsIs,
        intent: "Owns sampling-mode rules and round progression independently from pipeline execution.",
    },
    BoundaryMapping {
        boundary: "align",
        owner_module: "crate::align",
        primary_api: &["CompetitiveAligner", "CompetitivePreset"],
        transition: BoundaryTransition::PortedAsIs,
        intent: "Owns competitive alignment and reference/index loading for the current minimap2 backend.",
    },
    BoundaryMapping {
        boundary: "adjudicate",
        owner_module: "crate::adjudicate",
        primary_api: &["FragmentAdjudicator", "AdjudicationInputs", "AdjudicationResult"],
        transition: BoundaryTransition::PortedAsIs,
        intent: "Owns fragment-level host/virus/ambiguous classification from alignment evidence.",
    },
    BoundaryMapping {
        boundary: "aggregate",
        owner_module: "crate::aggregate",
        primary_api: &["CandidateAggregator", "AggregationLevel"],
        transition: BoundaryTransition::PortedAsIs,
        intent: "Owns candidate-level breadth/nonoverlap aggregation from adjudicated fragment evidence.",
    },
    BoundaryMapping {
        boundary: "decision",
        owner_module: "crate::decision",
        primary_api: &["DecisionEngine", "apply_negative_control", "ReleaseGate"],
        transition: BoundaryTransition::Adapted,
        intent: "Owns candidate decision semantics, negative-control comparison, and release-status gating as separate concerns.",
    },
    BoundaryMapping {
        boundary: "report",
        owner_module: "crate::report",
        primary_api: &["write_report_bundle", "ReportWriter"],
        transition: BoundaryTransition::Adapted,
        intent: "Owns report-bundle materialization from pipeline outputs and preserves report/audit contract shape.",
    },
    BoundaryMapping {
        boundary: "audit",
        owner_module: "crate::audit",
        primary_api: &["run_audit_verify", "AuditVerifyReport"],
        transition: BoundaryTransition::Adapted,
        intent: "Owns machine-checkable validation of report bundles against reference and calibration artifacts.",
    },
];

pub const LEGACY_BOUNDARY_SEAMS: &[BoundaryMapping] = &[
    BoundaryMapping {
        boundary: "pipeline input classification seam",
        owner_module: "crate::pipeline::InputSpec",
        primary_api: &["legacy internal enum removed in task 3"],
        transition: BoundaryTransition::PendingRemoval,
        intent: "Legacy input-shape ownership belonged in crate::io and is intentionally retired from pipeline orchestration.",
    },
    BoundaryMapping {
        boundary: "reference index selection seam",
        owner_module: "crate::align::discover_reference_input",
        primary_api: &["legacy bundle lookup retained until task 10"],
        transition: BoundaryTransition::PendingRemoval,
        intent: "Current aligner-side bundle probing is kept only until aligner lifecycle work lands in task 10.",
    },
];

pub fn boundary_mapping(name: &str) -> Option<&'static BoundaryMapping> {
    V2_BOUNDARY_MAP
        .iter()
        .chain(LEGACY_BOUNDARY_SEAMS.iter())
        .find(|mapping| mapping.boundary == name)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    #[test]
    fn v2_boundary_map_lists_every_expected_domain_boundary() {
        let boundaries = V2_BOUNDARY_MAP
            .iter()
            .map(|mapping| mapping.boundary)
            .collect::<Vec<_>>();

        assert_eq!(
            boundaries,
            vec![
                "cli",
                "reference bundle",
                "calibration profile",
                "screen pipeline",
                "i/o readers",
                "sampling",
                "align",
                "adjudicate",
                "aggregate",
                "decision",
                "report",
                "audit",
            ]
        );
    }

    #[test]
    fn boundary_names_are_unique_across_active_and_legacy_maps() {
        let mut seen = BTreeSet::new();

        for mapping in V2_BOUNDARY_MAP.iter().chain(LEGACY_BOUNDARY_SEAMS.iter()) {
            assert!(
                seen.insert(mapping.boundary),
                "duplicate boundary entry: {}",
                mapping.boundary
            );
        }
    }

    #[test]
    fn legacy_seams_are_marked_pending_removal() {
        assert!(LEGACY_BOUNDARY_SEAMS
            .iter()
            .all(|mapping| mapping.transition == BoundaryTransition::PendingRemoval));
    }
}
