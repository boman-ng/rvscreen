# ADR 0002: ambiguity resolution and group detection v4

Status: accepted
Date: 2026-05-07

## Decision

New rvscreen report bundles use `rvscreen.report_bundle.v4` for ambiguity-aware detection. The primary biological detection call is made at `group` level from `CandidateAggregator::with_level(AggregationLevel::VirusGroup)`. Accession rows remain supporting evidence and carry an explicit `accession_resolution` label.

## Rationale

A high `ambiguous_fragments / accepted_fragments` ratio is evidence about resolution and interpretability, not by itself proof that strong viral evidence is biologically absent. Therefore ambiguity is no longer a hard veto for strong group-level positive calls. It remains auditable through decision reasons and accession-resolution labels.

## Contract

- `results/detection_calls.group.tsv` is the canonical group detection artifact.
- `results/candidate_calls.accession.tsv` and legacy `candidate_calls.tsv` are accession evidence artifacts.
- The v4 call header adds `accession_resolution` before `decision`.
- Group rows use `aggregation_level=group` and `accession_resolution=not_applicable`.
- Accession rows use `accession_resolution=high`, `low_resolution`, or `unresolved`.
- Strong accepted-core, high-ambiguity accession evidence is not required to be accession-positive; its corresponding group call can be positive when abundance, absolute support, breadth, non-overlap, enrichment, evidence-strength, and background gates pass.
- Weak-core or poorly supported high-ambiguity evidence remains non-positive.
- Multi-look representative early stop for group positives requires continuity of the same `group`; high-confidence accession positives can still stop immediately.
- `summary/result_overview.json` uses `rvscreen.result_overview.v2` with `top_detections` and `top_accession_evidence`.

## Consequences

This separates organism/group detection from accession or strain resolution. A positive group detection is an algorithmic biological signal, not a clinical causality claim. A low-resolution accession row means the organism/group evidence is stronger than the exact accession assignment.

## Alternatives rejected

- Keep ambiguity as a hard positive veto: over-suppresses strong EBV/HHV-4-like signals.
- Remove ambiguity annotations entirely: loses specificity and explainability evidence.
- Make all strong high-ambiguity accession rows positive: conflates detection with accession resolution.
- Generate group calls only in summary code: report display must not create decision semantics after the decision step.
- Full ambiguity-source decomposition now: scientifically useful but unnecessary complexity for this first-pass rule-based fix.

## Future work

Calibrate group-detection and resolution thresholds on validation panels, and later split ambiguity into host-competitive, viral-multimap, repeat/low-complexity, quality, and background classes.
