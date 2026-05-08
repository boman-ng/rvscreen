# ADR 0001: report schema v3 correctness contract

Status: accepted
Date: 2026-05-07

## Decision

New rvscreen report bundles use `rvscreen.report_bundle.v3`. The v3 contract records corrected non-overlap support, coverage-block counts, explicit denominator provenance, accession-first multiplicity, no default negative-control background decision path, and structured release-gate blockers.

## Contract

- `nonoverlap_fragments` is computed with `max_nonoverlapping_intervals_v1`: per-contig half-open intervals, zero-length intervals excluded, sorted by `(end, start)`, greedy accept when `start >= last_accepted_end`.
- `coverage_interval_blocks` preserves the prior merged coverage-block count.
- The default candidate decision does not use negative-control background ratios. Missing controls are metadata only.
- `multiplicity_family` follows the actual aggregation level; production remains accession-first for this patch.
- Fractions use `total_sampled_fragments`; `unique_fraction` remains a legacy alias, not molecule-unique evidence.
- `audit/release_gate.json` persists `status`, `primary_blockers`, and `secondary_blockers` and audit verification recomputes them.

## Note on schema versioning

The v3 contract described above was never released independently of v4. ADR
0002 was accepted in the same change set, so report bundles produced by this
repository jump directly from v2 to v4. The v3 corrections (non-overlap
intervals, structured release-gate blockers, denominator provenance, and
removal of the negative-control background ratio from the default decision
path) are still part of the v4 contract.
