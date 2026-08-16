# FACTOR F0 preregistration

Date: 2026-07-26  
Role: development-only information audit; no predictive model and no confirmation claim.

## Question

Can three public kinase releases support the user's proposed
Functional-Anchor × Chemical-Token Orthogonal Residualization (FACTOR) mechanism without requiring
exact target–ligand overlap across releases?

This is a user-supplied route and does not consume one of the agent's at-most-three candidate slots.
F1 is forbidden unless every F0 gate below passes.

## Sources and firewalls

1. KIRHub 2026 wild-type matrix: one percent-inhibition environment. The existing audited target and
   ligand registries are used; the numerical matrix is not read.
2. Reinecke et al. 2024 Kinobeads apparent-Kd panel, CC BY 4.0: only the processed observation
   registry's target, ligand, scaffold, endpoint and document fields are read; `affinity` and
   `replicate_sd` are not read.
3. The `Christmann2016` provenance slice of the published Kinase200/Papyrus 05.6 release. Papyrus
   05.6 is CC BY-SA 4.0. Rows must have a single source, document and assay, `pchembl_value_N=1`,
   exact relation, and exactly one endpoint flag. The activity columns are not read.

Papyrus-ChEMBL31 is excluded because it can overlap the permanently quarantined ChEMBL confirmation
chain. Sharma2016 is excluded because its eligible slice contains only four targets. Theisen is
excluded because its released table lacks assay/document identifiers. ActFound is excluded because
its assay tasks do not provide a usable protein-condition path. None of these exclusions may be
reversed after seeing F0 results.

The existing ChEMBL confirmation partition is not opened. F0 is development-only even if it passes.
A future confirmation must be a new independent source.

## Label-blind primitive contract

- Protein adapter: the KLIFS-aligned 85-residue canonical pocket. An anchor token is
  `(aligned_position, residue_identity)` or `(aligned_position, residue_physicochemical_class)`.
  Target, family, group, document and source identifiers are prohibited from the token.
- Chemical adapter: non-hashed count-Morgan radius-1/2 local environments, BRICS fragments and
  localized RDKit pharmacophore families. Full-molecule fingerprints are prohibited from the
  interaction path. Bemis–Murcko scaffold is used only for the holdout audit.
- Rare-token stability: tokens appearing in only one entity are retained for novel-primitive
  accounting but cannot connect documents. Tokens present in more than 80% of documents cannot
  create a graph edge.
- For bounded computation, pair-overlap uses a deterministic label-blind subset of eight anchor
  tokens and eight carrier tokens per entity. Selection is by SHA-256 lexical order and is fixed
  before any outcomes.

## Holdout and coverage audit

For every observation, the held axes are its complete KLIFS family, Bemis–Murcko scaffold and source.
An anchor is covered only if it occurs in another family and another source. A carrier is covered
only if it occurs in another scaffold and another source. Token weights are inverse document/entity
frequency and are normalized within an entity.

Samples are separated into:

- `known_combination`: covered primitives and at least one sampled anchor×carrier pair seen in
  another source;
- `new_combination`: covered primitives but no sampled pair seen in another source;
- `novel_primitive`: anchor or carrier coverage below the frozen threshold.

The document primitive graph joins two environments only when they share at least two
non-ubiquitous anchor tokens and two non-ubiquitous carrier tokens. Endpoint-specific environments
remain distinct.

## Power and source weighting

The primary prospective inference unit is the unique
`KLIFS-family × Murcko-scaffold × source` cell. Sources receive equal total weight, and cells within
each source receive equal weight. The Kish effective sample size is computed from these frozen
weights. The normal paired-effect envelope is

`MDE80 = (z_0.975 + z_0.80) / sqrt(n_eff)`.

This is an optimistic label-blind design envelope, not an estimate of biological noise. Environment-
and source-level counts are reported separately and may block interpretation even when the analytic
MDE passes.

## Frozen F0 gates

All must pass:

1. at least 95% of observations map to a valid 85-residue pocket;
2. at least 99% map to a valid carrier set and Murcko scaffold;
3. weighted anchor coverage median at least 0.95 and 10th percentile at least 0.80;
4. weighted carrier coverage median at least 0.90 and 10th percentile at least 0.70;
5. at least 70% of observations are not `novel_primitive`, with new-combination and
   novel-primitive strata reported separately;
6. the largest document-graph component contains at least 80% of eligible environments and all
   three sources;
7. after frozen equal-source weighting, no source contributes over 40%, and every source contains
   at least 10 families and 20 scaffolds;
8. grouped MDE80 is at most 0.03;
9. at least 20 independent endpoint-specific environments and at least three independent sources
   remain.

Pass: `FACTOR_F0_PASS_AUTHORIZE_FIXED_FEATURE_F1`.  
Fail: `FACTOR_F0_FAIL_STOP_BEFORE_F1`.

No threshold, source, token definition, seed or endpoint rule may be changed after the audit.
