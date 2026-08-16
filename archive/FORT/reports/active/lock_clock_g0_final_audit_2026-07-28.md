# LOCK/CLOCK G0 final adversarial audit

**Date:** 2026-07-28

## Accepted artifacts

| artifact | SHA-256 |
| --- | --- |
| `reports/active/lock_clock_g0_label_free.json` | `3548f0a165168fed576b008a71c8473c772f6f5a34722a62bd70f7b7483e8df8` |
| `manifests/lock_clock_g0_reordering_preparation.json` | `9ffc54b6c70a1896204305f5ef4f9cf2bb9ef438b5cb719eec67c99712a4802e` |
| `reports/active/lock_clock_g0.json` | `555522b0956834b807972b002b7de79b01dba100a9680a7d4d449e6cb7f022e0` |
| `reports/active/lock_clock_g0_amendment_a1.md` | `de08390c8718d76393f3feb3942f834ac2813dd4665102815d3f19a851527aae` |

The accepted implementation bundle is
`39b49405a0d4fac3d63faa66011ec8866c9de4266575f68e140bcda948ace835`.
It binds the runner and the three eligibility/profile helpers used by the audit.

## A1 integrity result

The original formal sequence was invalidated because its artifact comparator accepted discrete type
changes such as boolean `True` versus integer `1`. A1 requires exact kind agreement, exact boolean and
integer equality, and permits the frozen `1e-12` tolerance only between two floating-point values.

The invalidated artifacts remain audit-only:

| artifact | SHA-256 |
| --- | --- |
| `reports/active/lock_clock_g0_invalidated_pre_a1_label_free.json` | `9367fa1c257fe0bd0960de912907aa25393d61ec034054e0b5c7e8c70f396372` |
| `manifests/lock_clock_g0_invalidated_pre_a1_reordering_preparation.json` | `bd4d467be5ba05479fafc43cb01c0a3246ea577a6ee30c628e89f84eb7b65e43` |
| `reports/active/lock_clock_g0_invalidated_pre_a1.json` | `73902075627c9562ac9d036a978fb18d5575e9aff577b2f9a0cd67d8ccfca44b` |

An independent structured comparison removed only each top-level `source` field. All three accepted
artifacts then matched their invalidated counterparts exactly under the A1 comparator. The rerun changed
integrity enforcement, not the scientific result.

## Gate recomputation

The four G0-L gates were independently recomputed from the accepted metrics and exactly matched the
recorded booleans:

| gate | accepted metric | threshold | result |
| --- | ---: | ---: | --- |
| normalized LOCK minimum eigenvalue | `2.07e-16` | `>= -1e-8` | pass |
| family/composition residual energy | `0.3851` | `>= 0.05` | pass |
| CKA with pooled ESM-2 | `0.1633` | `<= 0.95` | pass |
| within-family non-constant pair fraction | `0.9992` | `>= 0.80` | pass |

The separate low-dimensional claim failed: effective rank was `289.36`, and top-16 centered energy was
`0.2340` against a frozen `0.80` requirement. No larger rank was authorized.

G0-R returned `LOCK_G0_REORDERING_NOT_IDENTIFIED_STOP`. The primary contrast missed the frozen
materiality threshold, LOCK was inferior to exact aligned identity, and it failed to beat both
BLOSUM-label permutation and parameter-matched random PSD geometry. These are mechanism failures, not
mere power failures. The positive position-shuffle and wrong-target contrasts establish that alignment
and target identity matter; they do not identify BLOSUM substitution semantics.

## Claim boundary

1. LOCK/CLOCK effectiveness in local mutation landscapes is external literature evidence.
2. Substitution geometry is not identified in the current DTA graph.
3. Prediction of target-specific ligand reordering by that geometry is not established.
4. Strict dual-cold performance improvement is untested.

The conservation-weighted arm is a label-free positional prior, not structure-conditioned CLOCK.
KirHub is ligand-warm, kinase-only, and single-source. No GP, operator, scalar-gated correction,
matrix/tensor completion model, Transformer, Mamba, or other affinity predictor was trained.

## Final disposition

Only the prospective A0 reliability panel may proceed: 12 targets from at least 6 families, at least
16 shared scaffold-diverse ligands, two operationally independent sites, one pKi or pKd endpoint,
complete randomized factorial inclusion, retained inactive/censored outcomes, and all target, ligand,
assay, document, and provenance firewalls. It estimates reliability and variance; it is not powered for
a `0.03` predictive gain.

Final category: **3 - current data cannot identify the mechanism; new prospective measurement
conditions are required.**

`sealed_test_consumed=false`.

## Verification

An independent data/firewall agent recomputed all accepted hashes, the implementation bundle, the four
G0-L gates, the failed low-dimensional claim, the exact `0.0000257008` primary-threshold shortfall, and
the aligned-identity, BLOSUM-permutation, and random-PSD failures. It found no remaining blocker. The
three accepted artifacts were also independently confirmed to be type-strictly identical to the
invalidated pre-A1 artifacts after removing only their top-level `source` fields.

Using `D:\anaconda\envs\drug\python.exe`, the focused LOCK/CLOCK and pocket-oracle suite returned
`16 passed`, and the full repository suite returned `365 passed`. `git diff --check` reported no
whitespace error; it emitted only pre-existing CRLF-to-LF warnings for modified tracked files.
