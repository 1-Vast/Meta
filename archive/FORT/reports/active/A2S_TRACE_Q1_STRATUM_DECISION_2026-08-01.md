# A2S-TRACE Q1 decision — the support-information admission boundary

Date: 2026-08-01
Artifacts: `reports/active/a2s_trace_q1_stratum_2026-08-01.json`,
`reports/active/a2s_trace_q1_records_2026-08-01.parquet`,
`research/a2s_trace_stratum.py`, `tests/test_a2s_trace_stratum.py`
Roles opened: `fit`, `probe`. `locked` source role and the A2S recipient roster were never requested.
Device: NVIDIA GeForce RTX 4060 Laptop GPU (`D:\anaconda\envs\drug`).

**Decision: `INFORMATION_ADMITTED_IN_A_LOCAL_RELATION_STRATUM`.**
The programme's blocking contradiction C9 is resolved. It was not a disagreement between corpora.

---

## 1. What was measured

One corpus (ChEMBL-37 dualcold pKi TRAIN under the balanced v2 lock: 222 / 110 / 107 fit / probe /
locked homology components, zero target, homology, document and assay overlap across roles), one
frozen component-cross-fitted ridge base, one fixed Tanimoto KRR (ridge 0.1), one true residual
derangement, one norm-matched wrong-target arm, one paired component bootstrap (2,000 draws).

Two axes were varied and nothing else:

| Axis | Levels |
|---|---|
| **Support policy** | `random_within_target` (BindingDB-style), `scaffold_disjoint` (random support, queries scaffold- and connectivity-cold to it), `provenance_disjoint` (the frozen v2 policy: support and query disjoint on scaffold, connectivity, document and assay) |
| **Support→query relation stratum** | nearest support Tanimoto of the query, binned `<0.20`, `0.20–0.35`, `0.35–0.55`, `≥0.55`, plus the pooled `all` |

12,246 probe episodes over 3 episode seeds; aggregation fixed as
episode draws → seed/target mean → component mean → paired component bootstrap.

Admission rule, preregistered in the code: a cell is admitted when the KRR-minus-base CI lower bound
exceeds the 0.005 MDE **and** (for k ≥ 3) the correct-minus-deranged CI lower bound exceeds zero.

## 2. Result

**FACT — the gain lives in one stratum, in every policy.** KRR minus frozen base, target-macro CI,
paired component 95 % lower bound:

| policy | k | `<0.20` | `0.20–0.35` | `0.35–0.55` | `≥0.55` | pooled |
|---|---:|---:|---:|---:|---:|---:|
| `random_within_target` | 1 | −0.001 | −0.003 | −0.005 | −0.014 | +0.007 |
| `random_within_target` | 3 | −0.003 | −0.003 | −0.003 | **+0.036** | +0.023 |
| `random_within_target` | 5 | −0.007 | −0.006 | +0.006 | **+0.048** | +0.033 |
| `scaffold_disjoint` | 3 | −0.003 | −0.003 | −0.000 | **+0.023** | +0.015 |
| `scaffold_disjoint` | 5 | −0.000 | +0.001 | +0.012 | **+0.031** | +0.026 |
| `provenance_disjoint` | 3 | −0.012 | −0.006 | +0.004 | −0.020 | +0.001 |
| `provenance_disjoint` | 5 | −0.006 | −0.013 | −0.006 | −0.004 | +0.009 |

At `t ≥ 0.55`, k=5, `random_within_target`: absolute target-macro CI rises 0.5167 → 0.5801,
NDCG@10 0.7247 → 0.7824 (gain +0.054 [+0.041, +0.069]), and RMSE improves with a lower bound of
+0.496 pKi. The direction is identical on every metric.

**FACT — the controls agree, cell by cell.** In every admitted cell the correct-minus-deranged CI
lower bound is positive (+0.021 to +0.062) and the correct-minus-**norm-matched-wrong-target** lower
bound is positive (+0.019 to +0.057). In every null cell both are also null or negative. The
information and the assignment specificity switch on and off together.

**FACT — the episode-constant channel is exactly worth zero.** The level (mean-residual) estimator
returned a CI delta of exactly 0.0000 with a degenerate bootstrap in all 45 policy × k × stratum
cells. C1 is now a verified property of this harness rather than an assumption.

**FACT — the support-Gram inverse is load-bearing.** Nadaraya–Watson recovers about half of the KRR
gain (`t ≥ 0.55`, k=5, random policy: +0.024 vs +0.048 LCB).

## 3. Interpretation

**INFERENCE — C9 is resolved and neither prior result is overturned.** The `provenance_disjoint`
policy draws queries at mean nearest-Tanimoto 0.19–0.30, i.e. almost entirely inside the two measured
null bins. Its global null was a correct measurement of a stratum that contains no usable support
information. The BindingDB positive was a correct measurement of a stratum that does. The boundary is
the **support→query chemical relation**, not the dataset, the vendor, or the provenance closure.

**INFERENCE — this is local SAR transport, and must be called that.** The admitted stratum is
explicitly *not* scaffold-cold with respect to the support. A gain measured there is a claim about
transporting a measured residual to a near analogue on the same target, under a declared support
policy. It is not a claim about global dual-cold DTA, and the null strata must be reported beside it
every time.

**INFERENCE — k=1 behaves as theory predicts.** Every stratum-conditional k=1 cell is null; the
pooled k=1 cell is weakly positive because a single support compound still yields a query-dependent
weight `T(q,s)`, so ordering can move through between-query similarity variation alone.

## 4. Consequences that are now binding

1. Any mechanism phase runs **inside the admitted stratum** and is judged against **fixed Tanimoto
   KRR at equal support information**, not against the frozen base. Beating the frozen base there is
   free.
2. The `provenance_disjoint` construction is retired as a mechanism-training substrate. It is
   retained as the off-stratum non-inferiority control.
3. No number from this gate may be quoted without its policy and its stratum.
4. `probe` is a development role. Confirmation requires freezing the protocol and opening `locked`
   once. Recipient labels stay sealed.

## 5. What this gate does not establish

- It does not show that any *learned* object beats the fixed analytic smoother — that is Q2.
- It does not show that protein information is load-bearing; the estimator is ligand-only.
- It does not establish a practical drug-discovery effect size; no utility threshold was preregistered
  for the ranking gain, only a statistical MDE.
- The retained rows exclude provenance-crossing rows by design, so the estimand is the
  provenance-quarantined source pool, not all of ChEMBL-37 pKi.
