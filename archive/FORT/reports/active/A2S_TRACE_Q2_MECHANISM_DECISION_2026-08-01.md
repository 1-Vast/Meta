# A2S-TRACE Q2 decision — learned per-pair transport reliability

Date: 2026-08-01
Artifacts: `reports/active/a2s_trace_q2_mechanism_2026-08-01.json`,
`reports/active/a2s_trace_q2_records_2026-08-01.parquet`,
`reports/active/a2s_trace_headroom_2026-08-01.json`,
`research/a2s_trace.py`, `research/a2s_trace_headroom.py`, `tests/test_a2s_trace.py`
Preregistration: `A2S_TRACE_MECHANISM_EXPLORATION_PROMPT_2026-08-01.md`
Design: `A2S_TRACE_MECHANISM_ANALYSIS_2026-08-01.md`
Roles opened: `fit` (training), `probe` (development measurement). `locked` and the A2S recipient
roster were never requested. Device: RTX 4060 Laptop GPU, `D:\anaconda\envs\drug`.

**Decision: `POSITIVELY_CONTROLLED_NULL_LEARNED_TRANSPORT_NOT_ADMITTED`.**

One positive result survives and it is not the mechanism: a single global transport scale.

---

## 1. What was run

19,611 `fit` training episodes / 4,796 inner-validation episodes (OOF fold 4) / 4,029 `probe`
episodes, two declared support policies, k ∈ {1,3,5}, three model seeds, evaluated in the
Q1-admitted stratum (nearest support Tanimoto ≥ 0.55) with 74–76 homology components. Aggregation
and bootstrap exactly as in Q1. The KRR ridge was selected on `fit` (0.03), so the analytic baseline
is not handicapped.

The ladder is nested in parameterisation **and** in optimisation: a three-epoch warm-up trains only
the analytic scalars, so each learned rung starts from the converged restriction below it.

## 2. Result

**FACT — absolute target-macro CI / NDCG@10, admitted stratum, mean over three seeds.**

| rung | k=3 CI | k=3 NDCG@10 | k=5 CI | k=5 NDCG@10 |
|---|---:|---:|---:|---:|
| frozen base (R0) | 0.5090 | 0.7341 | 0.5114 | 0.7350 |
| level channel (R1) | 0.5090 | 0.7341 | 0.5114 | 0.7350 |
| NW smoother (R2) | 0.5361 | 0.7563 | 0.5414 | 0.7606 |
| CKA-NNLS | 0.5486 | 0.7694 | 0.5658 | 0.7844 |
| fixed Tanimoto KRR (R2b) | 0.5576 | 0.7762 | 0.5730 | 0.7907 |
| **+ global transport scale (R2c)** | **0.5650** | **0.7830** | **0.5830** | **0.8008** |
| static kernel mixture | 0.5649 | 0.7826 | 0.5829 | 0.8004 |
| + learned query gate (R3) | 0.5651 | 0.7829 | 0.5830 | 0.8009 |
| **+ learned pair reliability = TRACE (R4)** | 0.5649 | 0.7828 | 0.5828 | 0.8006 |
| modulation without gate (R5) | 0.5652 | 0.7829 | 0.5827 | 0.8007 |
| TRACE, protein zeroed (R6) | 0.5649 | 0.7828 | 0.5828 | 0.8006 |
| TRACE, scalar pair features only (R7) | 0.5649 | 0.7828 | 0.5828 | 0.8006 |

**FACT — paired component bootstrap, 95 % interval.**

| contrast | k=3 CI | k=5 CI |
|---|---|---|
| fixed KRR − base | +0.0535 [+0.0360, +0.0731] | +0.0680 [+0.0514, +0.0850] |
| **global scale gain (R2c − R2b)** | **+0.0104 [+0.0042, +0.0169]** | **+0.0089 [+0.0027, +0.0149]** |
| **M1: TRACE − R2c (the headline)** | **−0.00015 [−0.00063, +0.00030]** | **−0.00003 [−0.00059, +0.00053]** |
| M1b: same on NDCG@10 | −0.00017 [−0.00063, +0.00020] | +0.00001 [−0.00035, +0.00037] |
| gate gain (R3 − R2c) | −0.00016 [−0.00065, +0.00030] | +0.00004 [−0.00033, +0.00041] |
| low-capacity modulation (R7 − R2c) | −0.00015 [−0.00063, +0.00030] | −0.00003 [−0.00059, +0.00053] |
| TRACE − static mixture | −0.00002 [−0.00055, +0.00051] | +0.00002 [−0.00035, +0.00045] |
| TRACE − CKA-NNLS | +0.0230 [+0.0122, +0.0360] | +0.0169 [+0.0092, +0.0249] |
| **M4: correct − deranged** | +0.0896 [+0.0670, +0.1146] | +0.0907 [+0.0702, +0.1118] |
| **M5: correct − norm-matched wrong target** | +0.0602 [+0.0429, +0.0788] | +0.0831 [+0.0649, +0.1023] |
| **M3: residual-null − base** | +0.00000 [+0.00000, +0.00000] | +0.00000 [+0.00000, +0.00000] |

**Gate verdicts.** M2 (nesting) passes exactly — the R2/R2b restrictions reproduce NW and KRR to
floating point, asserted in `tests/test_a2s_trace.py`. M3 passes bitwise. M4 and M5 pass with wide
margins. **M1 and M1b fail.** M7: protein-zero costs 0.0000 CI, so the protein channel is not
load-bearing and nothing here may be called protein-conditioned DTA.

## 3. The null has power

**FACT — synthetic positive control.** The identical learner, the identical episodes, the identical
probe components, in a world where transport reliability really is pair-dependent (support residuals
transported only between molecules close on one held-out descriptor):

| k | base | fixed KRR = R2c | TRACE | oracle | TRACE − R2c | recovered |
|---:|---:|---:|---:|---:|---|---:|
| 3 | 0.7966 | 0.8227 | 0.8501 | 0.8846 | **+0.0262 [+0.0154, +0.0370]** | 42 % of the oracle gap |
| 5 | 0.8034 | 0.8247 | 0.8407 | 0.8830 | **+0.0159 [+0.0078, +0.0238]** | 27 % of the oracle gap |

**INFERENCE.** The pipeline detects a pair-reliability effect of +0.016 to +0.026 CI on unseen
components. On real data it measures −0.0001 with a 95 % upper bound of +0.0005. The null therefore
excludes an effect roughly **30× smaller** than the one the control recovers. This is a measurement,
not a failure to look.

## 4. Where the remaining headroom is, and why it is unreachable here

**FACT — hindsight oracles on the same probe episodes** (`a2s_trace_headroom_2026-08-01.json`;
every oracle reads the labels it is scored on, so these are ceilings, not methods):

| oracle (k=5, `random_within_target`) | absolute CI | headroom over fixed KRR |
|---|---:|---|
| fixed Tanimoto KRR | 0.5857 | — |
| best **episode-level scale** in hindsight | 0.6636 | +0.078 [+0.068, +0.088] |
| best **episode-level support subset** in hindsight | 0.6453 | +0.059 [+0.052, +0.067] |
| best **per-query support subset** in hindsight | 0.9487 | +0.353 [+0.327, +0.380] |

**INFERENCE — the action class is not the constraint.** A per-query binary reliability rule — exactly
TRACE's action class, restricted to on/off — can reach CI 0.93–0.95 with hindsight. The mechanism can
express near-perfect rankings. What is missing is any **label-free predictor of which choice is
right**. Combined with §3, the conclusion is specific:

> **In provenance-quarantined ChEMBL-37 pKi, the reliability of residual transport between a support
> compound and a near-analogue query is not predictable from label-free ligand chemistry or pooled
> protein sequence at the resolution these features provide.** It is not that the class is too weak,
> and not that the learner is too weak. Both were measured.

The mean hindsight-optimal episode scale is 2.2–2.5, while the single meta-learned global scale
converges to ≈1.5. The gap between them (+0.078 hindsight vs +0.009 achieved) is the size of the
prize for any future *episode-level* magnitude router — and TAMSK-family routers are exactly that
class. This is the most concrete remaining target the programme has.

## 5. The one positive, stated at its true size

**FACT.** A single meta-learned, target-independent scalar multiplying the transport improves the
admitted-stratum ranking over standard fixed Tanimoto KRR by +0.0104 [+0.0042, +0.0169] CI at k=3 and
+0.0089 [+0.0027, +0.0149] at k=5, with NDCG@10 moving in the same direction and RMSE not degraded.

**INFERENCE — what it means and what it is not.** The frozen base carries roughly twice the
within-episode prediction spread of the transport (SD 0.96 vs 0.46 at k=5) while ordering at chance,
so the base is over-weighted relative to the support evidence. The scale is a precision re-weighting
between two channels. It is *not* rank-null (it is a scale on a query-varying quantity, not a shift),
but it is also **not an adaptation mechanism**: it has zero target-specific parameters *and* zero
query dependence beyond what the kernel already supplies. It belongs in the baseline of every future
comparison in this programme, and any method that reports a gain over unscaled KRR without granting
the baseline this scalar is reporting a hyperparameter, not a mechanism.

**FACT worth keeping.** The convex ranking surrogate (RankNet logistic) will not find this scale: its
optimum sits far below the CI optimum because it keeps paying for confidently-wrong pairs. Only the
bounded smoothed-CI surrogate finds it. Any future ranking work on this substrate should use the
bounded surrogate.

## 6. Admissibility audit of the delivered mechanism

| condition | status |
|---|---|
| 1. Learned across source episodes | yes — amortised, no recipient-specific fitting |
| 2. Identifiable at k ≤ 5 | yes — **zero** target-specific quantities; C3 does not bind |
| 3. Query-dependent | yes — `m_qi` and `α_q` vary per pair and per query |
| 4. Structurally abstaining | yes — `r_S ≡ 0 ⇒ Δ ≡ 0` bitwise, proven and tested |
| 5. Bounded | yes — `\|Δ_q\| ≤ max_i \|r_i\|`, tested |
| 6. Nested-falsifiable | yes — R2b reproduces fixed KRR and R2 reproduces NW to floating point |
| 7. Not shortcut-driven | yes — weights never read a label; derangement and wrong-target contrasts are structural, and both are large |

All seven hold. The mechanism is admissible and it is null. That is the cleanest form this result
could take.

## 7. What must not be claimed

- Not "TRACE improves DTA." It ties the analytic bar to within ±0.0006 CI.
- Not "learned kernels do not work." The claim is bounded to this stratum, these features, this base
  and 74–76 components.
- Not "protein conditioning fails." Pooled ESM-2 through a 32-dim projection contributed nothing;
  residue-level or pocket-level conditioning was not tested here.
- Not a recipient result. `probe` is a development role and was inspected repeatedly during
  development. Confirmation requires freezing this protocol and opening `locked` once.

## 8. Recommended next action

The evidence points at one target, not at a seventh architecture:

1. **Episode-level magnitude, not pair-level reliability.** The hindsight ceiling for an episode-level
   scale is +0.078 CI against +0.009 achieved by a single global scalar. A *predictable* per-episode
   scale — from support-only, label-free quantities such as support spread, Gram conditioning and
   support-query similarity mass — is the only measured gap with real room in it. It is also exactly
   the TAMSK claim, now with a measured ceiling and a bar that already includes the global scale.
2. **Do not open `locked` for TRACE.** There is nothing to confirm.
3. Keep the Q1 stratum map as the programme's reusable instrument. Every future support-adaptation
   result must report its policy and its relation stratum.
