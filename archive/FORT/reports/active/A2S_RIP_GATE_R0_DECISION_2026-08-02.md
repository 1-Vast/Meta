# A2S-RIP Gate R0 — decision

Date: 2026-08-02
Artifacts: `reports/active/a2s_rip_gate_r0_2026-08-02.json`,
`reports/active/a2s_rip_gate_r0_records_2026-08-02.parquet`,
`research/a2s_rip_gate.py`, `tests/test_a2s_rip_gate.py`
Preregistration: `A2S_RIP_MECHANISM_DERIVATION_2026-08-02.md` §3.7
Roles opened: `fit`, `probe`. `locked` and the recipient roster never requested.

**Decision: `SELECTION_CEILING_REAL; NOT_REACHABLE_FROM_AN_OBSERVABLE_MARGIN`.
A2S-RIP is retracted as specified. Preregistered triggers P4 and P5 both fired.**

> **Note on the machine verdict string.** The artifact records
> `RIP_CEILING_ADMITTED`. That label was computed by the preregistered
> `decide()` function and refers to **R0a only** — the existence of a selection
> ceiling. It is not the mechanism verdict. The mechanism verdict is governed by
> P4/P5 in §3.7 of the derivation, and both fired. The function is left unchanged
> rather than re-fitted after seeing the result.

---

## 1. Gate R0a — the selection ceiling is real and large. **PASS**

Hindsight-selected subset of the empirical-Bayes head, minus applying the whole head
(probe, within-target scaffold-disjoint, 50 components):

| k | best coverage | ΔCI vs wholesale | 95 % interval |
|---:|---|---:|---|
| 3 | 0.5 | **+0.0621** | [+0.0542, +0.0705] |
| 5 | 0.4 | **+0.0748** | [+0.0626, +0.0886] |

Against the frozen base at k=5: oracle selection reaches **+0.0735 [+0.0645, +0.0834]**
at 40 % coverage, while the wholesale head is −0.0013. **FACT.** Applying a noisy head to
the right 40 % of compounds is worth more than the entire fully-supervised head measured
in G4 (+0.052). The premise of A2S-RIP is correct: the problem at k ≤ 5 is *where* to
apply the head, not whether the head exists.

**R0d(oracle) PASS** — the ceiling is not a magnitude artefact: oracle selection beats a
wholesale head rescaled to the same mean edit magnitude by +0.0577 (k=3) and +0.0665 (k=5),
both with lower bounds above +0.05.

## 2. Gate R0b — the observable margin barely discriminates. **TECHNICAL PASS, MATERIAL FAIL**

AUC for predicting whether a proposed edit points the right way, on unseen components:

| k | statistic | AUC | 95 % interval |
|---:|---|---:|---|
| 3 | evidence margin `z_q` | 0.5536 | [0.5273, 0.5779] |
| 3 | `\|δ_q\|` alone | 0.5507 | [0.5258, 0.5744] |
| 5 | evidence margin `z_q` | 0.5553 | [0.5289, 0.5796] |
| 5 | `\|δ_q\|` alone | 0.5514 | [0.5257, 0.5752] |

**FACT.** The margin clears chance, but by ~5 AUC points. **FACT.** The posterior
covariance — the device that made A2S-RIP a mechanism rather than a heuristic — adds
**0.004 AUC** over the edit magnitude alone. It is **not load-bearing**.

## 3. Gates R0c/R0d — the certification works; there is nothing worth certifying

**R0c PASS (validity).** A harm-rate threshold fitted on `fit` targets transfers to
`probe` essentially unchanged: at α = 0.40, fit harm 0.399 → probe harm 0.387 (k=5),
0.396 → 0.397 (k=3). Fisch-style cross-task threshold transfer *works on this substrate*.
But no coverage met α = 0.20 or 0.30 on any target — the estimator is never that accurate.
And the gain at the certified coverage is **+0.0001 [−0.0155, +0.0138]**.

**R0d FAIL on the implementable rule.** Margin-selected intervention minus a
magnitude-matched wholesale head: **+0.0007 [−0.0010, +0.0023]** at k=3, and at k=5 no
coverage beats the matched control at all.

**The decisive table.** CI gain over the frozen base, probe, k=5:

| coverage | oracle | **margin (implementable)** | random | harm rate (margin) |
|---|---:|---:|---:|---:|
| 0.2 | +0.0606 [+0.0536] | **−0.0002 [−0.0128]** | −0.0012 | 0.373 |
| 0.4 | **+0.0735 [+0.0645]** | **−0.0006 [−0.0170]** | +0.0008 | 0.399 |
| 0.6 | +0.0705 [+0.0614] | −0.0011 [−0.0188] | −0.0035 | 0.418 |
| 1.0 | −0.0013 [−0.0194] | −0.0013 [−0.0194] | −0.0013 | 0.445 |

**FACT.** The implementable rule is indistinguishable from random selection at every
coverage. **P4 fired.** Its residual advantage over wholesale is fully explained by the
reduction in edit magnitude. **P5 fired.**

## 4. Why this cannot be repaired inside the current representation

**INFERENCE.** The gap is not a tuning gap. Reaching the oracle curve requires a
statistic with AUC of order 0.75+; the two available observables give 0.55, and they give
the *same* 0.55, so combining them cannot close a 20-point gap. Selecting a better
conformity score, a learned `ψ_θ`, extra label-free features or a different `α` are all
moves inside a space whose best ingredient is worth 5 AUC points.

The binding difficulty is upstream: **at k ≤ 5 the empirical-Bayes head spreads its mass
over all 26 coordinates, so no per-compound quantity derived from it is sharp.** This is
the same finding as G2/G4 seen from a third angle, and it is a property of the
*representation*, not of the estimator or the decision layer.

**Per the standing instruction: this cannot be resolved within k ≤ 5 steps.** A new branch
is opened (§5), and the retained positive results are carried into it:

1. the selection ceiling (+0.074 CI at 40 % coverage) — the prize;
2. cross-task threshold transfer works (R0c) — the certification layer is sound and
   reusable the moment a sharp statistic exists;
3. the magnitude-matched and random-selection controls — now mandatory in all successors.

## 5. The alternative hypothesis, tested immediately

Before opening the branch, the competing structural hypothesis was tested on the same
probe splits: is the per-target head **sparse in coordinates**, even though G2 showed it is
not **low-rank in the source subspace**?

Truncating each target's own head to its top-`s` coordinates by magnitude, versus
truncating to the top-`s` source principal directions (52 targets, 50 components):

| budget `s` | **coordinate-truncated** | rank-truncated (G2's operation) |
|---:|---:|---:|
| 2 | −0.0083 [−0.0254, +0.0077] | −0.0051 [−0.0159, +0.0054] |
| 3 | −0.0008 [−0.0195, +0.0171] | −0.0019 [−0.0112, +0.0067] |
| 5 | +0.0128 [−0.0047, +0.0303] | −0.0047 [−0.0198, +0.0088] |
| **8** | **+0.0343 [+0.0152, +0.0523]** | +0.0085 [−0.0043, +0.0234] |
| 26 (full) | +0.0542 [+0.0315, +0.0761] | +0.0542 [+0.0315, +0.0761] |

**FACT.** Eight of 26 coordinates retain **63 %** of the full head's gain with a lower
bound well above the MDE, while eight principal directions retain **16 %** with an interval
crossing zero. Coordinate sparsity is materially better structure than low-rankness at a
matched budget.

**FACT.** Across 52 targets the top-weighted coordinate takes **20 distinct values** out of
26. Targets use *different* small coordinate sets.

> **INFERENCE — this reinterprets G2 rather than contradicting it.** A collection of
> sparse vectors with *different supports* has an approximately flat covariance spectrum.
> G2's flat spectrum (top-3 directions = 34.7 % of variance) is therefore **predicted by**
> the sparse hypothesis, not evidence against it. The programme mis-read a signature of
> heterogeneous sparsity as an absence of structure.

**FACT.** Effective sparsity in this generic descriptor basis is `s ≈ 8`, not 2–3.
Standard sparse-recovery arithmetic needs `O(s·log(d/s))` measurements: `8·log(26/8) ≈ 9–10`.
**That is precisely the measured knee of the G4 learning curve (k ≈ 10).** The theory
retrodicts the number.

**HYPOTHESIS carried into the branch.** To move the knee to `k ≈ 5`, do not build a better
estimator — build a **representation in which the per-target head is 2–3 sparse**. That is
a quantitative, falsifiable target, and §6 gives it a biological basis.

## 6. Status of each preregistered prediction

| # | Prediction | Outcome |
|---|---|---|
| P1 | RIP − wholesale LCB > 0.005 at k=3,5 | **FAIL** (−0.0006 at k=5) |
| P2 | RIP ≈ wholesale at k=20–40 | not reached (P1 failed first) |
| P3 | Realised harm ≤ certified α | **PASS** — the one component that worked |
| P4 | Random selection destroys the gain | **FIRED** — margin ≡ random |
| P5 | Magnitude-matched control destroys the gain ⇒ retract | **FIRED** — retracted |
| P6 | Risk–coverage monotone | **PASS** — harm rises monotonically with coverage, 0.373 → 0.445 |
| P7 | Effect of order 0.005–0.02, not 0.05 | correct for the implementable rule (≈0); the *ceiling* is 0.074 |
