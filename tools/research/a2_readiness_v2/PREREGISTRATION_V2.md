# Preregistration v2 — **SUPERSEDED for Stage P** (2026-08-16)

> Stage P's authoritative preregistration is
> `tools/research/stageP_cpc/PREREGISTRATION.md`. It resolves a contradiction
> in this document: §1 described Stage P as one over-driven arm while §7 listed
> five (`A0frozen`, `A0repro`, `CPCpos`, `CPC`, `CPCwrong`, `CPCrand`,
> `A3perm`). The prerequisite is **two newly trained arms** — `A0repro` and
> `CPCoverdrive` — at three matched seeds. The remaining arms belong to the
> later admission stage and must not be counted in the prerequisite's cost or
> gates.
>
> §§2-8 below remain the reference for the *admission* stage. §1 is superseded.

Frozen before any candidate is trained or observed. Supersedes
v1's preregistration draft (consolidated into `tools/research/a2_readiness/SUPERSEDED.md`).

This document does **not** authorize CPC. Phases 1-3 removed CPC's measured
premise (E3 withdrawn) and relocated its target (E4 → F10), and no evidence
now predicts that it helps. What it authorizes is one cheap, decisive
prerequisite — Stage P — whose outcome determines whether an admission-grade
CPC arm is worth running at all. §§2-6 specify that arm completely, so that if
Stage P passes, nothing is designed after seeing a result.

## 0. What must not be claimed, whatever happens

Fixed now, so that a positive result cannot quietly widen its own scope:

* **no pocket, contact, binding-site or "biologically localized" claim.** The
  protein path is exactly invariant to residue-slot order (F9, measured 2.4e-08
  pK). The attention reads an unordered bag of sequence-window summaries. There
  is no common-frame complex geometry for any of 17,717 DTA cells.
* **no SOTA or excellence claim.** Every number here is `meta_val` development
  evidence on the population that selects the configuration.
* **no "protein specificity" claim from an uncentered control.** Uncentered
  wrong-protein gaps measure the level (F7). Only centered contrasts speak to
  ordering.
* **no meta-learning claim from ligand-side transfer.** The one positive result
  of this cycle (`embed`'s protein-independent SAR direction, Δ`r` +0.2623) is
  ligand-side and must be reported as such.
* **no attribution of a gain to the architecture** when ligand overlap,
  external data, retrieval, calibration shift or evaluation selection could
  produce it.
* `meta_test` stays sealed. Opening it requires a written authorization
  recorded in the artifact (`QPSMPData.seal_record()`), and only for a
  candidate that has already passed every gate below.

## 1. Stage P — the prerequisite (this is what is authorized)

**Question.** Is protein-conditioned within-target ordering learnable *at all*
on this data, when the trunk is free and the objective explicitly demands it?

Phases 1-3 measured a trunk that was never asked for it. Two explanations
predict everything observed equally well: (a) the objective never asked, or
(b) 346 training targets at 9-21 ligands each, on sequence + 2D inputs, do not
contain enough within-target signal for any model to learn it. Stage P
separates them, and it is the cheapest experiment that can.

**Design.** One arm, matched to A0 in every respect (architecture, 1200 steps,
seeds, bank, budget) except the objective, which is replaced by the *most
favourable possible* version of the CPC hypothesis:

* `protein_contrast_loss_weight` raised to **2.0** (A0: 0.5) and computed on
  the **centered** prediction;
* fired on **every** episode including k=0 (see §3);
* the regression terms retained at A0's weights so the arm remains a DTA model.

This is deliberately an over-driven arm, not an admission candidate. Its only
job is to answer whether the signal exists. If an objective weighted 4× and
fired on every step cannot produce protein-conditioned ordering, a weight-0.5
version will not.

**Measurement.** The Phase 1 instrument, unchanged: the centered
wrong-protein contrast `r(interaction) − r(interaction | wrong P)` across all
five frozen donor strata, three seeds, component-paired bootstrap.

**Gates.**

| gate | requirement |
|---|---|
| **P1** (primary) | the centered contrast reaches **≥ +0.05** (the preregistered smallest effect of interest) with a positive component-paired lower bound, at the `nearest` stratum |
| **P2** | the effect does not come from donor destruction: correct-protein `r(interaction)` must not fall below A0's 0.2206 − 0.05 (§4) |
| **P3** | the protein-driven shift is **aligned** with truth: corr(shift, centered truth) ≥ +0.10, against A0's −0.014 and randinit's +0.033 |
| **P4** | the shift is **reproducible across seeds**: mean pairwise cosine of the shift vectors ≥ +0.30, against randinit's −0.003 |
| **P5** | k=0 MSE does not regress by more than 0.10 pK against A0's 2.1488 |

**Stopping rule.** If P1 fails, **the protein-conditioned-ordering family is
closed** and the measured boundary is recorded: on this protocol, this feature
set and this budget, within-target ordering is not protein-conditionable by an
objective that demands it. §§2-6 are then not run, and the next cycle moves to
the calibration lane (M0) or to the ligand-side direction in §7.

If P1 passes but P3 or P4 fails, the objective is producing arbitrary
protein-dependent movement — the randinit failure mode — and the family is
also closed, with that specific cause recorded.

**Cost.** One configuration × 3 seeds × 1200 steps, matched to A0's budget.
No new parameters. One extra forward per episode over A0 on the 240 k=0 steps.

## 2. Stage A3 — the admission-grade CPC arm (specified, not authorized)

Runs **only** if P1-P5 all pass. Nothing below may be changed after Stage P is
observed.

**Mechanism.** One training innovation, zero model change:

```text
p  = zero_shot(P_correct, L_q)          q = zero_shot(P_donor, L_q)
p̃  = p − mean_q(p)                      q̃ = q − mean_q(q)      ỹ = y − mean_q(y)
L_cpc = softplus( ( ‖p̃ − ỹ‖²/Q − ‖q̃ − ỹ‖²/Q ) / T )
```

Centering removes `protein_value(P)` exactly, because it is constant across the
queries of one target. Therefore `∂L_cpc/∂(protein_head) ≡ 0` and the level
branch cannot satisfy the term — verified on the real `SimilarityGrammarModel`
by 11 structural probes (`a2_readiness/tests/`), including the contrast showing
A0's uncentered form *does* reach `protein_head` (> 1e-9).

Weight **0.5**, replacing the existing `protein_contrast_loss_weight`. No other
loss term changes. Total added trainable parameters: **0**.

## 3. Episode coverage, warmup, update count and cost — resolved

The incumbent's protein contrast is inside `if adapt and support_size > 0`
(`train_qpsmp.py:774`), so despite supervising `full.zero_shot` — a k=0
quantity — it never fires on a k=0 episode. `support_size` cycles over
`(0,1,2,3,5)`, so it runs on **960 of 1200 steps**. Warmup does not gate it in
A0: `representation_warmup_fraction = 0.0` gives `warmup_steps = 0`, and the
loop starts at `step = 1`, so `phase_a` is never true.

**Decision: CPC fires on every episode, including k=0.** The term is a
zero-shot objective; gating it on the presence of support is an incidental
coupling, not a design. This makes the k=0 endpoint — the quantity the whole
cycle is about — supervised on every step rather than four fifths of them.

| quantity | A0 | Stage P / A3 |
|---|---:|---:|
| steps with the protein contrast | 960 / 1200 | **1200 / 1200** |
| effective updates carrying the term | 2,880 episodes | **3,600 episodes** (+25%) |
| extra wrong-protein forwards | 2,880 | 3,600 (+720, +25%) |
| added parameters | — | **0** |
| added peak memory | — | none (the extra forward is `torch.no_grad`-free but sequential, and A0 already performs it on 4/5 steps) |
| wall-time estimate | 1.00× | ~1.05-1.08× |

Warmup is set explicitly to `representation_warmup_fraction = 0.0` in both
stages, matching A0, so the term is never silently disabled.

## 4. The donor-destruction shortcut — a hard gate, not a statistic

`binding_contrastive_loss([correct, wrong], T)` is `softplus((correct − wrong)/T)`,
minimised by making `correct` small **or** `wrong` large (F8). A contrastive
protein objective can therefore "succeed" entirely by ruining the donor
prediction, which is not protein specificity — it is a learned aversion.

Every arm reports both sides separately:

| reported per arm | correct protein | wrong protein |
|---|---|---|
| centered MSE | ✅ | ✅ |
| within-target `r` | ✅ | ✅ |
| CI | ✅ | ✅ |
| Spearman | ✅ | ✅ |

and the gap is decomposed:

```text
gap = [ r_correct − r_correct(A0) ]  −  [ r_wrong − r_wrong(A0) ]
          \_______ improvement _______/    \_______ degradation _______/
```

**Hard gate G-SHORTCUT.** The arm fails if `r_correct` falls below A0's
`r_correct` by more than 0.05, regardless of how large the gap is. An
improvement in the gap that is more than half attributable to the degradation
term is recorded as a **failure**, not as a partial success.

## 5. A3perm — defined exactly

The label-permutation control, specified to the level where two people would
implement it identically.

For each episode, within the recipient target only:

1. permute the **query label vector** `y` by a fixed per-episode permutation
   drawn from `np.random.default_rng(stable_seed("A3perm", seed, target, draw))`;
2. **preserve** exactly: the target identity, the protein inputs, the query
   ligand multiset and its order, the support set and its labels, the query
   panel size, and the label **mean and standard deviation** within the target
   (a permutation preserves both by construction — no renormalisation is
   applied or needed);
3. **do not** permute across targets, across components, or between support and
   query;
4. reject the permutation and redraw if it is the identity, for panels of
   size ≥ 2.

Prediction: because the labels no longer correspond to the ligands, the
centered CPC term is being asked to order toward noise. `r_correct` must
collapse to the Phase 1 shuffled-label baseline (A0: −0.008), and the centered
protein contrast must not exceed +0.02. An arm that still shows a protein
contrast under A3perm is producing protein-dependent output that is independent
of the labels — the arbitrary-movement failure mode — and **fails**.

## 6. Donor matching — a frozen distance rule, so magnitude is not a confound

Both the wrong-protein and the random-protein controls draw from the **same**
frozen rule (`_donors.py::stratified_donors`), fixed before Stage P:

* candidates: `meta_val` targets in a **different** homology component (so both
  proteins are equally unseen, and the substitution is legal);
* geometry: pooled ESM embeddings whitened on **`meta_train` only**;
* strata: whitened-cosine quantiles `nearest / q25 / median / q75 / farthest`.

The `random-distance-matched` control is the `median` stratum, not a uniformly
drawn protein: a uniform draw would differ from the nearest donor in *distance*
as well as in identity, confounding the two. Reporting all five strata makes
the distance dependence itself a measurement — as in Phase 1, where the level
shift rose 0.215 → 0.342 pK across the strata while the ordering shift did not
move.

## 7. Arms

Every arm matched in budget, seeds, bank and evaluation.

| arm | what it is |
|---|---|
| `A0frozen` | the retained incumbent checkpoints |
| `A0repro` | the incumbent configuration **retrained here** — the matched control that absorbs retraining noise |
| `CPCpos` | positive-only centered shape control at matched loss weight: the centered term on the correct protein alone, no counterfactual. Separates "centering the objective helps" from "the *contrast* helps" |
| `CPC` | the candidate |
| `CPCwrong` | matched wrong-protein control (`nearest` stratum) |
| `CPCrand` | random distance-matched control (`median` stratum) |
| `A3perm` | ligand/label permutation control (§5) |
| `A2int` | **not run.** Phase 2 rejected the internal-representation A2-min on four representations; the conditional that would have authorized it is not met |

## 8. Statistics and stopping

* three fixed seeds (20260815/16/17), matched budget, nested k=0/1/2/3/5;
* component-level paired bootstrap over the 19 `meta_val` components, 9,999
  draws, seed 20260816; seeds averaged **within target** before resampling;
* a contrast is RESOLVED only if the interval excludes zero **and** the effect
  reaches the smallest effect of interest (0.05 in `r`); otherwise
  RESOLVED_NEGLIGIBLE, DECISIVE_NULL or UNDERPOWERED (`_frozen.py::verdict`);
* **the power problem is stated in advance**: same-configuration retraining
  moves aggregate k=0 `r` by 0.051 (R14 screening, A0frozen vs A0repro). Any
  claimed improvement below ~0.10 is inside two retraining draws and will be
  reported as unresolved regardless of its interval;
* a hard gate failure stops the family immediately, the negative result is
  preserved, and the next evidence-supported hypothesis is taken up.

## 9. Admission (unchanged from the standing contract)

Opening `meta_test` requires, in addition to every gate above: aligned
improvement in MSE **and** CI **and** Spearman, no calibration regression,
correct support-label dependence at k≥1, correct-protein dependence measured on
the **centered** contrast, consistency across all three seeds, and a written
authorization recorded in the artifact. Smoke tests are correctness checks and
are never performance evidence.
