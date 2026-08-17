# Stage P preregistration — authoritative (2026-08-16)

Frozen before any arm is trained or observed. **This document supersedes every
other statement of Stage P's arms, cost or gates**, including
`tools/research/a2_readiness_v2/PREREGISTRATION_V2.md` §1 and any cost figure
in `task.md` or `history.md`. Where they conflict, this file governs.

## Question

Can the current `ContactGrammar` architecture learn useful **protein-conditioned
within-target ordering** when an objective explicitly demands it?

No objective in R0-R14 ever asked. Every wrong-protein control in the project is
computed on *uncentered* error, which the additive `protein_value(P)` branch
satisfies on its own — it is constant across a target's queries, so a 0.215 pK
level shift separates correct from donor and the gradient is extinguished by
`protein_head` before reaching the ligand-varying path (DATAFLOW_AUDIT F6/F7).

This stage decides between two explanations that fit every measurement made so
far equally well:

* **(a)** the objective never asked, or
* **(b)** 346 training targets at 9-21 ligands each, on sequence + 2D inputs,
  do not contain enough within-target signal for any model to learn it.

## Arms — exactly two, three matched seeds

| arm | what it is |
|---|---|
| `A0repro` | the incumbent configuration, **retrained here**. Not the frozen R3R4 checkpoints: those are a diagnostic reference and cannot absorb this run's retraining variance |
| `CPCoverdrive` | identical model, seeds, data, budget and regression losses; the protein contrast computed on the **centered** prediction, fired on **every** episode, at weight **2.0** |

Seeds: 20260815, 20260816, 20260817. Six runs total.

`CPCpos`, `CPCwrong`, `CPCrand` and `A3perm` belong to the **later admission
stage** and are not trained, counted or reported here. Including them in a
prerequisite that only asks whether the signal is learnable at all would
conflate the question with its follow-up.

## The single training change

```text
p = zero_shot(P_correct, L_q)     q = zero_shot(P_donor, L_q)
p̃ = p − mean_q(p)   q̃ = q − mean_q(q)   ỹ = y − mean_q(y)
L_cpc = softplus( ( ‖p̃−ỹ‖²/Q − ‖q̃−ỹ‖²/Q ) / T ),   T = 0.1
```

Centering removes `protein_value(P)` exactly, so `∂L_cpc/∂(protein_head) ≡ 0`
at every parameter value — verified on the real `SimilarityGrammarModel` by
`tests/test_centered_contrast.py` (13 probes), alongside the contrast showing
the incumbent's uncentered form *is* satisfied by a level shift.

Implementation: `scripts/train_qpsmp.py::centered_protein_contrast`, selected by
`--protein-contrast-form centered`. The default remains `uncentered`, so every
recorded arm is bit-unchanged by the flag's existence.

**Why over-driven.** Weight 2.0 against the incumbent's 0.5, on every episode
rather than 4 of 5. This is not an admission candidate — it is the most
favourable version of the hypothesis. If a 4×-weighted objective firing on every
step cannot produce protein-conditioned ordering, a weight-0.5 version will not.

## Gradient routes

| route | expected |
|---|---|
| `protein_head` (the level branch) | **identically zero** from `L_cpc`; nonzero only from the regression terms |
| `interaction_head`, `contact_weight` | nonzero from `L_cpc` |
| `grammar` (attention), `protein_encoder` | nonzero from `L_cpc` via the donor forward |
| `ligand_encoder` | nonzero, but shared with the regression terms — reported to detect the ligand branch being starved |

Per-branch gradient norms are recorded at fixed steps for both arms.

## Donor construction

* training donors come from **`meta_train` only**, from a **different homology
  component**. Verified by construction: `QPSMPData.draw_episode(split, …)`
  builds its donor pool from `self.components[split]` and excludes the recipient's
  own component, so a `meta_train` episode cannot receive a `meta_val` donor;
* no `meta_val` protein or embedding enters training donor selection;
* evaluation donors are the frozen `meta_val` stratified rule
  (`_donors.stratified_donors`, whitened on `meta_train` only), used **only**
  for frozen evaluation.

## Metrics, reported separately for both protein conditions

For each arm, seed and k ∈ {0,1,2,3,5}:

1. correct-protein centered MSE, within-target `r`, CI, Spearman;
2. wrong-protein centered MSE, `r`, CI, Spearman;
3. **improvement term** `r_correct(CPC) − r_correct(A0repro)`;
4. **donor-degradation term** `r_wrong(A0repro) − r_wrong(CPC)`;
5. the fraction of the correct-minus-wrong gap contributed by each side;
6. alignment of the protein-induced shift with centered truth;
7. seed-to-seed cosine of the protein-induced shift vectors;
8. k=0 MSE / CI / Spearman and calibration;
9. per-branch gradient norms; parameter count; peak GPU memory; wall time.

Seeds are averaged **within target** before the component-paired bootstrap over
the 19 `meta_val` components (9,999 draws, seed 20260816).

## Gates

**Primary (P1).** `r_correct(CPC) − r_correct(A0repro)` at k=0 must have a
component-paired lower bound **above zero** and a mean of at least **+0.05**,
the preregistered smallest effect of interest.

**A larger correct-minus-wrong gap is not success.** If the gap widens while
`r_correct` does not improve, the arm has learned to damage the donor
prediction, which is an aversion, not specificity.

Secondary, all evaluated only if P1 passes:

| gate | requirement |
|---|---|
| P2 donor-degradation share | the degradation term contributes **less than half** of the gap |
| P3 alignment | corr(protein-induced shift, centered truth) ≥ **+0.10** (A0 measures −0.014; a random init +0.033) |
| P4 reproducibility | mean pairwise cosine of the shift vectors across seeds ≥ **+0.30** (a random init measures −0.003) |
| P5 no calibration regression | k=0 MSE not worse than `A0repro` by more than 0.10 pK, calibration not worse by more than 0.05 |

## Stop rules

* **P1 fails** → stop. Conclude exactly: *centered-objective training on the
  current `ContactGrammar` at this budget and protocol does not produce
  protein-conditioned within-target ordering.* Do **not** conclude that
  protein-conditioned architectures in general are impossible, and do not run
  the admission stage.
* **P1 passes, P3 or P4 fails** → stop. The objective is producing arbitrary
  protein-dependent movement — the random-initialisation failure mode measured
  in `NOISE_AND_LEAKAGE_AUDIT.md` §7 — and the family closes with that cause
  recorded.
* **P1 passes, P2 fails** → stop, recorded as donor destruction.
* No post-hoc weight, schedule or arm may be added to rescue a failed gate.

## Power, stated in advance

Same-configuration retraining moves aggregate k=0 `r` by **0.051** (R14
screening, A0frozen vs A0repro). The +0.05 minimum effect therefore sits at
about one retraining standard deviation, which is why P1 additionally requires a
component-paired lower bound above zero rather than a point estimate. Three
seeds is the budget; it is not generous, and an unresolved P1 will be reported
as unresolved rather than as a failure or a pass.

## Exact commands

```bash
# A0repro — the matched control (seed in {20260815,20260816,20260817})
conda run -n drug python -m scripts.train_qpsmp \
  --arch similarity_only --steps 1200 --seed <SEED> \
  --split-directory dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1 \
  --output report/meta_fewshot/stageP_cpc_20260816/A0repro_seed<SEED>

# CPCoverdrive — the candidate
conda run -n drug python -m scripts.train_qpsmp \
  --arch similarity_only --steps 1200 --seed <SEED> \
  --protein-contrast-form centered --protein-contrast-loss-weight 2.0 \
  --split-directory dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1 \
  --output report/meta_fewshot/stageP_cpc_20260816/CPCoverdrive_seed<SEED>

# frozen evaluation of both arms
conda run -n drug python -m tools.research.stageP_cpc.evaluate \
  --stage report/meta_fewshot/stageP_cpc_20260816 \
  --output tools/research/stageP_cpc/STAGE_P_meta_val.json
```

## Runtime

Pre-run calibration on a 40-step probe gave 0.70 s/step, implying ~14 min per
run. **The measured in-run rate is 1.14 s/step** (900 steps in 1,022 s),
because the calibration excluded the periodic validation passes. Corrected
estimate: **≈23 min per run, ≈2.3 h for six runs** sequential.

Recorded here rather than silently: the 1.5 h figure that appeared in `task.md`,
`history.md` and `CURRENT_MODEL_EVIDENCE.md` before the run is superseded by
this measurement. The centered form adds one wrong-protein forward on the 240
k=0 steps the incumbent skips (+20% of those forwards, ≈ +5% wall time).

## Standing constraints

Query labels are loss/metric targets only. `meta_val` is read once, after
`meta_train` selection is frozen. `meta_test` is not available to this stage:
its labels are used for no fitting, selection or reported metric, and the
process-isolation incident recorded in `GOVERNANCE_INCIDENT.md` remains open.
No result may be described as pocket-aware, contact-resolved or biologically
localized — the protein path is exactly invariant to residue-slot permutation.
