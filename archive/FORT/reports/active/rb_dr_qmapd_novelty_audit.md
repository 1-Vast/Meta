# RB-DR-QMAPD — novelty audit (2026-07-26)

Rao–Blackwellized Design-Robust Quotient-Space Completion-Marginalized Amortized Posterior
Distillation. This audit is written before any oracle gate is run and distinguishes what already exists
from the residual contribution. Per the spec, no component is claimed novel individually.

## Provenance note (honesty)

The task references three reports as source-of-truth: an HSBO audit, a Q-MAPD proposal, and a
Design-Robust Q-MAPD recommendation. **None of these files exists in this repository.** The authoritative
record here is `history.md` plus the existing reports (`hqgbma_*`, `hier_bayes_meta_*`, `panel_gate_*`,
`orrc_eb_blueprint_v2/v3`, `cfri_*`, `bm0/bm1/bm2`). "HSBO" and "Q-MAPD" appear to belong to a different
research lineage not present here. This audit and program therefore build on the actual record and do not
cite non-existent reports. Every method claim below is checked against `history.md`.

## Prior-art verification (all 2026, checked on the web)

| paper | id | what it already provides |
|---|---|---|
| AdaMBind — meta-learning task-adaptive DTA | Nat Commun 2026 | MAML few-shot DTA with easy-to-hard task scheduling (its cold-start gains are on random-similarity splits that ICLR'25 shows are inflated) |
| Neural Operator Processes (NOPs) | arXiv 2606.22946 | neural-process conditioning + operator decoding; **query-aligned attention**; probabilistic prediction from sparse/partial observations |
| Amortized Energy-Based Bayesian Inference | arXiv 2605.15407 | observation-conditioned transport map trained by an **energy-distance** objective; likelihood-free, no Jacobians |
| Flow Annealing Posterior Sampling (FAPS) | arXiv 2606.22346 | function-space posterior sampling from sparse/noisy observations, variable query discretisation |
| Conservative NPE via DRO (DRO-NPE) | arXiv 2605.28516 | **Wasserstein/entropic DRO** worst-case objective for conservative, better-calibrated amortized posteriors |
| Nearly Optimal Bayesian Inference for Structural Missingness | arXiv 2601.18500 | **posterior integration over missingness** (avoid plug-in imputation bias); near-Bayes-optimal under an SCM prior |
| Privileged Information Distillation for LMs (π-Distill) | arXiv 2602.04942 | teacher with **privileged information**, student acts without it at inference |
| Self-Distilled RLVR (RLSD) | arXiv 2604.03128 | privileged teacher causes **information leakage/instability** if it sets gradient *direction*; fix = teacher sets *magnitude*, real signal sets direction |

Mathematical foundations (older, permitted): conditional expectation and Rao–Blackwellization; law of
total probability / total covariance; latent-variable marginal likelihood; measurement-error
random-effects (REML/EM); strictly proper multivariate scores (energy, variogram); entropic
distributionally robust optimisation; paired grouped bootstrap; the weighted Hodge quotient (already in
`research/orrc_v2.py`, audited).

## What is NOT novel (do not claim)

Neural processes; sparse context-to-function operators; query-aligned attention (NOPs); function-space
posterior prediction (NOPs/FAPS); privileged-information distillation (π-Distill, RLSD); missing-variable
integration (Structural-Missingness, law of total probability); energy-distance training (Amortized
Energy-Based BI); distributionally robust optimisation (DRO-NPE); Hodge projection (this repo).
Rao–Blackwellization, the law of total covariance, and measurement-error random effects are classical.

## Residual contribution being tested (and its honest status)

The only candidate contribution is a **specifically constrained combination** aimed at recovering a
dual-cold, target-specific SAR *ranking* function from few support measurements:

1. the interaction estimand is defined in the **exact weighted-Hodge quotient** (train-only teacher),
   removing target/ligand main effects that every prior DTA distillation ignores;
2. a **completion-design-marginalized (Rao–Blackwellized) teacher law** `P̄_j = Σ_m w_jm P_jm` with the
   **full law-of-total-covariance** (within- + between-completion), so distillation matches a marginal,
   not a single imputed completion (grounds: Structural-Missingness, LTC);
3. **entropic design-robustness** across several *legal, label-blind* completion designs (grounds:
   DRO-NPE), so the student is not tuned to one arbitrary observation design;
4. a **pointwise-consistent** amortized set-to-function student whose per-query prediction is invariant
   to which unrelated queries share the batch (grounds: NOPs query-aligned attention), avoiding the
   transductive leakage that has repeatedly produced retracted gains here (task.md frozen negative #4);
5. the privileged teacher is a **bounded weighting/matching signal, never the sole gradient direction**
   (grounds: RLSD leakage finding), with a query-label anchor at coefficient 1.

**This combination is a hypothesis, not a result.** The repository record is a long series of real
signals (biological, pretraining, interaction, uncertainty, protein-conditioning) that did **not**
transfer to the cold-target estimand — most recently `HQGBMA_STAGE_D_FAIL_STOP` (frozen ESM predicts no
transferable interaction subspace) and `PROTEIN_CONDITIONED_PRIOR_NOT_LOAD_BEARING`. The novelty of the
combination confers **no presumption** of success. The decisive question is prior to the architecture:
does additional same-target evidence contain reproducible, scaffold-generalising SAR information beyond
the same teacher at k shots (Stage O1)? If not, none of the above is implemented. Falsification has
priority over completeness.

## Sources

AdaMBind (nature.com/articles/s41467-026-70554-5); NOPs (arxiv.org/abs/2606.22946); Amortized
Energy-Based BI (arxiv.org/abs/2605.15407); FAPS (arxiv.org/abs/2606.22346); DRO-NPE
(arxiv.org/abs/2605.28516); Structural Missingness (arxiv.org/abs/2601.18500); Privileged Information
Distillation (arxiv.org/abs/2602.04942); Self-Distilled RLVR (arxiv.org/abs/2604.03128).
