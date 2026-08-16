# HQ-GBMA Stage D — decision (2026-07-26)

Preregistration `reports/active/hqgbma_preregistration.md`; design `reports/active/hqgbma_design.md`.
Train-only, single seed 1729, CUDA (`drug`, RTX 4060). No development, confirmation, Davis, or sealed
label was read. Report: `reports/active/hqgbma_stage_d.json`.

## Question

Can a protein sequence predict the low-dimensional Grassmann subspace `range(P_t)` of the target's
identified interaction, generalising across held homology components — so that few-shot ranking
adaptation confined to that subspace is protein-specific?

## Construction (as preregistered)

Metz pKi dense panel, TRAIN cells only (112 targets / 101 components / 619 ligands, registry
`94da6bb5…173e`). Ambient space = 64-d ligand-feature PCA; exact projected interaction label
`r_int = M_X^W y` (Hodge projection audited: KKT/idempotence/LSMR-LSQR all converged). Each target's
OWN empirical interaction direction `β̂_t` (ridge over its own edges) and sandwich covariance `V_t`
(deliberately **not** `Bᵀu_t`). Protein map `R_t=f(ESM_t)` on Stiefel(64, r) via QR retraction, trained
by error-corrected containment; five frozen homology-component folds. Rank chosen by nested train-only
folds → `r=6` (inner-CV containment 0.073 / 0.076 / 0.220 for r=3/4/6).

## Result (error-corrected containment fraction, held components)

| arm | mean | grouped 95% CI |
|---|---:|---|
| **global** (shared non-protein subspace) | **0.340** | [0.163, 0.548] |
| protein(random) | 0.314 | [0.137, 0.562] |
| protein(shuffle) | 0.237 | [0.109, 0.396] |
| protein(pooled-only, non-pLM) | 0.201 | [−0.026, 0.442] |
| **protein(true ESM)** | **0.108** | [−0.143, 0.324] |

Paired contrasts (grouped component bootstrap, 78 components):

| contrast | mean | LCB95 | criterion |
|---|---:|---:|---|
| true − global | **−0.232** | −0.434 | fail (needed ≥0.02, LCB>0) |
| true − shuffle | −0.130 | −0.511 | fail (needed LCB>0) |
| true − random | −0.206 | −0.549 | fail |
| true − pooled-only | −0.094 | −0.196 | fail |
| stability (true − shuffle, top-1% ligands removed) | −0.124 | −0.281 | fail |

Numerics passed: Stiefel orthonormality max error 4.2e−7; projection audit converged; EB prior PD.
Signal exists: global containment 0.340 > 0, positive-signal fraction 0.777 (87/112 targets).

## Verdict

```text
HQGBMA_STAGE_D_FAIL_STOP
PROTEIN_DOES_NOT_PREDICT_A_TRANSFERABLE_INTERACTION_SUBSPACE
ARCHITECTURE_MECHANISM_RESULT_ONLY__NO_INDEPENDENT_PREDICTIVE_CLAIM
```

The central research question's second half is answered **no** on this substrate, decisively and with a
mechanism. Two facts stand together:

1. **The interaction is real and low-rank — but global.** A single shared subspace (no protein input)
   contains 34% of the noise-corrected energy of held targets' own empirical interaction directions,
   with LCB95 well above zero, at effective rank ~6. This is consistent with Gate PA (`p=0.000488`) and
   PD-M (3–8 explainable directions): there is a transferable low-rank interaction basis across
   kinases.
2. **Protein-conditioning the subspace is not merely unhelpful — it is harmful.** Making the subspace a
   function of ESM generalises *worse* than the shared subspace (−0.232, LCB −0.083) and worse than
   shuffled, random, and pooled-only protein. The `f(ESM)→Stiefel` map fits training-target idiosyncrasy
   that does not transfer across homology components; on held targets an arbitrary shared basis is
   strictly better than a protein-selected one.

This is the same signature as `BM1_RR_FAIL_STOP`, `PANEL_GATE_PC_FAIL_STOP`, and the 2026-07-26
hierarchical protein-conditioned covariance prior (`PROTEIN_CONDITIONED_PRIOR_NOT_LOAD_BEARING`), now
reproduced at the strongest formulation available: a hard Grassmann confinement with exact quotient
supervision and error-corrected teacher coefficients. Across four independent mechanisms —
covariance orientation, precision, unconstrained interaction head, and subspace selection — frozen
ESM-2 does not carry protein-specific SAR-reordering information that transfers dual-cold on this panel.

## What is (honestly) established, and the constructive alternative

* **Established**: an identifiable, low-rank, *shared* interaction subspace transfers across held
  homology components. The active ingredient for within-target reordering is the few-shot
  **support-conditioned** posterior operating in that shared subspace — not the protein.
* **Not established (and now unlikely on this substrate)**: any protein-conditioned subspace/covariance/
  interaction that beats its own destructive controls dual-cold.
* **Constructive route, if pursued (separate preregistration; NOT a rescue of this gate)**: a
  protein-*free* HQ-GBMA — the exact quotient interaction basis + a single global shared subspace +
  the two-posterior (calibration/ranking) meta-adapter of `model/grassmann_bayes.py`. This would be a
  few-shot ranking model whose value comes from support, not protein; it cannot and would not claim
  protein specificity. On spent development rows it is at best `ARCHITECTURE_MECHANISM_RESULT_ONLY`.

## Discipline

A single failed mandatory criterion is a fail-stop. No rank, width, epoch, loss weight, or threshold
was changed after the result (`r=6` was the *largest* grid value and it still failed; increasing it is
inadmissible). Stage E (protein-conditioned episodic meta-training) is **not** authorised. The reusable,
tested primitives (`model/grassmann_bayes.py`, exact quotient via `research/orrc_v2.py`) and the design
report are retained. No confirmation/Davis/sealed access; the Metz development rows remain spent and
unread here.
