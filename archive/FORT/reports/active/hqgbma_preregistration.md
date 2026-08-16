# HQ-GBMA — preregistration (2026-07-26)

Registered **before any Stage result is read**. Design: `reports/active/hqgbma_design.md`. This freezes
(a) the staged plan and (b) the decisive train-only **Stage-D protein→Grassmann-subspace gate**, which
answers the second half of the central research question: *can a protein sequence predict the
low-dimensional subspace in which the identified target–ligand interaction lives, generalising across
held homology components?* No development, confirmation, Davis, or sealed label is read at any point in
this preregistration's scope.

## Admissibility

The OPEN-S audit returned `NO_OPEN_POWERED_INDEPENDENT_PANEL` and the Metz development rows are spent.
Therefore only **train-only mechanism gates** are admissible now. Any later predictive score on spent
development rows is reported as `ARCHITECTURE_MECHANISM_RESULT_ONLY__NO_INDEPENDENT_PREDICTIVE_CLAIM`,
never as external validation. Firewalls retained: target id, UniProt accession, homology component,
ligand parent connectivity, Bemis–Murcko scaffold, Tanimoto<0.95, document, assay. Statistical unit =
homology component.

## Stage-C prerequisite (already substantially established)

The interaction must be identifiable in the exact quotient space (true pairing beats target and ligand
derangement). Gate PA (`p_adapt=0.000488`) and PD-M (3/8 feature-explainable directions, stable under
top-1% ligand removal) already establish this on the identical substrate and projector. Stage D
proceeds on that basis; the quotient projection audit (`project_block_w`: relative KKT `<1e-8`,
idempotence `<1e-7`, LSMR/LSQR `<1e-6`) is re-asserted in the Stage-D run.

## Stage-D substrate and construction (frozen)

* Metz dense kinase panel, TRAIN cells only, registry sha256
  `94da6bb5a59c2911672fde982530c8dd6a673c194b2b2d7b4638df7768c8173e`.
* Ambient space = ligand-feature space `R^m`, `m=64` (64-bin count-Morgan + 10 descriptors → 64-d
  centred PCA, `ligand_feature_table`). Target features = ESM-2 pooled → 32-d PCA (`panel_target_features`).
* Exact projected interaction label `r_int = M_X^W y` (`Substrate.residual`, audited).
* Per **training** target `t`: `β̂_t` = ridge (`ρ=1.0`) of `r_int` on ligand features over `t`'s own
  edges; `V_t` = sandwich covariance with `σ̂²_t` the per-target residual variance. `β̂_t` is the
  target's OWN empirical interaction direction (deliberately **not** `Bᵀu_t`).
* Five frozen homology-component folds (`component_folds`, seed 1729). The protein map `R_t=f(ESM_t)`
  (light net on the 32-d ESM PCA → Stiefel(64,r) via QR retraction) is trained on out-of-fold targets
  by the error-corrected containment loss and evaluated on the held fold's targets. One seed (1729).

## Arms (identical held targets, identical β̂_t/V_t, identical folds)

* **protein(true)** — `R_t=f(real ESM_t)`.
* **global** — one shared subspace `R_0` (top-`r` eigenvectors of `Σ_t (β̂_tβ̂_tᵀ − V_t)` on the
  training fold) applied to every held target; the "interaction is low-rank but not protein-specific"
  null.
* **protein(shuffle)** — ESM↔target assignment deranged within exposure-matched blocks among held
  targets.
* **protein(random)** — Gaussian ESM features at matched scale.
* **protein(pooled-only)** — ESM replaced by amino-acid composition + length (trivial protein
  descriptor); tests whether the pLM is load-bearing over a non-pLM protein representation.

## Metric

Error-corrected containment fraction (higher = better), per held target, averaged per homology
component:
```
c_ec(t) = [ β̂_tᵀ P_t β̂_t − tr(P_t V_t) ] / [ β̂_tᵀ β̂_t − tr(V_t) ] ,   P_t = R_t R_tᵀ .
```
Targets with non-positive noise-corrected signal energy `β̂_tᵀβ̂_t − tr(V_t) ≤ 0` are excluded (no
signal to contain). Paired contrasts use grouped component bootstrap (10,000 draws).

## Frozen pass criteria (a single failure ⇒ `HQGBMA_STAGE_D_FAIL_STOP`)

1. **rank** `r` selected by nested train-only folds lies in `{3,4,6}`.
2. **numerics**: `‖R_tᵀR_t − I‖∞ < 1e-5` for every held target; Hodge projection audit converged.
3. **protein beats the non-protein null**: `c_ec(protein_true) − c_ec(global)` paired mean `≥ 0.02`
   **and** grouped `LCB95 > 0`.
4. **protein specificity — shuffle**: `c_ec(true) − c_ec(shuffle)` grouped `LCB95 > 0`.
5. **protein specificity — random**: `c_ec(true) − c_ec(random)` grouped `LCB95 > 0`.
6. **pLM is load-bearing**: `c_ec(true) − c_ec(pooled_only)` grouped `LCB95 > 0`.
7. **signal exists**: `c_ec(global) > 0` and the fraction of held targets with positive noise-corrected
   signal energy `≥ 0.5` (guards against a vacuous pass).
8. **stability**: after removing the top-1% residual-energy ligands and re-fitting `β̂_t`,
   `c_ec(true) − c_ec(shuffle)` grouped `LCB95 > 0` still holds.

No threshold, rank, width, epoch, or loss weight is increased after a failed result. Passing Stage D
authorises **only** a review of the Stage-E episodic meta-training predictive diagnostic (which, on
spent development rows, is `ARCHITECTURE_MECHANISM_RESULT_ONLY`). It does not authorise multi-seed
runs, Hierarchical MoT, long training, or confirmation/Davis/sealed access. Three seeds are authorised
only after a one-seed Stage-D pass.
