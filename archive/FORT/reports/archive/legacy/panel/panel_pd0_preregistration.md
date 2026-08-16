# Gate PD0 preregistration - ORRC train-only recoverability and identifiability

Registered 2026-07-25 before any PD0 statistic was computed. Train-only on the dense panel
(`dataset/public/chembl_37/processed/panel_metz/`, registry sha256
`94da6bb5a59c2911672fde982530c8dd6a673c194b2b2d7b4638df7768c8173e`, 12,574 train cells, 112 targets,
101 homology components, 619 anchor ligands). No panel development cell, no panel confirmation cell
and no ChEMBL confirmation label is read. PB and PC keep their verdicts unchanged; nothing in PD0 may
reinterpret them.

## Estimand

With observed train edges `e=(t,l)` in `Omega`, labels `y` and the target/ligand incidence design
`X`, the interaction residual is `r = M_X y` where `M_X` is the exact Hodge projector onto the
orthogonal complement of the incidence column space. Frozen features define the edge design
`Z_e = u_t (x) v_l` (row-major reshape, so the bilinear form is `u_t^T B v_l`), and exact nuisance
orthogonality is imposed in **observed-edge space**:

```text
Z_perp = M_X Z ,   g_Omega = Z_perp vec(B)
```

Separately centering the target and ligand feature matrices is *not* used and is not claimed to be
orthogonality on this incomplete panel.

The fit is the convex problem

```text
min_{B,S}  0.5 || r - Z_perp vec(B) - P_Omega(A_t S A_l^T) ||^2
           + lam_B ||B||_*  + lam_L ||S||_*
```

`A_t` and `A_l` are orthonormal bases of the `D`-weighted orthogonal complements of the frozen target
and ligand feature spaces, with `D_t` and `D_l` the observed-cell counts per target and per ligand:
`A_t` spans `{x : U^T D_t x = 0}` and `A_l` spans `{x : V^T D_l x = 0}`. Because both bases are
orthonormal in the standard inner product, `||A_t S A_l^T||_* = ||S||_*`, so the latent penalty is an
exact nuclear norm and the constraint set is a linear subspace: the problem stays convex and the
proximal operator is exact. Missing cells are never filled with zero; `P_Omega` samples only observed
edges and its adjoint scatters only observed edges.

## Frozen features and bases

Target `u_t`: 32 centered PCA coordinates of the frozen pooled ESM-2 650M embedding, the identical
basis Gate PA used. Ligand `v_l`: 64 centered PCA coordinates of the 64-bin count-Morgan plus ten
physicochemical descriptor basis, the identical basis Gate PA used. No affinity label selects either
basis. `U = [1, u]` and `V = [1, v]` include an explicit constant so the latent complement cannot
reintroduce a target or ligand main effect.

## Frozen readout

The out-of-sample readout is the product reference with **training-only** feature means:

```text
g_hat(t,l) = (u_t - u_bar_train)^T B (v_l - v_bar_train)
```

`Z_perp` is the training-edge estimating equation only and is never used as a prediction rule. In
PD0 every held-out statistic uses the inner-training-fold means, never the held-out fold's.

## Frozen hyper-parameter selection

Nested validation over **training homology components** using the existing fixed fold map
(`research/panel_power.component_folds`, seed 1729, five folds), which is the same map Gate PB used.
The grid is scale-free: `lam_B in {0.50, 0.25, 0.10, 0.05} x lam_B_max` and
`lam_L in {1.00, 0.50, 0.25, 0.10} x lam_L_max`, where `lam_max` is the smallest penalty that makes
the corresponding block exactly zero (`lam_B_max = ||U^T diag(r) V||_op` on the projected design,
`lam_L_max = ||A_t^T R_Omega(r) A_l||_op`). `lam_L = 1.00 x lam_L_max` forces `L = 0` exactly and is
the no-latent member of the grid.

Selection criterion, fixed here: the mean over held-out components of the within-component Pearson
correlation between `r` and the frozen product readout, where the fit saw none of that component's
cells. Ties resolve toward the larger penalty. The grid is restricted in advance to candidates whose
full-train effective rank is at most 8, because rank <= 8 is the model class the ORRC contract
defines, not an outcome; if no candidate satisfies that, PD0 fails on the rank criterion.

## Fixed spectral rule

`effective_rank(B) = #{ i : sigma_i >= 0.05 * sigma_1 }`. The exact numerical rank
(`sigma_i > 1e-10`) is reported alongside but is not the gating quantity.

## Statistics and inference

For each of the five folds, ORRC is fitted on the other folds' components and the leading eight
singular directions `(sigma_i, a_i, c_i)` of that fold's `B` produce a held-out per-edge rank-1
prediction `g_hat_i`. The unit of inference is the homology component: for every component the
statistic is the within-component Pearson correlation between `r` and `g_hat_i` on that component's
held-out edges. This gives one value per component per direction, and each component is scored by a
fit that never saw it.

* **Component bootstrap**: 10,000 resamples of the 101 components; LCB95 of the mean.
* **Permutation nulls**: 4,096 deterministic graph-exposure-matched permutations of the target
  features and, separately, 4,096 of the ligand features, applied to the held-out edges with the
  fitted directions held fixed. Because the directions were fitted on other components, this null is
  exchangeable. Blocks are matched on label-free graph exposure (target: edges and distinct ligands;
  ligand: edges and distinct targets). `p = (1 + #{null >= observed}) / (n + 1)`.
* Edge resampling is never used.

## Criteria

| id | criterion | threshold |
|---|---|---|
| PD0-1 | exact projection audit | relative KKT `< 1e-8`, idempotence `< 1e-7`, LSMR/LSQR disagreement `< 1e-6` |
| PD0-2 | convex solver audit | primal feasibility `<= 1e-8` relative, nuclear-norm dual feasibility violation `<= 1e-3`, complementarity gap `<= 1e-3`, deterministic repeat agreement `<= 1e-6` |
| PD0-3 | effective transferable rank | in `[1, 8]` |
| PD0-4 | at least one direction `i <= 8` is feature-explainable | target-permutation `p <= 0.01` **and** ligand-permutation `p <= 0.01` **and** component-bootstrap `LCB95 > 0` |
| PD0-5 | stability | PD0-4 still holds for at least one direction after removing the top 1% residual-energy ligands |

Auxiliary-feature energy `||Z_perp vec(B)||^2 / ||r||^2` and latent-only energy
`||P_Omega(L)||^2 / ||r||^2` are both reported. A latent component without a stable
feature-explainable direction returns `ORRC_PD0_FAIL_STOP` and does not authorize prediction.

## What a pass authorizes

Only Gate PD1: one deterministic five-fold leave-homology-component-out development run on exactly
the Gate PB rows, preceded by its own power record that freezes `max(0.03, MDE80)` from
component-level **arm heterogeneity** in the already observed PB contrasts rather than same-arm
retraining noise. PD0 authorizes no development prediction by itself, no few-shot posterior, no
signed prior mean, no multi-seed run, no Hierarchical MoT, no long training and no confirmation
access.

---

## SUPERSEDED, 2026-07-25

This preregistration is withdrawn by the audit verdict `BLUEPRINT_REQUIRES_MAJOR_REVISION` and is
replaced by `reports/active/orrc_eb_blueprint_v2.md`. It is kept unedited for traceability.

Three of its specifications are now known to be wrong or inadmissible:

1. it imposed exact orthogonality through the unweighted projector `M_X` and constrained the latent
   block by full-grid two-sided feature orthogonality. Section C.1 of the revised blueprint gives an
   explicit 3x3-minus-one-cell counterexample: that constraint permits the latent block to align 29%
   with the interaction direction, so the `B`/`L` split was not identifiable;
2. it restricted the penalty grid to candidates of effective rank `<= 8`, which reimposes a
   non-convex rank constraint through the selection rule. Effective rank must be a post-hoc frozen
   spectral rule with permutation and component-bootstrap support;
3. its pass clause authorised Gate PD1 on the Gate PB development rows. ORRC-EB was designed after
   those rows were observed, so they cannot provide independent confirmation for this route. A
   predictive gate now requires a newly registered independent panel or separately approved untouched
   confirmation access (blueprint section F).

The run that this document authorised is recorded as void in
`reports/active/panel_gate_pd0_status.md`. No result from it is carried forward.
