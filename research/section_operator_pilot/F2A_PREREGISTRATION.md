# F2A preregistration: protein-anchored bioactivity-atlas section

Date frozen: 2026-08-08, after retiring F1B and before evaluating this branch
on either development panel.  KCGS numeric outcomes remain unread.

## Failure addressed

F1B's learned DeepSets head collapsed to an almost constant map.  It marginally
improved interaction prediction over its own support-free surface but ignored
support-label permutations and lost to a location-only update.  F2A therefore
removes the learned adapter and replaces the learned pair basis with a
source-atlas construction whose protein dependence is explicit.

## Construction

Let `R` be the double-centred PKIS1 activity matrix.  A source-only ligand KRR
maps a new ligand to a predicted activity profile over the PKIS1 kinases:

\[
    h(L)=k_L(L,L_s)(K_L+\lambda I)^{-1}R.
\]

For a new protein pocket `P`, three label-blind KLIFS/KiSSim kernels produce
normalized weights over the source kinases.  They respectively use the fixed
hinge, DFG/back-pocket, and front/solvent SiteAlign views:

\[
    v_a(P,L)=\sum_j w_{aj}(P)h_j(L),\quad a\in\{H,D,F\}.
\]

The section basis is

\[
    s_0=(v_H+v_D+v_F)/3,\quad
    s_1=v_H-v_D,\quad s_2=v_F-v_D.
\]

All three surfaces are centred over the unlabeled task ligand set.  Given an
unordered support set, the only adapted coordinates are the positive-ridge
posterior coefficients `(tau,u1,u2)` in

\[
    \hat r=s_0+\tau+u_1s_1+u_2s_2,
    \qquad d_{adapt}=3\leq k=5.
\]

This is a finite reference-atlas realization of a biological section: the
protein fixes the admissible profile neighbourhood and support labels locate a
low-dimensional section inside it.  A wrong protein changes both the atlas
neighbours and the query surfaces; it is not represented by a freely rescaled
protein factor.

## Frozen estimation protocol

- Ligand kernel: Tanimoto on the 1,024-bit radius-two Morgan fingerprint already
  present in the preregistered nuisance vector.  `lambda` is chosen by
  three-fold generic-Murcko-scaffold-cold PKIS1 profile reconstruction from
  `{0.01,0.1,1,10,100,1000}`.
- Protein kernels: mean squared SiteAlign-property distance after the fixed
  KiSSim soft subpocket masks.  Each RBF temperature is the median nonzero
  source-source distance; no activity label selects it.  To avoid dilution by
  biologically remote kinases, each view is normalized over its eight nearest
  source pockets (fixed before outcome evaluation).
- Section penalties are selected by three-fold simultaneous scaffold-cold and
  kinase-group-cold PKIS1 episodes.  Location penalties are
  `{1,10,100}` and tangent penalties are `{0.01,0.1,1,10}`.  Selection minimizes
  five-shot query MSE with deterministic D-optimal support on `[1,s1,s2]`.
- Raw transfer prediction adds the same source-only ligand and protein
  main-effect regressors used in F0/F0R.
- Primary supports contain five distinct scaffolds; all query compounds sharing
  a support scaffold are excluded.

## Controls and admission gate

The gate is identical to F1B.  On PKIS2 at `k=5` with D-optimal support, the
correct atlas section must have positive lower 95% target-bootstrap bounds for
raw MSE reduction versus support-free, location-only, uniform/zero protein,
nearest non-self protein, wrong-target support, and permuted support.  It must
also beat support-free and nearest-protein interaction prediction.  All raw
contrasts must have positive point estimates on Anastassiadis2011.

Random-support results and `k=20` are diagnostic.  Random and D-optimal raw MSE
will not be directly contrasted because their scaffold-excluded query sets can
differ.

## Stopping and integration

Failure retires this exact atlas.  Passing only licenses a later mapping of the
bounded section statistic into the frozen `F(z) -> B(z) -> K(beta)` path; it
does not itself validate a scalar bypass.  No production or frozen-theory file
is changed by this experiment.
