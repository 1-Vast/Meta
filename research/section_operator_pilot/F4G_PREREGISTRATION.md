# F4G preregistration: gauge-fixed biological atlas section

Date frozen: 2026-08-08, after F3K showed that enlarging the support set does
not make the flexible two-tangent section reliably identifiable, and before
evaluating this reduced law.  KCGS numeric outcomes remain unread.

## Structural repair

The failed sections allowed two support-fitted biological tangent directions.
Even with a joint protein/ligand atlas, those directions can fit support noise
and can partly compensate a wrong pocket.  F4G removes this gauge freedom.

The three KLIFS/KiSSim atlas views are averaged into one fixed, task-centred
biological interaction curve

\[
    s(P,L)=\{v_H(P,L)+v_D(P,L)+v_F(P,L)\}/3.
\]

The support set may fit only a location and a scalar amplitude around the
source prior `(tau,a)=(0,1)`:

\[
    \hat r(P,L)=\tau+a\,s(P,L),\qquad d_{adapt}=2\leq k=5.
\]

Both coordinates are obtained by a positive-ridge posterior.  A wrong protein
changes the entire curve `s(P,.)`; one amplitude cannot reconstruct an
arbitrary changed query shape.  Joint support-label permutation changes the
support covariance with `s`, while reordering support pairs leaves the solve
exactly invariant.

## Frozen protocol

- The ligand activity atlas, eight-nearest-source protein anchors, KiSSim soft
  subpocket masks, and additive main-effect regressors are identical to F2A.
- Ligand KRR regularization is chosen by PKIS1 scaffold-cold profile
  reconstruction.
- Location penalty is selected from `{1,10,100}` and amplitude penalty from
  `{0.1,1,10,100}` using five deterministic random five-shot episodes per
  target on three simultaneous scaffold-cold and kinase-group-cold PKIS1
  folds.  The amplitude prior is fixed at one.
- Primary external support policy is random selection of five distinct generic
  Murcko scaffolds over seeds `20260808..20260827`.  Queries sharing a support
  scaffold are excluded.
- Interaction curves are centred over the finite unlabeled task ligand
  universe before comparison with double-centred outcomes, so the location
  coordinate cannot affect the interaction endpoint.

## Gate

At full coverage, PKIS2 must show positive lower 95% target-bootstrap MSE
reduction versus support-free amplitude, location-only, uniform/zero protein,
nearest non-self protein, wrong-target support, and permuted support.  Centred
interaction prediction must beat support-free and nearest-protein controls.
The corresponding eight point estimates must be positive on
Anastassiadis2011.

Passing licenses the two-dimensional biological statistic for a later
law-valued operator integration.  It does not license a scalar bypass, ranking
claim, or edit to frozen/production files.
