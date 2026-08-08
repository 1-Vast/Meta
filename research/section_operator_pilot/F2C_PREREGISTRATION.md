# F2C preregistration: conformal identifiability domain for the atlas section

Date frozen: 2026-08-08, after retiring unconditional F2A and before running
this selective branch on a development panel.  KCGS numeric outcomes remain
unread.

## Mathematical purpose

F2A showed that the protein-anchored atlas is useful on average, but a nearest
non-self pocket can be observationally indistinguishable for some tasks.  The
frozen theory already permits a partial task domain and abstention; forcing a
point prediction outside the identifiable domain would contradict that
structure.

F2C augments the F2A coordinates with a source-calibrated certificate

\[
    z=(\tau,u_1,u_2,c),\qquad d_{adapt}=4\leq k=5,
\]

and defines the section only when `c>0`.  Abstention is a model output, not an
evaluation exclusion chosen from query outcomes.

## Support-only certificate

For each five-shot support set, compute the following without query labels:

1. leave-one-support-out errors of the correct atlas, support-free atlas,
   location-only update, uniform-protein atlas, and nearest-pocket atlas;
2. the corresponding four LOO error reductions relative to the correct atlas;
3. D-optimal condition ratio of `[1,s1,s2]`;
4. correct-versus-nearest surface separation over the unlabeled task ligands;
5. support residual spread, correct support-fit error, and section coefficient
   norm.

A fixed `StandardScaler + Ridge(alpha=10)` maps these features to the smallest
query MSE advantage of the correct section over all raw controls and over the
support-free/nearest interaction controls.  Training examples are generated
only from three simultaneous scaffold-cold and kinase-group-cold PKIS1 folds.
Ten deterministic random-support episodes per held target are used in each
fold. Predictions for calibration are cross-fitted by the held source fold.

Let `e = predicted_margin - observed_margin`.  The one-sided 80% split-conformal
quantile `q` of cross-fitted `e` is frozen.  The deployment certificate is

\[
    c=\widehat{margin}(S,Q,P)-q.
\]

The task is admitted iff `c>0`.  This is a marginal lower certificate under the
source episode exchangeability assumption; it is not claimed to be a
distribution-free conditional guarantee under panel shift.

## Fixed episode and endpoint protocol

- The F2A source-selected ligand ridge and section penalties are recomputed
  without development data.
- Primary support policy is random selection of five distinct generic Murcko
  scaffolds with seeds `20260808..20260827`.  F2A showed that D-optimal selection
  can change the excluded query population, so D-optimal remains diagnostic.
- Every arm is evaluated on exactly the same admitted episode/query cells.
- Interaction prediction is corrected here: the predicted residual curve is
  centred over the task's finite unlabeled ligand universe before comparison
  with the panel's double-centred outcome.  This removes the adapted constant
  `tau`; the previous per-query residual endpoint incorrectly retained it.
- Bootstrap resamples target clusters, retaining all seeds within a target.

## Admission gate

PKIS2 primary success requires:

1. admitted episode coverage at least 20%, spanning at least 30 target clusters;
2. positive lower 95% target-bootstrap bounds for raw MSE reduction versus
   support-free, location-only, uniform protein, nearest protein, wrong-target
   support, and permuted support;
3. positive lower 95% bounds for centred interaction MSE reduction versus
   support-free and nearest protein; and
4. exact support-order invariance.

Anastassiadis2011 must have at least 10 admitted target clusters, at least 10%
episode coverage, and positive point estimates for the same contrasts.  These
coverage floors prevent a vacuous certificate.

Passing licenses only the biological statistic and partial-domain certificate
for later `F(z) -> B(z) -> K(beta)` integration.  It does not license a scalar
bypass or any ranking claim.
