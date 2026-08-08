# F6P preregistration: component-identifiable orthogonal product section

Date frozen: 2026-08-08, after F5C and before evaluating this decomposition.
KCGS numeric outcomes remain unread.

## Resolution of the identifiability mismatch

Previous gates incorrectly required one support-adapted mean predictor to prove
both (i) ligand-coupled biological interaction and (ii) task-level location.
Those are different identifiable components.

F6P fixes the biological interaction at the source atlas curve `s(P,L)` and
adapts only a nuisance/task location:

\[
 \tau(S,P)=\Pi_{[-0.5,0.5]}
 \frac{\sum_{(L_i,r_i)\in S}\{r_i-s(P,L_i)\}}
      {k+\lambda_\tau},
 \qquad
 \hat r(P,L)=s(P,L)+\tau.
\]

Thus `d_adapt=1<=k=5`.  A wrong protein can correct only a constant and cannot
change its query interaction shape.  Joint support-pair reordering leaves the
map invariant.  Scrambling labels across the same support ligands also leaves
this *location-only* coordinate invariant by algebra; this is a declared null,
not evidence of failure.  Ligand-label coupling is not claimed for `tau`.

## Frozen protocol

- `s(P,L)` is the average of the hinge, DFG/back-pocket, and front/solvent
  eight-nearest-anchor atlas views from F2A, centred over the unlabeled task
  ligand universe.
- Ligand KRR regularization is selected by PKIS1 scaffold-cold profile
  reconstruction.  `lambda_tau` is selected from `{1,10,100}` on five random
  distinct-scaffold episodes per target in simultaneous PKIS1 scaffold-cold and
  kinase-group-cold folds.
- External evaluation uses 20 fixed random support seeds; queries sharing a
  support scaffold are excluded.

## Componentwise gate

Raw PKIS2 prediction must beat support-free, location-only, uniform protein,
nearest non-self protein, and wrong-target support with positive lower 95%
target-bootstrap bounds.  The five point estimates must be positive on
Anastassiadis2011.

The biological interaction curve must beat the nearest-protein curve with a
positive PKIS2 lower bound and positive Anastassiadis point estimate.  Its
interaction prediction is exactly equal to the support-free curve by
construction and is checked as an equality invariant.  Label scrambling is
also checked as an exact prediction invariant and is not included as an
improvement contrast.

Passing admits two separately identified facts: the protein-dependent atlas
interaction and the support-dependent location statistic.  Their later joint
effect must still enter the frozen law-valued `F(z) -> B(z) -> K(beta)` path;
no scalar bypass or ranking claim is licensed.
