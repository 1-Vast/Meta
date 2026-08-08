# F1B preregistration: section-conditioned bilinear mechanism operator

Date frozen: 2026-08-08, before implementation or inspection of KCGS numeric
outcomes.

## Why this branch exists

F0 and F0R established three facts: dense kinase panels contain a large
double-centred interaction component; few-shot support labels contain useful
target-specific information; and a rank-two SVD section can nevertheless make
the supplied protein pocket dispensable.  The last failure is structural.  In
a separable expansion

\[
    r(P,L)=\sum_q u_q\,a_q(P)b_q(L),
\]

free per-axis adaptation can absorb a wrong value of `a_q(P)`.  Protein
derangement is therefore not identifiable merely by increasing the quality of
the ridge solve or the support design.

F1B removes that gauge freedom.  It learns three *joint* pair surfaces

\[
    (s_0,s_1,s_2)=H\!\left(
      \{\phi_a(L)\odot\psi_a(P)\}_{a\in\{H,D,F\}}
    \right),
\]

where `H`, `D`, and `F` are fixed hinge, DFG/back-pocket, and front/solvent
biochemical views.  A support-permutation-invariant section network may adapt
only four bounded coordinates

\[
    z(S,P)=(\tau,u_1,u_2,c),\qquad d_{adapt}=4\leq k=5,
\]

and predicts the interaction residual by

\[
    \hat r(P,L)=s_0(P,L)+c\{\tau+u_1s_1(P,L)+u_2s_2(P,L)\}.
\]

Changing `P` changes every query surface.  Two tangent coordinates cannot in
general reconstruct that change over an entire query set.

This is an experimental bridge to the frozen theory, not a production change.
No file under `model/` or the frozen theory directory may be edited unless the
admission gate below passes.

## Fixed biological construction

1. Ligands use only RDKit-derived, centre-restricted Morgan fingerprints from
   the five already preregistered pharmacophore channels plus twelve ordinary
   molecular descriptors.  No compound identifier is a feature.
2. Proteins use the 85 aligned KLIFS pocket residues encoded by the public
   SiteAlign physicochemical table.  Residues are softly masked by their fixed
   KiSSim distance to the hinge, DFG, and front-pocket centres.  No target or
   family identifier is a feature.
3. The three ligand/protein views are fitted with source-only unsupervised PCA,
   followed by a learned bilinear Hadamard fusion.  The masks and view/channel
   assignments are fixed before labels are read by the training loop.
4. Supervision for the pair surfaces is the PKIS1 double-centred activity
   matrix.  This algebraically removes ligand and target main effects.  Raw
   predictions add separately fitted source-only ligand and protein main-effect
   regressors.

## Meta-learning protocol

- Source: PKIS1 only.
- Architectural checkpoint selection: a single deterministic, simultaneous
  Murcko-scaffold-cold and kinase-group-cold source split.  The full source
  model is then retrained for the selected number of steps.
- Episodes contain five support compounds from distinct generic Murcko
  scaffolds.  Queries exclude every support scaffold.
- The adapter is a DeepSets map: identical element encoder, mean and second
  moment aggregation, then a bounded four-coordinate head.  Thus support order
  cannot affect the result.
- Loss terms are query residual MSE, support-free residual MSE, and a
  double-difference rectangle loss.  Protein derangement is an evaluation
  control and is not an optimization margin.
- Primary support selection is deterministic D-optimal selection on
  `[1,s1,s2]`, one compound per scaffold.  Random support is reported
  separately, but its raw MSE is not compared directly with D-optimal MSE
  because their scaffold-excluded query sets differ.

## Frozen hyperparameters

- PCA dimensions: 24 per ligand view and 12 per protein view.
- Bilinear width: 16 per view; pair-head width: 32; DeepSets width: 32.
- Optimizer: AdamW, learning rate `1e-3`, weight decay `1e-4`, gradient norm
  clipped at `5`.
- Episode batch: 24; five supports and 24 queries per episode.
- Loss weights: episodic `1.0`, support-free `0.5`, double-difference `0.25`,
  coordinate penalty `1e-3`.
- Validation every 50 steps, maximum 1,500 steps, patience 8 validations.
- Training seed `20260808`; evaluation seeds `20260808..20260827` for random
  support.  Bootstrap uses 10,000 target-cluster draws.

No hyperparameter is selected from PKIS2 or Anastassiadis2011.

## Controls and admission gate

The primary gate is PKIS2, `k=5`, D-optimal support, with target-cluster paired
bootstrap confidence intervals.  The correct operator must have strictly
positive lower 95% confidence bounds for MSE reduction versus all of:

1. support-free pair prediction;
2. shrinkage-only location calibration;
3. zero/mean protein features with the same labels;
4. the nearest non-self pocket with the same labels;
5. support labels from a different target; and
6. permuted support labels.

The same six contrasts must have positive point estimates on the independent
Anastassiadis2011 development panel.  Correct-pocket interaction prediction
must additionally beat a zero interaction and nearest-pocket interaction on
PKIS2.  Support-order invariance is an exact test, not a statistical endpoint.

Only if all conditions pass is the biological coordinate admitted for a later
law-valued `F(z) -> B(z) -> K(beta)` integration.  A scalar affinity decoder is
diagnostic only and is never evidence that the frozen operator itself has been
implemented.

## Read firewall and stopping rule

- PKIS2 and Anastassiadis2011 are consumed development panels.
- KCGS numeric outcomes, DAVIS labels, and recipient labels remain unread.
- If the gate fails, F1B is retired.  Allowed next branches are a genuinely
  different joint pair basis, a different biologically specified observation
  law, or an explicit non-identifiability result.  Reweighting these same three
  surfaces against the failed development panels is forbidden.
