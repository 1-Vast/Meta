# Addendum ADD-2 to PREREGISTRATION.md (issued before any real-data training)

Date: 2026-08-20. Base prereg SHA a7b17e8a..., ADD-1 SHA aa8d06af....
Reason: the synthetic instrument qualification exposed an evaluation-definition
defect in the base prereg, on planted data, BEFORE any real-data run in this
stage. Frozen base rules are never edited; this addendum supersedes the
affected rules.

## Defect

The base prereg made "residual R2" the primary metric (targets
r = c - m_hat). On residual targets the analytic ligand-pattern prior A0 is
ill-posed: it predicts exactly the component the residualization removed,
producing large negative R2 (-3.0 on the planted instrument), so any Delta
R2 over A0 is inflated and meaningless. Checkpoint selection on residual
val loss additionally selected a near-init model because val residual loss
is dominated by pair-memorization overfit after few epochs.

## Correction (frozen)

1. PRIMARY metric everywhere (smoke gates, ladder, SPB verdict): R2 of
   predictions against the FULL centered target c_p(L), pooled over test
   cells. The A0 prior is well-posed on this target by construction.
2. Checkpoint selection: validation weighted MSE against the FULL centered
   target c (the deployable quantity), patience 150, cosine schedule,
   lr 3e-3, epochs 900 (AM-5 already frozen).
3. Residual-target R2 is retained as a SECONDARY diagnostic only, never
   used for gates or selection.
4. SPB SUPPORTED rule (i) becomes: A5 full-centered R2 > A0 full-centered
   R2 with paired parent bootstrap 90% CI excluding 0, evaluated per seed
   and pooled over seeds; remaining sub-rules unchanged.
5. Phase-4 smoke gate (b) unchanged in wording but computed on the full
   centered target; gate (c) unchanged.

## Also frozen here (AM-4, AM-5 rollup)

- AM-4: router alpha = mean-pool skip + ligand-conditioned routed deviation
  (strictly superset of the A1 bilinear; enables the ladder comparison).
- AM-5: lr 3e-3, cosine T_max=900, epochs 900, patience 150 - adopted after
  the instrument diagnosed undertraining at lr 1e-3 / 300 epochs.

The instrument qualification itself is re-run under the corrected protocol;
its verdict gate (recovery >= 0.25 Delta R2 over A0 on the FULL centered
target, 90% parent bootstrap CI excluding 0) supersedes the base wording.
