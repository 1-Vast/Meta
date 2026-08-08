# E0 Evidence Consolidation And Synthetic Identifiability

Updated: 2026-08-07

Decision: `MAP_REALIZATION_OR_OPTIMIZATION_DEFECT`. No repair, real affinity training, DAVIS
access or typed-interaction work is authorized.

## Consolidated Finding

The teacher oracle attains holdout CI `1.00000`;
its deranged-protein oracle CI is `0.51250`,
for partner delta `+0.48750`. The exact
8 x 6 x 5 sufficient statistic reconstructs the direct teacher with maximum
error `2.19e-07`.

The frozen MAP checkpoint reaches train correct CI
`0.85313` and holdout correct
CI `0.68553`. Its holdout
partner delta is `+0.03618`.
The identifiable loss is therefore at T2 -> T3. Existing artifacts do not
separate hypothesis-class realization from optimization because no epoch or
gradient trace was persisted.

One of the eight holdout derangements has local sequence identity >=40%.
This is a control-contract violation, but it does not explain the aggregate
failure: that task retains strong teacher and MAP correct-versus-deranged
contrast. The holdout is also a lexical sample covering only
`1.40%` of fold-4 tasks and
`11.27%` of fold-4 closure components;
it is insufficient for an architecture-wide generalization claim.

## Boundary Table

| Boundary | Observable | Verdict |
|---|---:|---|
| Teacher attainable | oracle correct CI / partner delta | PASS: 1.00000 / +0.48750 |
| Derangement changes teacher | oracle correct - deranged CI | PASS: +0.48750 |
| Teacher sufficient statistics | max reconstruction error | PASS: 2.19e-07 |
| Frozen geometry retention | teacher defined from and exactly reconstructed by P1B geometry | PASS_BY_CONSTRUCTION |
| MAP realization | holdout correct / deranged CI | FAIL: 0.68553 / 0.64934 |
| Optimization convergence | training/gradient trace | NOT AUDITABLE |
| Holdout diversity | tasks / closures / proteins | INSUFFICIENT: 8 / 8 / 8 |
| Derangement contract | pairs at or above 40% local identity | PARTIAL VIOLATION: 1 / 8 |
| Corpus consistency | rows before / after task floor | PASS: 152934 / 152737 |

## Provenance

The 197-row transition is exact: 22 nonempty tasks fell below
20 valid connectivity compounds after enforcing the frozen P1B input contract.
Another 12 governed tasks had zero
model-valid rows, so the governed-to-model task transition is
3817 ->
3783.
Reconstructed and materialized activity-ID sets match:
`true`.
No affinity value field was materialized.

## Stop Rule

The historical synthetic Gate remains failed. E0S is diagnostic only. A future
repair must be separately registered and cannot reinterpret this audit as an
E0-S or typed-interaction authorization.
