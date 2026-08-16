# K-LBP v2 R1 decision

## Prior result invalidation

The first R1 result produced on **2026-07-28** is **INVALIDATED**. It materialized the affinity
column before dropping `y` and used non-registered strike-loss and missingness handling. Its verdicts
and `klbp_r1_coordinates.npz` must not be consumed.

## Current audit

- `det_proxy_card`: `R1_COORDINATE_SURVIVES`
- `klifs_pocket_composition`: `R1_NON_AUTHORIZED_DIAGNOSTIC_ONLY`

Authorized surviving coordinates: `det_proxy_card`.

`klifs_pocket_composition` is a non-authorizing diagnostic only. It is not the dependency-gated
structure-token or geometry coordinate registered in `task.md` Part 9.

## Claim boundary

R1 reads no affinity label and authorizes no affinity, training, or predictive claim. Survival means
only that the audited coordinate did not cross the frozen taxonomy, ESM-redundancy, or popularity stop.
