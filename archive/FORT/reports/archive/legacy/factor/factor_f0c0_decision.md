# FACTOR-C F0-C0 decision

Date: 2026-07-26  
Verdict: `FACTOR_F0C0_REPRESENTATION_UNIDENTIFIED_STOP`

Authoritative output: `reports/active/factor_f0c0.json`, SHA-256
`F1AA1DEC48AF7226FF2534A0C20700491C206995454930AEEC7FE01189F88CF8`.
Preregistration SHA-256:
`D2EE11EF74DD73EF76ACA48007EE0D56689BEC7E04F0928D7FFF6690C01E5AAB`.

No activity/inhibition value, quarantined ChEMBL confirmation row or sealed test was read.

## Result

| Gate family | Result | Decision |
|---|---:|---|
| external source-balanced coverage median / q10 | 0.644 / 0.421 | fail vs 0.90 / 0.70 |
| inner KIRHub-held calibration median / q10 | 0.836 / 0.591 | fail |
| inner Christmann-held calibration median / q10 | 0.852 / 0.551 | fail |
| inner Reinecke-held calibration median / q10 | 0.766 / 0.530 | fail |
| true-minus-decoy grouped LCB95 | +0.435 | pass |
| role, bond and motif-attachment reconstruction | all margins pass in all folds | pass |
| effective-rank ratios | 0.340–0.362 | pass vs 0.25 |
| pair/motif calibration false coverage | 0.046–0.053 | approximately at gate |
| atom calibration false coverage | 0.294–0.334 | fail |

Strict corpus firewalls are intact in every fold: zero held-vs-atlas parent-connectivity overlap, zero
held-vs-atlas Murcko-scaffold overlap, and zero inner train-vs-validation scaffold overlap.

## Why this is not a real-support verdict

The fixed atom representation keeps only 13 nonconstant robust-scaled dimensions. For common neutral
atoms, permuting formal charge and incident-bond summaries frequently leaves the vector unchanged.
Those unchanged records were counted as chemistry-broken decoys, so the atom false-coverage gate is
not a valid semantic anti-cheat test. At the same time, inner scaffold-OOD coverage fails, showing
that the fixed rules do not yet provide an adequate continuous ontology even before external-source
transfer.

The large and stable true–decoy gap, strong structural reconstruction and noncollapsed effective rank
show that the representation is not empty. They do not rescue the invalid bandwidth calibration.
Consequently the low external 0.644/0.421 coverage cannot be attributed confidently either to true
public-data chemical insufficiency or solely to discrete aliasing.

## Consequence

F1-C is locked. The user-specified F0-C1 route is now scientifically permissible because C0 failed
at representation/calibration rather than establishing a real-support failure. F0-C1 must be
separately preregistered, use only unlabeled inner-train molecular graphs, enforce the same strict
corpus firewall, and use decoys that are verified to change chemistry carrier-by-carrier. No external
coverage result may select its architecture or bandwidth.

The 96 pair-carrier cap is a low-cost audit bound, not a final-model constraint. If the mechanism is
later proven, full pair expansion and formal training should run on a higher-memory machine with a
compute-matched replication.

