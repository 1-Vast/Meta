# LOCK/CLOCK G0 pre-A1 invalidation

The first G0-L, preparation, and G0-R artifacts are invalid and are retained only for audit:

| artifact | SHA-256 |
| --- | --- |
| `lock_clock_g0_invalidated_pre_a1_label_free.json` | `9367fa1c257fe0bd0960de912907aa25393d61ec034054e0b5c7e8c70f396372` |
| `lock_clock_g0_invalidated_pre_a1_reordering_preparation.json` | `bd4d467be5ba05479fafc43cb01c0a3246ea577a6ee30c628e89f84eb7b65e43` |
| `lock_clock_g0_invalidated_pre_a1.json` | `73902075627c9562ac9d036a978fb18d5575e9aff577b2f9a0cd67d8ccfca44b` |

An independent audit completed after these files were written and showed that `same_report()` did not
enforce the amendment's exact comparison for discrete fields. Python equality accepted `1 == True`,
and the numeric branch accepted an integer and a nearby float. A hand-edited G0-L or preparation
artifact could therefore have passed a recomputation check despite a discrete type change.

There is no evidence that the three files were edited or that their numerical calculations changed.
Nevertheless, their integrity gate did not satisfy the pre-result contract, so none is interpreted.
Amendment A1 changes only typed artifact comparison and its regression tests. It changes no source,
coordinate, randomization, seed, fold, estimator, threshold, bootstrap, MDE, or verdict rule.

