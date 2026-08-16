# Verification of the eleven binding audit findings

Numerical authority: `AUDIT_VERIFICATION.json`, produced by
`scripts/audit_stage10_claims.py`. Every value is recomputed from the retained
Stage 9/10 artifacts, the three accepted `similarity_only` checkpoints and the
governed corpus. No number is taken from a narrative report.

**All eleven findings are confirmed.** The original artifacts are unchanged; only
the claims made about them are corrected.

| # | finding | verdict | recomputed evidence |
|---|---|---|---|
| 1 | Stage C's 12.3% is a development result, not confirmed | accepted | follows from 2 |
| 2 | `beta`, source and `w` selected on the same data used for inference | confirmed | `stage10_retrieval_prior.py` sweeps 3 sources x 5 weights on `meta_val`; `BOOTSTRAP_meta_val.json` resamples that same population |
| 3 | 44 targets, **10** components (not 11) | confirmed | 44 targets, 10 components in `BLEND_meta_val.rows.jsonl` |
| 4 | only 6 of 10 components improve | confirmed | 6 positive, 4 negative (-0.155 to -0.012) |
| 5 | 305/624 query cells have the exact ligand in `meta_train` | confirmed | 305/624 = 48.88% |
| 6 | on ligand-disjoint targets the gain is small and unresolved | confirmed | 12 targets, 5 components: **+0.050 [-0.074, +0.175]** |
| 7 | Stage A's 25.5% is transductive | confirmed | `shape - shape.mean() + level.mean()`, both means over the query panel |
| 8 | the prior is an offline evaluator | confirmed | no `ARCHITECTURES` entry, no checkpoint tensor, absent from `evaluate_qpsmp.py` |
| 9 | `w=0` uses `beta=8`, not the learned scale | confirmed | learned 7.9749 / 7.9849 / 7.9897; `w=0` is *approximately*, not exactly, the checkpoint |
| 10 | raw pooled ESM cosine is highly compressed | confirmed | see below |
| 11 | `meta_test` is consumed | accepted | policy |

## The two findings that change what happens next

**Finding 6 is the decisive one.** Restricting the Stage C blend to the 12
targets whose query ligands never appear in `meta_train`:

| population | targets | components | k=0 MSE reduction | 95% component interval |
|---|---:|---:|---:|---|
| all Stage C targets | 44 | 10 | +0.198 | **[+0.016, +0.405]** |
| exact-ligand-disjoint | 12 | 5 | +0.050 | **[-0.074, +0.175]** |

The headline effect is roughly four times smaller and no longer resolved once
exact ligand recall is removed. The 12.3% figure is therefore **development
evidence about a population that is half exact-recall**, not evidence of
chemical-neighbourhood transfer.

**Finding 10 reopens protein representation.** Cosine similarity of raw pooled
ESM vectors, `meta_val` targets against 399 `meta_train` targets:

| representation | mean | std | 1st-99th pct | mean spread across the top 16 |
|---|---:|---:|---|---:|
| raw pooled ESM | 0.897 | 0.050 | 0.766 - 0.979 | **0.024** |
| train-mean-centred | -0.006 | 0.320 | -0.584 - 0.855 | **0.238** |
| train-whitened | 0.000 | 0.062 | -0.115 - 0.202 | **0.420** |

Every raw similarity lies in a band of width 0.21 around 0.90, and the 16
nearest training targets are separated by 0.024. A `softmax(16 * sim)` over that
band is nearly uniform, so Stage A's `protein_neighbor_esm` was close to a
global mean **by construction**. Centring on a `meta_train`-only mean widens the
usable spread tenfold. Stage A's conclusion "protein retrieval is weak" is
therefore correct **only for raw pooled cosine** and does not generalise to all
protein representations. This is exactly the new evidence Stage R2 must test.

It does **not** reopen Mac-Diff locality, conformer routing, PBCNet2.0 or
Cartesian equivariance: those were rejected on structural-input and multi-seed
training evidence, not on retrieval sharpness.

## Evidence vocabulary used from here on

| grade | meaning |
|---|---|
| **exploratory** | measured once; no preregistration; may use oracle or transductive quantities |
| **development** | preregistered gates, but hyperparameters selected on the same population used for inference |
| **conditional** | resolved, but conditional on a stated unresampled factor (trained seeds) or an uncontrolled population property (ligand overlap) |
| **confirmed** | selection separated from inference by an outer fold, preregistered gates, positive component-level lower bound on a population never used for selection |

Regrading the two headline results:

* Stage A composed retrieval, 25.5% -> **exploratory, transductive upper bound**;
* Stage C blend, 12.3% k=0 -> **development, conditional on 48.9% exact ligand
  overlap**; the exact-free effect is +0.050 [-0.074, +0.175], **unresolved**.

No result in this project is currently **confirmed** under this vocabulary.
