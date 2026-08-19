# Implementation amendments — 2x2 oracle-covered subset diagnostic (2026-08-19)

The frozen preregistration (ee844b2b...) is unchanged: split, seeds,
total update count, loss weights, sampling STREAMS namespace, and gates
are untouched. The following implementation repairs are recorded here
because they change training semantics or reported quantities.

## A1. Absolute-cell sampling key now includes the minibatch index

train2x2.py, joint cells: the L_abs batch is now drawn with
rng_key(cell, "abs", ep, mb) — the contrast-minibatch index mb enters
the SHA-256 keyed stream. Previous launches re-permuted with the SAME
key inside every minibatch and took the first 512 cells, so the same
<=512 absolute cells were consumed ~12x per epoch, inflating the
effective weight of L_abs and biasing the joint-vs-centered comparison.
Fix: unique keyed selection per minibatch (each minibatch sees its own
permutation slice). Split/seed/update count/loss weights unchanged.
(Note: the archived Stage-1 train1a.py used the analogous per-epoch
pattern, but Stage-1 artifacts are frozen and were NOT regenerated.)

## A2. Point estimates are OBSERVED quantities, not bootstrap means

train2x2.py effect_boot: the reported point estimate is now
observed_pair_mean_effect (the frozen primary estimand: pair-mean R2
difference) and observed_parent_mean_effect (parent-mean of within-
parent means), with bootstrap_ci {lo2.5, hi97.5, draws 2000, cluster
parent} reported alongside. The bootstrap mean is explicitly recorded
as bootstrap_mean_not_a_point_estimate and is NOT used for status. The
frozen status rule (established iff lo2.5 > 0 and |observed point| >=
0.05; absent iff |point| < 0.02) now applies to the observed pair-mean
point.

## A3. Leave-one-parent-out sign stability covers all five effects

rep_main_joint, rep_main_centered, obj_main_klifs, obj_main_esm, and
interaction all carry leave_one_parent_out_sign_stable; a parent whose
removal empties a cell's rows is marked "nan" (cannot occur on the
matched data where all four cells share the same six test parents).

## A4. Reporting-stage shape alignment (var_dec)

Ls tensor is aligned to len(rows) in the KLIFS variance-decomposition
block; training objective untouched (reporting only).

## Launch ledger (audit trail)

- launch 1 (job r00): crashed in klifs_joint L_abs torch-vs-numpy
  indexing (pre-epoch-1). Recorded in commands.jsonl.
- launch 2 (job 13j): crashed post-training of klifs_joint in var_dec
  (44 vs 64). Recorded.
- launch 3 (job 4gq): all four cells trained; crashed in effect_boot
  leave-one-parent-out (string-vs-int index bug, A3). Cell aggregates
  in the job log are PREVIEW ONLY and were not reused.
- launch 4 (post A1-A3): full four-cell rerun from scratch.
