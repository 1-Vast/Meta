# Stage 2 preregistration: relative-transport discriminator

Written after arm A produced data and **before** arms B, C and D finished.
Arm A's numbers are quoted below as the comparator; the gates are stated in
terms of it.

## Arms

Matched seed 20260812, 800 steps, 4 episodes/step, hidden 384 / embed 192 /
48 contact types / 5 ligand layers, cosine schedule, lr 6e-4, backbone scale
1.0, binding weight 1.0. The trunk, encoders and zero-shot endpoint are
identical in every arm; only the transport differs.

| arm | `--arch` | `--difference-loss-weight` | isolates |
|---|---|---:|---|
| A | `grammar` | 0.0 | Stage 4 comparator (bounded rescaling gate) |
| B | `relative` | 0.0 | the signed difference operator alone |
| C | `relative` | 0.5 | + explicit per-(query, support) difference supervision |
| D | `relative_noreliability` | 0.5 | C minus the leave-one-out label-consistency credit |

## Evaluation

**`meta_val`, all 44 eligible episodes.** The consumed `meta_test` split is not
used for any Stage 2 or Stage 3 decision. Checkpoint selection also used
`meta_val` (2 targets per component), so these numbers are optimistically
biased; the bias is identical across arms, which is what a discriminator needs.
Ranking metrics come from `scripts/evaluate_arms_ranking.py`, run uniformly on
all four saved checkpoints after the chain completes.

## Arm A comparator, `meta_val`, 44 episodes

| k | full | zero-shot | level / sar-cut | permuted |
|---|---:|---:|---:|---:|
| 0 | 1.857 | 1.857 | 1.857 | — |
| 1 | 1.415 | 1.857 | 1.509 | 2.886 |
| 2 | 1.198 | 1.857 | 1.303 | 1.361 |
| 3 | 1.201 | 1.857 | 1.269 | 1.311 |
| 5 | 1.184 | 1.857 | 1.202 | 1.266 |

Note that on `meta_val` the permutation control is correctly signed for arm A at
every k, unlike on the consumed `meta_test` bank where it inverted at k=2, 3 and
5. Split-level behaviour differs, and no Stage 2 conclusion may be transferred
to `meta_test` without re-measurement.

## Gates

| id | requirement |
|---|---|
| S1 | `full` beats `level_only` at every k in {1,2,3,5} by a margin at least as large as arm A's (0.094 / 0.105 / 0.068 / 0.018) |
| S2 | permuted-support MSE exceeds `full` at every k in {1,2,3,5} |
| S3 | **CI and Spearman of `full` are not below those of `level_only` at any k** |
| S4 | k=0 within 0.05 of arm A (the trunk is unchanged; larger drift means the transport is disturbing trunk training) |
| S5 | `full` MSE not worse than arm A at any k in {1,2,3,5} |

S3 is decisive: it is the exact gate the Stage 4 admission run failed, and the
whole reason this mechanism was designed. A candidate that improves MSE while
failing S3 is rejected, because that is the shrink-toward-the-level pathology
the previous transport already exhibited.

## Continuation rule

Promote at most one relative arm to Stage 3. Prefer the simplest arm that
passes: B over C over D. If S3 fails for every relative arm, the signed
difference operator is rejected, the source change is reverted, and the
`grammar` trunk stands as the retained candidate with its level-calibration
attribution unchanged.
