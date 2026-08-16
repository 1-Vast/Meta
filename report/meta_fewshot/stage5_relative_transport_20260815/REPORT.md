# Stage 5 result: signed relative transport is rejected, and why

Numerical authority: `GEOMETRY_COVERAGE_AUDIT.json`, `VAL_*.json`,
`RANKING_meta_val.json`, `RANKING_stage4_on_meta_val.json`, `arm_*/RESULT.json`.
Gates were fixed in `STAGE2_PREREGISTRATION.md` before arms B, C and D produced
data.

## 1. Geometry: closed on data, not on preference

| quantity | value |
|---|---:|
| governed holo complexes (raw mmCIF on disk) | 14,906 |
| DTA targets with an exact holo sequence | 15 / 499 |
| DTA targets with a containment-level holo sequence | 110 / 499 |
| DTA ligands sharing a holo SMILES | 84 / 9,880 |
| **DTA cells with a common-frame protein-ligand complex** | **0 / 17,717** |

The 15 sequence-matched targets carry structures bound to *different* ligands.
`pilot20k_structure_supervision_v2` and `r0b_exact_geometry_v3` store contact
maps and distance bins — invariant summaries, **not coordinates**.

Every Cartesian and equivariant family (PBCNet2.0, TensorNet, PaiNN, MACE,
Equiformer, E2Former, SE(3)-EGNN) fails on the same input constraint, so the
choice among them is moot here. `model/cartesian.py` remains verified by
`tests/test_cartesian.py` (O(3) with reflection, translation, symmetric
traceless rank-2, permutation, padding, gradients, cross-sample edge rejection)
and unused. The active model raises on coordinate inputs rather than accepting
them silently. Reopening this needs **new data**, not a new architecture.

## 2. The mechanism that was tested

PBCNet2.0's Siamese *relative* idea with the geometry removed, plus AdaMBind's
task-difficulty idea as a leave-one-out label-consistency credit, with no MAML,
no inner loop and no test-time gradients:

```text
f = f0(q) + s(n) * sum_k w_qk * [ r_k + delta(P, L_k -> Lq) ]
delta(a->b) = m(e_a,e_b) - m(e_b,e_a)          exactly antisymmetric
w_qk        = softmax_k( tau <key_q,key_k> + c_k ),  c_k = -|r_k - LOO_k| / kappa
```

17 Stage 1 gates pass, including exact antisymmetry, `delta(a->a) = 0`, k=0
identity, support-permutation invariance, query equivariance, label-permutation
sensitivity, padding invariance, level-only nesting, signed-effect recovery at
k=1, private-mechanism rejection, and no dead trainable tensor. The antisymmetric
readout bias was removed because it cancels identically and is unidentifiable.

**The mechanism is algebraically sound. It does not work on real data.**

## 3. Stage 2 result, `meta_val`, 44 episodes, matched seed and budget

| arm | k0 | k1 | k2 | k3 | k5 |
|---|---:|---:|---:|---:|---:|
| A `grammar` | 1.857 | 1.415 | 1.198 | 1.201 | 1.184 |
| B `relative` | **1.599** | **1.312** | **1.168** | 1.199 | **1.143** |
| C `relative` + difference loss | 1.753 | 1.438 | 1.228 | 1.251 | 1.177 |
| D `relative` no reliability + diff loss | 1.749 | 1.405 | 1.208 | 1.182 | 1.146 |

`full` minus `level_only` (positive means the query-specific channel helps):

| arm | k=1 | k=2 | k=3 | k=5 | permutation gap at k=2 |
|---|---:|---:|---:|---:|---:|
| A `grammar` | **+0.095** | **+0.105** | **+0.067** | **+0.017** | **+0.163** |
| B `relative` | +0.001 | -0.002 | -0.046 | -0.019 | +0.000 |
| C `relative`+diff | -0.000 | -0.000 | -0.052 | -0.018 | +0.000 |
| D no reliability | -0.008 | -0.003 | -0.003 | -0.001 | +0.000 |

**Gate outcome: S1 fails for every relative arm, and S2 fails at k>=2.** The
difference operator collapses to `delta ~ 0` with flat weights — exactly the
designed safe floor. The permutation gap is **identically zero** at k>=2, which
is the algebraic signature of uniform weights and a null operator: with flat `w`
and `delta = 0` the transport is `mean(r)`, which is permutation invariant.

Neither explicit per-(query, support) difference supervision (arm C) nor
removing the reliability credit (arm D) changed this. Per the preregistered
continuation rule, the signed relative transport is **rejected**.

Arm B has the best absolute MSE of the four, but its `full` equals its
`level_only`: that improvement is a trunk effect from a quieter transport
gradient, not a mechanism effect. It is not evidence for the mechanism.

## 4. Why both mechanisms failed — the useful finding

Ranking metrics on the *same* `meta_val` episodes, `full` against `level_only`:

| checkpoint | budget | k=2 CI full / level | k=3 CI full / level | k=5 CI full / level |
|---|---:|---|---|---|
| Stage 5 arm A | 800 steps | 0.5757 / 0.5700 | 0.5912 / 0.5700 | 0.5826 / 0.5700 |
| Stage 4 seed 20260812 | 2000 steps | 0.5539 / 0.5715 | 0.5606 / 0.5715 | 0.5621 / 0.5715 |
| Stage 4 seed 20260813 | 2000 steps | 0.5674 / 0.5601 | 0.5416 / 0.5601 | 0.5496 / 0.5601 |
| Stage 4 seed 20260814 | 2000 steps | 0.5993 / 0.5914 | 0.5828 / 0.5914 | 0.5861 / 0.5914 |

At 800 steps the gate **improves** the concordance index at k=2, 3 and 5. At
2000 steps it **degrades** it in 9 of 12 (seed, k) cells, on the same split, with
the same architecture. The earlier Stage 4 conclusion was therefore not a
`meta_test` peculiarity: it is what more optimization does to this mechanism.

The mechanistic explanation, supported by both stages:

* the transport is trained on squared error;
* the MSE-optimal use of `k` noisy support residuals is shrinkage toward their
  mean, because a per-query reweighting of noise adds variance;
* a level shift is constant across queries and therefore cannot change ranking;
* so gradient descent buys MSE by increasing shrinkage and pays for it in
  within-target discrimination, and it buys **more** of that trade the longer it
  trains.

The signed difference operator does not escape this. It is also trained on
squared error, and its MSE optimum is `delta = 0` — the same shrinkage solution
expressed in a different parameterisation. That is why it went inert rather than
going wrong.

## 5. What this implies for the next hypothesis

The blocker is the **objective**, not the operator. A query-specific channel
trained to minimise squared error on cold-target episodes will converge to level
calibration whatever its functional form. The next candidate must make ranking
the primary training signal for the transport — the existing
`pairwise_ranking_loss` carries weight 0.5 against a dominant MSE term and is
applied to the whole prediction rather than to the query-specific component.

Two further constraints established here:

* adding supervision to a failed mechanism did not rescue it (arm C was worse
  than arm B at k=0, 1, 2, 3);
* the label-consistency credit produced no label binding when the operator it
  depends on is null; it cannot be evaluated independently of a working
  difference operator.

## 6. Disposition

* `model/relative_grammar.py`, `tests/test_relative_grammar_synthetic.py` and
  these results are retained as rejected-candidate evidence. The module is
  opt-in (`--arch relative`), the trainer default remains `bpsf`, and the
  `grammar` comparator is byte-identical to Stage 3/4.
* No admission, performance or mechanism claim is made for the relative
  transport.
* No Cartesian equivariance, atomic-level recognition or binding-mode claim is
  made anywhere: the input does not support one.
