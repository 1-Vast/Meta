# Active experiment protocol

## Status

```text
S7/L2B Phase 1 B5 ................. COMPLETE, DEVELOPMENT PASS 6/6
Phase 2A teacher audit ............ COMPLETE
S2R synthetic direct-W witness .... COMPLETE, PASS
S3R real structural direct-W ...... COMPLETE, FAIL AT R1
S4R-A ligand representation audit . COMPLETE, INFORMATIVE
S4R graph-aware ligand repair ..... COMPLETE, FAIL AT R1
S5D estimand/collapse diagnostic .. COMPLETE, MECHANISM FALSIFIED, E1-E3 FAIL
C0 untouched corpus and closure ... COMPLETE, ALL GATES PASS
C1 exact-coupling information ..... COMPLETE, FAIL AT C1a
C2 correspondence router .......... NOT PREREGISTERED, NOT TRAINED
heldout-B / R6 / affinity ......... NOT OPENED
active training stage ............. NONE
```

## S2R

The gauge-free direct matrix estimator passed three calibration seeds and one
sealed synthetic seed. Sealed held-out component-macro `AP_bidir = 0.6620`.
This established trainability only; it did not establish biology.

## S3R

Frozen inputs: ESM2 residue states and mean-pooled 41-D ligand atom features.
Trainable object: one `1280 x 41` direct matrix, unit Frobenius norm. Training:
210 fixed Adam updates, hierarchical pair/construct/component aggregation, no
hyperparameter or seed selection. All five Gates failed; terminal verdict
`REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED`. Module participation and
deterministic replay passed, so the failure was scoped to the measurement
basis.

## S4R-A

Label-blind representation audit, zero residue label reads. The mean-pooled
41-D basis has pair-difference effective rank `6.183`, 687 distinct ligand
graphs share a bit-identical vector, and `85.2%` of the difference-norm
variance is heavy-atom-count difference. All six audited Morgan candidates
cleared every A-gate; the capacity-parsimony rule selected radius 1, `d = 128`.
Verdict `GRAPH_LIGAND_REPRESENTATION_AVAILABLE_AND_INFORMATIVE`.

## S4R

One axis changed: the ligand statistic. Baseline `g_base` is the 41-D atom
marginal mean; candidate `g_graph` is the per-heavy-atom count of the 128 most
frequent train radius-1 Morgan environments. Identical protein branch, gauge,
estimator, loss, sampler, split, seeds, control maps and 210-update stream —
verified by three cross-stage anchors, including a bit-exact reproduction of
the S3R candidate by the `baseline41` arm.

Primary panel: 46,818 pairs, 112 closure components, one common mask.

```text
candidate     0.046856      baseline41    0.035880
b5diff        0.031582      foreign       0.046212
context       0.027357      ligand_only   0.025472
permuted      0.036293      zero_W        0.025472
chem_shuffle  0.051322      (non-gating)
```

```text
R1  candidate - chance       +0.021384  [LCB +0.016064]  FAIL (< +0.05)
R2  candidate - B5           +0.015273  [LCB +0.008173]  FAIL (< +0.03)
R3  candidate - foreign      +0.000644  [LCB -0.009226]  FAIL (< +0.03)
R3b candidate - ligand-only  +0.021384  [LCB +0.016064]  FAIL (< +0.03)
R4  candidate - context      +0.019498  [LCB +0.013599]  FAIL (< +0.03)
R5  candidate - permuted     +0.010563  [LCB +0.003880]  FAIL (< +0.05)
C1  candidate - baseline41   +0.010976  [LCB +0.004939]  non-gating
C2  baseline41 - chance      +0.010408  [LCB +0.006920]  non-gating, = S3R R1
```

Terminal verdict: `REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED`.

The representation change was real — the above-chance gain doubled and the
candidate now beats its capacity-matched permuted-label learner, which S3R did
not. It is nevertheless not ligand-conditioned: foreign ligand pairs are almost
free and the chemistry shuffle scores higher, so the learned direction is a
construct-level residue-change prior. `W`'s leading singular values are nearly
flat (`0.209, 0.189, 0.187, 0.177, 0.170, ...`), consistent with a ligand
argument that barely steers the residue direction.

Module participation and deterministic replay passed. Heldout-B was neither
created nor read. R6 was not opened. Affinity reads were zero. No production
code or frozen mathematics was modified.

## S5D

No training, zero new parameters, frozen S4R checkpoints reused. Registered
mechanism for R3: the estimator collapses ligand differences onto about one
residue direction per protein. **Falsified.**

```text
rho_dg    data-side upper bound       0.4550
rho_graph candidate residue fields    0.4793   excess over data  0.0138
rho_base  baseline41 residue fields   0.5758
true-vs-foreign field cosine          0.4487   over 46,817 pairs
rule required median >= 0.80 and excess >= 0.10
```

The estimator steers on the ligand, and the graph statistic produces more
diverse fields than the baseline. D2 then measured the symmetric-difference
conditional estimand, which cancels pocket membership exactly: 40,157 eligible
pairs, 107 components, median 7 changed residues, median gain fraction 0.50.

```text
candidate 0.655030   foreign 0.655470   chance 0.643744
baseline41 0.638830  permuted 0.628586
E1 candidate - chance    +0.011285  [LCB -0.007749]  FAIL (< +0.05)
E2 candidate - foreign   -0.000440  [LCB -0.021814]  FAIL (< +0.03)
E3 candidate - permuted  +0.026444  [LCB -0.002977]  FAIL (< +0.03)
E4 candidate - baseline  +0.016199  [LCB -0.005947]  non-gating
E5 baseline41 - chance   -0.004914  [LCB -0.023844]  non-gating
```

Terminal verdict: `LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED`; the D2 Gates fail
independently. Ligand information arrives at the residue field intact and
points somewhere biologically wrong under both estimands. Heldout-A has now
been consumed three times and no fourth estimand variant on it is permitted.

## C0/C1

Audit only, zero parameters, on a corpus no stage had touched: 24,874 exposed
PDB ids excluded, 2,836 untouched local mmCIF entries, 2,039 admissible systems
and 1,862 scored after the CCD-scaffold rule. P1B semantics are respected —
`contact_prob(i,s)` is "any residue in slot `s` contacts atom `i`", never
additive mass, and multiple residues in a slot may contact one atom.

The registered `M4` mapping rule failed its own fail-closed check at `23/40`.
P1B's sequence comes from BioLiP column 20, not the mmCIF entity sequence.
Amendment 01 corrected the rule to the true P1B path before any statistic was
read; the check then passed `60/60`. The C1 run started under the rejected
mapping was stopped and discarded unread.

```text
G0a components                 496      >= 60     PASS
G0b largest fraction        0.0811     <= 0.25    PASS
G0c minimum detectable eff 0.00453     <= 0.05    PASS   (null sigma 0.04053)
```

The union closure gave 89 components but exceeded the giant-component cap, so
the registered DataSAIL-style two-dimensional fallback was used. The 3-mer
prefilter was measured, not trusted: 3,037 true identity edges, 0 missed.

```text
within-slot AP empirical               0.985611
within-slot AP fixed-degree rewire     0.953959
within-slot AP atom shuffle            0.985611   exact no-op, degenerate
within-slot AP geometry shuffle        0.993948   shuffling makes it EASIER
complete-edge AP, additive marginal    0.753449
C1a empirical - rewire  +0.031652 [LCB +0.029690]  needs +0.05   FAIL
C1b 162,276 positive units / 100,563 checkerboards               PASS
C1c replicate jaccard 0.830 over 17 cross-entry pairs            PASS
```

Terminal verdict: `EXACT_EDGE_COUPLING_NOT_SUPPORTED_BY_TEACHER`. The empirical
AP of `0.9856` leaves only `0.0144` of headroom above a pure contact-degree
predictor, so the `+0.05` margin is unreachable in principle on this statistic.
At `6.0 A` a slot holds about three sequence-adjacent and therefore spatially
adjacent residues, and there is almost nothing left to deconvolve. The C2
router was not preregistered and not trained.
