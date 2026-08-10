# Preregistration — P1R2B-PHASE2B-S4R

## Single-axis graph-aware ligand representation repair of the real structural direct-W stage

Stage identifier: `P1R2B-PHASE2B-S4R_GRAPH_AWARE_LIGAND_DIRECT_W`

Written 2026-08-10, after `PHASE2B_S4R_REPRESENTATION_AUDIT.json`
(`GRAPH_LIGAND_REPRESENTATION_AVAILABLE_AND_INFORMATIVE`, commit `b6327d4`) and
before any S4R runner, checkpoint, prediction, arm score or Gate value exists.

## 1. Authorization chain

```text
Phase 2A  LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING
S2R       BINARY_ORDINAL_IDENTIFIABILITY_REPAIRED      sealed synthetic 0.6620
S3R       REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED  candidate 0.035880
S4R-A     GRAPH_LIGAND_REPRESENTATION_AVAILABLE_AND_INFORMATIVE
```

S3R passed every numerical, participation, firewall and replay check, so its
failure is not an optimizer, runtime or implementation failure and not proof of
an untrainable estimator. Phase 2A had already shown the labels carry
same-construct ligand conditionality. The unresolved boundary is the
measurement basis.

S4R-A then measured that basis. On 39,435 label-blind heldout-A ligand-graph
pairs the mean-pooled 41-D representation has pair-difference effective rank
`6.183` out of a numerical rank of 33, 687 distinct ligand graphs share a
bit-identical vector, and `85.2%` of the variance of the difference norm is
explained by the heavy-atom-count difference alone. The selected graph-aware
statistic reaches effective rank `20.927`, retains `92.2%` of the baseline
difference energy in its own linear span, and places `35.5%` of its own
difference energy outside anything the baseline difference can linearly
express — `35.4%` after also removing heavy-atom-count and log-count
differences.

## 2. Question

Exactly one axis changes:

> Replacing the global mean of 41-D atom-local features by one frozen
> graph-aware 2D ligand statistic, and changing nothing else, does the real
> ligand-conditioned residue direction become identifiable under closure shift?

This is not a search for a better score. A failure closes the pose-free
representation repair route.

## 3. The single changed axis — frozen

Baseline, exactly as in S3R:

```text
g_base(L) = mean over heavy atoms of the 41-D atom-local one-hot descriptor
D_base    = 41
```

Candidate, selected by the S4R-A capacity-parsimony rule and now immutable:

```text
encoder      rdkit.Chem.rdFingerprintGenerator.GetMorganGenerator(radius=1)
identifiers  GetSparseCountFingerprint(mol).GetNonzeroElements(), unfolded 32-bit
vocabulary   the 128 identifiers with the most distinct TRAIN-split ligand
             graphs, ties by ascending identifier
             sha256 a200a4b986af1850fdb1d244f2e002c9b5ae707a114d8a3635053edb215ed877
g_graph(L)[j] = (count of vocabulary environment j in L) / (heavy atoms of L)
D_graph    = 128
```

Both statistics are an arithmetic mean over heavy atoms of a per-atom
descriptor, both are invariant to atom permutation, and both are normalized by
heavy-atom count. The only difference is that the per-atom descriptor is
graph-aware: a radius-1 Morgan identifier is a hash of an atom's own invariants
together with the multiset of its bonded neighbours, so it encodes bond
connectivity and local functional-group arrangement that an atom marginal
cannot.

The vocabulary is derived from train-split ligand structures only. It reads no
label. Heldout-A environments outside the vocabulary contribute zero; the
audit measured heldout-A coverage `0.999975`.

No graph network is introduced, nothing is trained inside the encoder, and the
encoder is byte-frozen by the SHA-256 above.

## 4. Everything else is held fixed

Byte-identical to S3R:

- frozen ESM2-650M residue states `H_P`, never trained;
- protein nuisance basis `Q_P = onb{1, b^P}` over the frozen `b_prior`,
  float64 Gram-Schmidt, tolerance `1e-10`;
- one direct matrix `W`, no bias, projected to `||W||_F = 1` after every step;
- the S1R all-residue bidirectional pairwise logistic loss on the RMS-normalized
  score;
- pair -> construct -> closure-component -> batch aggregation, and the same
  hierarchical sampler `C_MAX=2`, `P_MAX=8`, `BATCH_COMPONENTS=16`;
- the frozen closure split `make_split`, its train / heldout-A definition and
  its component and ligand-graph firewall;
- Adam, `lr=1e-3`, no weight decay, gradient clip `5.0`, exactly 210 updates;
- seeds `SEED_PARAM=20260901`, `SEED_SAMPLER=20260902`, `SEED_BOOT=20260903`;
- the frozen control map file `phase2b/control_maps.json`, sha256
  `e187a5f00f0b66328877bacd93b22471fe607e382e811f2674ecfc4a9dec9c33`;
- R1-R5 controls and margins, unchanged in both definition and value.

No learning rate, budget, sampler, margin, threshold, seed, loss or
representation search is permitted at any point.

### 4.1 Stream identity contract

The training stream is a deterministic function of the train pair list, the
construct-to-component map and `EPOCHS`, none of which this stage changes. The
runner must therefore reproduce the S3R stream exactly:

```text
semantic_sha256 == 4bc68d54884437ded999fbb5f8fc8997b47b456f2bdde4e0fbafd9df3dcdc3ef
updates          == 210
```

Any mismatch is `S4R_CONTRACT_OR_LABEL_FIREWALL_FAIL_CLOSED`.

### 4.2 Baseline replication contract

The `baseline41` arm is the S3R candidate re-executed inside S4R. With an
identical stream, seed, split, loss and 41-D representation it must reproduce
the S3R result exactly:

```text
|macro AP_bidir(baseline41) - 0.03588006089257408| <= 1e-12
```

This is a cross-stage determinism anchor, not a Gate. A mismatch means some
surface other than the ligand representation moved and is
`S4R_CONTRACT_OR_LABEL_FIREWALL_FAIL_CLOSED`.

## 5. Trainable object and capacity

```text
W_graph in R^(1280 x 128)   163,840 parameters
W_base  in R^(1280 x 41)     52,480 parameters
```

The candidate uses `3.125x` the baseline's parameters, which is the smallest
increase among the six audited representations that cleared every information
gate. Capacity is not free, and it is not controlled by argument: R5 trains an
identical `1280 x 128` learner on permuted labels through the identical stream,
so any score reachable by capacity alone appears there.

## 6. Data, splits and firewall

- Primary development panel: the S3R heldout-A eligible pair set, which is
  frozen by `build_pairs` — unordered, within one exact construct,
  scaffold-distinct, different ligand graph, non-empty residue symmetric
  difference. Expected census `46,818` pairs and `112` closure components, and
  train `226,765` pairs over `554` components. Any deviation fails closed.
- **Heldout-B is not created, not written and not read by this stage.** Only
  `train` and `heldoutA` residue label views exist.
- Heldout-A has already been consumed by S3R. It is development evidence only.
  A pass here authorizes a separately registered independent structural
  confirmation and nothing else.
- The training process may open `train_residue_masks.json.gz` and no other
  label view. The scoring process may open `heldoutA_residue_masks.json.gz` and
  no other. Both are enforced by an opened-path log that fails closed.
- Zero reads of affinity values, ChEMBL, BindingDB, DAVIS, KIBA, recipient,
  metaval, R6 or biological `z`.

## 7. Arms — frozen

Trained, each through the identical 210-update stream:

```text
candidate    g_graph, true train labels
repeat       g_graph, true train labels, independent process, determinism
permuted     g_graph, frozen within-construct derangement of the train labels
baseline41   g_base,  true train labels
```

Inference-only, all evaluated on exactly the same pair, construct and component
mask as the candidate:

```text
b5diff        frozen B5 residue prior differential, no ligand branch
foreign       candidate W with the frozen foreign-ligand pair map
context       candidate W with within-amino-acid-type residue context shuffling
ligand_only   candidate W with each protein's residue states replaced by their
              residue-mean, so the score carries ligand information but no
              residue-resolved protein context
zero_W        the zero matrix, the analytic chance arm
chem_shuffle  candidate W with the frozen within-construct ligand derangement,
              secondary diagnostic, non-gating
```

`ligand_only` is analytically degenerate: a residue-constant score lies inside
`span{1} subset span{Q_P}`, so the gauge projection annihilates it exactly.
Registering it is a structural proof that no ligand-only shortcut can survive
the estimand, and its AP must equal the chance arm within `5e-3`.

## 8. Gates — frozen

Every Gate is a paired closure-component bootstrap, `10,000` resamples, seed
`20260903`, one-sided 95% lower bound. A Gate passes only when the observed
delta reaches its margin **and** the lower bound is above zero.

| Gate | contrast | margin |
|---|---|---:|
| R1 | candidate - chance | +0.05 |
| R2 | candidate - frozen B5 differential | +0.03 |
| R3 | candidate - foreign ligand pair | +0.03 |
| R3b | candidate - ligand-only | +0.03 |
| R4 | candidate - residue-context corruption | +0.03 |
| R5 | candidate - trained permuted-label learner | +0.05 |

R1-R5 are the S3R Gates with identical contrasts and identical margins. R3b is
the ligand-only control required by this stage's authorization; it joins the
shortcut family of R3 and R4 and takes that family's margin.

Reported with a paired bootstrap and explicitly **non-gating**:

```text
C1  candidate - baseline41       the direct value of the changed axis
C2  baseline41 - chance          the S3R effect, re-measured in this stage
```

C1 is not a Gate because the stage's authorization fixes R1-R5, and because a
candidate that beat the baseline while failing R1 would still not identify a
residue direction. C1 is nevertheless the quantity that determines what the
next research step may be, and it is reported with its full interval.

## 9. Module participation and replay — preconditions

Interpreted **before** any biological Gate:

1. gradient of `W` finite and non-zero at every update;
2. relative `W` movement `>= 0.05`;
3. `| ||W||_F - 1 | <= 1e-5` at every step and at load;
4. heldout raw score variance `>= 1e-8`;
5. `zero_W` macro AP equal to the analytic chance macro within `5e-3`;
6. `ligand_only` macro AP equal to the chance macro within `5e-3`;
7. `context` and `foreign` macro AP strictly below the candidate;
8. candidate, repeat, permuted and baseline41 share one stream semantic hash;
9. repeat reproduces `W`, every raw prediction and every per-pair metric to
   `<= 1e-7`, and the two prediction files have equal SHA-256;
10. every arm scores the identical pair set with the identical
    pair-to-construct-to-component map, verified by one common-mask SHA-256.

Any failure is `GRAPH_LIGAND_REPRESENTATION_TRAINING_FAILED` and the biological
Gates are not interpreted.

## 10. Terminal verdicts

Exactly one, by earliest failed boundary:

```text
S4R_CONTRACT_OR_LABEL_FIREWALL_FAIL_CLOSED
GRAPH_LIGAND_REPRESENTATION_TRAINING_FAILED
REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED          R1 or R2 fails
GRAPH_LIGAND_STATISTIC_SHORTCUT_DEPENDENT            R3, R3b, R4 or R5 fails
GRAPH_AWARE_RESIDUE_DIRECTION_IDENTIFIED_IN_DEVELOPMENT
```

`GRAPH_LIGAND_REPRESENTATION_NOT_INFORMATIVE` was already resolved against by
S4R-A and is not reachable here.

## 11. Stopping rules

- One run. No seed, budget, capacity, representation or threshold may be
  changed after any S4R metric is read, and no failed Gate may be rescued.
- `REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED` closes the pose-free
  representation repair route. It does not authorize attention, a larger PLM,
  a parallel branch, a second protein encoder, geometry, pose, typed channels,
  affinity supervision, knowledge graphs, PU learning or few-shot adaptation.
- `GRAPH_AWARE_RESIDUE_DIRECTION_IDENTIFIED_IN_DEVELOPMENT` authorizes exactly
  one thing: writing a separate preregistration for an independent structural
  confirmation on a panel that is not heldout-A. It does not open heldout-B in
  this stage, does not open affinity, does not admit anything to `z` and does
  not authorize production integration.
- Heldout-B, R6 and affinity remain closed under every outcome of this stage.

## 12. Mathematical and biological boundary

The graph-aware statistic is an upstream biological measurement. It may not
enter `z` because it improves a residue AP. The operator

```text
A(F, z) = K(B(z) F(z))
```

is unchanged, and so are CSMO, Band, the positive ridge, the simplex and the
fixed mesh. This stage identifies at most a bounded binary residue-ranking
direction under closure shift. It does not identify residue-atom coupling,
interaction energy, affinity, selectivity, a few-shot section or a validated
end-to-end DTA model.

If the candidate ever passes independent confirmation, the remaining sequence
is unchanged:

```text
independent structural confirmation
  -> source-affinity correct > ligand-only and correct > wrong-protein
  -> support-rank / coverage identifiability
  -> bounded biological z admission
  -> unchanged probability-law operator
```
