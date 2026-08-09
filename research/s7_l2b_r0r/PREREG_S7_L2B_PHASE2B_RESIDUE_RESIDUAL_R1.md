# Preregistration — S7/L2B Phase 2B **R1**

## One ligand-conditioned residue residual over a frozen protein-only pocket prior

Stage identifier: `S7_L2B_PHASE2B_LIGAND_CONDITIONED_RESIDUE_RESIDUAL_R1`

Written: 2026-08-10. Repository commit at registration: `0bd1702`.

**Supersedes** `PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL.md`
(SHA-256 `ae6d1a0186bb37af86f3b6eb98c513bce7e67a8745aaf5a3811ce5c9b98ab477`),
which is retained byte-identical and marked
`SUPERSEDED_BEFORE_EXECUTION_DESIGN_DEFECT`. That document was never executed and
produced no result; the defects are enumerated in `PHASE2B_DESIGN_AUDIT.md`
(D1–D11).

**Authorizing evidence.** Phase 2A verdict
`LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING`,
`report/s7_l2b_r0r/PHASE2A_VERDICT.json`, evidence closure committed as
`0bd1702`.

This registration is frozen and committed **before** any Phase 2B
implementation file is written and before any Phase 2B number exists. The
Phase 2A chronology is unaffected by that commit and is not improved by it.

---

## 1. Question

Does the existing **frozen** sequence + 2-D ligand representation support a
ligand-conditioned residue-level correction beyond the generic protein-pocket
prior, on unseen protein components and unseen ligand graphs?

This is a structural auxiliary experiment. It is not an affinity experiment, it
reads no affinity value, and it changes nothing in `model/`, `scripts/`,
production `z`, CSMO, Band, mesh or `A(F,z) = K(B(z)F(z))`.

## 2. Frozen protein-only prior

Taken directly from the frozen B5 checkpoint
(`dataset/processed/s7_l2b_r0r/sealed_preds_b5/B5_checkpoint.pt`, whose SHA-256
is recorded in `P1_B5_GATE.json`), using only its protein branch:

```text
r_r(P)     = GELU( W_h ( LayerNorm( h_r ) ) )
b_r^P(P)   = b + alpha * w_pi( r_r(P) )
```

`h_r` is the frozen ESM2-650M final-layer residue state already extracted and
hashed. The expression contains **no atom and no ligand term**, so ligand
independence is structural, not asserted.

Frozen properties, all verified before training:

- computed for **every** train and held-out construct;
- materialised to disk as **float64** and hashed;
- gauge: raw logit units, no centering, no rescaling, no clipping;
- B5 remains in `eval()` mode with `requires_grad=False` on every parameter.

A label-free diagnostic reports the Pearson and Spearman correlation between
`b^P` and the Phase 2A per-complex additive residue marginal `alpha_r`. **No
model choice, threshold, gate or arm may depend on that diagnostic.** It exists
only to document how far the superseded definition was from this one.

## 3. The only trainable object

```text
g(L)              = mean over the existing deterministic 41-D atom features
delta_raw_r(P,L)  = sum_{k=1..8} ( U h_r )_k * ( V g(L) )_k
```

| item | frozen value |
|---|---|
| rank `K` | **8**, exactly |
| `U` | `R^{8 x 1280}`, no bias |
| `V` | `R^{8 x 41}`, no bias |
| new parameters | **10,568**, and no others |
| initialisation | Xavier uniform, torch seed `20260901` |
| dtype | float32 forward, float64 projection |

Frozen and not trained: ESM2, the B5 checkpoint, the 41-D atom features, the
residue features, the split, the closure map. Forbidden and absent: any
additional PLM, GNN, attention stack, geometry branch, typed-interaction branch,
affinity head, PU loss, knowledge graph, parallel SSL module, ligand encoder
replacement, or any further trainable tensor.

No hyperparameter search. No seed sweep. No model selection.

## 4. Residual gauge

Per protein `P`, build an orthonormal nuisance basis by modified Gram–Schmidt in
**float64** over the columns `[ 1 , b^P(P) ]`, in that order, dropping any column
whose residual norm after orthogonalisation falls below `1e-10` relative to its
original norm. This handles a constant or numerically degenerate `b^P(P)`.

```text
delta(P,L) = ( I - Q_P Q_P^T ) delta_raw(P,L)
```

`Q_P` depends only on `P`, so it commutes with the same-protein difference and
cannot carry ligand information.

Fail-closed requirement, checked on every evaluated protein:

```text
|| Q_P^T delta || / ( 1e-30 + || delta || )  <=  1e-8
```

Violation is a contract error terminating the stage, not a reported number.

## 5. One consistent score

```text
s_r(P,L)            = b_r^P(P) + delta_r(P,L)
Delta s_r(P;La,Lb)  = s_r(P,La) - s_r(P,Lb)  ==  delta_r(P,La) - delta_r(P,Lb)
```

Because `b^P` is ligand-independent it must cancel **exactly**. This is verified
numerically to `<= 1e-12` in absolute value on a fixed sample of held-out pairs,
fail-closed.

## 6. Eligible pairs

An eligible pair is an **unordered** pair `{L_a, L_b}` of records sharing one
exact construct (`seq_key`, so residue indices are identical by construction),
with different `graph_key`, both Murcko scaffolds non-empty and different, and
non-empty symmetric difference `R_a Δ R_b`.

`(L_a, L_b)` and `(L_b, L_a)` are **one** dependent observation. Gain and loss
are scored inside that one pair, never as two.

Pairs with `|R_a Δ R_b| = 0` are excluded by this frozen rule and **counted**.

Expected census, to be **verified and reported exactly, not assumed**:

```text
train      ~766 constructs, ~554 closure components, ~226,765 eligible pairs
held-out A ~175 constructs, ~112 closure components, ~46,818 eligible pairs
```

Any material departure from these figures is reported as an exclusion census; it
does not by itself terminate the stage.

## 7. Metric — over all aligned residues

For each eligible pair, over **all** `L` residues of the construct:

```text
AP_gain    score  Delta s      positives  R_a \ R_b
AP_loss    score -Delta s      positives  R_b \ R_a
AP_change  score |Delta s|     positives  R_a symmetric-difference R_b
```

The **primary** statistic is the bidirectional differential AP,
`AP_bidir = ( AP_gain + AP_loss ) / 2`, computed inside the unordered pair.

Chance level per pair is the exact expectation of AP under a constant score
(one tied block), computed with the same closed-form tie-aware estimator used in
Phase 2A (amendments 01/02); `chance_bidir` is the mean of the gain and loss
chance levels.

Aggregation order, frozen:

```text
residue -> unordered ligand pair -> construct -> protein closure component
        -> component macro
```

Inference unit is the protein closure component. Residues, pairs and constructs
are never inference units. Intervals: one-sided 95% lower bounds from
**10,000** paired closure-component bootstrap resamples, seed `20260903`.

The symmetric-difference-restricted AP of the superseded design is reported only
as a **secondary conditional sign diagnostic** and may never be quoted as
deployment evidence.

## 8. Training contract

| item | frozen value |
|---|---|
| optimizer | AdamW |
| learning rate | `1e-3` |
| weight decay | `1e-4` |
| epochs | 6 |
| gradient clipping | 5.0 |
| checkpoint | final epoch only |
| model selection | none |
| early stopping | forbidden |
| parameter seed | `20260901` |
| pair-sampling seed | `20260902` |
| bootstrap seed | `20260903` |
| control-map seed | `20260904` |
| synthetic seed | `20260905` |

**Loss.** For one pair, over all `L` residues, with
`p_r = sigmoid( Delta s_r )`:

```text
L_pair = mean over the non-empty groups of
           mean_{r in R_a \ R_b}            BCE( p_r , 1.0 )
           mean_{r in R_b \ R_a}            BCE( p_r , 0.0 )
           mean_{r not in R_a xor R_b}      BCE( p_r , 0.5 )
```

The soft target `0.5` on unchanged residues is the statement `Delta s_r = 0`
there. Group means make the objective invariant to pocket size and to the
overwhelming majority of unchanged residues, and they match the all-residue
evaluation.

**Hierarchical balance.** The objective aggregates
`pair -> construct -> component -> batch`, so a construct with 318 ligands
cannot outweigh one with 2.

**Deterministic hierarchical sampler.** Per epoch, every training closure
component is visited; within a component at most `2` constructs are sampled;
within a construct at most `8` eligible pairs are sampled; batches contain `16`
components. Selection uses seed `20260902 + epoch`. The selected
`(epoch, component, construct, pair)` identifiers are **materialised and
hashed**.

## 9. Controls — materialised and hashed before training

**`R3` two-ligand foreign-pair control.** Both ligands are replaced. For each
eligible held-out pair `{L_a, L_b}`, a foreign pair `{L'_a, L'_b}` is drawn from
the **training** ligand pool such that all four `graph_key`s are distinct, all
four scaffolds are distinct, `L'_a` and `L'_b` differ from each other, and each
is matched to its counterpart by nearest heavy-atom count with the pooled-feature
Euclidean distance as tiebreak. No fixed points. Selection uses no label and no
score. One map, shared by every compared arm, seed `20260904`.

**`R4` residue-context control.** Frozen contextual ESM2 residue states are
shuffled **within the same amino-acid type within the same protein sequence**.
Residue composition and the per-type state distribution are preserved; the
contextual assignment is destroyed. The corrupted `h` is used consistently for
`b^P`, `Q_P` and `delta`. A length-matched, closure-disjoint, score-blind
wrong-protein arm is reported as secondary only.

**`R5` trained permutation control.** The identical `U`/`V` architecture,
optimizer, sampler, budget and seeds, trained on one frozen within-construct
ligand-label derangement. The derangement is a single frozen mapping (no fixed
points where the construct size permits; constructs of size 1 are not
permutable and are excluded and counted), seed `20260904`, persisted and hashed,
and it is **not changed after any metric is seen**.

**`RANDOM_FEATURE_NULL`.** The identical architecture at initialisation,
untrained. A numerical sanity check only. It is **not** called capacity-matched.

**Chemistry-shuffle diagnostic.** The same frozen derangement applied at
inference time, reported as a secondary arm.

## 10. Gates — frozen

All primary effects use construct-balanced closure-component macro inference
with one-sided 95% component-bootstrap lower bounds. No gate may be lowered,
replaced or reinterpreted after any result is seen.

| id | contrast | margin | requirement |
|---|---|---:|---|
| `R1` | candidate `AP_bidir` − per-pair chance | ≥ 0.05 | LCB95 > 0 |
| `R2` | candidate − frozen B5 differential | ≥ 0.03 | LCB95 > 0 |
| `R3` | candidate − two-ligand foreign-pair control | ≥ 0.03 | LCB95 > 0 |
| `R4` | candidate − residue-context corruption | ≥ 0.03 | LCB95 > 0 |
| `R5` | candidate − trained permuted-label learner | ≥ 0.05 | LCB95 > 0 |
| `R6` | `G_2B` pair AP − sealed B5 pair AP | ≥ −0.005 | LCB95 ≥ −0.005 |

**Frozen B5 differential baseline (`R2`).** `Delta s_{B5,r} = m_r(P,L_a) −
m_r(P,L_b)` where `m_r(P,L)` is the mean over the atom axis of the **sealed** B5
pair logits for that record. The per-record atom-mean vectors are precomputed
from the sealed float16 memmap and hashed. (Analytically the ligand-dependent
part of this is `(P r_r) · (Q gbar(L))`, a rank-32 bilinear differential; the
candidate is rank 8, so `R2` is not a capacity comparison in the candidate's
favour.)

**`R6` auxiliary pair score.** `G_2B_ra(P,L) = G_B5_ra(P,L) + delta_r(P,L)`,
`delta_r` broadcast along the atom axis, scored with the identical tie-aware
estimator and the identical evaluation mask as Phase 2A, paired by component
against the sealed B5 component table. `R6` is a **non-inferiority** gate;
`LCB95 > 0` is **not** required.

The scaffold-strict held-out-B analysis is **secondary**. It must not show a
sign reversal, and if it is underpowered that is reported rather than
suppressed.

## 11. Replicate reproducibility reference

Renamed from "replicate oracle ceiling". It is **not** a mathematical ceiling — a
model that denoises annotation error may legitimately exceed the agreement
between two noisy annotations of the same system.

It is computed only on the preregistered matched subset where the replicate
structure exists (same construct, same `graph_key`, different `pdb_id`), the
subset size is stated, and it is not extrapolated to the full held-out panel.
**It cannot determine any PASS.**

## 12. Module-participation audit — fail-closed, before any gate is read

1. nonzero `U` and `V` gradient norms throughout training;
2. relative parameter movement `||theta_final − theta_init||_F / ||theta_init||_F`
   ≥ **0.05** for each of `U` and `V`;
3. non-degenerate activation variance of `U h_r` and `V g(L)` on held-out data
   (each ≥ `1e-8`);
4. zeroing `U h_r` collapses `AP_bidir` to within `0.005` of chance;
5. zeroing `V g(L)` collapses `AP_bidir` to within `0.005` of chance;
6. residue-context shuffle degrades `AP_bidir`;
7. ligand shuffle degrades `AP_bidir`;
8. same-seed execution reproduces the checkpoint hash and the predictions
   bit-identically;
9. synthetic recovery passes (section 13).

Only **output-level** `delta` claims are permitted. The bilinear factorisation
carries a rotation gauge — `U -> R U`, `V -> R V` leaves `delta` unchanged for
any orthogonal `R` — so individual `U`/`V` channels are not interpretable and no
claim about them will be made.

## 13. Preflight and synthetic trainability

Before any real-label metric is opened:

1. hash every input, split, cache, checkpoint and control map;
2. verify every required ESM2 `seq_key` is present;
3. verify train/held-out closure-component overlap is exactly zero;
4. verify held-out ligand `graph_key` overlap with training is exactly zero;
5. report held-out-B scaffold overlap;
6. prove `g(L)` is invariant to atom permutation (`<= 1e-12`);
7. prove swapping ligand order flips the sign of `Delta s` exactly;
8. prove the hierarchical aggregation is invariant to pair duplication;
9. exercise the constant-`b` and degenerate-`b` projection fallback;
10. verify float64 projection orthogonality against the `1e-8` tolerance;
11. verify gradients reach `U` and `V` through the projection;
12. **recover a synthetic rank-8 projected differential teacher**;
13. verify same-seed checkpoint and prediction determinism;
14. audit file access and prove affinity-value reads remain exactly zero.

**Synthetic teacher.** Draw `U*`, `V*` with seed `20260905`, form the projected
`delta*`, and on each training pair define the synthetic gain set as the top
`m = 8` residues by `Delta delta*` and the synthetic loss set as the bottom `m`.
Train the identical pipeline on those synthetic labels and evaluate `AP_bidir`
on held-out synthetic pairs. **Required: `AP_bidir >= 0.50`.** The teacher lies
exactly in the hypothesis class, so failure indicts the optimizer, projection or
aggregation, never the biology.

No hyperparameter may be tuned against the synthetic holdout after this contract
is frozen.

Failure of any item terminates the stage with
`PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED`.

## 14. Terminal verdicts — earliest-failure rule, exactly one

```text
1. PHASE2B_CONTRACT_OR_ARTIFACT_FAIL_CLOSED
     a required input, split, baseline, map, hash or projection contract failed

2. PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED
     the exact model/projection/optimizer pipeline failed its synthetic or
     numerical control; no biological conclusion is permitted

3. PHASE2B_MINIMAL_RESIDUAL_NOT_IDENTIFIED
     R1 or R2 failed. Rejects ONLY frozen ESM2 + pooled 41-D atom features +
     rank-8 projected bilinear residue residual. It does not close every
     sequence + 2-D model class.

4. PHASE2B_SHORTCUT_DEPENDENCE
     R1 and R2 passed but R3, R4 or R5 failed. Reject the candidate.

5. PHASE2B_RESIDUE_DIFFERENTIAL_IDENTIFIED_BUT_B5_INTEGRATION_FAILED
     R1-R5 passed, R6 failed.

6. STRUCTURAL_LIGAND_CONDITIONED_RESIDUE_STATISTIC_IDENTIFIED_IN_DEVELOPMENT
     all gates and the module-participation audit passed
```

Verdict 6 authorizes **only** a separately preregistered sealed structural
confirmation. It does not authorize affinity, few-shot adaptation or `z`
admission.

## 15. Interpretation limits

Held-out A has already been read during Phase 1, Phase 2A and this design. It is
**development evidence, not independent confirmation**.

A Phase 2B PASS may establish only *ligand-conditioned residue localization in
development*. It may not establish exact residue–atom coupling, physical
interaction energy, affinity direction, selectivity, off-target prediction,
few-shot adaptation, a biological coordinate admitted into `z`, or a validated
end-to-end DTA model.

## 16. Theory boundary

`A(F,z) = K(B(z)F(z))` is unchanged. Phase 2B's BCE/AUPRC objective is an
upstream structural auxiliary experiment; the frozen theory provides **no**
guarantee for pairwise, listwise, ranking or AP objectives, and none is claimed.

No variable-length residue vector may enter `z`. Any future `z` bridge must
separately define a frozen finite-dimensional bounded statistic, a measurable
query-label-free `z(S,Q,gamma)`, a fixed compact metric domain `Z`, support rank,
conditioning and query coverage, abstention or zero adaptation outside
support-observable directions, and the continuous affinity-band loss with
positive ridge that the frozen theory requires.

## 17. Registered expectations

1. `AP_bidir` is expected to be **modest in absolute terms**: Phase 2A measured
   that only ~44% of the alternative-ligand mask difference is
   ligand-attributable, the rest sitting at replicate noise level.
2. `R2` and `R5` are expected to be the discriminating gates. `R2` is
   demanding because the frozen B5 differential is already a rank-32 bilinear
   in the same frozen states.
3. A **large** `AP_bidir` should first be suspected of construct leakage — the
   same PDB entry or a near-duplicate ligand on both sides of a pair — and
   `PHASE2A_CONSTRUCT_GROUPS.json` is the first artifact to inspect.
