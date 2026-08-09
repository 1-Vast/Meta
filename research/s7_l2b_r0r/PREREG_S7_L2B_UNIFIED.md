# PREREG S7_L2B_UNIFIED — exact-residue localization

Registered: 2026-08-09, **before any model exists and before any AP is computed**.
Stage: `S7_L2B_R0R_RECONSTRUCT_DATA_BASELINE_AND_CONTRACT` → modelling phase.
Everything below is fixed by this document. Amendments require a new commit with
a stated reason; in-place edits are forbidden.

All inputs to this document are **label-blind**: edge counts, identities and
topology only. No model has been fitted, no AP computed, no affinity value read.

## 1. Estimand and its scope

> Can an exact-residue localizer, from protein sequence and a 2D ligand graph,
> predict which (residue, ligand-heavy-atom) pairs are in non-covalent contact,
> **beyond a ligand-only model and beyond mismatched-partner controls**, on
> proteins it has never seen?

**Scope, stated up front:** generalisation is over **proteins**. The
ligand-disjoint strata additionally test generalisation to unseen ligands, but
the inference units are protein components, so the primary claim is
protein-level.

### 1.1 Closure — methodological correction, registered

R0R-2 measured that union-**merging** ligand identity into the inference
partition chains transitively through promiscuous cofactors (ATP/ADP/NAD/heme
class) and yields a **93.19 %** giant component over 14,589 complexes, making
component inference impossible. Ablation localised the cause: protein closure
alone gives 1,994 components with a 3.11 % largest; adding exact ligand graph
alone drives it to 76.22 %.

**Adopted construction.** Inference units are **protein closure components**
(exact PDB, exact sequence, UniProt, 40 % alignment-verified homology). Ligand
closure is enforced as a **disjointness filter** between train and held-out
sets, not as a merge relation. This achieves protein *and* ligand leakage
control without transitive chaining.

This is a change of construction, not a relaxation: the held-out sets below are
*more* strictly separated from training than a merged partition would have
allowed, because a merged partition could not have been evaluated at all.

## 2. Data, frozen

Source: MONN at `f2b62ccf49c18a9502aa0eb0d582c6e0735ef200`, all six file
SHA-256 verified, non-commercial use, never redistributed. Neither affinity TSV
is ever opened. Edge corpus reproduced deterministically and bit-identically
(`MONN_LOCAL_REPRODUCTION_AUDIT.json`).

Split seed **20260810**, target held-out fraction 0.20, whole components only.

| Set | Complexes | Components | Largest frac | Ligand graphs | Positive edges |
|---|---:|---:|---:|---:|---:|
| train | 9,759 | 1,376 | 0.0377 | 7,540 | 151,066 |
| **held-out A (primary)** — protein- **and** exact-ligand-graph-disjoint | **2,415** | **196** | **0.1520** | 2,157 | 36,046 |
| held-out B (secondary, strict) — additionally scaffold-disjoint | 1,881 | 160 | 0.1648 | 1,691 | 29,073 |
| held-out C (diagnostic) — protein-disjoint only | 2,979 | 293 | 0.1480 | 2,462 | 44,732 |

**Confirmation cohort — SEALED, not to be scored in development.**

| Set | Complexes | Components | Largest frac | Positive edges |
|---|---:|---:|---:|---:|
| protein-disjoint | 710 | 325 | 0.0197 | 3,383 |
| **+ exact-ligand-graph-disjoint (primary)** | **486** | **225** | **0.0288** | **2,323** |
| + scaffold-disjoint (strict) | 347 | 170 | 0.0403 | 1,603 |

**Confirmation admissibility criteria, frozen here:** ≥200 inference components
and largest component ≤25 %. Both hold for the primary confirmation stratum
(225 components, 2.88 %). These are set from bootstrap-power reasoning — a
component bootstrap needs on the order of 10² independent units and no dominant
unit — not from what the data happens to yield. Publication/time closure
(R0R-3) is an **additional** requirement before confirmation may be opened; it
is not yet built, so confirmation remains closed regardless of development
outcome.

## 3. Primary target and metric

Binary residue × ligand-heavy-atom contact.

Evaluation scores the **complete valid matrix** for every held-out complex.
Training may subsample negatives; evaluation may not.

```text
1. average precision per complex over its complete matrix
2. mean complex AP within each protein closure component
3. equal-weight macro average across components
4. paired whole-component bootstrap, 2,000 replicates, seed 20260811
```

Rows/pairs/complexes are never treated as IID. AP implementation, dtype
(float64 for scoring) and deterministic tie handling are pinned; ties are broken
by a fixed lexicographic (residue_index, atom_slot) order so AP cannot depend on
row order. A complex with zero positives is excluded from AP and recorded.

## 4. Arms — identical rows, masks, targets, budget

| Arm | Residue input | Ligand input | Role |
|---|---|---|---|
| `B0` | — | — | per-complex prevalence constant |
| `BL` | **none** | atom features | ligand-only shortcut |
| `B4` | explicit non-PLM local sequence features | atom features | matched low-capacity baseline |
| `B5` | frozen ESM2 per-residue embeddings | atom features | the candidate |
| `BP` | residue features from a **different** protein | atom features | wrong-protein control |
| `BX` | residue features | atom features from a **different** ligand (different exact graph and different scaffold) | wrong-ligand control |
| `BM` | residue features with residue **order shuffled** within the complex | atom features | motif/positional control |

`B4` and `B5` differ **only** in the residue representation. Same atom
representation, same head, same rank, same negative sampler, same epochs,
optimiser, learning rate, weight decay, seed, and same evaluation mask.

Every contrast uses one common analysis mask: the intersection of validity over
all arms. An arm-mask mismatch is a contract failure.

Forbidden as inputs anywhere: target ID, PDB ID, UniProt ID, assay ID, CCD code,
dataset membership, and any global ligand or protein identifier.

## 5. Head — fixed

```text
r_i = GELU(W_h LN(x_i) + b_h)          r_i in R^128
pi_i  = sigmoid(w_pi  . r_i + b_pi)     residue prior
rho_j = sigmoid(w_rho . a_j + b_rho)    atom propensity
s_ij  = b + alpha logit(pi_i) + beta logit(rho_j) + (P r_i) . (Q a_j)
```

Interaction rank **32**. Six epochs. AdamW, lr 1e-3, weight decay 1e-4, batch by
complex, seed 20260812. Deterministic negative sampler: exactly **6 unique true
negatives per positive**, drawn from residue-local, atom-local and uniform
categories with uniform backfill; never a positive, never a duplicate; a
negative-sampling manifest is written.

**Parsimony, registered:** only the **binary** channel is modelled. Typed
channels are NOT added at this stage. All seven PLIP channels were measured as
prevalent enough to be evaluable, but constraint 9 forbids adding them until an
experiment shows the binary channel is information-limited in a way typing would
fix. Adding them now would be capacity without evidence.

## 6. Gates — frozen thresholds, one-sided LCB95 > 0 throughout

Primary, on held-out A, macro-AP over components:

```text
G1  B4 - B0  >= 0.02      the baseline must beat prevalence at all
G2  B4 - BL  >= 0.02      protein information must add over ligand-only
G3  B4 - BP  >= 0.02      the CORRECT protein must be required
G4  B4 - BM  >= 0.02      residue order/position must be load-bearing
G5  B4 - BX  >= 0.02      the CORRECT ligand must be required
```

`B5` is authorised **only if** `B4` establishes the pipeline is sound, and `B5`
must then satisfy the same five Gates **plus**:

```text
G6  B5 - B4  >= 0.02      the PLM must earn its parameters
```

Held-out B (scaffold-strict) is a required robustness report; a sign reversal
there is a stop condition.

**No threshold may be changed after any AP is observed.**

## 7. Registered escalation rule for B5

`B5` (ESM2-650M, ~2.5 GB download, frozen weights) is added **only** if the
evidence indicates missing protein information, defined in advance as:

> `B4 − BL` fails G2, **or** `B4 − BL` passes but `B4` macro-AP remains below
> 0.10, indicating the residue representation is the binding constraint.

If `B4` already satisfies every Gate with a large margin, escalating to a 650M
PLM is capacity without evidence and is **not** authorised by this document.

## 8. Terminal verdicts

```text
S7L2B_EVALUATOR_CONTRACT_FAIL_CLOSED
S7L2B_TRAINABILITY_FAIL          (synthetic control not recovered)
S7L2B_LIGAND_ONLY_SHORTCUT       (G2 fails: no protein information)
S7L2B_PARTNER_CONTROL_FAIL       (G3/G4/G5 fail)
S7L2B_BASELINE_ESTABLISHED_PLM_NOT_INDICATED
S7L2B_BASELINE_ESTABLISHED_PLM_INDICATED
S7L2B_PLM_BELOW_GATE
S7L2B_EXACT_RESIDUE_LOCALIZATION_DEVELOPMENT_PASS
```

Only the last authorises opening the sealed confirmation cohort, and only after
R0R-3 publication/time closure is built. Nothing here authorises affinity,
few-shot adaptation, or any change to `K(B(z)F(z))`, CSMO, Band, simplex,
positive ridge, mesh, `theory/`, `model/` or production `z`.

## 9. Trainability control

Before believing any negative result, a synthetic control fits a **known**
function of the frozen inputs under the identical pipeline. Failure to recover it
is `S7L2B_TRAINABILITY_FAIL` and is an optimisation defect, **not** evidence
about biology.

## 10. Module participation

For every trainable branch: gradient norm, parameter update norm, activation
variance, saturation fraction, branch ablation, gradient blocking, input
shuffle, held-out contribution, and incremental value over the simpler arm.
A branch that receives gradients but adds no held-out value over `BL` or `B4` in
the contrast it exists to explain is **removed**, not reported as a success.
