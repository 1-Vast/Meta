# Preregistration — P1R2B-PHASE2B-S4R-A

## Frozen ligand representation availability and information audit

Stage identifier: `P1R2B-PHASE2B-S4R_A_LIGAND_REPRESENTATION_AUDIT`

Written 2026-08-10, after `PHASE2B_S3R_GATE.json`
(`REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED`) and before any S4R training
run, checkpoint, prediction or evaluation metric exists.

This document registers **only** the label-blind representation audit. It does
not authorize training. Training requires the separate S4R preregistration,
which may only be written after this audit has produced a persisted artifact.

## 1. Question

S3R failed on the basis

```text
frozen ESM2 residue states  x  mean-pooled 41-D ligand atom features
```

The leading hypothesis is that the global mean over the 41-D atom features
destroys ligand graph topology, functional-group arrangement, atom-local
identity in context, and scaffold-specific substructure, so no linear map
`W` can express the residue direction the labels require.

Before spending a training budget, decide by measurement alone:

> Does a frozen graph-aware 2D ligand statistic preserve information that the
> mean-pooled 41-D representation provably cannot express?

If no member of the declared candidate family does, the terminal verdict is
`GRAPH_LIGAND_REPRESENTATION_NOT_INFORMATIVE` and no model is trained.

## 2. What the estimator can actually see

The S3R/S2R estimand is

```text
s(P, La, Lb) = (I - Q_P Q_P^T) H_P W [ g(La) - g(Lb) ],   W in R^(1280 x D_g)
```

`H_P` and `Q_P` are protein-side and are held fixed by this repair. For a fixed
protein the residue direction is `H_P (W dg)` with `dg = g(La) - g(Lb)`. The
reachable set of residue directions is therefore exactly the image of the
**ligand difference vectors** under one linear map.

Every audit statistic below is consequently computed on `dg`, not on `g`. A
representation that is rich in `g` but degenerate in `dg` cannot help this
estimator.

## 3. Inputs and firewall

Label-blind. This audit opens:

- the frozen S7/L2B corpus assembly `s7_dataset.build()` (structures, sequences,
  scaffolds, atom names, closure components);
- the frozen MONN `mol_dict` RDKit molecules;
- the frozen 41-D `atom_features` and `g_of` mean.

This audit does **not** open, and must record zero reads of:

- any MONN residue edge mask or residue label view;
- heldout-B;
- ChEMBL, BindingDB, DAVIS, KIBA, recipient or metaval affinity values;
- any S3R score, checkpoint or Gate value.

Because `build_pairs` filters on the residue symmetric difference, the audit may
**not** use the S3R eligible pair list. It uses the strictly label-blind
superset

```text
same seq_key, different graph_key, both Murcko scaffolds present and distinct
```

deduplicated to unordered ligand-graph pairs, computed separately inside the
train split and inside the heldout-A split. All A-gates below are evaluated on
the heldout-A superset, since generalization under closure shift is the failure
being diagnosed.

## 4. Candidate family — frozen

The repository already contains a frozen graph-aware 2D ligand encoder:
`rdkit.Chem.rdFingerprintGenerator.GetMorganGenerator`, used by
`pa2_teacher_conditionality.py` for Phase 2A scaffold distinctness. No new
trainable graph network is introduced, and none is permitted by this audit.

Morgan/ECFP environment identifiers are the canonical Weisfeiler-Lehman
relabelling of a molecular graph: the identifier of atom `i` at radius `r` is a
hash of its own atom invariants together with the multiset of its bonded
neighbours' identifiers at radius `r-1`. They therefore encode bond
connectivity and functional-group arrangement, which an atom-marginal mean
cannot.

The declared candidate family is

```text
radius r in {1, 2}
vocabulary size d in {128, 256, 512}
```

with, for ligand `L` having `n(L)` heavy atoms,

```text
g_graph(L)[j] = count of vocabulary environment j in L   /   n(L)
```

The vocabulary is the top-`d` Morgan environment identifiers ranked by the
number of **distinct train-split ligand graphs** containing them, ties broken
by ascending identifier. Vocabulary construction reads train ligand structures
only; it reads no label and no heldout-A structure.

`GetMorganGenerator` is instantiated with `fpSize = 2**20` and the count
fingerprint's non-zero elements are used directly, so vocabulary entries are
exact environment identifiers and no folding collision occurs.

### 4.1 Why the pooling is unchanged

The baseline is

```text
g_base(L) = mean over heavy atoms of a 41-D atom-local one-hot descriptor.
```

The candidate is

```text
g_graph(L) = mean over heavy atoms of a d-D atom-centred neighbourhood one-hot
             descriptor, truncated to a frozen train-derived vocabulary.
```

Both are an arithmetic mean over heavy atoms of a per-atom descriptor, both are
invariant to atom permutation, and both are size-normalized the same way. The
**only** axis that changes is whether the per-atom descriptor is atom-local or
graph-aware. Nothing else in the experiment may change.

## 5. Statistics — frozen

Let `B` be the heldout-A pair-difference matrix of the 41-D baseline and `M`
the heldout-A pair-difference matrix of a candidate, both over the same
label-blind pair superset and the same pair order.

1. **Inventory.** Distinct ligand graphs, scaffolds, constructs and closure
   components, per split; distinct Morgan environments per radius, and how many
   heldout-A environments are absent from the train vocabulary.
2. **Exact collapse.** Number of groups of distinct ligand graphs whose
   representation vectors are bit-identical, and the number of graphs covered.
3. **Near collapse across distinct scaffolds.** Fraction of scaffold-distinct
   within-construct heldout-A pairs whose representation cosine similarity
   exceeds `0.999`, and the fraction with `||dg|| = 0`.
4. **Effective rank.** `exp` of the Shannon entropy of the normalized squared
   singular spectrum of the mean-centred matrix, for `g` over all distinct
   ligand graphs and for the heldout-A difference matrix. Numerical rank and
   the ratio `sigma_1 / sigma_k` conditioning at `k = min(20, D_g)` are also
   reported.
5. **Incremental information beyond the 41-D mean.**
   `INC = ||M - [1 B] beta||_F^2 / ||M||_F^2` with `beta` the least-squares
   solution. This is the fraction of candidate difference energy that no linear
   function of the baseline difference can express.
6. **Retention of the 41-D mean.** `RET = ||B - [1 M] gamma||_F^2 / ||B||_F^2`.
   Small `RET` means the candidate is an enrichment of the baseline rather than
   a swap of one lossy view for another.
7. **Size and composition marginals.** `R^2` of `||dg||` on `|dn_atoms|`, for
   baseline and candidate; and the size-and-baseline-orthogonal incremental
   information `INC_perp`, regressing `M` on
   `[1, B, dn_atoms, d log n_atoms]`.

## 6. A-gates — frozen thresholds

Every threshold is expressed against the mean-pooled baseline measured on the
same pairs, because the baseline is the object being repaired. A candidate
`(r, d)` is **admissible** when all of the following hold on heldout-A.

```text
A1  effective rank of the candidate difference matrix
      >= 3 x effective rank of the baseline difference matrix
A2  INC  >= 0.25          at least a quarter of the candidate difference
                          energy is unreachable from the baseline difference
A3  RET  <= 0.10          the candidate retains at least 90% of the baseline
                          difference energy in its own linear span
A4  coverage >= 0.99      at least 99% of heldout-A label-blind pairs have a
                          non-zero candidate difference
```

`INC_perp`, conditioning, exact collapse and the size `R^2` are **reported and
non-gating**; they enter the S4R preregistration as declared confounds, not as
admission criteria.

## 7. Selection rule — frozen

If no `(r, d)` in the declared grid is admissible, terminate with
`GRAPH_LIGAND_REPRESENTATION_NOT_INFORMATIVE`.

Otherwise select the admissible pair minimizing `(d, r)` lexicographically.

This is a **capacity-parsimony** rule, not a performance rule. The standing
project constraint forbids adding capacity, and the trainable object scales as
`1280 x D_g`. Among representations that provably clear every information
gate, the smallest is therefore the correct single-axis test. No evaluation
metric, no label, and no training result may influence this choice, and the
rule may not be re-run after any S4R metric is read.

## 8. Disclosure of prior unregistered exploration

Before this document was written, an unregistered scratchpad script computed
the inventory of section 5.1 and the section 5.4-5.7 statistics over the
section 4 grid. That exploration was strictly label-blind. It is disclosed
here because it means the A2, A3 and A4 numeric values in section 6 were chosen
by an author who had already seen candidate values, and A1 was expressed as a
baseline multiple for the same reason.

The mitigations are: every threshold is anchored to the measured baseline
rather than to a candidate value; the selection rule of section 7 is
capacity-parsimony and cannot be steered toward a better-scoring candidate; and
the registered audit recomputes every statistic from source so the disclosed
exploration contributes no number to any artifact. Nothing in that exploration
touched a residue label, an affinity value or an S3R metric.

## 9. Deliverables

```text
report/s7_l2b_r0r/PHASE2B_S4R_REPRESENTATION_AUDIT.json
report/s7_l2b_r0r/PHASE2B_S4R_REPRESENTATION_AUDIT.md
```

The JSON records this document's SHA-256, the executing commit, the full grid
of statistics, every A-gate outcome, the selected `(r, d)`, the frozen
vocabulary SHA-256, the parameter count implied for `W`, and
`residue_label_reads = 0`, `affinity_value_reads = 0`.

## 10. Terminal verdicts of this audit

Exactly one:

```text
S4R_AUDIT_CONTRACT_FAIL_CLOSED
GRAPH_LIGAND_REPRESENTATION_NOT_INFORMATIVE
GRAPH_LIGAND_REPRESENTATION_AVAILABLE_AND_INFORMATIVE
```

Only the last authorizes writing the S4R training preregistration, and it does
so for the selected `(r, d)` alone.

## 11. Boundary

This audit measures a ligand representation. It does not measure biology, does
not admit any statistic to `z`, and does not modify

```text
A(F, z) = K(B(z) F(z)).
```

A representation that clears every A-gate is still only an upstream biological
measurement whose usefulness is unknown until the registered S4R Gates are run.
