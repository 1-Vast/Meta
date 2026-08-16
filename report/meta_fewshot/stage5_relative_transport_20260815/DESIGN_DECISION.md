# Architecture audit and design decision

Date: 2026-08-15. Written before any Stage 2 arm produced data.

## 1. What the local data actually supports

`GEOMETRY_COVERAGE_AUDIT.json` (reproduce with
`scripts/audit_geometry_coverage.py`):

| quantity | value |
|---|---:|
| governed holo complexes (`pilot20k_holo_governed_v2`) | 14,906 |
| unique receptor sequences in that corpus | 14,906 |
| unique ligand chemotypes in that corpus | 2,869 |
| raw mmCIF coordinates present on disk | yes |
| BindingDB Ki deployment cells | 17,717 |
| DTA targets with an **exact** holo sequence | 15 / 499 |
| DTA targets with a containment-level holo sequence | 110 / 499 |
| DTA ligands sharing a holo SMILES | 84 / 9,880 |
| **DTA cells with a common-frame protein-ligand complex** | **0** |
| DTA targets with a common-frame complex | 0 |

`pilot20k_structure_supervision_v2` materializes `contact`, `distance_bin`,
`atom_mask`, `residue_mask` at 128x128 — **rotation- and translation-invariant
summaries, not coordinates**. `r0b_exact_geometry_v3` likewise stores
`distance_angstrom` and `distance_bin`, not positions.

The decisive number is zero. Not 3%, not 1%: **no BindingDB deployment pair has
a solved complex.** The 15 exact-sequence targets have holo structures bound to
*different* ligands (HEM and similar CCD chemotypes), which is not the pair the
DTA task asks about.

`model/cartesian.py` is correct — `tests/test_cartesian.py` verifies O(3)
equivariance including reflection, translation invariance, symmetric-traceless
rank-2 structure, node permutation, padding safety, finite gradients without
geometry, and rejection of cross-sample edges. It has no legal input here.

## 2. Candidate comparison

| design | requires | applicable to this task? | verdict |
|---|---|---|---|
| Cartesian rank-2 joint protein-ligand encoder (PBCNet2.0 / TensorNet style) | common-frame complex coordinates | **no** — 0/17,717 cells | rejected on data, not on merit |
| PaiNN / MACE / Equiformer / E2Former / SE(3)-EGNN on the complex | same | **no** — identical blocker | rejected |
| Separate equivariant protein and ligand encoders, invariant fusion | independent frames | ligand conformers are not in the bank; protein coordinates exist for 15/499 targets; fusing independently framed structures is exactly the manufactured-geometry failure mode the contract forbids | rejected |
| Structural auxiliary supervision (contact/distance co-training on the holo corpus) | invariant labels only — available | legal, but the holo corpus shares ~0 ligands and 110/499 targets at containment level, and it adds a second data stream plus an extra loss on top of an unfixed mechanism | deferred, documented |
| Geometry-optional hybrid with no-coordinate fallback | — | the geometry branch would never activate on a single DTA episode; it buys nothing measurable | rejected as ceremonial |
| **Signed reference-query relative transport (PBCNet2.0's Siamese idea with geometry removed)** | ligand graph + protein sequence only — available for every cell | **yes** | **selected** |

No architecture can conjure a complex that does not exist. Since every
equivariant family fails on the same input constraint, the choice among them is
moot for this task, and selecting one would force either manufactured geometry
or a branch that never fires.

## 3. What is borrowed, and from where

* **PBCNet2.0** ([bioRxiv 2025.06.04.657800](https://www.biorxiv.org/content/10.1101/2025.06.04.657800v1),
  [Nat Chem Biol](https://www.nature.com/articles/s41589-026-02241-x)) — the
  **Siamese reference-query relative formulation**: predict the *signed
  difference* between a reference and a query rather than an absolute
  correction, with sign-flip consistency. Borrowed. Its Cartesian rank-2
  message passing over 8.6M complex pairs is **not** borrowed: it needs
  complexes.
* **AdaMBind** ([Nat Commun](https://www.nature.com/articles/s41467-026-70554-5),
  [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13102954/)) — **target as task**
  and **task-difficulty weighting**. Borrowed as a leave-one-out label-consistency
  credit on each support observation. Its MAML inner loop, test-time gradient
  updates and easy-to-hard scheduling are **not** borrowed: the contract
  excludes inner loops and multi-stage recipes.
* **MetaDTA / attentive neural processes**
  ([openreview](https://openreview.net/pdf?id=yzlif16IASM)) — support-attention
  used as a weight distribution over support affinities. Already present in the
  Stage 4 model; the admission run showed it is **not sufficient alone**, which
  is what motivates the difference operator.
* **TensorNet** ([arXiv 2306.06482](https://arxiv.org/pdf/2306.06482)) — read for
  the rank-2 Cartesian decomposition and its efficiency argument; not used, for
  the reason above.

## 4. Selected mechanism

```text
f(P,Lq,S) = f0(P,Lq) + s(n) * sum_k w_qk * [ r_k + delta(P, L_k -> Lq) ]
r_k         = y_k - f0(P,L_k)                       (label locked, detached)
s(n)        = n / (n + lambda),   s(0) = 0
delta(a->b) = m(e_a,e_b) - m(e_b,e_a)               (exactly antisymmetric)
w_qk        = softmax_k( tau <key_q, key_k> + c_k )
c_k         = -|r_k - LOO_k| / kappa                (label-consistency credit)
```

`e` is the interaction embedding already produced by the grammar trunk, so the
trunk, encoders and zero-shot endpoint are **unchanged** and `--arch grammar`
remains an exact comparator.

### Why this addresses the measured Stage 4 failures

| Stage 4 failure | mechanism response |
|---|---|
| `full` beat `sar_cut` with no positive bootstrap lower bound; the gate could only rescale a residual toward the level | `delta` adds a signed, query-dependent shift with no bound tying it to the level, so the correction is not a shrinkage |
| CI fell 0.647 -> 0.571-0.610, Spearman 0.372 -> 0.169-0.257 | `delta` is trained on *differences*, the quantity that determines ranking; and at `delta == 0` the model returns exactly the previous safe level calibration, so the floor is neutral rather than harmful |
| permuted support was **better** at k=2,3 — support identity unused | `c_k` measures whether each support label agrees with the structural difference operator applied to the other supports. Under permutation that agreement collapses, so the permuted arm is worse **by construction**, not by luck. This is inactive at k=1 (no leave-one-out) and active exactly where the control failed |
| k=1 was scalar | `f = f0(q) + s*(r_1 + delta(1->q))` is query-specific at k=1 |

### Algebraic properties (all tested)

`delta(a->a) = 0`; `delta(a->b) = -delta(b->a)`; transport affine in the support
labels; `n=0` returns `f0` exactly; `delta = 0` with flat weights returns the
shrunken support mean exactly; support permutation invariant; query permutation
equivariant; label permutation **not** invariant; padding invariant; every
trainable tensor receives gradient. The readout bias was removed because it
cancels identically in an antisymmetric construction and would be an
unidentifiable parameter.

## 5. Honest limits declared up front

* No Cartesian equivariance, atomic 3D recognition, or binding-mode claim is
  made or testable on this task; the model refuses coordinate inputs with an
  explicit error naming the coverage audit.
* Using the holo corpus would require either a second training stream
  (deferred) or docking/co-folding to invent complexes (forbidden here).
* Selection and Stage 2/3 decisions use **meta_val**. The historical meta_test
  split has been consumed by prior architecture search and is reported only at
  the final stage.
