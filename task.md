# Current task

## Non-negotiable objective

Build and validate a **trainable few-shot drug-target affinity predictor for
unseen targets**. Targets are meta-learning tasks. The model must use large
open datasets to learn transferable biological knowledge, then adapt from
`k=1/2/3/5` support affinities without reading query labels.

Success requires all three layers:

1. **Biology:** a partner-specific protein-ligand measurement that adds
   affinity information beyond ligand-only and wrong-protein controls.
2. **Meta-learning:** source target episodes learn a shared low-dimensional
   mechanism basis; support labels estimate only the new target's identifiable
   section.
3. **Mathematics:** the admitted bounded statistic enters the unchanged
   `A(F,z)=K(B(z)F(z))` operator with explicit rank, conditioning, query
   coverage and abstention.

## Minimal model contract

```text
phi(P,L) in R^288                 frozen audited biology candidate
U in R^(288 x d), d <= 5          source/meta-learned task subspace
m(P,L) = U^T phi(P,L)
f0 = f_L(L, endpoint) + w0^T m(P,L)
a_t = positive-ridge support section in row(M_support)
y_hat(P_t,L_q) = f0 + a_t^T m(P_t,L_q)
```

`U` and `w0` are learned through target-wise support/query episodes. This is
the trainable meta-learning core. The closed-form section replaces a free
MAML inner network; it is smaller, deterministic and aligned with the support
identifiability requirement. The Wan et al. AdaMBind paper is a methodological
reference for target-as-task episodic learning and task adaptivity, not a
requirement to reproduce all of its modules.

## Open-data roles

- BindingDB curated Ki/Kd and Klaeger Kdapp: endpoint-specific quantitative
  constraints; never merge endpoint scales.
- Kinobeads, PKIS and PKIS2: within-panel ordinal/profile pretraining, not
  absolute Ki/Kd calibration.
- PDSP: non-kinase development stratum after panel/provenance census.
- Davis, KIBA, recipient and future time-split data: closed confirmation only.

Training data may be large and dependent. Population claims still require
document, assay, protein-homology, ligand-scaffold and publication-time
controls. Optimization authorization and scientific-admission authorization
are separate Gates.

## Phase 0 outcome (2026-08-10)

The episodic stage was attempted and stopped at its first precondition.

```text
FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE
```

Label-blind census, zero affinity label reads, Ki only. Leakage is zero across
target, ligand, scaffold, document and protein-homology-40. Source supply is
ample (442 targets; 220 usable at `k=5`). The evaluation split is the binding
constraint: 68 targets, `24/19/18/16` usable at `k=1/2/3/5` and `24/18/9/8`
scaffold-disjoint, against a declared minimum of 30 and a declared power
ceiling `MDE_d <= 0.600` versus the observed `0.622`.

No model was preregistered, none was trained, and no threshold was moved. The
hypothesis of target-specific coefficient heterogeneity is **untested**, not
refuted.

A label-blind follow-up localized the cause and **rejected the data-absent
explanation**. Recounting the same governed projection under few-shot rules
only — single chain, Ki, `>= k+3` ligands — and dropping the cycle-positive
quotient requirement yields 25,072 rows and 910 targets, of which
584/499/459/**394** support `k=1/2/3/5` and 218 have `>=8` ligands across `>=2`
documents. The quotient-shaped corpus kept 12,457 cells and 236 `k=5` targets.

```text
EVALUATION_PANEL_LIMITED_BY_ESTIMAND_MISMATCHED_CORPUS
```

The cycle-positive filter belongs to the crossed-rectangle estimand, which needs
closed rectangles; few-shot needs only per-target ligand depth. `MDE_d` falls
from `0.622` at 16 targets to `0.249` at 100.

Not yet established: independent **component** depth after protein-40, scaffold
and document closure. Target depth is demonstrated; component depth is not, and
the 85.86% giant-component pathology may persist.

## Immediate next stage (revised)

Preregister a **few-shot-shaped corpus rebuild**, then re-run the unchanged
Phase 0 gate:

1. rebuild from the same governed projection under few-shot admission rules;
2. apply the *same* strict conflict closure — document, protein-40 homology and
   Murcko scaffold union — with no relaxation;
3. report component counts and giant-component share **before** freezing a
   split, failing closed if component depth is still insufficient;
4. re-run Phase 0 with thresholds unchanged (`>=30` evaluation targets at
   `k=5`, `MDE_d <= 0.600`).

No model is authorized until that gate passes.

## Stage after the corpus rebuild passes

Preregister one experiment that changes only the failed sharing assumption:

```text
shared w FAIL
  -> source-learn U, w0 across target episodes
  -> evaluate unseen-target k=1/2/3/5 sections
```

The experiment must compare support-free, zero-section, correct support,
foreign support and permuted-support controls. It must report target-macro
CI/Spearman/RMSE increments, support rank, conditioning, query coverage and
abstentions. It may not lower the failed shared-linear Gate or use confirmation
labels for model selection.

## Complexity boundary

Do not add a new PLM, GNN branch, cross-attention stack, knowledge graph, pose
module, typed energy head or learned support encoder. A new component is
allowed only after a registered ablation proves the current minimal model lacks
the required information. No high-dimensional pair tensor enters `z`.
