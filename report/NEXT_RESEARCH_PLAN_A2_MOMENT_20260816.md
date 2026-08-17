# Next Research Plan: Protein-Conditioned SAR Moment Update (A2 Moment Meta)

Status: **preregistered design; not yet implemented or trained**. This document
supersedes scattered follow-up ideas. An executor must first read `task.md`,
`report/CURRENT_MODEL_EVIDENCE.md`, and this document in full.

## 1. Why this is the only active model family

R7-R14 closed three axes: replacing the ranking loss again, unconstrained
query-specific attention, and larger shape MLPs. The current evidence is that
A0 still has weak zero-shot ordering, while fixed Morgan/Tanimoto transport can
use support labels reliably for k>=2. Seven learned query-specific channels were
deployment-inert or damaged calibration. The next question is therefore not
whether more capacity helps, but whether frozen A0 representations contain a
**low-dimensional, protein-conditioned SAR coordinate system that at most five
labels can identify**.

The research contribution is deliberately concentrated in two claims:

1. **Model innovation: Protein-Conditioned SAR Moment Update.** Support labels
   form a low-rank moment statistic without a closed-form solver, inner loop, or
   deployment gradient. The update is linear in labels, invariant to support
   order, and exactly zero at k=0.
2. **Training innovation: Correlation-Preserving Counterfactual Meta-Training.**
   This is enabled only after the model mechanism passes a real-data admission
   gate. It must preserve absolute MSE and within-target ordering while learning
   falsifiable label binding from a counterfactual registry independent of the
   evaluation controls.

The biological claim must remain narrow. The coordinate describes empirical
within-target SAR directions induced by ligand changes. It is not an atomic
contact map, binding pose, or free-energy decomposition. None of the current
17,717 DTA cells has a legal common-frame protein-ligand pose. Cartesian and
PBCNet2.0 ideas therefore provide relative-modeling inspiration only and are not
part of this stage.

## 2. The only permitted A2-min operator

For the frozen A0 complex representation `e0(P,L)`, first center it with 32
fixed, label-blind meta-train ligand anchors. A trainable projection `A_phi`
then produces `z(P,L) in R^R`:

```text
r_i = stopgrad(y_i - f0(P,L_i))
c_S = (1/k) sum_i r_i z(P,L_i)
eta(k) = eta_inf * k / (k + lambda)
delta(P,L_q,S) = eta(k) <c_S, z(P,L_q)>
f = f0 + delta
```

`A_phi` and two positive shrinkage scalars are the only new trainable quantities
in A2-min. Do not add attention, a hypernetwork, a solver, ridge regression, a
pseudoinverse, query labels, externally retrieved labels, or additional data.

## 3. Staged execution and hard gates

### S0: governance, power, and algebra

- Freeze the datasets, A0 checkpoint, episode bank, Morgan/Tanimoto comparator,
  seeds, and evaluation code.
- Select hyperparameters only on meta-train component folds. Meta-val permits
  one development confirmation after the choice is frozen. Keep double-cold
  meta-test excluded (logical exclusion after parsing; see the governance incident).
- Compute the number of components and targets, minimum detectable effect, and
  bootstrap stability for every k before training.
- Structural tests must prove k=0 identity, support permutation invariance,
  query permutation equivariance, odd correction under label-sign reversal,
  label-to-query Jacobian rank no greater than `min(Q,k,R)`, padding invariance,
  and the absence of every query-label path.
- Synthetic gates must cover recovery of anisotropic shared low-rank SAR,
  abstention on private mechanisms, monotonically worsening performance as
  label corruption increases, and the absence of stable gains under random or
  wrong proteins.

Any structural-gate failure stops the family before BindingDB training.

### S1: frozen-representation discriminator (A2-min)

Train A2-min only, with an MSE-primary objective and one ordinary
forward/backward pass. The first search is restricted to `R in {4,8,16}`, chosen
inside meta-train component folds. Compare matched-budget arms:

- B0: frozen A0;
- B1: A0 plus scalar support-level calibration;
- B2: fixed Morgan/Tanimoto residual transport;
- B3: A2-min;
- B4: ligand-only A2 with the protein condition removed;
- B5: random frozen coordinates;
- B6: wrong-protein A2;
- B7: matched-wrong labels and support-label permutation.

B3 must produce non-scalar query spread at k=1 and beat B1. At k=2/3/5 it must
be no worse than B2 in MSE, CI, and Spearman, with positive component-paired
bootstrap evidence. B4/B5/B6/B7 must materially weaken the benefit. A target
mean shift, a gain confined to repeated/high-similarity ligands, or an effect
unchanged by the wrong protein does not establish protein-conditioned SAR.

### S2: the single permitted representation repair

S2 is authorized only if B3 beats B1 and random coordinates, fails to beat B2,
and the protein counterfactual proves that real protein-conditioned signal is
present. Replace the linear projection with one normalized `96 -> 32 -> 16`
MLP. Keep every other contract, arm, seed, and budget unchanged. If it still
fails to beat B2, close the A2 family. Do not add attention or capacity.

### S3: central training innovation

Only after S1 or S2 passes, compare:

- C0: MSE only;
- C1: MSE plus a low-weight, valid-sample-masked correlation-preservation term;
- C2: C1 plus counterfactual meta-training.

Training counterfactuals must come from a registry independent of evaluation.
Permuted labels must preserve the support mean. At k=1, matched-wrong support
must match residual magnitude. The wrong-protein arm must replace only the
protein condition and must not mix in donor level or donor labels. C2 must
simultaneously lower clean MSE, improve CI/Spearman, and enlarge the gap between
correct and corrupted support. Otherwise the training innovation fails.
Activity-cliff weighting is allowed only as a later, isolated ablation and must
not be introduced with C2.

### S4: optional attention

Only after C2 passes, and only if error analysis identifies support-evidence
allocation as the remaining bottleneck, one low-capacity, label-blind attention
module may estimate effective k. Label values must still enter only linearly
through the moment statistic. Attention is not a core claim. Delete it if it has
no identifiable k=1 effect, reduces to Tanimoto, or perturbs k=0.

## 4. Promotion, stopping, and reporting

Each stage starts with a short bug-finding smoke test and then uses matched
budgets, multiple fixed seeds, and full nested-k evaluation. Short-run numbers
are never performance evidence. A final candidate must satisfy all of the
following:

- no significant k=0 regression; genuine non-scalar adaptation at k=1; aligned
  k=2/3/5 improvements over A0, scalar level, and Tanimoto;
- MSE, CI, Spearman, and activity-cliff sign improve without trading one metric
  against another;
- component bootstrap conclusions are stable and all three seeds point in the
  same direction;
- label permutation, matched-wrong support, wrong protein, ligand-only input,
  and random coordinates fail as predicted;
- parameter count, gradient coverage, peak memory, and wall time are recorded,
  with no dead trainable parameters;
- meta-test is opened once only after the candidate and all development choices
  are frozen and explicit authorization has been recorded.

A failed gate is a valid terminal result. Retain only `PREREGISTRATION.md`,
`RESULT.json`, `REPORT.md`, necessary prediction rows, and a loadable admitted
checkpoint for each stage. Delete duplicate smokes, logs, and failed checkpoints
after consolidating their verdict. Synchronize `history.md`, `task.md`,
`report/CURRENT_MODEL_EVIDENCE.md`, and `report/EVIDENCE_LEDGER.md` after every
decision.

## 5. Relationship to M0/MSA

M0 is an independent protein-side calibration diagnostic and must not be mixed
with A2 in one stage or attribution claim. A passing M0 result authorizes one
subsequent protein-feature ablation; it does not prove meta-learning. A passing
A2 result does not prove that MSA features help. Keep both evidence chains
separate.
