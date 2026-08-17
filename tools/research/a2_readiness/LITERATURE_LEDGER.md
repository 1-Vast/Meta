# Literature ledger

Each entry separates **prior art** (what the source establishes), **local
adaptation** (what must change for this protocol), and **new contribution**
(what would remain ours, and therefore must be falsified on its own).
Literature is inspiration and constraint, never evidence that a local
mechanism will work.

## Directly load-bearing

### Ignored counterfactual conditioning

*Prior art.* The conditional-generative literature names the exact failure
mode measured here: a conditional model is free to ignore its conditioning by
finding a solution with `p(x|c) = p(x)`, and **counterfactual conditioning can
be ignored even when observational conditioning is not**. Enforcing it needs
an explicit training signal rather than architectural conditioning
([Unifying Causal Representation Learning with the Invariance Principle,
arXiv:2409.02772](https://arxiv.org/abs/2409.02772); and the counterfactual-
invariance line, [Learning Counterfactually Invariant Predictors](https://files.sri.inf.ethz.ch/wfvml23/papers/paper_12.pdf),
[Counterfactual Invariance to Spurious Correlations](https://openreview.net/pdf?id=BdKxQp0iBi8)).

*Local adaptation.* Our conditioning variable is the protein, and the
observation is that conditioning is honoured in the **level** and ignored in
the **ordering** — a partial, decomposable version of the phenomenon that the
general framing does not distinguish. The `level`/`ordering` split is exact
here because `protein_value` is constant within a target, so we can enforce
conditioning on precisely the component that is being ignored.

*New contribution.* Measuring the split (300:1 level-to-ordering response),
showing that the ignored component is recoverable (110× at random init), and
enforcing conditioning on the centered component only.

### Counterfactual invariance *and sensitivity* in protein binding

*Prior art.* [Counterfactual Peptide Editing for Causal TCR–pMHC Binding
Inference (arXiv:2604.13256)](https://arxiv.org/html/2604.13256v1) pairs an
invariance loss (predictions must not move under conservative substitutions)
with a **contrastive loss that amplifies sensitivity at key residues**. The
principle borrowed is narrow: a counterfactual training signal can be used to
*increase* a model's responsiveness to a biological input, not only to
suppress spurious responsiveness.

*Local adaptation.* We have no residue-level binding labels and no
common-frame structure, so we cannot define "key residues". Our counterfactual
is a whole-protein substitution by a similarity-matched donor from a different
homology component, which the R5 contract already constructs.

*New contribution.* Applying the sensitivity direction to a *centered*
prediction so that the level branch is algebraically excluded from satisfying
it.

### Shortcut learning in drug–target prediction

*Prior art.* [AI-Bind (Nature Communications 2023)](https://www.nature.com/articles/s41467-023-37572-z)
showed state-of-the-art DTI models fail on novel structures because they
exploit the topology of the protein–ligand bipartite network rather than node
features. [Do Protein-Ligand Models Learn Binding Sites or Just Binding
Likelihood? (arXiv:2605.24045)](https://arxiv.org/pdf/2605.24045) asks the
same question of the protein side.

*Local adaptation.* The double-cold split already removes the topological
shortcut this literature targets (zero exact-ligand, scaffold, component and
document overlap). Our measured shortcut is a *different and finer* one that
survives that split: the protein is used for the target level and ignored for
within-target ordering.

*New contribution.* Naming and quantifying a shortcut that persists after the
standard fix, with a level/ordering decomposition and a positive control.

### Regression-compatible ranking, and its limit

*Prior art.* [RCR, arXiv:2211.01494](https://arxiv.org/abs/2211.01494)
formalises the incompatibility between pointwise regression and listwise
ranking optima.

*Local adaptation and outcome.* Implemented and tested in R14. The alignment
identity held exactly, and the term proved **inert** at this model's operating
point (1.7% of the regression gradient). Recorded as a closed axis. Listed
here because it is the reason the loss-form direction is not revisited.

### Neural collapse under regression

*Prior art.* [The Prevalence of Neural Collapse in Neural Multivariate
Regression, arXiv:2409.04180](https://arxiv.org/abs/2409.04180): under MSE
with weight decay, last-layer features collapse onto a subspace whose
dimension matches the target dimensionality; the collapse is driven by the
regularisation.

*Local adaptation.* Our target is one-dimensional, the worst case. This is a
plausible partial account of the readout collapse measured in E4, but it
predicts feature collapse generally, whereas we measure a *selective*
collapse — the level survives, the ordering does not. The selectivity needs
the objective-level explanation (F6), not the geometric one alone.

## Surveyed, informing the comparison, not adopted

| direction | source | why not adopted |
|---|---|---|
| FEAT set-to-set embedding adaptation | [arXiv:1812.03664](https://arxiv.org/abs/1812.03664) | support-set mechanism; inert at the measured k=0 bottleneck and degenerate at k=1 |
| Set Transformer | Lee et al.; [overview](https://www.emergentmind.com/topics/set-transformer) | same; ≤5 elements makes set self-attention near-degenerate |
| Matching networks / metric-based few-shot | classic line | the incumbent's Tanimoto transport is already a fixed-metric matching network |
| CNP / Attentive NP | [arXiv:2210.09211](https://arxiv.org/abs/2210.09211) | legal in shape (no test-time gradient) but a k≥1 mechanism; would be the eighth query-specific channel |
| ActFound pairwise meta-learning | [Nature Mach. Intell. 2024](https://www.nature.com/articles/s42256-024-00876-w) | uses external ChEMBL assays; would be external data, and it reports strong results *without* protein information, which cuts against the protein-conditioned thesis |
| SQRL similarity-quantised difference learning | [arXiv:2501.09103](https://arxiv.org/abs/2501.09103) | subsumed by the retained Tanimoto comparator |
| GraphCliff / cliff-aware contrastive | [arXiv:2511.03170](https://arxiv.org/pdf/2511.03170) | R9 measured the cliff weight as a net negative for ranking |
| Optimal-transport protein aggregation | [Bioinformatics Advances 2025](https://academic.oup.com/bioinformaticsadvances/article/5/1/vbaf060/8088230) | protein-representation interventions were rejected in Stage 8 on 3/3 seeds; and E4 shows the protein is already read — the deficit is in the readout, not the encoder |
| Cartesian / equivariant encoders | PBCNet2.0, TensorNet, MACE, Equiformer | 0/17,717 legal common-frame inputs |

## Honest statement of novelty

The proposed training signal is a **composition of two mechanisms this
repository already contains** — within-target centering and a protein
counterfactual contrast — motivated by a measurement showing that composing
them is exactly what is missing. The machine-learning idea (counterfactual
sensitivity training) is prior art. What is new is the *diagnosis* that
identifies where to apply it, and the algebraic argument that centering makes
the level branch unable to satisfy it. The contribution should be claimed at
that level and no higher.
