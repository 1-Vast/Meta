# Method-ladder closure map (2026-08-18, repaired 2026-08-18 post-completion)

The method-ladder cycle listed eight named families. The 2026-08-17/18 research
cycle superseded the ladder; this map records, for each family, what was
**actually implemented and tested**, and how far the resulting verdict reaches.

The first version of this file marked six families "closed by measurement".
That was an overstatement: for most of them the named method's **defining
operator was never implemented**. What was measured was a nearby experiment
that shares a motivation. The post-completion review
(`report/POST_COMPLETION_REVIEW_20260818.md`) requires the distinction, so this
table now carries an explicit evidence class per family.

## Evidence classes

| class | meaning | what the verdict licenses |
|---|---|---|
| **direct** | the named method's defining operator was implemented and trained on the governed protocol | a claim about that method, at this scale, protocol and budget |
| **partial** | the defining mechanism was implemented in an adapted form, with a named component of the original absent | a claim about the implemented variant only; the original remains untested |
| **proxy** | a different operator addressing the same motivation was tested; the named method was never instantiated | a claim about the motivation being unrewarded *on the tested operator* — **not** a falsification of the named method |
| **not instantiated** | no implementation exists | nothing |

No family in this table is in the **direct** class. That is the honest state of
the ladder: it was superseded by the stage programme before any of its named
methods was built.

## The eight families

| # | family | what was actually implemented | measured successor | evidence class | verdict, at its true scope |
|---|---|---|---|---|---|
| 1 | multimodal representation collapse + basis reallocation | contrastive coembedding (InfoNCE) and regression alignment, trained 3 seeds; collapse metric measured 0.99859 -> 0.9908. **Basis reallocation was not implemented** — no dimension/capacity budget was reallocated between modalities | Stage K/K2 | **partial** | K-REG gave the first all-k resolved MSE improvement across 3 seeds; the k=0 centered gain did not survive pooling -> NOT CONFIRMED, nothing promoted. Collapse is reduced but not removed. Basis reallocation remains untested |
| 2 | Gradient Blending / OGM | orthogonal level/shape routing (Stage E) and a decoupled frozen-feature head (Stage Q). **OGM's on-the-fly gradient modulation and Gradient Blending's held-out overfitting-to-generalization weights were not implemented** | Stage E, Stage Q | **proxy** | proxy negative; direct method not instantiated. Two routing/decoupling operators failed the ranking gate. This does not falsify OGM or Gradient Blending, whose defining operator is a per-modality *gradient rescaling*, not a routing constraint |
| 3 | Disentangled Gradient Learning | the Stage E routing ablations, i.e. the same operator as family 2. **No gradient-disentanglement objective was implemented** | Stage E routing ablations | **proxy** | proxy negative; direct method not instantiated. Same evidence and same limit as family 2 |
| 4 | attention MIL / Set Transformer / adaptive pooling | a panel-set level head over `cat(protein summary, panel mean, panel max)` — **fixed order-invariant mean/max pooling followed by an MLP**. No gated attention-MIL weights, no induced set-attention block, no pooling-by-multihead-attention | Stage E/J/L/Q panel and level heads | **proxy** | proxy negative; direct method not instantiated. Learned *panel-conditioned level* improves k=0 calibration and degrades k>=1 ranking in four compositions. Whether a learned attention pooling would behave differently is untested |
| 5 | DrugBAN-style bilinear interaction | a learned pairwise edge MLP over (query embed, support embed) **ligand-ligand** pairs, with fixed Tanimoto as an additive anchor. **DrugBAN's bilinear attention over drug-substructure x protein-residue pairs was not implemented**; the pairing is a different one | Stage F | **proxy** | proxy negative; direct method not instantiated. The tested ligand-ligand kernel is inert against fixed Tanimoto — the fifth learned-kernel family to fail. DrugBAN's drug-protein bilinear map is untested here |
| 6 | FS-CAP-style episodic scale | an assay-aware zero-shot level head (journal/panel/protein covariates) plus a paired cross-target level alignment term. **FS-CAP's defining mechanism — encoding support compounds conditioned on the target and predicting the query from a support compound-activity context — was not implemented** | Stage J | **proxy** | proxy negative; direct method not instantiated. The paired term added nothing measurable and the assay-aware head was rejected on ranking. FS-CAP's support-conditioned encoder is untested |
| 7 | AdaMBind-style task valuation / label-noise robustness | task valuation **was** implemented and tested: the Stage A `A2` selector scored candidate tasks by post-adaptation query loss and support/query gradient cosine, inside a real inner/outer-loop trainer. **Label-noise robustness was not implemented** | Stage A/B | **partial** | the selector is REJECTED on measured cause, not noise: it chose tasks with gradient cosine +0.9897 against a +0.6555 population mean, i.e. redundancy rather than informativeness. The framework was NOT PROMISING then REJECTED. The label-noise half of the family is untested |
| 8 | MMP-cliff transformation learning | pairwise signed-gap supervision over embed pairs (Stage F), and the Stage L2 measurement of a protein-independent signed-SAR direction (r +0.270). **Pairs were not matched-molecular-pair identified**, so the "MMP" constraint of the family is absent | Stage F, Stage L2 | **partial** | the implemented pairwise operator is inert against fixed Tanimoto and cliff sign never improved by a resolved amount. A genuinely MMP-constrained transformation learner is untested |

## What this map does and does not close

**Does:** record that the ladder is no longer an open work item — every family
has a measured successor experiment with a recorded verdict, and no family is
waiting on a decision.

**Does not:** falsify OGM, Gradient Blending, Disentangled Gradient Learning,
Set Transformer, attention MIL, DrugBAN or FS-CAP. Their defining operators
were never instantiated in this repository. Reopening any of them is a new
experiment, not a repeat, and the correct standard for reopening is the one in
the review: state the new function class or information source that makes the
test materially different from the proxy that failed.

The ladder harness (`tools/research/method_ladder/_shared/`) is retained as
tooling.
