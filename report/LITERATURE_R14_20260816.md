# Literature survey for the R14 cycle (2026-08-16)

Read as inspiration, not authority. Every item below is separated into
**prior art** (what the source establishes), **local adaptation** (what would
have to change for this protocol, this data, these constraints) and **new
contribution** (what would remain genuinely ours and therefore has to be
falsified on its own). Nothing here is evidence about MetaSieve until it is
measured on `bindingdb_ki_double_cold_v1`.

The survey is organised by the four **measured** unresolved causes from
`report/BOUNDARY_20260816.md`, not by field, because the point is to find
mechanisms for causes we already localized.

---

## Cause 1 — the zero-shot calibration/ordering Pareto conflict

The single most-repeated finding of R3R4-R13: every arm that improves the
shape term pays for it in CI, and every arm that improves calibration loses
shape. The premise was measured, not assumed — under ordinary training the
level and shape objectives have gradient cosine **-0.334** on the interaction
trunk and **-0.532** on the ligand encoder (R3R4). The project's response was
**gradient routing** (send each term only to the module that owns it). R11
falsified routing as a general fix: "the level/shape routing trades
calibration for shape on every architecture tested."

### Prior art

**Regression Compatible Ranking (RCR)** — [Bai et al., "Regression
Compatible Listwise Objectives for Calibrated Ranking with Binary
Relevance", arXiv:2211.01494](https://arxiv.org/abs/2211.01494) (CIKM 2023) —
is the closest thing in any field to a *diagnosis* of this conflict rather
than a workaround. Its argument, in our notation:

* a pointwise regression loss is minimised when `link(s_i) -> y_i` — a
  per-item, scale-bearing target;
* a listwise softmax cross-entropy is minimised when
  `exp(s_i)/Σ_j exp(s_j) -> y_i/Σ_j y_j` — a normalized, scale-free target;
* these are **different global minima**, so the two terms fight for the whole
  of training. The conflict is a property of the loss pair, not of the
  architecture or the optimizer.

Their fix is to replace `exp(·)` in the listwise normalizer with the *same
link the pointwise loss uses*, giving a generalized `ListCE(T, s, y) =
-(1/C) Σ_i y_i log[T(s_i)/Σ_j T(s_j)]`. With `T = σ` alongside a sigmoid
pointwise loss, the pointwise optimum `σ(s_i) -> E[y_i]` **automatically**
satisfies the listwise optimum, because `σ(s_i)/Σ_j σ(s_j) -> P_i/Σ_j P_j`
follows from it. The two components become mutually aligned instead of
competing. Deployed in YouTube's production ranking.

Adjacent: [JRC, "Joint Optimization of Ranking and Calibration with
Contextualized Hybrid Model", arXiv:2208.06164](https://arxiv.org/pdf/2208.06164)
and ["Learning to Rank when Grades Matter", arXiv:2306.08650](https://arxiv.org/html/2306.08650)
both formalise the same tension; JRC's remedy is a hybrid two-logit model
rather than an alignment proof, and it is closer to a workaround.

### Local adaptation required

RCR is derived for **binary relevance** with a sigmoid link and a listwise
group = one query's document list. Our labels are **continuous pK** (roughly
4-11) and our group is **one protein's query panel within one episode**. The
translation is not free:

* `T` must be a link whose range is positive on our label support so the
  normalizer is well defined and the alignment identity survives. In pK space
  the labels are already strictly positive, so `T = identity` on a
  positive-shifted score is the natural choice, but the shift is a modelling
  decision that changes the effective weighting of low-affinity ligands, and
  that has to be measured, not assumed;
* our panels are 16 queries, far smaller than a web ranking list, so the
  normalizer's variance is a real concern;
* our regression term is MSE, not sigmoid cross-entropy, so the alignment
  proposition has to be re-derived for the squared-error optimum
  `s_i -> y_i` rather than `σ(s_i) -> E[y_i]`. This is a genuine derivation
  step, and if it does not hold, this direction dies before implementation.

### New contribution, if it survives

Testing **alignment instead of routing**. The project has already measured
that the conflict exists and that routing around it fails on every
architecture. RCR says the conflict is removable at its source. Nobody has
applied that to cold-target DTA, and — more to the point — nobody has applied
it to a case where the *routing alternative has already been falsified with
three seeds*. That negative result is what makes this a sharp experiment
rather than another loss-weight sweep: R11 is the control the literature does
not have.

Explicitly **not** claimed: that this is novel machine learning. It is a
known objective construction applied to a place where we have unusually good
evidence that the problem it solves is the binding one.

---

## Cause 2 — collapse of the direct shape readout under regression training

R13, measured: the MLP shape branch "collapses to near-zero spread under the
shape variance term" (shape std 0.00-0.06). Stage 9, measured: the trained
zero-shot endpoint "has essentially no within-target ligand discrimination"
— re-centring it on the true target mean gives 0.7403 against a flat
constant's 0.7430, and CI 0.525 against a 0.500 coin flip. R9, measured: the
endpoint spread across one episode's queries is 0.087-0.186 pK against a
0.93 pK label spread.

### Prior art

The phenomenon has a name outside chemistry. [Andriopoulos et al., "The
Prevalence of Neural Collapse in Neural Multivariate Regression",
arXiv:2409.04180](https://arxiv.org/abs/2409.04180) establishes **Neural
Regression Collapse**: under MSE with weight decay, last-layer features
collapse onto a low-dimensional subspace whose dimension matches the target
dimensionality, and the collapse is *driven by the regularization* — with
regularization parameters at zero there is no collapse. Our target is
one-dimensional (a scalar pK), which is the worst case: the theory predicts
collapse to a **one-dimensional** feature subspace, i.e. exactly the
"target-level constant plus nothing" degeneracy the project has measured five
times.

The self-supervised literature has the standard counter-mechanism:
[VICReg, arXiv:2105.04906](https://arxiv.org/abs/2105.04906) adds an explicit
hinge on the per-dimension standard deviation, and follow-ups
([orthogonality regularization, arXiv:2411.00392](https://arxiv.org/pdf/2411.00392))
attack the same dimensional collapse. Important caveat for us: **this is what
R13 already tried and it failed.** R13's `shape_variance` term *is* a
variance regularizer, and the branch collapsed anyway. A variance term in the
loss is a soft pressure that the optimizer can trade away; it is not a
guarantee.

On inductive bias: [bilinear MLPs (ICLR 2025)](https://proceedings.iclr.cc/paper_files/paper/2025/file/7504142a20a3e1fe9dd7de42f475828c-Paper-Conference.pdf)
show bilinear forms parameterize pairwise interactions with a low-rank
structure that ordinary ReLU MLPs do not recover. This is consistent with our
own measurement: R13 recorded that with the variance term removed the MLP
branch survives "but still underperforms the bilinear readout on this task".
The bilinear readout was never the problem; its *supervision* was.

Conditioning: [FiLM, feature-wise linear
modulation](https://distill.pub/2018/feature-wise-transformations/) is the
standard way to let one input control the scale of another representation at
negligible parameter cost.

### Local adaptation required

The lesson we should take is **not** "add another variance loss" — that is
R13, already falsified. It is that a *soft* spread penalty cannot beat a
collapse the objective's own optimum is driving. If spread is required, it
should be **architecturally guaranteed** rather than penalized: normalize the
shape score by a statistic of a fixed, learned, protein-independent anchor
set (inductive — the anchors are model parameters, never the query panel),
and let a FiLM-style protein-conditioned scalar set the amplitude.

The anchor-normalization is *not* free of risk: dividing by an anchor spread
that itself shrinks reproduces the collapse one level down, so any such
design needs a gate that measures the anchor spread directly.

### New contribution, if it survives

A shape readout whose **within-target spread cannot be optimized to zero**,
combined with a protein-conditioned amplitude that makes "how steep is this
target's SAR" a learned, protein-specific quantity rather than a global
constant. That second half is the part that also speaks to Cause 3, which is
why these two are being considered as one mechanism rather than two.

---

## Cause 3 — insufficient protein-specific k=0 signal

Measured: wrong-protein gaps cross zero at k=0 for every arm. R0: raw pooled
ESM cosine spans a band of width 0.21 around 0.90 with a **0.024** spread
across the nearest 16 training targets, so `softmax(16*sim)` was near-uniform
by construction; train-only centring widens that spread to 0.238.

### Prior art

The 2025-2026 protein-LM literature converges on the same diagnosis:
**mean-pooling residue embeddings destroys functional specificity.**
[Aggregating residue-level protein language model embeddings with optimal
transport (Bioinformatics Advances 2025)](https://academic.oup.com/bioinformaticsadvances/article/5/1/vbaf060/8088230)
and [Context-Aware Protein Representations Using PLMs and Optimal Transport
(bioRxiv 2026)](https://www.biorxiv.org/content/10.64898/2026.01.24.701517v1.full)
both treat residue embeddings as a *distribution* and compare distributions
(sliced Wasserstein) rather than their means. [Isotropy and Geometry of
Pretrained Protein LMs, arXiv:2510.10655](https://arxiv.org/pdf/2510.10655)
gives the geometric reason the raw cosines are compressed into a narrow band.

### Local adaptation required — and a standing prohibition

This is the direction where the project's own record most constrains us.
Stage 8 (Mac-Diff sequence-window locality prior) **passed 14 structural
gates and still regressed k=0 in 3/3 seeds**. The R0 audit reopened "protein
representation for k=0" only in the narrow sense that the *retrieval*
conclusion was scope-limited to raw pooled cosine — and then immediately
noted that adding the sharper train-centred ESM retriever *did not help*,
"which lowers the prior on protein-representation interventions again."

So: an optimal-transport protein aggregation is a legitimate idea in the
literature and a **poor bet here**, on our own multi-seed evidence. It is
recorded as surveyed and **not selected**. If protein specificity is to
improve in R14, the cheap route is the FiLM-style protein-conditioned
amplitude in Cause 2 — which uses the protein representation we already have
to control a quantity the model currently cannot express at all — not a new
protein encoder.

---

## Cause 4 — the absence of useful k=1 query-specific adaptation

Seven query-specific channels across three families, all deployment-inert
while disturbing calibration, under both MSE-primary and ranking-primary
objectives.

### Prior art

[ActFound (Nature Machine Intelligence 2024)](https://www.nature.com/articles/s42256-024-00876-w)
is the strongest result in this space: pairwise *within-assay* relative
learning plus meta-learning across 35,644 ChEMBL assays, reported to beat
PBCNet "even without using any target protein information." That last clause
is the important one for us, and it cuts against the project's thesis rather
than for it.
[SQRL, arXiv:2501.09103](https://arxiv.org/abs/2501.09103) reformulates
activity prediction as relative-difference learning between structurally
similar pairs, with similarity used to select and weight pairs.
[GraphCliff, arXiv:2511.03170](https://arxiv.org/pdf/2511.03170) and
[activity-cliff-informed contrastive learning](https://pmc.ncbi.nlm.nih.gov/articles/PMC11643338/)
add cliff-awareness as an inductive bias.
On the adaptation side, [Conditional Neural Processes for Molecules,
arXiv:2210.09211](https://arxiv.org/abs/2210.09211) and the TabPFN/TabFM
in-context line show that a frozen set-conditioned model can absorb a small
labelled context in one forward pass with no gradient update — which is the
legal shape of adaptation under our constraints.

### Local adaptation required

Every one of these would be, in our terms, an eighth query-specific channel.
The mandate is explicit that they must not be revived "without a materially
different, falsifiable mechanism", and none of the above supplies one: they
supply *better similarity metrics and better pair selection*, and the project
already has a validated fixed Morgan/Tanimoto residual weighting at k>=2
(Stages 6/7) that does exactly that job.

### Decision

**Cause 4 is not addressed in R14.** The retained few-shot behavior stays the
Stage 6/7 Tanimoto-weighted residual shift at k>=2 and the scalar level shift
at k=1. R14's admission target still requires that adding adaptation must not
degrade k=0, so k>=1 remains measured — but no new k=1 mechanism is proposed,
because the evidence does not support one and the mandate forbids a
cosmetic one.

This is a deliberate narrowing, and it is the reason R14 can keep to two core
contributions.

---

## Surveyed and explicitly not selected

| direction | source | why not |
|---|---|---|
| optimal-transport / distributional protein aggregation | Bioinformatics Advances 2025; bioRxiv 2026 | Stage 8 rejected a protein-representation prior on 3/3 seeds; R0 lowered the prior again |
| CNP / TabPFN-style in-context adaptation | arXiv:2210.09211; TabPFN/TabFM | an eighth query-specific channel; no materially different mechanism |
| ActFound-style cross-assay pairwise meta-learning | Nature Mach. Intell. 2024 | uses external bioactivity corpora; would be external data, not architecture |
| SQRL similarity-quantized difference learning | arXiv:2501.09103 | subsumed by the retained Tanimoto weighting; would be an eighth channel |
| activity-cliff contrastive / cliff-weighted losses | PMC11643338; GraphCliff | R9 measured the cliff weight as a *net negative* for ranking; dose response already run |
| curriculum / difficulty reweighting | curriculum-learning line | a reweighting axis; R9 already showed pair reweighting is not the lever |
| Cartesian / equivariant encoders | PBCNet2.0, TensorNet, PaiNN, MACE, Equiformer | 0 of 17,717 cells have a common-frame complex; closed on data |
| retrieval priors | Stage 10 | falsified by R0; named baseline only |

---

## What this survey selects

Exactly two, mutually reinforcing, one model and one training:

* **A (model)** — a shape readout with an architecturally guaranteed
  within-target spread and a protein-conditioned amplitude, replacing a soft
  variance penalty that R13 measured to be insufficient. Addresses Causes 2
  and 3.
* **B (training)** — a **regression-compatible** within-target ranking
  objective whose optimum coincides with the calibration optimum, replacing
  gradient routing that R11 measured to fail on every architecture. Addresses
  Cause 1.

They reinforce: B removes the reason the level branch had to be routed away
from the trunk, which is what broke calibration in R11; A supplies the
within-target resolution that B's ranking term needs something to push on.
Neither is decoration — the falsification plan in the R14 preregistration
requires each to be shown necessary by a matched ablation.

**No implementation follows from this document.** The premises are tested
first by no-training and synthetic diagnostics (`report/R14_DIAGNOSTICS_*`),
and the two derivations that could kill the design before any GPU time — the
squared-error alignment identity, and whether anchor normalization merely
moves the collapse — are the first things measured.
