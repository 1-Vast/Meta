# Meta-learning and protein-ligand identifiability review

Date: 2026-07-31

> **Decision update:** use [`innovation_gate_decision_2026-07-31.md`](innovation_gate_decision_2026-07-31.md) for the corrected strict v2 audit and IDG-RBP result. The evidence review below remains valid, but its pre-v2 rectangle counts are historical.

## Scope and sources

The supplied paper was checked against the open full text, its source-data record,
and the authors' public implementation before any code change:

- Wan et al., **AdaMBind**, Nature Communications 17, 3734 (2026), DOI
  [`10.1038/s41467-026-70554-5`](https://doi.org/10.1038/s41467-026-70554-5),
  PMC [`PMC13102954`](https://pmc.ncbi.nlm.nih.gov/articles/PMC13102954/).
- AdaMBind source repository at commit `01a169a6d62fba0d6c003f47bfba539e55f5b344`:
  [`Moohyun-w/AdaMBind`](https://github.com/Moohyun-w/AdaMBind).
- AdaMBind source data, Figshare DOI
  [`10.6084/m9.figshare.30963823.v1`](https://doi.org/10.6084/m9.figshare.30963823.v1).
- Finn et al., **Model-Agnostic Meta-Learning**, ICML 2017,
  [`arXiv:1703.03400`](https://arxiv.org/abs/1703.03400).
- Oreshkin et al., **TADAM: Task dependent Adaptive Metric for Few-Shot
  Learning**, NeurIPS 2018,
  [`arXiv:1805.10123`](https://arxiv.org/abs/1805.10123).
- Patacchiola et al., **Learning to Learn, a Bayesian Perspective** / Deep
  Kernel Transfer, NeurIPS 2020,
  [`10.48550/arXiv.2008.05414`](https://doi.org/10.48550/arXiv.2008.05414).

## What AdaMBind actually changes

AdaMBind uses a conventional MAML-style DTA predictor (a molecular GNN and a
protein 1D CNN) and adds an adaptive task sampler. The public `Scheduler` scores
candidate tasks from query loss and support/query gradient cosine similarity,
then samples tasks with a categorical policy. The policy is trained from
validation query loss. This is an allocation policy over meta-training tasks;
it is not a new protein-ligand interaction operator or a source of new labels.

The paper reports two task splits: random and a CD-HIT 40%-sequence-identity
novel-target split. Each novel target still receives either 5 or 40 labeled
same-target support pairs. The public `DataSplit.py` randomly divides rows within
each target into support and query. The reported protocol does not close ligand
scaffolds/chemical neighbours, assay or document identity, endpoint, or source
provenance between support and query. Therefore its gains answer a different
estimand from FORT's strict dual-cold task.

The paper's own caveat is also material: it uses a simple graph encoder without
explicit bond-edge attributes and presents the scheduler as an adaptation aid.
None of these details establish that the protein branch changes a ligand's
relative affinity ordering under a crossed-target design.

## How the prior art maps to FORT

| Mechanism | What it can legitimately provide here | Why it is not yet a primary fix |
| --- | --- | --- |
| MAML inner/outer updates | A capacity-matched gradient baseline for k=5 | Existing FORT MAML did not beat calibration; gradients cannot recover absent protein information |
| TADAM task-conditioned metric | A task-dependent rescaling of an already meaningful embedding | It assumes a transferable task metric; wrong-protein invariance means that assumption is unverified |
| AdaMBind easy-to-hard scheduler | A secondary task-allocation ablation after the estimand passes | Scheduler changes which labels are seen during training, not the information in a fixed crossed panel; query-label policy training must remain training-only |
| Deep Kernel Transfer / Bayesian posterior | Exact support-conditioned inference and uncertainty | Current finite-rank Bayesian route already failed ligand-only and gradient controls |
| Contrastive meta-learning | A possible representation objective | Contrastive positives/negatives must have a biological identity and provenance contract; arbitrary target similarity would encode shortcuts |

The architecture consequence is intentionally conservative. A TADAM-like metric
or AdaMBind scheduler is admissible only as an equal-budget secondary branch after
the protein-conditioned residual probe shows a positive correct-protein increment.
It cannot be used to rescue a failed protein path.

## Current evidence and decision rule

The existing AnchorDelta retraining improved aggregate ranking but was unchanged
under wrong-protein replacement. The present TRAIN-only audit now reports exact
ligand rectangles, reversal rates, repeat noise, and homology-pair bootstrap units
without using development or sealed labels. The decisive gate is not the raw
rectangle count: it is source/assay comparability, component concentration, and a
correct-protein advantage over protein-free and homology-matched wrong-protein
controls.

If the audit fails strict comparability or the low-capacity residual probe fails,
the correct action is to stop protein-conditioned architecture expansion. If it
passes, the next implementation is a rank-4/8 bilinear residual probe with exact
antisymmetry and nested target/homology-cold splits; only then may a scheduler
ablation be added. No result from AdaMBind's random or 40%-identity benchmark is
used as evidence for the strict FORT claim.
