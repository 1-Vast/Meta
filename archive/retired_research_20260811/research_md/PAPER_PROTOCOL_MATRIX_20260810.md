# Paper-locked protocol matrix

Date: 2026-08-10

This document separates literature reproduction, the MetaSieve main benchmark,
and strict confirmation. Papers are methodological references, not a template
to copy wholesale. No missing paper detail may be silently replaced by a local
choice, and no borrowed convention is exempt from a MetaSieve-specific reason.

## Lane A: AdaMBind reproduction

Primary source: https://www.nature.com/articles/s41467-026-70554-5

| Item | Frozen paper-described setting |
|---|---|
| task | one protein target per meta-task |
| datasets | BindingDB, KIBA and Davis as released by the authors |
| random task split | 8:1:1 train/validation/test |
| novel-target split | CD-HIT clusters at 40% sequence identity, then 8:1:1 |
| headline shots | 5 and 40 |
| support sweep | 5, 10, 20, 30 and 40 |
| repeats | 5 independent runs |
| metrics | MSE, concordance index, R2, Spearman and Pearson |
| optimizer/loss | Adam and MSE; exact remaining values must come from the release |

The paper uses MAML, an adaptive task sampler and label-noise augmentation.
MetaSieve may compare with this lane, but cannot rename its closed-form ridge
as an AdaMBind reproduction.

## Release audit blockers

Official code: https://github.com/Moohyun-w/AdaMBind

Versioned archive DOI: https://doi.org/10.5281/zenodo.18595084

The locally audited public snapshot is commit
`01a169a6d62fba0d6c003f47bfba539e55f5b344`; its acquisition manifest is
`dataset/raw/adambind_public_01a169a6/acquisition_manifest.json`.

- The Zenodo files are access-restricted as of this audit, so the versioned
  archive cannot currently be downloaded or hashed locally.
- Public `train.py` defaults to `nums=10`, not the paper's 5/40 headline.
- Fixed split files are used only when `data/{dataset}_{1,2,3}.txt` exist;
  otherwise the code creates a random 8:1:1 split.
- The public repository does not expose the described CD-HIT split manifests
  or their construction command.
- The public training code combines the initial train and validation indices
  before resampling meta-train/meta-validation tasks. This conflicts with a
  literal independent validation reading of the paper.
- The paper describes 10 outer epochs of 100 meta-updates, while public code
  defaults to five outer iterations.
- The paper describes a three-layer protein 1D CNN, while the public model path
  contains one Conv1d call.
- Public training and evaluation paths index different label columns in places;
  the tensor schema must be tested before any result can be trusted.
- No license file was present in the pinned public repository tree. Downloaded
  files therefore remain isolated and are not redistributed by this project.

Until the versioned archive resolves these items, the only honest label is
`paper-described reimplementation`, not `exact reproduction`. An exact-code
lane and a corrected-protocol lane must be reported separately.

## Lane B: MetaSieve literature-aligned main benchmark

This lane borrows only conventions that improve comparability or data quality:

- protein target as the meta-task and target-level 8:1:1 splits from AdaMBind;
- a 40% protein-similarity novel-target profile as the headline generalization
  setting;
- single-protein, valid small-molecule, measurement-type separation and
  within-assay replicate median principles from CARA;
- support/query episodes and per-target metric distributions;
- RDKit molecular graphs as a conventional ligand input.

It deliberately does not copy AdaMBind's MAML, adaptive sampler, one-hot/CNN
protein encoder, mixed benchmark endpoints or scalar head. It does not copy
CARA's assay-as-task formulation, 50-shot support, VS/LO classifier, or
ChEMBL-specific field layout. MetaSieve keeps exact Ki as its primary estimand,
uses its biological interaction representation, the `d<=5` identifiable
Meta-Section, and eventually the separately gated CSMO law output.

Every borrowed item must have a source and applicability statement. Every
deviation must identify whether it is required by the estimand, the core
innovation, or a repository governance constraint, and must have a matched
control when it can affect performance.

### Open assay-to-target contract

CARA aggregates replicates within a ChEMBL assay, while AdaMBind groups all
observations of a protein into one task. BindingDB can contain the same
target-ligand pair under multiple document/protocol panels. These rules do not
compose automatically. Before materialization, MetaSieve must preregister one
of the following and test its sensitivity:

- retain panel-specific observations but keep every ligand wholly on one side
  of support/query;
- take a median within panel, then an equal-panel median for the target-ligand
  task label;
- restrict the main benchmark to an official processed paper release whose
  pair-level aggregation is already defined.

No option may be selected by the resulting target count or model score.

## Lane C: MetaSieve strict confirmation

This is governed by `PREREG_CORE_META_SECTION.md`, not AdaMBind or CARA:

| Item | Frozen MetaSieve setting |
|---|---|
| endpoint | exact, positive, uncensored Ki only; pK=9-log10(nM) |
| split | document/scaffold/protein-40 dependency closure |
| shots | 1, 2, 3 and 5 nested supports |
| support controls | correct, zero, foreign and permuted labels |
| biology controls | full, wrong-protein and capacity-matched ligand-only |
| inference unit | target, with dependency-component sensitivity |
| core family | frozen 288D T-BASIS to d<=5 positive-ridge section |

This lane is a leakage-resistant stress test, not the sole training benchmark
and not a reproduction of AdaMBind, CARA, MetaDTA, FS-CAP or R2-D2. Those
papers justify comparison arms and protocol precedents only.

## Reference policy

The intended division is not a numerical 80/20 copying rule. It is an ownership
rule:

| Surface | Default owner |
|---|---|
| baseline reproduction | the cited paper and its pinned release |
| generic chemistry parsing | established RDKit/paper convention |
| main task and endpoint | MetaSieve scientific estimand |
| protein/interaction representation | MetaSieve innovation |
| support adaptation | MetaSieve identifiable Meta-Section |
| probabilistic output | MetaSieve CSMO, only after its Gates |
| strict confirmation | MetaSieve leakage and inference governance |

Similarity to a paper is not evidence by itself; novelty is concentrated in
the scientific mechanism and adaptation contract, not manufactured through
arbitrary preprocessing differences.

## Method boundaries

- MetaDTA is an attentive neural-process DTA model, not MAML. Its support sizes
  and pIC50 task differ from the MetaSieve Ki k<=5 estimand.
- R2-D2 validates differentiable closed-form ridge in few-shot image
  classification; applying it to DTA is a hypothesis, not published DTA proof.
- The final FS-CAP paper must be cited separately from its older preprint;
  preprint dataset/result numbers cannot be presented as final-paper results.
