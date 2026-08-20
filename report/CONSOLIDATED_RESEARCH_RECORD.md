# MetaSieve Consolidated Research Record

**Authority date:** 2026-08-19  
**Scope:** BindingDB-Ki cold-target DTA development, few-shot adaptation, and
the separately governed protein-conditioned interaction investigation.

This is the single narrative summary for the root `report/` directory. It
replaces duplicated root-level plans, literature reviews, handoffs, boundary
summaries, and completion narratives. Numerical claims remain auditable in
`EVIDENCE_LEDGER.md` and in the leaf stage artifacts; this file is not a
replacement for frozen preregistrations or machine-readable results.

## 1. Current Decision

The project has not produced a reproducible model with MSE <= 1.00 pK^2 on
the strict double-cold BindingDB-Ki protocol. The current evidence supports a
useful ligand-side few-shot signal, but does not establish a transferable
protein-conditioned interaction signal. These are separate conclusions:

- performance value is not rejected merely because the protein-mechanism claim
  is unresolved;
- protein-conditioned claims require a valid positive control and a
  cross-component, label-bound result;
- `meta_test` remains governed as a held-out split and must not be used for
  tuning.

The active work is the Q2d synthetic qualification ladder. Its purpose is to
decide whether the measurement and learning harness can recover a deliberately
identifiable interaction before any biological interpretation. While Q2d is
running, no production `model/` or `scripts/` change is authorized.

## 2. What Has Been Established

### 2.1 Positive and reusable findings

1. Fixed Morgan/Tanimoto support weighting recovers a real within-target
   ligand-neighborhood signal and improves few-shot ranking at `k >= 2` in
   matched within-checkpoint tests. It is a comparator and possible practical
   baseline, not evidence of protein conditioning.
2. Pairwise/ranking-oriented training produced the strongest observed
   within-target shape and activity-cliff signal, but calibration and ranking
   compete under the strict protocol. The result is not a uniformly superior
   production model.
3. Local ESM windows around a validated mutation/pocket coordinate carry more
   mutation-selective information than globally pooled ESM. The original
   distance-ratio instrument was invalid; the corrected X0c qualification
   preserves the capability conclusion only where its new definition passes.
4. The repository has strong data-contract, split, permutation, censoring,
   and audit tests. Failures found in training and evaluation have generally
   been reproducible implementation or estimand failures, not silently
   ignored.

### 2.2 Rejected or unresolved claims

- Output-side query-specific adapters, reliability gates, and several
  HyperSAR/TERM variants were repeatedly deployment-inert or harmful.
- A protein embedding used as a per-target key does not demonstrate
  transferable protein chemistry. Shuffled-protein controls reproduced the
  apparent gain in the prior protein-conditioned SAR experiment.
- Cartesian/PBCNet2/TensorNet-style complex geometry is not a valid active
  path for the BindingDB main bank: there are no usable common-frame
  protein-ligand complexes in the governed cells.
- Cross-target SAR fields built from unconstrained ID-space factors fail cold
  ligand generalization. This is a harness diagnosis, not evidence that all
  biological interaction fields are absent.
- The strict protein-conditioned interaction task is currently unresolved,
  not biologically falsified. A valid external positive control and an
  identifiable synthetic qualification remain prerequisites.

## 3. Evidence Timeline

### R0-R8: few-shot and transport search

The early cycles tested level calibration, SAR transport, support-conditioned
gates, retrieval, ranking losses, and chemistry-aware controls. The stable
lesson was that target-level calibration is large, while many learned
query-specific mechanisms either collapse to uniform weighting or perturb
calibration without improving ranking. Fixed ligand similarity was the first
reproducible positive transport comparator.

### R9-R15: boundary, literature, and programme consolidation

The boundary cycles separated practical performance from mechanism attribution.
They showed that absolute cross-assay pK calibration is a major nuisance;
document and assay effects cannot be inferred from protein and ligand inputs
alone. Literature review therefore shifted the programme toward paired
contrasts, local pocket representations, censoring-aware estimands, and
separate practical and mechanism tracks.

### Stage X0-X0c: measurement qualification

The variant-coordinate layer was corrected and externally checked. Historical
numbering, construct offsets, old-residue agreement, stable hashing, KLIFS
coverage, and restricted-data handling are now explicit. The original I2
distance-ratio test was invalid because WT and mutant windows were extracted
at different coordinates and the edit token had a zero cross-protein
denominator. X0c preserves the corrected mapping evidence but records the
one-hot planted-signal harness failure: oracle recovery is possible, while the
tested representation/training path is not yet sufficient.

### Q2d: current qualification ladder

Q2d-1 demonstrated that an unconstrained random protein-row x ligand-column
factor is not cold-generalizable. Q2d-1b/1c/1d progressively test
feature-conditioned factors, span restrictions, oracle mappings, and runner
dataflow. A failed rung must close the corresponding interpretation; it must
not be repaired by changing the frozen gate after observing results.

## 4. Why the Main Metric Is Difficult

The strict protocol combines protein-component cold splitting, ligand/scaffold
coldness, `k=0..5`, absent common-frame geometry, heterogeneous assay/document
calibration, and absolute pK MSE. This is a valid stress test but not the only
practical DTA setting. A useful future programme must report two tracks:

1. **Practical track:** protein-cold targets with `k=5/10/20/40`, realistic
   chemical-series support/query relationships, and explicit assay metadata
   where available.
2. **Mechanism track:** strict double-cold, paired contrasts, shuffled and
   matched-wrong controls, cluster bootstrap, and frozen positive controls.

Failure on the mechanism track blocks the protein-interaction claim; it does
not erase a valid ligand-only practical result.

## 5. Method and Biology Implications

The most defensible biological route is not to copy a large 3-D complex model
without poses. It is to use validated pocket-local sequence/structure priors,
matched WT-mutant or ligand-identical contrasts, and a model whose interaction
component is directly supervised and separately identifiable from protein and
ligand main effects. Relevant literature families include kinase PCM,
mutation-to-DeltaDelta prediction, local protein-language-model variant
scoring, paired metric learning, set/attention meta-learning, and
censoring-aware bioactivity modelling. These are inspirations and baselines;
they do not substitute for an in-protocol positive control.

The next biological qualification must use a same-platform, ligand-identical
WT/variant panel where possible. Percent remaining activity is not silently
treated as pK; ATP/Km and assay context remain explicit covariates. Any external
dataset must be audited for censoring, construct identity, mutation numbering,
and license restrictions before model training.

## 6. Current Gates And Next Actions

1. Finish the frozen Q2d qualification ladder and record PASS/FAIL with all
   negative controls.
2. If and only if the feature-conditioned synthetic positive control passes,
   run the corrected local-ESM representation arm under a new preregistration.
3. Build the legal Saifudeen/kinase panel census before training; count usable
   paired rows after censoring and construct matching.
4. Run B1/B2/C in order: same-study positive control, pocket-local attribution,
   then cold-variant generalization. Do not connect the module to production
   DTA before these gates pass.
5. In parallel, benchmark practical few-shot baselines per dataset without
   merging incompatible assay labels. Report RMSE/MSE, CI, Spearman, ranking
   and component-level bootstrap intervals.

## 7. Evidence And File Policy

- `EVIDENCE_LEDGER.md` is the compact numerical index.
- `CURRENT_MODEL_EVIDENCE.md` is the current baseline evidence.
- `CORE_TASK1_UNRESOLVED_TERMINAL_20260817.{md,json}` and the completion
  manifest are retained as machine-auditable terminal records.
- `report/mechanism/`, `report/meta_fewshot/`, and `tools/research/` contain
  leaf preregistrations, results, tests, and hashes; they are not deleted as
  part of narrative cleanup.
- Restricted raw datasets stay outside Git. Manifests, hashes, and semantic
  audits are retained.

## 8. Superseded Root Narratives

The following root-level narrative files were merged into this record and are
removed from the active root after references are updated:

`AGENT_HANDOFF_A2_MOMENT.md`, `BOUNDARY_20260816.md`,
`BOUNDARY_20260817_NIGHT.md`, `COMPLETION_STATEMENT_20260818.md`,
`COMPLETION_STATEMENT_CORE_TASK1_20260817.md`,
`core_task_1_protein_conditioned_signal.md`, `FINAL_STATE_20260818.md`,
`LITERATURE_R14_20260816.md`, `LITERATURE_RESEARCH_SYNTHESIS_20260815.md`,
`NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md`,
`POST_COMPLETION_REVIEW_20260818.md`,
`protein_conditioned_signal_investigation.md`.

The recent Stage X qualification documents and the latest programme review
remain separate until the active Q2d work is complete; they are current
inputs, not redundant historical summaries.
