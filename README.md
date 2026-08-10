# MetaSieve-DTA

Mechanism-first few-shot drug-target affinity research with a frozen
probability-law operator.

## Current status

```text
MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED
GEOMETRY_AND_PAIR_COMPATIBILITY_IDENTIFIED
OPEN_BINDINGDB_QUOTIENT_TRAINING_PIPELINE_EXECUTABLE
CQ_R1_DEVELOPMENT_INTERACTION_OBSERVED
CQ_TBASIS_LINEAR_AFFINITY_WITNESS_NOT_OBSERVED
TARGET_COEFFICIENT_HETEROGENEITY_NOT_YET_TESTED
K_SHOT_SECTION_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_FEWSHOT_DTA_MODEL
```

MetaSieve remains a trainable deep-learning bioinformatics research system.
Frozen ESM2 and ligand encoders supply biological states; trainable interaction
statistics are admitted only after partner, nuisance, affinity and transfer
controls. The mathematical core remains the frozen constrained probability-law
operator:

```text
A(F,z) = K(B(z)F(z))
```

The open-data path now runs end to end in development. BindingDB Articles 202608
produced 12,457 governed Ki cells in 320 quotient-positive panels. Strict
document, protein-homology and ligand-scaffold closure leaves 31 conflict
components, but the largest contains 85.86% of cells. This is enough for source
optimization and not enough for a population claim.

The first real open-affinity training run fitted one panel-balanced linear
response on the structurally validated 288D radial T-BASIS. It explained only
`0.000709` of development quotient variance. Correct-pair performance did not
beat zero interaction, foreign ligand or deranged protein with a positive 95%
lower bound. Therefore one population-shared radial affinity direction is not
identified.

The intended endpoint is still unseen-target few-shot DTA. The next distinct
hypothesis is a source-learned target-coefficient subspace with `d<=5`, using
dense open profiling panels for ordinal training and BindingDB/Klaeger for
endpoint-specific quantitative constraints. Evaluation must be target-held-out
and scaffold-held-out at `k=1/2/3/5`; adaptation is restricted to the support
row space and must report rank, conditioning, query coverage and abstention.
This stage is not yet preregistered.

## Repository boundaries

- `theory/FINAL_FROZEN_THEORY/`: authoritative mathematics.
- `model/`: passed mathematical, encoder and geometry primitives; no validated
  assembled few-shot DTA pipeline.
- `scripts/`: passed data, sealing, structure, geometry and governance tools.
- `research/`: preregistered or executed research stages.
- `report/`: machine Gates, current status and evidence summaries.
- `history.md`: chronological failure and decision ledger.

Production `model/`, biological `z`, CSMO, Band and the law operator were not
changed by the BindingDB experiment. Davis, KIBA, recipient labels and external
confirmation remain closed.

## Read first

1. `report/CURRENT_RESEARCH_STATUS.md`
2. `report/crossed_interaction/OPEN_DATA_TRAINING_AND_FEWSHOT_ROUTE.md`
3. `report/crossed_interaction/cq_r2_tbasis_linear/weights.report.json`
4. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
5. `task.md`
6. `history.md`

## Verification

```powershell
conda run -n drug python -m pytest -q
```

Large third-party releases, embedding banks and caches are not redistributed;
see `DATA_AVAILABILITY.md`. Current consolidated regression: **224 passed**.
