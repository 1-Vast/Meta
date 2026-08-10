# MetaSieve current research status

Updated: 2026-08-10.

## Current verdict

```text
OPEN_BINDINGDB_QUOTIENT_TRAINING_PIPELINE_EXECUTABLE
CQ_R1_DEVELOPMENT_INTERACTION_OBSERVED
CQ_TBASIS_LINEAR_AFFINITY_WITNESS_NOT_OBSERVED
TARGET_COEFFICIENT_HETEROGENEITY_NOT_YET_TESTED
K_SHOT_SECTION_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_FEWSHOT_DTA_MODEL
```

BindingDB Articles 202608 now supports governed development training: 12,457
Ki cells, 320 cycle-positive panels, 31 strict conflict components, train
quotient rank 6,608 and development quotient rank 220. The largest component
share is 0.8586, so this is optimization data, not independent population
evidence.

The first real training run fitted one shared 288D T-BASIS linear response. It
explained 0.000709 of development quotient variance; every correct-versus-null
or partner control interval includes zero, and deranged protein is slightly
better at the point estimate. The shared radial direction is therefore not
identified.

## Next scientific question

The only high-value minimal route is target-coefficient heterogeneity on the
same frozen basis, learned from dense open profiling panels and evaluated as
target-held-out `k=1/2/3/5` episodes. Use Kinobeads/PKIS/PKIS2 for ordinal
pretraining, BindingDB/Klaeger for endpoint-specific quantitative constraints,
and PDSP for a non-kinase development stratum. Do not mix modalities or add a
larger representation before this test.

Any few-shot correction must be restricted to the support row space and report
rank, conditioning, query coverage and abstention. Raw pair maps and arbitrary
latents cannot enter `z`. The frozen operator `A(F,z)=K(B(z)F(z))` is unchanged.

## Current boundaries

- BindingDB training evidence is development-only; strict closure is dominated
  by one component.
- The 288D feature generator is structurally validated, not affinity-admitted.
- Davis, KIBA, recipient labels, heldout-B and external confirmation remain
  closed.
- No production `model/` or `scripts/` interface was changed.
- Historical S3R/S4R/S5D/C0/C1/X1A results remain in the evidence ledger,
  `history.md` and Git; they are not current execution authority.

## Read first

1. `report/crossed_interaction/OPEN_DATA_TRAINING_AND_FEWSHOT_ROUTE.md`
2. `report/crossed_interaction/cq_r2_tbasis_linear/weights.report.json`
3. `dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus/manifest.json`
4. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
5. `history.md`
