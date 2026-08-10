# MetaSieve-DTA

Mechanism-first few-shot drug-target affinity research with a frozen
probability-law operator.

## Current status

```text
MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED
GEOMETRY_AND_PAIR_COMPATIBILITY_IDENTIFIED
EXACT_RESIDUE_LOCALISATION_IDENTIFIED_IN_DEVELOPMENT
TEACHER_LIGAND_CONDITIONALITY_IDENTIFIED
S2R_SYNTHETIC_BINARY_ORDINAL_ESTIMATOR_PASS
S4R_A_LIGAND_MEAN_POOLING_COLLAPSE_MEASURED
S4R_GRAPH_AWARE_INCREMENT_REAL_BUT_NOT_LIGAND_SPECIFIC
REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED
S5D_LIGAND_STEERING_PRESENT_BUT_BIOLOGICALLY_MISDIRECTED
EXACT_EDGE_COUPLING_NOT_SUPPORTED_BY_TEACHER
AFFINITY_ENERGETICS_NOT_IDENTIFIED
K_SHOT_SECTION_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

MetaSieve remains a trainable deep-learning bioinformatics research system:
frozen protein and ligand encoders supply biological measurements, and small
trainable interaction heads are admitted only after falsification controls. Its
main mathematical object remains the frozen constrained probability-law
operator `A(F,z)=K(B(z)F(z))`.

The latest real structural experiment did not admit a new biological statistic.
S4R replaced only the mean-pooled 41-D ligand statistic with a frozen radius-1
Morgan per-heavy-atom statistic and held every other surface byte-identical to
S3R. The candidate scored `AP_bidir = 0.046856` against chance `0.025472`,
doubling S3R's above-chance gain and beating its capacity-matched
permuted-label learner. It still failed the registered `+0.05` R1 margin, and a
foreign ligand pair cost only `+0.000644`, so the recovered signal is a
construct-level residue-change prior rather than ligand-conditioned residue
selection. The pose-free ligand representation repair route is closed.

S5D then trained nothing and asked why. It registered the obvious mechanism —
that the estimator collapses ligand differences onto one residue direction per
protein — and falsified it: the true-versus-foreign residue field cosine is
`0.4487`, so the estimator plainly does steer on the ligand. An estimand that
cancels pocket membership exactly then found nothing either. Ligand information
arrives at full strength and points somewhere biologically wrong; the
conditional estimand route is closed too.

C0/C1 then tested the last standing hypothesis — correspondence — audit-only on
1,862 systems from raw RCSB coordinates that no stage had touched, excluding
24,874 exposed PDB ids. Every C0 admissibility Gate passed. C1 did not: exact
within-slot AP is `0.9856` against a `0.9540` degree-preserving null, leaving
only `0.0144` of headroom above a pure contact-degree predictor. Exact
atom-residue correspondence is very nearly a function of residue contact
degree, so the geometry-gated router was closed before any training.

## Repository boundaries

- `theory/FINAL_FROZEN_THEORY/`: authoritative mathematics.
- `model/`: passed mathematical, encoder and geometry primitives; no validated
  assembled DTA pipeline.
- `scripts/`: passed data, sealing, structure, geometry and governance tools.
- `research/`: preregistered or executed research stages, including S2R/S3R/S4R.
- `report/`: machine Gates, current status and evidence summaries.
- `history.md`: chronological failure and decision ledger.

No active training stage is authorized. Affinity values, heldout-B, R6,
few-shot adaptation, biological `z`, CSMO/Band and the frozen operator remain
untouched.

## Read first

1. `report/CURRENT_RESEARCH_STATUS.md`
2. `report/correspondence_router/C0_C1_EVIDENCE_CONSOLIDATION.md`
3. `report/s7_l2b_r0r/PHASE2B_S5D_EVIDENCE_CONSOLIDATION.md`
3. `report/s7_l2b_r0r/PHASE2B_S4R_EVIDENCE_CONSOLIDATION.md`
4. `report/s7_l2b_r0r/PHASE2B_S4R_REPRESENTATION_AUDIT.md`
5. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
6. `task.md`
7. `experiment.md`
8. `history.md`

## Verification

```powershell
conda run -n drug python -m pytest -q
```

Large third-party releases, embedding banks and caches are not redistributed;
see `DATA_AVAILABILITY.md`. Current consolidated regression: **193 passed**.
