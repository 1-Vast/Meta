# UBSE-G0P same-scaffold panel topology decision

Date: 2026-07-29  
Decision: `STOP_UBSE_SAME_SCAFFOLD_PANEL_TOPOLOGY_INADEQUATE`

## Result

The frozen label-blind same-scaffold panel gate fails two of seven criteria:

| Metric | Result | Gate |
| --- | ---: | --- |
| Panels / pair contrasts | 1,612 / 4,691 | pass |
| Exact targets / homology components | 993 / 494 | pass |
| Scaffolds / PubMed IDs | 965 / 904 | pass |
| Conflict components | 384 | scale part passes |
| Largest conflict component | 28.0397% | fail (`<=20%`) |
| Resource ceiling | 494 | pass (`>=423`) |
| Deterministic conflict-free packing | 452 | pass (`>=88`) |
| Largest homology panel share | 6.7618% | fail (`<=5%`) |
| Largest scaffold / PubMed share | 1.8610% / 4.0323% | pass |
| Frozen audit panels | 88 | pass |
| Residual closed training panels / contrasts | 1,524 / 4,551 | pass |
| Exact ChEMBL-TRAIN accession support | 25.277% | pass |

The audit read only identity and provenance metadata. It did not load a
binding-residue field, affinity field/value, coordinate, protected feature or
label, or sealed outcome.

## Interpretation

The result does not fail for lack of panels or conflict-free resources.
Unlike the exact-difference and recurring-edit routes, it has more than the
423-unit predictive resource floor and leaves a large train substrate after
closing an 88-panel audit set. Its binding defect is a single overrepresented
homology resource that also helps connect a 28% transitive conflict component.

G0P itself remains a formal STOP; its thresholds cannot be relaxed and its
result cannot authorize a student.

One bounded correction is scientifically admissible because it removes,
rather than synthesizes or reweights, the concentrated nuisance block:

- remove every homology component whose pre-removal panel share exceeds the
  unchanged 5% cap;
- recompute every unchanged G0P gate once;
- forbid any further threshold, scaffold, PubMed, or component-specific
  pruning selected from the new result.

This correction must be separately frozen as G0PB before execution. Failure
closes the current binding-residue-list student route. Passing may authorize
only an additive-null-controlled, affinity-blind student preregistration.

## Artifacts

- Preregistration:
  `reports/active/ubse_g0p_same_scaffold_panel_preregistration_2026-07-29.md`
- Result: `reports/active/ubse_g0p_seed1729.json`
- Panel manifest:
  `dataset/public/biolip2/processed/ubse_g0p_panels.parquet`
- Implementation: `research/ubse_g0p.py`
- Tests: `tests/test_ubse_g0p.py`
