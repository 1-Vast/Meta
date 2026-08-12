# MetaSieve research organization - 2026-08-12

## Scope

This pause audit covers the current SAR-delta repair chain after F-155. It does
not delete evidence files. It separates admitted mechanisms, failed mechanisms,
reproducibility assets and future cleanup candidates.

## Identified signals

```text
F-152  ChEMBL source-only SAR-delta signal passed its local Gate
F-153  BindingDB SAR-delta pair bridge beat zero-delta on development components
```

These signals are admitted only for their tested scopes. F-156 later showed
that they do not yet identify target-specific conditioning or a UniPert-style
target x chemical-transformation bridge. They do not authorize V1 integration
or biological claims.

## Failed mechanisms to retain

```text
F-149/F-150  ChEMBL affinity teachers failed closed
F-154        scalar SAR-delta neighbor potential failed closed
F-155        panel edge-distribution SAR-delta observable failed closed
F-156        bilinear target x ligand-delta attribution failed closed
```

F-154 and F-155 are especially important negative controls. They show that a
successful pair-delta model is not enough for the original panel CQ Gate unless
the observable remains target-partner specific after the additive quotient.
F-156 is the attribution correction: the F-153 predictive signal is real, but
not yet attributable to target-conditioned chemical-transformation response.

## Reproducibility assets to keep

```text
history.md
research/source_affinity/train_chembl_assay_sardelta.py
research/crossed_interaction/train_bindingdb_sardelta_cq_bridge.py
research/crossed_interaction/train_sardelta_potential_cq_observable.py
research/crossed_interaction/train_sardelta_edge_cq_observable.py
research/crossed_interaction/train_bindingdb_sardelta_attribution.py
tests/test_train_chembl_assay_sardelta.py
tests/test_train_bindingdb_sardelta_cq_bridge.py
tests/test_train_sardelta_potential_cq_observable.py
tests/test_train_sardelta_edge_cq_observable.py
tests/test_train_bindingdb_sardelta_attribution.py
report/source_affinity/
report/crossed_interaction/
```

## Cleanup candidates

No destructive deletion was performed in this pass. Future cleanup should first
review concrete paths and only delete or archive files that meet one of these
conditions:

```text
1. duplicate smoke output fully summarized by a retained RESULT.json and history.md
2. obsolete checkpoint or prediction artifact not referenced by a Gate record
3. pre-F-152 exploratory artifact superseded by an archived report
```

The current failed Gate scripts and result JSON files should remain in place
until a separate archive commit or path-level deletion review is done.

## Next scientific constraint

Do not continue with another scalar potential, unordered distribution summary,
or additive concat `[protein; ligand_delta]` model. If the research resumes,
the next falsifiable axis should first repair the attribution failure with a
true pair-indexed protein x chemical-transformation model. The admission
criterion must include:

```text
interaction beats ligand-delta
interaction beats wrong-target
interaction beats shuffled-target
antisymmetry holds for pair deltas
component bootstrap LCB95 > 0
support/query isolation preserved
```
