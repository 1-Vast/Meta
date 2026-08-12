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
F-158        UniPert-inspired pair-score difference failed closed
```

F-154 and F-155 are especially important negative controls. They show that a
successful pair-delta model is not enough for the original panel CQ Gate unless
the observable remains target-partner specific after the additive quotient.
F-156 is the attribution correction: the F-153 predictive signal is real, but
not yet attributable to target-conditioned chemical-transformation response.
F-157 symmetry Gate 0 passed and removed the target-main ordering shortcut.
F-158 then failed the registered U1 stop tree, so this UniPert-inspired
representation lane should stop rather than scale model capacity.

## Reproducibility assets to keep

```text
history.md
research/source_affinity/train_chembl_assay_sardelta.py
research/crossed_interaction/train_bindingdb_sardelta_cq_bridge.py
research/crossed_interaction/train_sardelta_potential_cq_observable.py
research/crossed_interaction/train_sardelta_edge_cq_observable.py
research/crossed_interaction/train_bindingdb_sardelta_attribution.py
research/crossed_interaction/audit_bindingdb_sardelta_symmetry.py
research/crossed_interaction/train_bindingdb_pair_score_difference.py
tests/test_train_chembl_assay_sardelta.py
tests/test_train_bindingdb_sardelta_cq_bridge.py
tests/test_train_sardelta_potential_cq_observable.py
tests/test_train_sardelta_edge_cq_observable.py
tests/test_train_bindingdb_sardelta_attribution.py
tests/test_audit_bindingdb_sardelta_symmetry.py
tests/test_train_bindingdb_pair_score_difference.py
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
additive concat `[protein; ligand_delta]` model, or larger UniPert-inspired
pair-score-difference encoder on the same SAR-delta supervision. F-157 passed
the data symmetry audit; F-158 failed U1 and triggers the stop tree. If the
research resumes, the next falsifiable axis should change the supervision
estimand, for example a governed 2x2 rectangle interaction target, rather than
only changing representation capacity. The admission criterion must include:

```text
interaction beats ligand-delta
interaction beats wrong-target
interaction beats shuffled-target
antisymmetry holds for pair deltas
rectangle/additive main effects are explicitly controlled if rectangles are used
component bootstrap LCB95 > 0
support/query isolation preserved
```
