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
F-159        PC-SAR source-only free target-specific oracle failed closed
F-161        exact transformation-ID rectangle low-rank oracle failed closed
```

F-154 and F-155 are especially important negative controls. They show that a
successful pair-delta model is not enough for the original panel CQ Gate unless
the observable remains target-partner specific after the additive quotient.
F-156 is the attribution correction: the F-153 predictive signal is real, but
not yet attributable to target-conditioned chemical-transformation response.
F-157 symmetry Gate 0 passed and removed the target-main ordering shortcut.
F-158 then failed the registered U1 stop tree, so this UniPert-inspired
representation lane should stop rather than scale model capacity.
F-159 followed the PC-SAR Meta-Kernel document's first stop-tree Gate and found
no free target-specific low-dimensional SAR headroom on the current main-v0 /
T-BASIS feature proxy.
F-160 is the current positive turn: BindingDB complete 2x2 rectangles show
strong observed-label double-difference magnitude after additive protein and
ligand effects cancel. Under the theory integration this is not yet latent
non-additivity because replicate/noise correction is missing. F-161 then showed
that exact ligand-pair transformation identity cannot transfer because
train/development share zero transformation IDs. F-162 replaced exact IDs with
a deployment-computable no-intercept quotient descriptor, but the hand-built
amino-acid composition x Morgan/physchem family failed G2 against zero,
ligand-only, wrong-protein and shuffled-protein controls. F-163 replaced only
the protein side with frozen ESM2 slot-region means. It beat zero but still
failed ligand-only, wrong-protein and shuffled-protein controls, so it
identified weak PLM-associated rectangle predictability rather than an admitted
cold target-specific quotient.

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
research/crossed_interaction/audit_bindingdb_rectangle_interaction.py
research/crossed_interaction/train_bindingdb_rectangle_lowrank.py
research/meta_fewshot/pcsar_oracle_gate.py
research/meta_fewshot/train_qpsmp_core.py
tests/test_train_chembl_assay_sardelta.py
tests/test_train_bindingdb_sardelta_cq_bridge.py
tests/test_train_sardelta_potential_cq_observable.py
tests/test_train_sardelta_edge_cq_observable.py
tests/test_train_bindingdb_sardelta_attribution.py
tests/test_audit_bindingdb_sardelta_symmetry.py
tests/test_train_bindingdb_pair_score_difference.py
tests/test_audit_bindingdb_rectangle_interaction.py
tests/test_train_bindingdb_rectangle_lowrank.py
tests/test_pcsar_oracle_gate.py
tests/test_train_qpsmp_core.py
report/source_affinity/
report/crossed_interaction/
report/meta_fewshot/qpsmp_core_smoke_separated_20260812/
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
research resumes, do not start full PC-SAR Meta-Kernel V2 on the current
main-v0/T-BASIS feature family because F-159 failed its source-only oracle
headroom Gate. F-162 and F-163 have now run two protein-side descriptor
replacements and G2 remains unadmitted, so G3a and G3b remain unauthorized. The
model-core contract has nevertheless been corrected to QPSMP: one retained
scalar potential must produce the zero-shot endpoint, centered support state,
delta/rectangle quotients and query certificate. The next falsifiable axis
should stay inside G2 and address protein specificity rather than ridge tuning
or endpoint integration. A train-only or externally supervised localizer over
frozen PLM slots is admissible only if it feeds the retained QPSMP scalar path
and preserves noise/replicate-aware or otherwise G2-compliant controls, a
deployment-computable chemical transformation descriptor, the no-intercept
quotient form, and component-bootstrap evaluation. The admission criterion must
include:

```text
interaction beats ligand-delta
interaction beats wrong-target
interaction beats shuffled-target
antisymmetry holds for pair deltas
rectangle/additive main effects are explicitly controlled if rectangles are used
observed-label magnitude is not overstated as latent non-additivity
train/development transformation representation has shared support
component bootstrap LCB95 > 0
support/query isolation preserved
quotients are derived from the retained scalar path, not an independent head
centered support state respects rank <= k-1
query certificate fields are emitted and calibrated before selective ranking
```

F-164 has now trained the corrected QPSMP endpoint scalar path with ligand-only
baseline features separated from crossed interaction features. It failed closed
against the matched level baseline, zero-support and foreign-support controls,
and it is worse than the previous V1 targeted-repair RMSE at every support size.
The only positive contrast is correct query against wrong query, which is
protein-identity sensitivity without useful endpoint calibration. Therefore the
next repair should not tune QPSMP ridge/task dimension on the same features; it
should return to the protein-specific G2 representation/localization problem.
