# PREREG — R2 identifiability-regime repair

> **POST-EXECUTION STATUS — INVALIDATED BY ITS OWN REGISTERED E1 FALSIFIER.**
> This file is preserved as the historical preregistration and is not silently
> rewritten.  Execution produced `gauge_ratio=1.0306` on meta-val and `1.0987`
> on meta-test, both above the §6 falsifier threshold `>1`.  The registered
> H0-regime therefore failed.  E0 supports a narrower descriptive conclusion:
> the average section gain is predominantly target-level calibration.  It does
> not imply that representation or data are theoretically unable to matter at
> `d=2,k=5`.  RFMS training is not authorized; its original wrong/wrong
> guarantee also fails without a nonzero quotient-exposure certificate.
>
> Implementation deviations discovered after execution: the original E0/E1
> code did not compute the required `protein_group_40` cluster macro; the
> original E2 loaded outcome-bearing cell rows and used a sparse additive fit as
> a capacity claim; the original E3 did not support matrix-shaped panels.  The
> repaired v2 scripts and final resolution are documented in
> `report/meta_fewshot/R2_MULTI_AGENT_RESOLUTION.md`.  These repairs are
> descriptive and open no Gate.

Registered: 2026-08-11. Registered before any R2 script was executed.
Supersedes nothing; amends the "next admissible action" clause of
`BIOLOGICAL_GAUGE_AUDIT_REPORT.md`.

## 0. Why this preregistration exists

The A0/A1 audit closed with `data acquisition, not architecture search`. That
conclusion assumed the v0 estimator *could* have used partner information and
did not. This preregistration records the opposite hypothesis, derived from the
estimator itself rather than from the data:

> **H0-regime.** At the deployed operating point (`d=2`, `k=5`, ridge `1.0`) the
> v0 predictor is, up to an `O(ridge)` term, invariant under episode-wise
> `G ∈ GL(d)` acting on the section coordinates. The partner can therefore
> contribute only through the support-weight vector `u ∈ R^k`, and the support
> fit re-identifies all `d` coefficients on its own. Wrong/wrong recovery is a
> property of the estimator, not evidence about biology.

If H0-regime holds, **no representation repair and no additional data can rescue
partner specificity at this operating point**, and A1's negative result is
uninformative about whether biological information exists.

## 1. Declared roles and label discipline

| Experiment | Opens labels? | Split | Declared role |
|---|---|---|---|
| E0 intercept null | yes (already-consumed) | `meta_val` primary, `meta_test` descriptive | **Descriptive diagnosis of a published result.** Opens no Gate. |
| E1 sufficient statistic | **no** | `meta_val` | Label-free estimator property. |
| E2 T-BASIS decomposition | **no** | whole corpus | Label-free representation property. |
| E3 crossed census | design only | local panels | Design census; no model fit, no held-out label read. |
| E4 RFMS development | yes | source + `meta_val` only | Development. **Not confirmation.** |
| E5 confirmation | yes | **fresh supply, see §5** | Confirmatory. Opened once. |

E1–E3 open no affinity label, so they cannot be an outcome-dependent expansion
of the hypothesis family. This is the reason A2's closure does not bind them.

## 2. E0 — the missing null (decisive, cheapest)

v0's battery contains no per-target intercept arm. Register the comparison:

```
intercept null:  yhat_q = pop(L_q) + mean_i( y_i - pop(L_i) )
full correct:    yhat_q = pop(L_q) + m_q · ridge(M, rho)
```

Preregistered decision rule, target-macro **and** cluster-macro:

- `MSE(intercept) - MSE(full) <= 0` → **`META_SECTION_EFFECT_IS_CALIBRATION`**.
  Every biological contrast in `MAIN_V0_RESULT.json` is then measuring
  calibration noise and must be withdrawn as evidence about biology.
- `0 < MSE(intercept) - MSE(full)` but `< 0.10 × (MSE(zero) - MSE(full))` →
  **`META_SECTION_EFFECT_PREDOMINANTLY_CALIBRATION`**.
- otherwise → `META_SECTION_EFFECT_NOT_EXPLAINED_BY_CALIBRATION`.

Prior expectation from the published table: `full_permuted` (2.047) already
retains 98.1% of the gain over `population_d0` (8.711 → 1.916), so the second
or first branch is expected. Permutation preserves the support label multiset,
so this is not a new hypothesis — it is a reading of the registered v0 result.

## 3. E1 — the sufficient statistic (label-free)

For each episode compute `u = m_q (MᵀM + λI)^{-1} Mᵀ ∈ R^k` under correct and
deranged partner features. Report, target-macro:

- `calibration_share = ‖Π_1 u‖ / (‖Π_1 u‖ + ‖Π_1⊥ u‖)`;
- `gauge_ratio = E‖u_correct − u_wrong‖ / E‖u_correct − 1/k‖`.

Decision rules: `calibration_share ≥ 0.9` →
`SUPPORT_WEIGHTS_ARE_ESSENTIALLY_UNIFORM`; `gauge_ratio ≤ 0.5` →
`PROTEIN_CHANNEL_NEAR_GAUGE_EQUIVALENT`. Either outcome corroborates H0-regime
without opening a label.

## 4. E2 / E3 — where the information went, and what is already on disk

E2 decomposes the frozen 288D T-BASIS into protein main, ligand main and
interaction variance shares, measures fixed-ligand partner dispersion, and
measures whether the protein main effect separates within-cluster homologs.
Registered thresholds: `interaction_share < 0.10` →
`TBASIS_IS_EFFECTIVELY_ADDITIVE`; `partner_dispersion_fraction < 0.10` →
`TBASIS_IS_LIGAND_DOMINATED`; `within/between α-distance > 0.7` →
`TBASIS_CANNOT_RESOLVE_HOMOLOGS`.

E3 counts, per local panel, rows, **interaction degrees of freedom**
`rows − (n_P + n_L − 1)` per connected block, and dependency components after
scaffold/cluster/assay closure. Registered rule: a panel qualifies as a
**training** supply at `interaction_df ≥ 5000` and as a **confirmation** supply
only at `dependency_components ≥ 30` with `largest_component_share ≤ 0.5`.

## 5. E4 / E5 — RFMS and its fresh confirmation supply

Mechanism, gates and required controls are specified in
`r2_reserved_fiber_section.py::describe()`. Additional registered constraints:

- `d ≤ 5`, `1 ≤ d_support < d`, `d_c = d − d_support ≥ 1`;
- the ligand frame `ψ` takes **no protein argument** (checked by construction);
- `c0(P)` writes only into the reserved block (checked by construction);
- report `c0_between_target_variance`; a collapse to a constant is a **FAIL**,
  because under collapse the wrong/wrong Gate would pass vacuously.

**Fresh confirmation supply.** `meta_test` (50 targets, 1,934 cells) and
`meta_val` (37 targets, nine clusters) are consumed. E5 must be run on a supply
disjoint from both, drawn from `BLK-BDB-PANELS` (85 assay-matched panels, 129
targets, 70 mmseqs40 clusters) and `PDSP-CORE` (non-kinase), with complete
mmseqs40-cluster and scaffold-component hold-out. The supply is frozen and
hashed before any E5 label is read. Kinase-only confirmation is not admissible
because requirement 5 of the research question is cross-dataset generality.

## 6. What would falsify the R2 programme

1. E0 shows the intercept null is clearly worse than `full_correct` at both
   target- and cluster-macro → H0-regime is wrong, the v0 section is doing
   ligand-specific work, and the diagnosis must be rewritten.
2. E1 shows `gauge_ratio > 1` → the partner channel is not gauge-like and the
   wrong/wrong recovery needs a different explanation.
3. E2 shows `interaction_share > 0.3` and homolog separation → T-BASIS is not
   the bottleneck and A1's negative is about probe power, not representation.
4. E3 shows `interaction_df < 5000` across all local panels → the crossed
   estimand is genuinely unavailable locally and `EXISTING_DATA_REQUIRES_NEW_
   ESTIMAND_CONSTRUCTION` stands.
5. RFMS with `d_c ≥ 1` still shows wrong/wrong recovery → implementation defect
   or `c0` collapse; not evidence for or against the mechanism until resolved.

## 7. Stopping rule

R2 stops, and no production migration is proposed, unless RFMS passes **all**
of: meta effect, support specificity, `correct < wrong/wrong`, a within-target
discrimination improvement, and cluster/component-macro independence, on the
E5 fresh supply. Loosening any of these to rescue RFMS is prohibited.
