# K-LBP v2 R4 — one-seed train-only mechanism pilot on Metz (preregistration)

**Frozen:** 2026-07-27, **before any R4 statistic was computed.**
**Program:** `task.md` Part 9 (K-LBP v2). **Stage:** R4. **Gating:** mechanism pilot only —
`R4_COORDINATE_LOAD_BEARING_TRAIN_ONLY` is **hypothesis-generating** and authorizes no predictive
claim, no R5, and no architecture comparison. **Prerequisites (hard):** R3 estimator certified
(`R3_ESTIMATOR_CERTIFIED`); at least one coordinate survived R1; for any LLM-card arm, R2 eligible.
**Design rationale:** `reports/active/model_blueprint_reconstruction_2026-07-27.md` §5, §10; task.md
§9.3–§9.9.

---

## 1. The single question

> **On the one substrate where the interaction is identified (Metz train), does the rank-1
> scalar-gated correction along an audited mechanism coordinate measure a nonzero, control-surviving
> `γ̂` — i.e., does the coordinate carry target-specific interaction information beyond the
> shared-global basis?**

Model (task.md §9.3, coefficient scale):

```text
beta_t = w_bar + gamma * (a^T k_t) * c + eps_t,   eps_t ~ N(0, V_t)
||a|| = ||c|| = 1, gamma >= 0,   gamma = 0  ==>  shared-global (nested null)
```

This is a **mechanism estimand** (a parameter), not a ratio and not a containment mean (both retired,
task.md §2.8). Within-target ranking predictions of the fitted model are reported as a non-gating
secondary view only.

## 2. What is read

* ChEMBL-37 Metz dense pKi panel, **TRAIN cells only** (PARC M0 eligible set, 111 targets / ~100
  components), via the Stage-D/PARC-M0 machinery unchanged: exact projected residual `M_X^W y`
  (`research.panel_pdm.Substrate`), per-target ridge coefficients and sandwich covariances
  (`research.hqgbma_stage_d.target_coefficients`), frozen 5-fold component folds
  (`research.panel_power.component_folds`).
* Coordinates: every R1 survivor; every R2-eligible card arm. As of this preregistration the candidate
  list is: `det_proxy_card`, `klifs_pocket_composition`, `card_named` (pending R2), `card_deidentified`
  (pending R2). Arms that failed their prerequisite stage are simply absent — they are not re-run here.
* No development, confirmation, Davis, or sealed label. `sealed_test_consumed=false`.

## 3. Fitting protocol (frozen; the R3-certified estimator, variant E)

1. **Global fit:** the R3 alternating-GLS estimator (variant E, errors-in-variables correction) on all
   eligible train targets, per coordinate arm. Report `γ̂`, `a`, `c`, and the top-loading coordinate
   fields.
2. **Cross-fitted destruction evaluation:** for each of the 5 frozen component folds, fit
   `(w̄, γ, a, c)` on out-of-fold targets, then score each held target's **weighted fit improvement**
   `δ_t = ||β̂_t − w̄||²_{V_t⁻¹} − ||β̂_t − w̄ − γ(aᵀk_t)c||²_{V_t⁻¹}`, aggregated to per-component means.
   Held targets never contribute to their own fold's fit.
3. **Inference:** per-component `δ` vector → component bootstrap (10,000 draws, seed 1729) for
   contrasts; sign test and Wilcoxon signed-rank on `δ` as the primary mechanism inference
   (task.md §2.8); `γ̂` CI from component bootstrap over refits.
4. **Rank-2 sensitivity (non-gating):** refit with a rank-2 correction `(γ₁(a₁ᵀk_t)c₁ + γ₂(a₂ᵀk_t)c₂)`;
   report `γ̂₂` and the rank-1→2 δ increment. Diagnostic only (R3 S4 informs its interpretation).

## 4. Arms (frozen)

| arm | coordinate | role |
| --- | --- | --- |
| per surviving coordinate | as listed in §2 | the candidates |
| `wrong_target` | candidate coordinate under the exposure-matched derangement (the PARC M0 permutation, reused) | **C3** |
| `random_coordinate` | Gaussian at matched dimension and column marginals | **C4 (pipeline-free)** |
| `taxonomy_centroid` | KLIFS group one-hot (same fitting protocol, same rank-1 form) | **the TR-0 arm — must be beaten** |
| `field_ablated` | candidate coordinate with its top-loading field (from the global fit) removed | **C5 attribution** |
| `shared_global` | γ ≡ 0 | nested null (reference for δ; not a fitted arm) |

## 5. Frozen gates (unit of inference: homology component; all intervals paired component bootstraps)

| gate | contrast / statistic | requirement |
| --- | --- | --- |
| **G-R4-1** | `γ̂` (global fit, variant E) | > 0 with component-bootstrap LCB95 > 0 for ≥ 1 surviving coordinate |
| **G-R4-2** | `δ(coordinate) − δ(wrong_target)` | LCB95 > 0 (C3) |
| **G-R4-3** | `δ(coordinate) − δ(random_coordinate)` | LCB95 > 0 (C4) |
| **G-R4-4** | `δ(coordinate) − δ(taxonomy_centroid)` | LCB95 > 0 (TR-0 lesson) |
| **G-R4-5** | `δ(coordinate) − δ(field_ablated)` | mean > 0 (attribution; LCB reported, mean-gated) |
| primary inference | sign / signed-rank on per-component `δ` | reported with the gates |

## 6. Frozen verdict rule

```text
prerequisite missing                  -> stage not run (recorded)
G-R4-1..G-R4-4 all pass (>=1 coord)  -> R4_COORDINATE_LOAD_BEARING_TRAIN_ONLY
                                          (hypothesis-generating; unlocks nothing by itself;
                                           R5 remains blocked by task.md §2.7)
any of G-R4-1..G-R4-4 fails          -> R4_COORDINATE_NOT_LOAD_BEARING_STOP
                                          (no re-run with other seeds, widths, ranks, or thresholds)
```

Multiple surviving coordinates are evaluated independently; a head-to-head paired contrast
(`δ(coord_A) − δ(coord_B)`) is reported as the registered coordinate-selection statistic.

## 7. Declared expected outcome (stated before running)

The honest prior, from the ledger: the pocket/chemistry coordinate shows a **positive δ** (PARC M0's
G3 was directionally positive; C2 measured pocket−ESM +0.0556 in a ligand-warm design), and the
genuinely uncertain gate is **G-R4-4** — the TR-0 precedent (own-group identity unnecessary) predicts
the taxonomy centroid retains most of the δ. If G-R4-4 fails while G-R4-1..3 pass, the correct reading
is that the coordinate carries real within-kinome signal that is **not resolvable beyond coarse
taxonomy** — a substantive negative for the mechanism-coordinate programme, not a procedural stop.

## 8. Prohibited rescues

No threshold, seed, rank, width, step-count, or fold change after a result. No arm added after a
result. No development/confirmation/sealed label. No containment statistic anywhere. A failed gate may
not be re-read as underpowered without a pre-registered power computation from the stored per-component
values. R4 authorizes no predictive claim, no R5, no architecture comparison.

## 9. Artifacts

```text
research/klbp_r4_pilot.py                     runner, deterministic seed 1729
reports/active/klbp_r4.json                   machine-readable result, parses with allow_nan=False
reports/active/klbp_r4_decision.md            verdict + what was NOT shown
tests/test_klbp_r4.py                         fold isolation proof, δ identity (train vs held),
                                              derangement fixed points, V_t⁻¹ weighting correctness
```
