# PARC M0 — coordinate-misspecification audit (preregistration)

**Frozen:** 2026-07-27, **before any M0 statistic was computed.**
**Program:** `task.md` Part 8 (PARC). **Stage:** M0. **Gating:** non-gating for any predictive claim;
gating only for whether PARC M1 may be attempted.
**Design rationale:** `reports/active/biological_knowledge_integration_review_2026-07-27.md`.

---

## 1. The single question

Six structurally independent parameterisations of the interaction have failed on the *same* protein
coordinate (pooled/frozen ESM-2), and HQ-GBMA Stage D showed the ESM-conditioned Grassmann map is
**worse than a protein-free shared basis** (`true − global = −0.232`, `LCB −0.083`).

> **Is the protein uninformative for the target-specific interaction direction, or was the coordinate
> wrong?**

M0 answers this by re-running Stage D with **exactly one changed input** — the protein coordinate — and
nothing else. Estimator, folds, seed, rank, ridge, optimiser, step count, ambient space, error-corrected
containment objective and bootstrap are all carried over unchanged.

## 2. Substrate and what is read

* ChEMBL-37 Metz dense pKi panel, **TRAIN cells only**, via `research.panel_gate_pa.load_panel_train()`.
* No panel development, panel confirmation, Davis, ChEMBL development/confirmation or sealed label is
  read. `sealed_test_consumed=false`.
* KLIFS v3.2 local snapshot `dataset/public/klifs_2026_07_22/raw/kinase_information.json.gz` supplies the
  85-residue aligned pocket per UniProt accession. It contains **no affinity information**.
* `dataset/public/chembl_37/processed/dualcold/target_sequences.json` supplies full sequences for the
  pocket-set-shuffle control. UniProt-derived, no activity.

## 3. Eligibility (frozen rule, matched across all arms)

A panel train target is **eligible** iff it has (i) a frozen ESM-2 entry, (ii) a *Human* KLIFS record with
a pocket string of **exactly 85 residues**, and (iii) a full sequence of length ≥ 85.
**Every arm is scored on exactly this eligible set**, so all contrasts are matched. Excluded targets are
enumerated in the report. Exclusion is by protein annotation only and is label-blind.

## 4. Arms (all frozen protein coordinates, PCA-whitened to d = 32)

| arm | coordinate | role |
| --- | --- | --- |
| `shared_global` | none — top-`r` eigenvectors of the noise-corrected scatter of out-of-fold `beta_hat` | **primary null** (`gamma = 0`); the arm that won Stage D |
| `esm_pooled` | frozen ESM-2 650M pooled, 1280 → 32 | incumbent; reproduces Stage D `protein_true` |
| `parc_pocket` | **PARC** — one-hot of the 85 aligned KLIFS pocket residues (85 × 20 = 1700) → 32 | the candidate coordinate |
| `parc_random_positions` | one-hot of **85 uniformly random positions from the same protein's own sequence** (1700 → 32) | **A2**, the decisive pocket control |
| `parc_pocket_composition` | 20-d amino-acid composition of the pocket (aligned position destroyed) | **A3 substitute**; isolates *positional* information |
| `parc_wrong_target` | `parc_pocket` under an exposure-matched target derangement | **A4** |
| `esm_wrong_target` | `esm_pooled` under the same derangement | A4 for the incumbent |
| `random_features` | Gaussian at matched scale and dimension | **A5** |

`parc_pocket_composition` carries at most 20 degrees of freedom by construction and is therefore
**reported non-gating** — it is an information-content ablation, not a capacity-matched control.

**A3 (Foldseek/3Di structure-token shuffle) is NOT run in M0** and is not claimed. No 3Di tooling exists
in the `drug` env and `dataset/structure/alphafold` was removed in the 2026-07-25 cleanup. The KLIFS
alignment is itself structure-derived, so M0 tests the *positional/pocket* half of PARC only. The 3Di half
remains an unaudited E-gate dependency.

## 5. Carried-over constants (unchanged from Stage D; none re-tuned)

```text
rank r          = 6          (Stage D's nested-CV selection, carried over; NO re-selection on the new
                              arm, so PARC receives no selection advantage)
feature dim     = 32         (Stage D panel_target_features dimensions=32)
ridge           = 1.0
map steps       = 400        lr = 5e-2       hidden = 64      weight decay = 1e-5
folds           = frozen component_folds(DualCold.panel()), 5 folds
seed            = 1729
MIN_ADVANTAGE   = 0.02       (Stage D MIN_PROTEIN_ADVANTAGE, carried over unchanged)
bootstrap       = research.panel_stats.bootstrap_lcb, component as the unit of inference
```

## 6. Frozen gates

Unit of inference is the homology component. All intervals are paired component bootstraps.

| gate | contrast | requirement |
| --- | --- | --- |
| **G0** | estimator sensitivity (§7) | on synthetic betas generated **from** the PARC coordinate, `parc − global` mean ≥ 0.02 and LCB95 > 0 |
| **G1** | `parc_pocket − esm_pooled` | LCB95 > 0 — *is the coordinate better than the incumbent?* |
| **G2** | `parc_pocket − shared_global` | mean ≥ 0.02 **and** LCB95 > 0 — *does it beat the protein-free basis?* |
| **G3** | `parc_pocket − parc_random_positions` | LCB95 > 0 — *is it the pocket, or any 85 residues?* |
| **G4** | `parc_pocket − parc_wrong_target` | LCB95 > 0 |
| **G5** | `parc_pocket − random_features` | LCB95 > 0 |

Reported, non-gating: `esm_pooled − shared_global` (Stage D reproduction), `esm_pooled − esm_wrong_target`,
`parc_pocket − parc_pocket_composition`, per-arm containment summaries, orthonormality error, projection
audit, positive-signal fraction.

## 7. G0 — the estimator-sensitivity positive control

A null result is only interpretable if the estimator *could* have detected coordinate dependence. Before
any real contrast is read, synthetic coefficients are generated in which the subspace is a **known
function of the PARC coordinate**:

```text
B_t   = orth( A_0 + sum_k u_t[k] * A_k )[:, :r]        A_* fixed Gaussian, seed 1729
beta_t^syn = B_t c_t + xi_t,   c_t ~ N(0, I_r),   xi_t ~ N(0, sigma^2 I_m)
sigma^2 = mean_t trace(V_t) / m      (noise matched to the real per-target estimation noise)
V_t^syn = sigma^2 I_m
```

The identical CV then scores `parc_pocket` and `shared_global` on these synthetic coefficients. `esm_pooled`
is scored too as a specificity check (it should sit near `shared_global`, since the synthetic subspace is a
function of PARC and not of ESM).

**If G0 fails the whole run returns `PARC_M0_ESTIMATOR_INSENSITIVE_NO_DECISION` and no scientific claim is
made in either direction.**

## 8. Frozen verdict rule

```text
G0 fails                          -> PARC_M0_ESTIMATOR_INSENSITIVE_NO_DECISION
G0 and G1..G5 all pass            -> PARC_M0_COORDINATE_LOAD_BEARING          (M1 authorised)
G0 and G1 pass, G2 fails          -> PARC_M0_COORDINATE_IMPROVES_ESM_ONLY_STOP
G0 passes and G1 fails            -> PARC_M0_COORDINATE_NOT_LOAD_BEARING_STOP
otherwise                         -> PARC_M0_COORDINATE_NOT_LOAD_BEARING_STOP
```

`PARC_M0_COORDINATE_IMPROVES_ESM_ONLY_STOP` is registered here as a distinct outcome because it is the
**declared expected result** (§9) and it must not be reportable as either a pass or a plain failure. It
stops M1 but it is a substantive, quantitative finding: it would localise the residual failure to the
coefficient map rather than to the coordinate.

## 9. Declared expected outcome (stated before running)

Most likely: **G1 passes, G2 fails** — PARC beats pooled ESM (C2/ASPIRE-P0 measured `pocket − ESM =
+0.0556 [+0.0403, +0.0714]` in a different, ligand-warm design) but still does not beat the shared global
basis. G3 is the genuinely uncertain gate: if `parc_random_positions` matches `parc_pocket`, then any
apparent PARC gain is a *sequence-composition* effect and not a pocket effect, and the correct reading is
that the pocket restriction is inert.

## 10. Prohibited rescues

No threshold may be changed after a result. No rank, width, step count, dimension, ridge, seed or fold may
be re-tuned after a result. No additional arm may be added after a result. No development, confirmation,
Davis or sealed label may be read. A G0 failure may not be repaired by weakening the synthetic signal.
M0 authorises no predictive claim, no F1–F4 work and no Mamba comparison under any outcome.

## 11. Artifacts

```text
research/parc_m0.py                       runner, deterministic seed 1729
reports/active/parc_m0.json               machine-readable result, parses with allow_nan=False
reports/active/parc_m0_decision.md        verdict + what was NOT shown
tests/test_parc_m0.py                     batched-loss equivalence, eligibility, control construction
```
