# K-LBP v2 R1 — external-information audit, affinity-free (preregistration)

**Frozen:** 2026-07-27, **before any R1 statistic was computed.**
**Implementation-audit amendment A1 (2026-07-28, before the first valid R1 run):** an initial runner
execution materialized the shared label-bearing panel frame before dropping `y`; its JSON/NPZ are
invalidated and are not a scientific result. Before rerunning, the implementation was corrected to
project only non-affinity columns and apply the TRAIN filter in the parquet reader. The following
clarifications resolve discrepancies found before any valid statistic was accepted: (i) explicit
missing indicators are appended to the coordinate; the zero placeholder in a missing value slot is
interpretable only jointly with that indicator and is not a standalone imputation; (ii) the T1
strike-loss denominator contains every non-constant expanded coordinate column, including the planted
taxonomy tripwire, while constant value/indicator columns are excluded; (iii) ridge alpha is selected
separately for every output column inside each outer training fold; (iv) the KLIFS amino-acid
composition K2 remains a reported R1 diagnostic but cannot authorize R4 because it is not Part 9's
dependency-gated structure-token/geometry coordinate; (v) the opening sentence about bootstrap
intervals applies if an interval is reported, but R1's frozen gates are point thresholds and report no
inferential interval. Thresholds, seed, and the K1 decision rule are unchanged. The runner now also
writes the registered decision file and coordinate-distance diagnostic. This amendment precedes the
first valid R1 result and may not be revisited after it.
**Program:** `task.md` Part 9 (K-LBP v2). **Stage:** R1. **Gating:** for coordinate survival only; a
surviving coordinate earns the right to enter R4. R1 authorizes **no** affinity claim, no training, and
no predictive statement of any kind.
**Design rationale:** `reports/active/model_blueprint_reconstruction_2026-07-27.md` §4, §10.

---

## 1. The single question

> **Does any candidate mechanism coordinate carry target-intrinsic information that is NOT already
> carried by (i) KLIFS family/group taxonomy, (ii) pooled ESM-2, or (iii) target popularity — without
> touching any affinity label?**

The precedent to fear is `TR0_PREMISE_FAIL_STOP`: a real coarse-kinase effect (+0.0400) proved not
resolvable to the KLIFS group (own-group-cold +0.0091 [−0.0076,+0.0256]). A coordinate largely
determined by taxonomy, ESM, or study depth would reproduce that null in R4 while passing every other
control. R1 exists to kill such coordinates **before** any affinity model is fitted.

## 2. What is read (no affinity labels anywhere)

* `dataset/public/klifs_2026_07_22/raw/kinase_information.json.gz` — KLIFS v3.2 local snapshot:
  85-residue aligned pocket strings, family/group assignment, subpocket and conformation annotations.
  **No affinity information.**
* `dataset/public/chembl_37/processed/dualcold/target_sequences.json` — full target sequences.
  UniProt-derived, no activity.
* Frozen ESM-2 cache `ESM_CACHE` (as in `research/panel_gate_pa.py`) — pooled 1280-d vectors.
  No affinity supervision.
* Target eligibility: exactly the PARC M0 eligible set (Human KLIFS 85-residue pocket + sequence ≥ 85),
  recomputed label-blind from annotations; **111/112 Metz panel train targets** expected. Excluded
  targets are enumerated.
* Train-panel **cell counts and ligand counts per target** (popularity proxies only — counts, never
  values). No pKi/pKd value is read, computed, or displayed anywhere in R1.
* **Nothing is downloaded; no API is called; no LLM card exists in R1.** The LLM-compiled card and the
  sequence-only de-identified card are R2 outputs; the R1 analyses below are pre-registered to be
  re-run on those cards as "R1-part-2" once R2 delivers them, with identical thresholds.

## 3. Coordinates built (deterministic, label-free)

| id | name | content | dimension |
| --- | --- | --- | --- |
| K1 | `det_proxy_card` | Deterministic proxy Mechanism Card: schema-subset fields computed in code from KLIFS + sequence only — `site_polarity_class` (charged-residue fraction in pocket, ord 0–2), `hbond_donor_availability` / `hbond_acceptor_availability` (donor/acceptor residue fractions, ord 0–2), `buried_hydrophobic_subpocket` (hydrophobic fraction ≥ frozen cut), `charged_residue_in_site` (dominant charge class, one-hot 3), `metal_coordination_in_site` (H/C/D/E at frozen KLIFS positions), `covalent_targetable_nucleophile` (Cys in pocket), `enzyme_or_receptor_class` (kinase→fixed class; included as a negative-control field expected to be struck by §4 T1), `documented_conformational_states` (KLIFS DFG/αC annotation coverage, ord 0–3), `induced_fit_reported` / `cryptic_pocket_reported` (KLIFS conformation-annotation presence), plus per-field missing indicators. Fields with no deterministic source (`site_sequence_conservation`, `cofactor_requirement`, `catalytic_mechanism_class`, `local_disorder_near_site`, `allosteric_site_documented`) **abstain** with missing indicators. | ≤ 24 + missingness |
| K2 | `klifs_pocket_composition` | 20-d amino-acid composition of the 85-residue pocket (the PARC M0 `parc_pocket_composition` arm) plus 4 chemistry aggregates (charged/polar/hydrophobic/aromatic fractions) | 24 |
| K3 | `esm_pooled` | pooled ESM-2 1280 → 32 PCA-whitened (the incumbent coordinate; carried as the redundancy reference, not a candidate) | 32 |

Coordinates K1/K2 are standardized (column z-scores over the eligible targets) after construction.
All construction code is deterministic; the construction seed is 1729 and the constructed matrices are
SHA-256 hashed into the R1 JSON.

## 4. Frozen audits and thresholds (all affinity-free)

Unit of inference: the eligible target (n ≈ 111). All intervals are target bootstraps (10,000 draws,
seed 1729). Three stop families, evaluated per coordinate K1, K2 (and later, per R2 card):

**T1 — taxonomy (C6).** For each coordinate column `x`, fit ridge regression from KLIFS **group**
one-hot (and separately KLIFS **family** one-hot) with 5-fold target-level CV; record CV R².
*Threshold:* any column with group-CV R² ≥ **0.50** (or family-CV R² ≥ 0.50) is **struck from the
coordinate** (task.md §9.4 hard schema rule 2). A coordinate losing ≥ 50% of its non-missing columns
to striking is `R1_COORDINATE_IS_TAXONOMY_STOP`.

**T2 — ESM redundancy.** For each surviving column, fit ridge from pooled ESM-2 (1280-d) with 5-fold
target-level CV (nested alpha selection inside training folds only). *Threshold:* median column CV R²
≥ **0.70** → `R1_COORDINATE_REDUNDANT_WITH_ESM_STOP`. Reported per-column; the coordinate-level
decision reads the median.

**T3 — popularity (C7 proxy).** Spearman between (i) each coordinate column and (ii) two popularity
proxies: train cell count per target, distinct ligand count per target. *Threshold:* median |ρ| ≥
**0.50** → `R1_COORDINATE_IS_POPULARITY_STOP`. Also recorded: coordinate fill rate vs popularity
(missingness itself must not be popularity-driven; |ρ| ≥ 0.50 is reported as a warning that forces
abstention-as-missing rather than imputation downstream — which is already the schema rule).

**R1-part-2 (after R2, identical thresholds):** T1–T3 re-run on the LLM named card, the
sequence-only de-identified card, and the deterministic-proxy card as the mandatory non-LLM baseline
arm (task.md §9.4 hard rule 7).

## 5. Reported, non-gating

Pairwise coordinate distance vs pairwise pocket-sequence identity (does the coordinate resolve beyond
global homology?); column-level R² tables; struck-column lists; eligibility audit; missingness map;
K1-vs-K2 agreement (do the proxy card and pocket chemistry agree where both are defined?).

## 6. Frozen verdict rule (per coordinate)

```text
T1 strike-loss >= 50% of columns   -> R1_COORDINATE_IS_TAXONOMY_STOP
T2 median column R² >= 0.70        -> R1_COORDINATE_REDUNDANT_WITH_ESM_STOP
T3 median |rho| >= 0.50            -> R1_COORDINATE_IS_POPULARITY_STOP
none of the above                  -> R1_COORDINATE_SURVIVES (eligible for R4, subject to R2 for cards)
```

Survival is **not** evidence that the coordinate predicts affinity (R1 reads no label). It is evidence
only that the coordinate is not explainable as taxonomy, ESM, or popularity.

## 7. Declared expected outcome (stated before running)

Most likely: K2 (pocket composition) survives T1 partially — chemistry aggregates like charged fraction
are expected to be partially taxonomy-predictable within kinases — and survives T2/T3; K1 loses its
class-label column to T1 striking (by design, the negative control) but retains site-chemistry columns.
The genuinely uncertain audit is T1 for the site-chemistry columns: if KLIFS group predicts them at
R² ≥ 0.5, the within-kinome mechanism-coordinate programme collapses to the TR-0 null.

## 8. Prohibited rescues

No threshold change after a result. No column re-weighting in place of striking. No imputation of
missing fields. No affinity label may be read to "check" an R1 decision. No new coordinate may be
added after results. R1 authorizes no training and no predictive claim.

## 9. Artifacts

```text
research/klbp_r1_information.py               runner, deterministic seed 1729
reports/active/klbp_r1.json                   machine-readable result, parses with allow_nan=False
reports/active/klbp_r1_decision.md            verdict per coordinate + what was NOT shown
tests/test_klbp_r1.py                         eligibility, coordinate construction, strike logic, no-label guard
```
