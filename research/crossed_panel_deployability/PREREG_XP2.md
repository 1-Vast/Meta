# Preregistration XP2 — Deployability Of The Crossed-Panel Interaction Section

Registered: 2026-08-08, before any XP2 arm, model, threshold or metric was scored.
Supersedes nothing. `PREREG_XP1.md` and every XP1 verdict stand unchanged.

## 0. Why this registration exists

XP1 established that a crossed kinase panel contains a large, reproducible,
low-rank protein-by-ligand interaction `gamma_ij`, and that a support-identified
section recovers part of it under strict kinase-group closure while no zero-shot
protein representation does. XP1 is **not production-admissible** because:

1. its strongest arm used `k = 16` support, against a frozen requirement of `k <= 5`;
2. ligands were reused across train and test, so the ligand loading `u(L)` was a
   lookup table and never had to be computed from chemistry;
3. protein closure was enforced, ligand-scaffold closure was not;
4. the latent factors carry a rotation/gauge ambiguity, so no coordinate can be
   given a biological name;
5. it produced a point predictor with no section radius, coverage certificate,
   realizability flag or abstention;
6. it measured interaction reconstruction, not calibration of the emitted
   probability law `z -> F(z) -> B(z)F(z) -> K(beta)`;
7. the PDSP replication fell below XP1's own non-negligibility floor.

XP2 tests whether the section survives all of these simultaneously.

## 1. Hypothesis under test

> A low-dimensional target-ligand interaction section can be identified from at
> most five support measurements, transferred to unseen protein groups **and**
> unseen ligand scaffolds, and constructed from deployable protein/ligand
> features without ID or lookup shortcuts.

Formally, whether

```
Gamma(P, L) ~= < u(L), v(P, S) >
```

can satisfy: rank `d <= 5`; support `k <= 5`; `u(L)` computed from ligand
structure only; `v(P,S)` from the permitted protein representation and support
only; validity for unseen protein groups and unseen ligand scaffolds; query
direction covered by the support-identified section; correct support beating
zero, foreign and permuted support; and positive value beyond ligand-only
prediction.

## 2. Frozen data releases

| ID | file | SHA-256 | license / provenance |
|---|---|---|---|
| `P-METZ-SRC` | `metz.xls` (Metz et al. 2011 Nat Chem Biol Table S1, the journal supplement itself) | `81731c4004823bd45fa3898e25d6491d799dfd0e0486fcc8c9c821f9419dd591` | publisher supplementary data; mirrored at `polinavino/kinase-selectivity-definitions@8ab79cae` |
| `P-KLAEGER-SRC` | `aan4368_Table_S2.xlsx` (Klaeger et al. 2017 Science Table S2) | `d28b91e62e78e5e011b60da27672875621fef5cdabbea793ac9cce4b98db2c32` | publisher supplementary data; same mirror |
| `P-KLAEGER` | `klaeger_matrix.csv` | `cdf66c7d4e7c1e3a35aeb6995abbfdaf15be80f3e07715524b2bb4449d871010` | derived matrix, verified against `P-KLAEGER-SRC` in XP2-A |
| `P-PDSP` | `KiDatabase.csv` | `45c9a18ac30f1fad350d1dde186bc1f226c5a75d474ca50f50713852a5637ac6` | NIMH PDSP Ki database, free public resource |
| `A-KLIFS` | `klifs_kinase_information_human.json` | recorded in `dataset/raw/crossed_panels/protein_annotation/acquisition_manifest.json` | KLIFS, academic open access |

**XP2 promotes `metz.xls` from provenance-only to primary source.** The derived
`metz_matrix.csv` is demoted to a cross-check. Reason: the supplement carries
`Canonical_Smiles`, per-cell censoring notation, and 3,858 compounds, all of which
the derived matrix discards.

Prohibited and not accessed: DAVIS, recipient labels, ChEMBL37 affinity values,
PKIS2, Anastassiadis. Label-read counters are asserted in every artifact.

## 3. Endpoint semantics and censoring (frozen)

The supplement encodes three cell states, verified in XP2-A:

| state | encoding | count over 3,858 x 172 | XP2 treatment |
|---|---|---|---|
| measured | numeric `pKi` | 103,118 | **analysis set** |
| left-censored | string `"< x"`, 50 distinct thresholds `4.0 … 6.2` | 154,175 | excluded from fitting and scoring; thresholds retained for the XP2 censoring control |
| not tested | blank | 405,482 | excluded |

XP1's single-floor description of this panel is corrected here. Excluding
censored cells is conditioning on the outcome; the destructive control in §9
uses the **actual per-cell thresholds** rather than one nominal floor.

## 4. Panel construction (frozen before scoring)

`BLK-METZ-XP2`, built from `P-METZ-SRC` by iterating to a fixed point:

1. keep compounds with non-null `Canonical_Smiles` that RDKit parses;
2. keep compounds with `>= 10` measured kinase cells;
3. keep kinases with `>= 50` measured compound cells;
4. keep kinases mappable to KLIFS (group + 85-residue pocket), using the frozen
   `SYMBOL_ALIASES` table from XP1.

Expected scale from the label-blind census: ~927 compounds x ~150 kinases,
~33,000 measured cells, density ~0.24, 519 Bemis-Murcko scaffolds. The realised
index and its SHA-256 are recorded in the panel manifest.

## 5. Closures (frozen)

- **Protein closure:** KLIFS **group** (the strict kinase split). Family and
  pocket-identity closures are reported as secondary only.
- **Ligand closure:** Bemis-Murcko scaffold, then single-linkage merge of
  scaffolds whose ECFP4 Tanimoto similarity is `>= 0.5`, giving **scaffold
  components**. Exact-compound exclusion is implied and additionally asserted.
- **Double held-out (XP2-D):** for fold `f`, hold out protein-group set `G_f`
  **and** scaffold-component set `S_f`. Test cells are `G_f x S_f`. Training
  cells are `(not G_f) x (not S_f)`. Support for a test protein is drawn from
  that protein's cells on **training scaffolds only**. Nothing from `S_f` enters
  training, feature fitting, hyperparameter selection or support.

## 6. Model family (frozen; capacity-controlled, smallest first)

Stage 1 (training block only): `mu + alpha_i + beta_j` by least squares on
measured cells.
Stage 2: rank-`d` factorisation of the training residual, `d <= 5`.
Stage 3: ligand landing `chi: chemistry -> u`, fitted on training ligands only.
Stage 4: support solve for `(b_j, v_j)` on `k <= 5` support cells with ridge.

Ligand feature arms, in the order they may be introduced:

| arm | features | capacity |
|---|---|---|
| `L-DESC` | 10 deterministic RDKit descriptors (MW, cLogP, TPSA, HBD, HBA, rotatable bonds, rings, aromatic rings, FractionCSP3, heavy atoms) | linear ridge |
| `L-ECFP` | Morgan r=2, 1024 bits, frozen | linear ridge (dual form) |
| `L-CHEMBERTA` | frozen `DeepChem/ChemBERTa-77M-MLM` mean-pooled embedding | linear ridge |
| `L-KRR` | ECFP4 Tanimoto kernel | kernel ridge — **only** if a linear arm already shows non-zero loading transfer |

Controls: `L-RANDOM` (random Gaussian ligand features), `L-MEAN` (mean training
loading for every ligand).

Protein arms for XP2-E: support-only; ESM-2 t30 full sequence; aligned KLIFS
85-mer pocket; protein+support; protein+ligand+support.

All preprocessing (standardisation, feature selection, PCA if any) and all
hyperparameter selection use **training cells only**, inside nested
leave-one-closure-component-out CV.

## 7. Registered arms and controls (XP2-C/D)

| arm | definition |
|---|---|
| `P0` | population `mu` |
| `LIG` | `mu + alpha_hat(L)` with `alpha_hat` **predicted from chemistry** (unseen ligands have no fitted `alpha`) |
| `ADD` | `LIG + b_j(S)` — additive plus support intercept |
| `SEC` | `ADD + <uhat(L), v_j(S)>` — the section, `uhat` from chemistry |
| `SEC-U*` | as `SEC` but support ligands use their fitted `u` (realistic reference-compound variant) |
| `ZERO` | `ADD` with `v_j = 0` (zero adaptation) |
| `PERM` | support labels permuted within the protein |
| `FOREIGN` | `v_j` solved on a different held-out protein's support |
| `RANDCORR` | ligand-matched random correction of matched norm |
| `ORACLE` | `v_j` solved on all of the protein's measured cells (ceiling) |

## 8. Metrics, units of independence, and intervals

Primary metric: **MSE on double-held-out measured cells**, in log units, and
`R2_gamma = 1 - SSE(arm)/SSE(ADD)`.

Registered contrasts:

- `Delta_deploy = MSE(LIG) - MSE(SEC)` — value over ligand-only chemistry;
- `Delta_interaction = MSE(ADD) - MSE(SEC)`;
- `Delta_specific_zero = MSE(ZERO) - MSE(SEC)`;
- `Delta_specific_foreign = MSE(FOREIGN) - MSE(SEC)`;
- `Delta_specific_perm = MSE(PERM) - MSE(SEC)`.

Independence units are **two-way**: protein closure components and scaffold
components. Intervals are the **wider** of (i) a cluster bootstrap over protein
components and (ii) a cluster bootstrap over scaffold components, 2,000
resamples each. Reporting only the narrower interval is prohibited.

Support size is evaluated at exactly `k = 1,2,3,4,5`. `k = 5` is primary. Rank
`d = 3` is primary; `d in {1,2,3,5}` reported. Seeds `{0,1,2,3,4}`.

## 9. Destructive controls

1. **Censoring control.** A synthetic panel additive by construction with the
   real `mu, alpha, beta`, the real residual sd, and the **real per-cell
   censoring thresholds** applied cell-by-cell. Falsification: if this
   reproduces a materially positive `Delta_specific_foreign`, XP2's result is a
   selection artefact and is withdrawn.
2. **Ligand-ID ablation.** Replace `uhat(L)` by a permutation of ligand
   identities. Must destroy the effect.
3. **Scaffold-leak probe.** Assert zero intersection between training and test
   scaffold components, and zero shared exact compounds.
4. **Random-feature and mean-loading controls** as in §6.

## 10. Gate (frozen)

`DEPLOYABLE_SECTION_STATISTIC_IDENTIFIED` requires **all** of:

1. XP1 evidence reproduced from immutable artifacts (XP2-A);
2. ligand loading transferable to unseen scaffolds: held-out loading `R2 > 0`
   with CI lower bound `> 0`, and beating `L-RANDOM` and `L-MEAN`;
3. at `k <= 5`: `R2_gamma >= 0.05` with CI lower bound `> 0.02`;
4. simultaneous protein-group and ligand-scaffold generalisation (XP2-D), same
   thresholds, on a test set opened once;
5. `Delta_specific_zero`, `Delta_specific_foreign`, `Delta_specific_perm` all
   with CI lower bound `> 0`;
6. `Delta_deploy` CI lower bound `> 0`;
7. external replication (XP2-F) reproduces the sign with CI lower bound `> 0`;
8. the statistic passes the frozen-theory interface audit (XP2-G).

Any failure yields the corresponding bounded verdict from the allowed list and
**stops** promotion. Thresholds may not be changed after any test set is opened.
A numerically positive but negligible effect (`R2_gamma < 0.05`) is a failure.

## 11. Allowed terminal verdicts

```text
XP1_EVIDENCE_NOT_REPRODUCIBLE
CROSSED_INTERACTION_REPRODUCED
LIGAND_SIDE_DEPLOYMENT_REPRESENTATION_FAILED
K_LE_5_SECTION_NOT_IDENTIFIED
PANEL_LOCAL_LOW_RANK_META_LEARNING
BIOLOGICAL_LANDING_NOT_IDENTIFIED
DOUBLE_HELD_OUT_SECTION_IDENTIFIED
EXTERNAL_REPLICATION_FAILED
DEPLOYABLE_SECTION_STATISTIC_IDENTIFIED
```

No XP2 outcome authorises affinity-energetics claims, biological `z` admission,
DAVIS evaluation, production integration, modification of the frozen theory,
CSMO or Band, or an end-to-end DTA claim.

## 12. External replication (XP2-F), frozen in advance

Primary external panel: `P-KLAEGER` (Klaeger 2017 kinobeads, `pKd_app`),
independent laboratory, independent measurement technology, disjoint compound
set from Metz. Structures resolved from drug names against an open chemical
registry; resolution is label-blind and audited.

Two separated conclusions, never merged:

1. **direction transfer** — every source parameter frozen, no external fitting;
2. **basis transfer** — only a preregistered scalar affine calibration
   (`a + b * prediction`) fitted on external **training** components.

Secondary external panel: `P-PDSP` (different protein class), reported but not
used for any decision.
