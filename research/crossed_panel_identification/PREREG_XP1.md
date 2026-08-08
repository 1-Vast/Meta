# Preregistration XP1 — Crossed-Panel Identification Of The Protein-By-Ligand Affinity Interaction

Registered: 2026-08-08
Registered before any contrast, model fit, arm score or Gate metric was computed.
Only label-blind structural census, endpoint semantics and censoring statistics
were computed before this document (repository precedent: `sigma_assay`
estimated before scoring).

## 0. Why this registration exists

The repository's affinity chain stopped at
`AFFINITY_INCREMENT_NOT_IDENTIFIED` / `BIOLOGICAL_Z_NOT_ADMITTED`. Two
structural facts from the existing evidence determine this design:

1. `E-AFF-R0` proved the historical readout (within-task concordance) is
   **exactly** invariant to per-task affinity location, so the historical
   negatives are scoped to within-task ranking only.
2. `E-AFF-X0` / `X0-FEAS` proved the governed ChEMBL37 corpus cannot supply
   independent crossed units, because document-keyed panels and homology-document
   closure are produced by opposite kinds of study. Acquiring a *panel* corpus —
   one lab, one assay, one complete ligand x protein rectangle — is the only
   construction that produces crossing.

XP1 therefore moves to externally acquired **complete crossed panels** and to a
**location-sensitive** readout (RMSE / variance explained in log units), not a
rank readout.

## 1. Hypotheses under test

- **H1 (representation failure).** A crossed panel contains identifiable,
  reproducible protein-by-ligand interaction signal `gamma_ij`, but no available
  biological representation of the protein predicts it for unseen proteins.
- **H2 (data identifiability failure).** `gamma_ij` is not reproducible above
  measurement noise, or is not separable from ligand and protein main effects, in
  any accessible panel. Then no architecture can recover it.
- **H3 (objective/parameterization failure).** `gamma_ij` is reproducible *and*
  predictable from available biological features, but only under a
  parameterization the current MetaSieve interface does not use (e.g. the signal
  is a low-rank section identified from the support set rather than a
  protein-feature-conditioned surface, or it is destroyed by the readout).

These are mutually exclusive at the level of the decision, and the arm structure
in section 5 is designed so that each hypothesis predicts a different arm
ordering (section 8).

## 2. Estimand

For panel `P`, ligand `i`, protein `j`, in log-affinity units,

```
y_ij = mu + alpha_i + beta_j + gamma_ij + eps_ij
sum_i alpha_i = 0, sum_j beta_j = 0, row and column sums of gamma = 0
```

The estimand is **`gamma_ij`** — the double-centered interaction. Not `alpha_i`
(ligand potency / promiscuity), not `beta_j` (protein main effect /
druggability), not `y_ij`.

Two derived estimands are also registered:

- **XP1-A (existence/reproducibility).** Is `gamma` reproducible above noise, and
  is its protein-side geometry a stable property of the protein?
- **XP1-B (transfer).** Does knowledge of an unseen protein — from its features
  (zero-shot) or from `k` labelled support ligands (few-shot) — predict
  `gamma_ij` on held-out cells of that protein?

## 3. Data releases (frozen)

| ID | Source | Endpoint | Shape | SHA-256 |
|---|---|---|---|---|
| `P-METZ` | Metz et al. 2011 Nat Chem Biol kinase panel, mirror `polinavino/kinase-selectivity-definitions@8ab79cae` file `metz_matrix.csv` | `pKi` (single lab, single assay family) | 704 compounds x 172 kinases, complete | `abe1e3c580478775a352ec5ee78ca565d4c863f0e3e642fdb21d956d8f9d4375` |
| `P-METZ-SRC` | same, `metz.xls` original supplement | provenance only | — | `81731c4004823bd45fa3898e25d6491d799dfd0e0486fcc8c9c821f9419dd591` |
| `P-KLAEGER` | Klaeger et al. 2017 Science kinobeads, `klaeger_matrix.csv` | `pKd_app` | 222 drugs x 343 kinases, complete | `cdf66c7d4e7c1e3a35aeb6995abbfdaf15be80f3e07715524b2bb4449d871010` |
| `P-KLAEGER-SRC` | same, `aan4368_Table_S2.xlsx` original supplement | provenance only | — | `d28b91e62e78e5e011b60da27672875621fef5cdabbea793ac9cce4b98db2c32` |
| `P-PDSP` | NIMH PDSP Ki database full CSV export, `https://pdsp.unc.edu/databases/kiDownload/download.php` | `Ki` (nM) -> `pKi` | 98,678 rows, 274 targets, 16,698 ligands | `45c9a18ac30f1fad350d1dde186bc1f226c5a75d474ca50f50713852a5637ac6` |
| `A-KLIFS` | KLIFS `kinase_information?species=HUMAN` | 85-residue aligned pocket, family, group, UniProt | 555 kinases, 521 with 85-mer | `479e8c3ad148f737...` (recorded in acquisition manifest) |

**Explicitly excluded by repository governance and not downloaded:** DAVIS
(`davis_affinity.csv`, `davis_proteins.csv`), Anastassiadis
(`anastassiadis_matrix.csv`). PKIS2 and Anastassiadis are consumed development
panels and are not used. No ChEMBL37 affinity value is read by XP1; the ChEMBL
`X1`/`X2` prohibition is untouched.

## 4. Endpoint semantics, censoring and analysis blocks (frozen)

- `P-METZ` is left-censored at `pKi = 4.0`; 59.16% of all cells sit exactly at the
  floor. Values are rounded to 0.1 log units (64 distinct levels).
- `P-KLAEGER` is floored at `pKd_app = 5.0`; 93.27% of cells sit at the floor.
- `P-PDSP` carries an explicit censoring flag `ki Note`: 34,300 rows `>` and 55
  rows `<`. All flagged rows are dropped from continuous analysis.

Frozen analysis blocks:

- **`BLK-METZ-60`** (primary). Greedy peel of `P-METZ` to uncensored density
  `>= 0.60`, deleting the row or column whose deletion most raises density,
  breaking ties toward columns. Result fixed at **704 compounds x 82 kinases**,
  57,728 cells, 34,764 uncensored (density 0.602). Censored cells are excluded
  from fitting and from every metric; a registered sensitivity arm repeats the
  primary contrast on `BLK-METZ-70` (704 x 41, density 0.702).
- **`BLK-PDSP-H`** (secondary). Human rows, `ki Note` null, `ki Val > 0`,
  `pKi = 9 - log10(Ki_nM)`. Cells are `(Unigene, Ligand ID)`; replicate rows in a
  cell are retained for the noise-ceiling estimate and averaged for the transfer
  arms.
- **`BLK-KLAEGER-B`** (tertiary, replication only). Binarised hit matrix
  `1[pKd_app > 5.0]`, restricted to kinases with `>= 5` hits and drugs with
  `>= 5` hits.

## 5. Arms

All arms predict the same held-out cells and are scored by the same metric. All
arms except `A0` and `A1` receive **the true protein's** support-derived scalar
offset, so that the correct-vs-wrong contrast isolates `gamma` and cannot be won
by a protein-level potency shift.

| Arm | Definition | Role |
|---|---|---|
| `A0` | `yhat = mu_hat` | population control |
| `A1` | `yhat = mu_hat + alpha_hat_i` | **ligand-only** `f(l)` (control A) |
| `A2` | `yhat = mu_hat + alpha_hat_i + beta_hat_j(S_j)` | additive + support location; *the* baseline for interaction (protein-ID shortcut control E) |
| `A3` | `A2 + gamma_hat_ij(protein features of j)` | **zero-shot correct protein** `f(l,p)` (control B) |
| `A4` | `A2 + u_i . v_hat_j(S_j)` , rank `r` | **few-shot interaction section**; `v_hat_j` from support residuals |
| `A5` | `A3` with features of a deranged held-out protein `pi(j) != j` | **wrong-protein control** (control C) |
| `A6` | `A4` with support residuals taken from a deranged held-out protein | **permuted-support control** (control D) |
| `A7` | `A2 + gamma_hat_ij` from a **random orthogonal** protein embedding | feature-content null |
| `AO1` | `A2 + ` best rank-`r` reconstruction of the held-out protein's own `gamma` column fitted on **all** its cells | oracle ceiling at rank `r` |
| `AO2` | `A2 + ` the held-out protein's own observed `gamma` | oracle ceiling (upper bound, diagnostic) |

`mu_hat`, `alpha_hat`, `u_i` and every kernel/regression parameter are estimated
on **training proteins only**. `beta_hat_j` and `v_hat_j` use only the `k`
support cells of protein `j`. Test cells are disjoint from support cells.

## 6. Splits, closure and leakage constraints

- **Unit of independence: the protein closure component.** Primary closure =
  KLIFS **family**; strict closure = KLIFS **group**; numeric closure =
  single-linkage clusters at KLIFS 85-mer pocket identity `>= 0.60`. Whole
  closure components are held out together. Reported for all three.
- **Cross-validation:** 5-fold over closure components, all folds reported.
- **Support/test disjointness:** within a held-out protein, `k` support ligands
  are drawn at random (seeded), the remaining ligands are test. `k` is registered
  at `{4, 8, 16, 32, 64}` and reported as a curve; `k = 16` is the primary.
- **Ligand main effect:** `alpha_hat_i` is fitted on training proteins only, so
  no held-out protein value enters `alpha_hat`.
- **Assay shortcut (control H):** `BLK-METZ-60` is a single laboratory, single
  assay panel, so assay identity is constant by construction and cannot act as a
  latent label. In `BLK-PDSP-H` the radioligand (`Hotligand`) and tissue
  (`source`) are recorded and a radioligand-stratified sensitivity arm is run.
- **Ligand/scaffold shortcut (control F):** `P-METZ` publishes no structures, so
  unseen-ligand generalization is **not** claimed on Metz; ligands are shared
  across train and test proteins by design (this is the meta-learning regime).
  A scaffold-disjoint arm is run only on `BLK-PDSP-H`, where SMILES exist.
- **Seeds:** 5 seeds `{0,1,2,3,4}` for support draws and derangements. Every
  reported quantity carries a 95% CI from a **cluster bootstrap over held-out
  closure components** (2000 resamples), not over cells.

## 7. Metrics

Primary metric: **`R2_gamma`** — variance of the held-out interaction residual
explained, relative to arm `A2`:

```
R2_gamma(arm) = 1 - SSE(arm) / SSE(A2)
```

computed on held-out uncensored test cells, pooled over folds; and
**`Delta RMSE`** in log units. Reported contrasts:

- `Delta_protein = RMSE(A1) - RMSE(A3)`  (protein information vs ligand-only)
- `Delta_interaction = RMSE(A2) - RMSE(A3)` and `RMSE(A2) - RMSE(A4)`
  (interaction information beyond additive protein location)
- `Delta_specific = RMSE(A5) - RMSE(A3)` and `RMSE(A6) - RMSE(A4)`
  (correct vs wrong protein / wrong support)

A rank readout is reported only as a secondary descriptive, because `E-AFF-R0`
established that it is blind to the location channel.

## 8. Predictions that separate the hypotheses

| Observable | H1 | H2 | H3 |
|---|---|---|---|
| XP1-A: `gamma` reproducibility across disjoint compound halves / independent platform | high | ~0 | high |
| `AO1`/`AO2` oracle `R2_gamma` on held-out proteins | high | ~0 | high |
| `A3` (zero-shot protein features) `R2_gamma` | ~0 | ~0 | ~0 or low |
| `A4` (few-shot support section) `R2_gamma` | ~0 | ~0 | **high** |
| `Delta_specific` for `A3`/`A4` | ~0 | ~0 | `> 0` |

H2 is distinguished by the **oracle** collapsing. H1 and H3 both keep the oracle
high and are separated by whether the *support-identified section* (`A4`)
recovers it while the *feature-conditioned surface* (`A3`) does not.

## 9. Gate (success / failure criteria), frozen

A candidate biological statistic is **admitted** only if, on `BLK-METZ-60` under
**group-level** closure (the strict split), all of:

1. `R2_gamma >= 0.05` with cluster-bootstrap 95% CI lower bound `> 0.02`;
2. `Delta_interaction > 0` with 95% CI lower bound `> 0` (correct protein beats
   the additive-plus-support baseline `A2`);
3. `Delta_specific > 0` with 95% CI lower bound `> 0` (correct protein beats the
   deranged-protein control);
4. `Delta_protein > 0` with 95% CI lower bound `> 0` (correct protein beats
   ligand-only `A1`);
5. the effect survives the feature-content null `A7` (CI of `A3 - A7` excludes 0);
6. the statistic's dimension is compatible with section 10.

Failure criteria: any of (1)-(5) violated, or the effect present only under the
permissive pocket-identity closure and absent under group closure.

**A numerically positive but experimentally negligible effect is a failure.**
`R2_gamma < 0.05` is registered as negligible.

## 10. Theory-imposed constraints on admissible `z`

From `theory/FINAL_FROZEN_THEORY`:

- `z = z(S,Q,gamma)` must be computable from the **observable support set, the
  query, and the declared specification** — at both meta-training and inference.
  A statistic requiring the query label is inadmissible.
- `Z` is a **compact metric domain represented by a finite union of compact
  cubes**, and the hypothesis class is a multilinear sieve on a mesh of `Z` with
  `nu_N ~ r_N^{-dim Z}` nodes and `D_N = (m+1) nu_N` parameters. The
  generalization term is `Gamma_N = O(sqrt(D_N log(Lambda N)/N))`. Therefore
  **`dim(Z)` is the binding resource**: every added biological coordinate costs
  exponentially in mesh nodes. An admissible biological `z` must be **bounded and
  low-dimensional** (single digits), not a 288-D or 1280-D embedding.
- `B(z)` depends on `z` only through the **finite** context map `kappa(z)`; all
  continuous `z`-dependence flows through `F(z) in Delta_m` mixing fixed anchors.
  So a biological coordinate can only move the emitted law inside the convex hull
  of `m+1` fixed bands.
- `(S-CONT)` requires the conditional base risk `L_0(z, beta)` to be uniformly
  continuous in `z`. This is a *regularity* requirement, not sufficiency: the
  theory optimizes against `E[L | zeta = z]`, so a `z` that discards
  affinity-relevant information is still mathematically legal — it simply moves
  the achievable risk floor. **Informativeness is therefore not admission;
  admission requires the empirical Gate in section 9.**

Consequently XP1 reports, for every candidate statistic, the **rank of the
interaction section** actually needed (`r` in `A4`, and the singular spectrum of
`gamma`), because that number *is* the required `dim(Z)` contribution.

## 11. Decisions that follow

| Outcome | Decision |
|---|---|
| XP1-A reproducibility ~0 and oracle ~0 | `DATA_IDENTIFIABILITY_FAILURE` for this panel class; stop expanding architecture; the requirement becomes a better panel, not a better model |
| Oracle high, `A3` and `A4` both ~0 | `REPRESENTATION_FAILURE`; redesign the biological representation |
| Oracle high, `A4` passes Gate, `A3` fails | `OBJECTIVE_OR_PARAMETERIZATION_FAILURE`; the interaction is a support-identified low-dimensional section, not a feature-conditioned surface; redesign the interface between the biological statistic and the frozen operator |
| `A3` (or `A4`) passes the **full** Gate incl. group closure and `Delta_specific` | candidate `BIOLOGICAL_Z_ADMITTED` **pending** an independent-source replication and a sealed novel-target transfer Gate, per `task.md` |
| Mixed / underpowered | `EVIDENCE_INSUFFICIENT` |

Nothing in XP1 authorizes DAVIS access, recipient-label access, ChEMBL `X1`/`X2`,
modification of the frozen operator, `model/` promotion, or an end-to-end DTA
claim. XP1 is a diagnostic on external panels.
