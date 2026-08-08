# Preregistration XP4 — Many-Panel Interaction Identification

Registered: 2026-08-08, before any XP4 arm, model, threshold or metric was
scored. `PREREG_XP1.md`, `PREREG_XP2.md` and every prior verdict stand unchanged.

## 0. Why this design, and why the previous ones could not work

Two structural obstructions have been established, and they are different:

- **`E-AFF-X0-FEAS`:** in a ChEMBL-style corpus, document-keyed panels and
  homology-document closure are produced by opposite kinds of study, so crossing
  and independence cannot coexist under a *cross-document rectangle* unit.
- **XP1/XP2:** a single large kinase panel has plenty of ligand-scaffold
  components (258) but its protein-side independence is capped by kinase
  taxonomy at **8** Manning groups. Every protein-side interval in XP1/XP2 rests
  on 8 clusters.

XP4 uses a third design that removes both: **many small, internally complete,
mutually independent crossed panels.** Crossing is *within* document; independence
is *across* documents. The label-blind census (`XP3_PANEL_GEOMETRY.json`) finds
**87** such panels after excluding every consumed and prohibited source, with a
median geometry of 3 targets x 29 ligands at 94% density.

This also fixes an estimand defect that Metz always had: because each panel is one
paper with one protocol, within-panel centring removes assay heterogeneity
exactly, and no cross-panel absolute comparison is ever required.

## 1. Hypothesis

> A deployment-observable biological basis `x(P,L)`, computed from protein
> sequence and ligand structure alone, predicts the **within-panel protein-by-
> ligand interaction residual** on held-out panels and held-out ligand scaffolds,
> beating capacity-matched random features and every pairing control.

This is the `B0` rung of the information ladder, evaluated for the first time with
an adequate number of independent protein-side units. `B1`-`B3` are only
justified if `B0` shows signal that a richer representation could sharpen, or if
`B0` fails in a way a richer representation could specifically repair.

## 2. Data release (frozen)

| ID | file | SHA-256 | licence |
|---|---|---|---|
| `P-BDB-ART` | `BindingDB_BindingDB_Articles_202608_tsv.zip` | `2529b1c572aa7b29...` (full value in the acquisition manifest) | CC BY 3.0 |
| `P-BDB-ASSAY` | `BindingDB_Assays_202608_tsv.zip` | `e6ce48748603b669...` | CC BY 3.0 |

Excluded by PMID before any value is read, and verified to contribute **0 rows**:
Metz 2011 `21572424`, Klaeger 2017 `29191878`, Davis 2011 `22037378`
(prohibited), Karaman 2008 `18183025`, Anastassiadis 2011 `21949673`,
PKIS2 `28767711`. No ChEMBL37 affinity value is read. DAVIS and recipient label
reads remain `0`.

## 3. Panel construction (frozen)

`BLK-BDB-PANELS`:

1. rows with a `Ki (nM)` value present, single-chain targets only
   (`Number of Protein Chains <= 1`), `Homo sapiens`;
2. non-null `Ligand SMILES` that RDKit parses, non-null UniProt ID and target
   sequence, non-null PMID;
3. drop the six consumed/prohibited PMIDs;
4. keep documents with `>= 2` distinct targets and `>= 20` distinct ligands;
5. within a document, a cell is the `(target, ligand)` pair; replicate rows are
   averaged in `pKi = 9 - log10(Ki_nM)`.

Expected from the census: 87 panels, 8,413 cells. The realised index and its
SHA-256 are recorded in the panel manifest.

## 4. Estimand

Within panel `d`, for target `j` and ligand `i`:

```
y_dij = mu_d + alpha_di + beta_dj + gamma_dij + eps_dij
```

Main effects are fitted **within each panel**, so panel-level assay offset,
panel-level potency scale and target druggability are all removed by
construction and cannot be exploited. The estimand is `gamma_dij`.

## 5. Closures (frozen)

- **Panel closure (primary independence unit):** whole documents are held out.
- **Protein closure across panels:** UniProt accessions are additionally merged
  by MMseqs2 clustering at 40% sequence identity, 50% coverage; whole clusters
  are held out with their panels. A panel is assigned to the fold of its
  most frequent protein cluster; any panel touching a held-out cluster is held out.
- **Ligand-scaffold closure:** Bemis-Murcko scaffolds merged by single-linkage
  ECFP4 Tanimoto `>= 0.5`; whole scaffold components are held out.
- **Double held-out test:** test cells are (held-out panels) x (held-out scaffold
  components), simultaneously.

## 6. Arms (frozen)

All arms predict the same held-out `gamma` cells.

| arm | definition | role |
|---|---|---|
| `Z0` | `gamma_hat = 0` | the estimand's own null; within-panel centring already removes both main effects |
| `BILIN` | `gamma_hat = phi(P)^T W psi(L)`, `W` low-rank `r <= 8`, fitted across training panels | **the hypothesis** |
| `RAND-P` | protein features replaced by capacity-matched random vectors | protein-content null |
| `RAND-L` | ligand features replaced by capacity-matched random vectors | ligand-content null |
| `RAND-BOTH` | both replaced | joint null |
| `PERM-PAIR` | protein-ligand correspondence permuted within the training panels | pairing null |
| `FOREIGN-P` | at test time, a foreign protein's features are substituted | specificity control |
| `ORACLE-R` | best rank-`r` reconstruction of the held-out panel's own `gamma` | ceiling |

`phi(P)` = frozen ESM-2 t30 mean-pooled embedding. `psi(L)` = Morgan r=2, 1024
bit. Both are frozen encoders; only the low-rank `W` is fitted. Ridge strength
and rank are selected by nested panel-grouped CV **inside the training panels
only**.

## 7. Metrics and inference

Primary: `R2_gamma = 1 - SSE(arm)/SSE(Z0)` on held-out cells, and paired MSE
contrasts. The independence unit is the **panel**; intervals are cluster
bootstraps over panels (2,000 resamples), reported alongside a scaffold-component
clustering, and the **wider** interval is used for every decision.

Registered contrasts:

- `Delta_bio = MSE(RAND-BOTH) - MSE(BILIN)` — biology over capacity-matched noise;
- `Delta_protein = MSE(RAND-P) - MSE(BILIN)` — protein content specifically;
- `Delta_ligand = MSE(RAND-L) - MSE(BILIN)`;
- `Delta_pair = MSE(PERM-PAIR) - MSE(BILIN)`;
- `Delta_specific = MSE(FOREIGN-P) - MSE(BILIN)`;
- `Delta_null = MSE(Z0) - MSE(BILIN)`.

## 8. Gate (frozen, identical floors to XP2 for comparability)

`INTERACTION_OBSERVABLE_FROM_DEPLOYMENT_INPUTS` requires **all** of:

1. `R2_gamma >= 0.05` with cluster-bootstrap 95% CI lower bound `> 0.02`;
2. `Delta_bio` CI lower bound `> 0`;
3. `Delta_protein` CI lower bound `> 0`;
4. `Delta_pair` CI lower bound `> 0`;
5. `Delta_specific` CI lower bound `> 0`;
6. all of the above under the **double** (panel x scaffold) held-out split,
   opened once.

A numerically positive but negligible effect (`R2_gamma < 0.05`) is a failure.
Thresholds may not be changed after the test set is opened. One registered
configuration is run; no architecture or hyperparameter sweep beyond the
nested-CV selection declared in §6.

## 9. What each outcome means

| outcome | verdict | next |
|---|---|---|
| Gate passes | `INTERACTION_OBSERVABLE_FROM_DEPLOYMENT_INPUTS` | proceed to the `k <= 5` section stage on this design |
| `ORACLE-R` high, `BILIN` ~ `RAND-BOTH` | `DEPLOYMENT_INPUTS_INSUFFICIENT` at rung `B0` | justifies `B1`/`B2` only if the missing information is nameable |
| `ORACLE-R` also ~0 | interaction is not estimable above noise in this design | `PUBLIC_DATA_INSUFFICIENT_FOR_DEPLOYABLE_INTERACTION_IDENTIFICATION` |

No XP4 outcome authorises DAVIS access, `model/` promotion, modification of the
frozen theory, CSMO, Band, `K` or the mesh, or any end-to-end DTA claim.
