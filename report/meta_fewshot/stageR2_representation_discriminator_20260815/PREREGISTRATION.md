# Stage R2 preregistration: frozen representation discriminator

No training, no model change. Run entirely under the **double-cold** protocol
(`bindingdb_ki_double_cold_v1`): development population `meta_val`, 41 targets,
19 components, 0 exact-ligand overlap, 81.6% of ligands below Tanimoto 0.40 to
every training ligand. `meta_test` is not touched.

## Why representations are not selected on MSE

The user constraint is explicit and Stage R0 explains it: ordinary
overlap-heavy `meta_val` MSE rewards recall. Under the double-cold protocol
recall is impossible, but MSE is still dominated by target-level calibration
(Stage 9: 59%), which a ligand representation cannot fix. Selection therefore
uses four criteria, with continuity primary.

## Candidates

Only representations that can be computed offline here are eligible.

**Ligand**

| arm | source | external data |
|---|---|---|
| `morgan` *(incumbent)* | Morgan r=2, 1024 bits, Tanimoto | none |
| `chemberta` | `DeepChem/ChemBERTa-77M-MLM`, mean-pooled, cosine | **yes** — masked-LM pretraining on public SMILES; no affinity labels |
| `rdkit_desc` | 8 physicochemical descriptors, `meta_train`-standardised, cosine | none |

**Protein**

| arm | source | external data |
|---|---|---|
| `esm_pooled` *(incumbent)* | cached ESM2-t30 pooled, cosine | ESM2 pretraining (already in the pipeline) |
| `esm_centered` | the same vectors minus the `meta_train` mean | none beyond the above |
| `esm_whitened` | `meta_train` covariance whitening | none beyond the above |
| `kmer3` | 3-mer profile cosine | none |

**Not run, and why.** The current `GINE` ligand encoder is only meaningful
inside a trained checkpoint, and every existing checkpoint was trained on the
*old* `meta_train`, which contains ligands that are now double-cold `meta_val`
ligands. Using it to select a representation would import that contamination, so
it is excluded. `facebook/esm2_t33_650M_UR50D` has an incomplete local cache and
does not load offline. MoLFormer, GraphMVP, PMMR and TM-Vec are not available
offline. These are recorded as untested, not as rejected.

## Criteria, fixed now

Aggregation is equal-component then equal-target throughout.

| # | criterion | definition | direction |
|---|---|---|---|
| **L1** | affinity continuity | per development target, Spearman between within-target pair similarity and `-abs(y_i - y_j)`. Every pair here is exact-free by construction | higher |
| **L2** | low-similarity performance | MSE of a fixed top-10 `meta_train` neighbour-mean predictor, restricted to the `< 0.40` tier | lower |
| **L3** | within-target ranking | CI and Spearman of that predictor against the truth | higher |
| **P1** | target-level calibration | MSE of a top-16 `meta_train` protein-neighbour target-mean predictor | lower |
| **P2** | protein specificity | P1 under a shuffled similarity vector minus P1 under the correct one | positive |

The neighbour counts (10 for ligands, 16 for proteins) are **fixed a priori and
not tuned**, so no representation can win by having its bandwidth optimised. A
top-1 and top-25 robustness column is reported but does not enter the decision.

## Decision rule

* Ligand: take the arm with the best **L1** provided its **L2** is not worse
  than the incumbent's. Otherwise retain `morgan`.
* Protein: take the arm with the best **P1** provided its **P2 > 0**. Otherwise
  retain `esm_pooled`.
* At most one of each is carried into Stage R3. If neither improves, the
  existing representations are retained and that is the result.

Any external-data arm that wins must have its pretraining corpus disclosed
separately from architecture gains in every downstream report.
