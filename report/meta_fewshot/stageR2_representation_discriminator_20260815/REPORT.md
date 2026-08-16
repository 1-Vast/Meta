# Stage R2: frozen representation discriminator

Numerical authority: `DISCRIMINATOR_meta_val.json`. Criteria fixed in
`PREREGISTRATION.md` before any result. Population: double-cold `meta_val`,
41 targets, 19 components, 1,411 cells, of which 1,141 are in the `< 0.40`
low-similarity tier. No training. `meta_test` untouched.

## Ligand representations

Neighbour count fixed at 10 a priori; 1 and 25 reported for robustness only.

| arm | **L1 continuity** | L2 MSE (< 0.40 tier) | L2 MSE (all) | L3 CI | L3 Spearman | external data |
|---|---:|---:|---:|---:|---:|---|
| **`morgan`** *(incumbent)* | **0.2690** [0.189, 0.345] | 2.5533 | 2.4648 | 0.5270 | 0.0637 | none |
| `chemberta` | 0.2604 [0.177, 0.342] | **2.2595** | **2.1502** | **0.5459** | **0.1288** | yes |
| `rdkit_desc` | 0.2298 [0.167, 0.293] | 2.2464 | 2.2453 | 0.5329 | 0.0973 | none |

**Decision: retain `morgan`.** The preregistered rule takes the best L1, and
Morgan wins it.

The honest qualification is that this is a **close and unresolved** call.
Morgan's and ChemBERTa's continuity intervals overlap almost completely
(0.269 [0.189,0.345] against 0.260 [0.177,0.342]), and ChemBERTa is better on
every other column — 11.5% lower MSE in the low-similarity tier and double the
Spearman. Had L2 been preregistered as primary, ChemBERTa would have won. The
rule was fixed first and is followed; the disagreement between continuity and
kNN accuracy is recorded as the reason to revisit ChemBERTa if the ligand branch
later proves to be the binding constraint. Retaining Morgan also keeps the
architecture claim free of external pretraining data, which has independent
value.

## Protein representations

| arm | **P1 level MSE** | P2 specificity | P2 interval | top-16 similarity spread |
|---|---:|---:|---|---:|
| `esm_pooled` *(incumbent)* | 2.5707 | +0.4648 | [-0.177, +1.154] | **0.0222** |
| `esm_centered` | 3.1069 | **-0.4916** | [-1.069, +0.110] | 0.2286 |
| **`esm_whitened`** | **2.4863** | +0.4870 | [-0.076, +1.045] | **0.3098** |
| `kmer3` | 2.7614 | +0.3202 | [-0.083, +0.694] | 0.0515 |

**Decision: select `esm_whitened`** under the preregistered rule (best P1, P2
point estimate positive).

Two things this measured that matter more than the 3.3% P1 gain:

* **The audit's compression finding is confirmed and repaired.** The spread of
  cosine similarity across the 16 nearest training targets is 0.022 for raw
  pooled ESM and 0.310 after `meta_train` whitening — a fourteen-fold widening.
  Any softmax over the raw band is nearly uniform by construction.
* **Centring alone is actively harmful** (P1 3.11, specificity **-0.49**):
  removing the mean without rescaling the axes amplifies low-variance
  directions. Only the full whitening helps. This falsifies the natural reading
  of the audit finding, which was that centring would be enough.

**No protein arm has resolved specificity** — every P2 interval crosses zero on
19 components. So the correct statement is that whitening produces a *sharper*
protein similarity, not a *demonstrably protein-specific* one.

## What is carried into Stage R3, and what is not

* **Ligand: unchanged.** Morgan/Tanimoto remains the chemical similarity, and
  the learned GINE remains the ligand encoder.
* **Protein: unchanged as a trunk input.** The trunk's protein projection is a
  learned linear map, so it can already express a whitening; feeding whitened
  vectors would not add capacity, and P1's advantage is small and unresolved.
* **`esm_whitened` is adopted wherever an explicit protein similarity is
  computed.** After Stage R0 killed retrieval, the live consumer is the
  falsification machinery: Stage R4 selects **similarity-matched wrong-protein
  donors** with it instead of random cross-component donors. That makes the
  wrong-protein control strictly harder to pass, which is the direction a
  control should move.

This is a negative result for representation change and is reported as one. The
bottleneck identified in Stages 9, R0 and here is not which frozen vectors are
used; it is that the trained model has no within-target ordering to begin with —
CI 0.527 for a Morgan kNN on this split, against the trained endpoint's 0.525 on
the old one.

## Not run

`facebook/esm2_t33_650M_UR50D` (incomplete local cache, does not load offline);
MoLFormer, GraphMVP, PMMR, TM-Vec (unavailable offline). The trained `GINE`
encoder was **excluded on contamination grounds**: every existing checkpoint was
trained on the old `meta_train`, which contains ligands that are now double-cold
`meta_val` ligands. These are untested, not rejected.
