# Stage W W1 preregistration — KIBA-only local interaction representation

Frozen **2026-08-17, before any W1 split metric or training/evaluation metric
was read.** W0 is complete: Davis failed and is closed for this surface; KIBA
passed and is the only W1 dataset. This document cannot be changed after the
first W1 split statistic is read.

## 0. Goal

On KIBA only, test whether a local protein × ligand-transformation interaction
representation can predict the crossed double difference

    D(family, p1, p2) = delta_y(p1, family) - delta_y(p2, family)

on protein components held out from training, where `family` is the frozen W0
soft family (murcko core + attachment environment + pharmacophore change
category). A positive result must be a **candidate surface result**, not a
final claim; independent replication remains mandatory.

## 1. Data and splits (label-blind)

* Rows: `dataset/raw/dta/kiba.tab` (SHA-256 recorded in W0 artifact).
* Targets: unique target sequences; components: the W0 CD-HIT40 clusters
  recomputed deterministically and saved in `W1_DATA.json`.
* Soft-family observations and target effects as in W0:
  `delta_y(target, family)` = median of same-target observations in that
  family; `D` rows = unordered target pairs within one family.
* **Component split.** Components sorted lexicographically; stable seed
  `20260821` selects 24 components (~19%) as `heldout`; the rest are `fit`.
  Every component is entirely on one side. Split is structure-only.
* Primary evaluation surface `heldout_repeated`: `D` rows whose both targets
  are heldout components and whose family also occurs in fit components.
* Secondary surface `heldout_cold`: heldout `D` rows whose family does **not**
  occur in fit components.
* Training rows: cross-component `D` rows whose both targets are fit
  components.
* `fit_unsampled`: a frozen 10% of fit cross-component D rows, never trained
  on, used only for the in-distribution shortcut diagnostic.

### Frozen W1 split admission (before any training)

1. `heldout_repeated` rows >= **500**;
2. heldout repeated components >= **10**;
3. repeated families >= **50**.

If any fails: stop, write the negative W1 report, train nothing.

## 2. Representation contract

* **Protein branch.** For each unique KIBA target, ESM-2 150M residue/region
  tokens are computed locally (128 ordered slots, 640-d, mask included) and
  recorded with model/hash provenance. The branch keeps ordered slot tokens;
  **no pooled protein summary, no target embedding, no target index, no
  component id, no assay id** reaches the operator.
* **Multiple latent pocket states.** Eight learned pocket queries cross-attend
  over the ordered protein slots; each pocket state is therefore a
  protein-conditioned region summary, not one global vector.
* **Ligand pharmacophore tokens.** The transformation context (core + R_a +
  R_b) is tokenised per heavy atom with pharmacophore features (element,
  aromatic, in-ring, HBD/HBA, formal charge, degree, hybridization). A
  two-layer self-attention encoder refines these tokens; tokens are preserved
  until the interaction readout.
* **Local interaction slots.** Each refined ligand token cross-attends over the
  eight pocket states, producing interaction slots. A learned interaction query
  aggregates the slots into one scalar `R(tau, p)`; this is a late readout,
  after the local cross-attention, and is not an early global pooling.
* `D_hat(tau,p1,p2) = R(tau,p1) - R(tau,p2)`. Identity, protein-pair
  antisymmetry and cycle consistency hold for every parameter setting.
* No Cartesian/3D geometry, no conformer, no pocket-structure supervision is
  used. "Pocket state" is a latent region summary and may never be called a
  biologically localized pocket.

## 3. Arms (same budget, identical init and batch order)

| arm | description |
|---|---|
| A `A_zero` | constant response; `D_hat` identically 0 |
| B `B_global` | edit token + global ESM pooled protein summary (negative reference) |
| C `C_local` | candidate local interaction operator |
| D `D_local_shuffled` | C architecture, trained on stable cross-component shuffled proteins |
| E `E_mean_tokens` | C architecture, target-independent mean protein slots (structurally `D_hat=0`) |
| F `F_label_shuffled` | C architecture, correct protein, within-family permuted `D` labels |

Wrong protein is evaluation only and drawn from the recipient's own population
(fit for fit rows, heldout for heldout rows), different CD-HIT40 component,
most similar admissible by ESM pooled cosine. Shuffled protein substitution is
a stable cross-component permutation within the same population. Residue-slot
permutation and a capacity-matched random protein representation are evaluated
as substitutions on the trained C arm.

Training: ordinary AdamW forward/backward; no ridge, no pseudoinverse, no
closed form, no correct-vs-wrong loss. Seed `20260821`; steps `3000`; batch
`256`; lr `3e-4`; weight decay `1e-4`; cosine schedule; Huber delta `1.0`;
gradient clip `5.0`; row sampling weight `1/sqrt(family_degree)`. No
checkpoint selection: fixed final parameters.

## 4. Frozen W1 success conditions

All on `heldout_repeated` unless stated. Paired contrasts use identical rows.
Two-way cluster bootstrap: components and families resampled with replacement,
row multiplicity = mean endpoint-component draw count × family draw count;
2,000 draws, seed `20260820`; effective independent units reported.

1. `C_local` vs `A_zero` Pearson >= **+0.05** AND vs `B_global` Pearson
   >= **+0.05** AND vs `E_mean_tokens` Pearson >= **+0.05** (constant arms
   have Pearson 0 by convention);
2. all three contrasts have bootstrap 95% lower bounds **> 0**;
3. C correct-input minus C shuffled-input Pearson >= **+0.05**, lower bound > 0;
4. C correct-input minus C matched-wrong-input Pearson >= **+0.05**, lower
   bound > 0;
5. C vs D arm: MSE delta < 0 AND Spearman delta > 0 AND CI delta > 0 AND sign
   accuracy delta > 0;
6. C minus F label-shuffled Pearson >= **+0.05**, lower bound > 0;
7. protein-induced shift (C correct minus C shuffled) Pearson with truth
   >= **+0.10**;
8. on `fit_unsampled`, C correct minus C shuffled Pearson >= **+0.05**, lower
   bound > 0;
9. leave-one-out influence: no single family or component accounts for more
   than **50%** of the C-vs-shuffled effect;
10. on `heldout_cold` (if >=100 rows), C correct minus C shuffled Pearson
    >= **0**; if <100 rows the surface is `not_evaluable` and the route cannot
    pass.

Residue-permuted and capacity-matched-random protein substitutions must be
reported; they are falsification controls and must not beat correct protein.

A single seed may reject; it may not confirm. If all ten pass, three fixed
seeds (`20260821/22/23`) run all six arms.

## 5. Stop rules

* Split admission fails -> stop, negative report, no training.
* Single-seed gates fail -> close **this KIBA soft-family local operator**; do
  not claim all local interaction representations are impossible.
* Single seed passes but multi-seed fails -> NOT CONFIRMED, nothing promoted.
* Multi-seed passes -> write the next-stage meta-learning proposal; do not
  integrate into `model/` or `scripts/`.

## 6. Verification

Structural tests must pass before training: antisymmetry/identity/cycle;
gradient coverage; no target-ID bypass; no label in splits; no cross-split
leakage; deterministic banks across `PYTHONHASHSEED`; provenance for the ESM
bank. Commands, environment, raw rows and JSON artifacts are retained in this
directory.
