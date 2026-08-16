# UBSE-G1 centered residue-contact student preregistration

Date: 2026-07-29  
Status: frozen before binding-residue load, ESM inference, or student training

## Claim under test

G0R established that BioLiP binding-residue lists contain a modest
cross-publication ligand-connectivity signal. G0PB established a feasible
same-target, same-PubMed, same-scaffold panel topology with an untouched
88-panel homology/scaffold/PubMed-cold audit set.

G1 tests the missing deployability claim:

> Can frozen sequence and 2D-ligand covariates predict *where two ligands on
> the same target and scaffold differ in residue contact*, beyond an exact
> additive target-plus-ligand null?

This is an affinity-blind Stage-1 semantic gate. It is not an affinity model,
coordinate teacher, or causal perturbation model.

## Frozen manifests and firewall

- Panel manifest:
  `dataset/public/biolip2/processed/ubse_g0pb_panels.parquet`
- Required manifest SHA-256:
  `4fea01e332eb3c60e41d76d5062d33cc95b13bc2e96b01df226532f78fe1b371`
- BioLiP source:
  `dataset/public/biolip2/processed/closed_registry.parquet`
- Required source SHA-256:
  `7905e4edf88073f564baa4b2d4fb50d496432bc4e15e97cccbfa0766b1b0638d`
- Allowed source columns:
  `target_key`, `sequence`, `pubmed`, `scaffold`, `conn`,
  `binding_residues_reindexed`
- Forbidden:
  every affinity field/value, coordinates, assay outcomes,
  development/confirmation features or labels, and sealed outcomes.

The G0PB `audit` rows are immutable and cannot be used for normalization,
early stopping, hyperparameter selection, or model selection.

## Frozen fit/validation split

Starting only from G0PB `train_candidate` panels:

1. run the G0P deterministic conflict-free packing with seed 1730;
2. take the first 64 panels as validation;
3. remove from fit every panel sharing validation homology component,
   scaffold, or PubMed;
4. keep the G0PB 88-panel audit unchanged.

This split is completed before reading contact labels. Training requires at
least 800 fit panels and 50 valid validation panels after label parsing.

## Frozen contact target

Parse every whitespace token in `binding_residues_reindexed` as one amino-acid
letter plus a 1-based integer sequence position. Reject an entire panel if any
listed ligand has:

- an empty or malformed residue list;
- a position outside its exact sequence;
- an unparseable 2D ligand; or
- a duplicate panel ligand after canonicalization.

For panel \(p\), ligand \(l\), and residue \(i\), define:

\[
Y_{pli}=1\{i\text{ is listed as a binding residue}\}.
\]

The within-panel centered target is:

\[
\widetilde Y_{pli}
=Y_{pli}-\frac{1}{|L_p|}\sum_{l'\in L_p}Y_{pl'i}.
\]

A contrast residue has nonconstant \(Y_{pli}\) across the panel's ligands.
At least 70 audit panels, 50 validation panels, and 800 fit panels must have
one or more contrast residues; otherwise G1 stops before GPU training.

## Frozen covariates

### Protein

- Model: `facebook/esm2_t6_8M_UR50D`
- Revision:
  `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`
- Local cached weights only; final hidden state; backbone frozen.
- Maximum window: 1,000 residues with 100-residue overlap.
- Overlapping residue embeddings are averaged.
- CUDA autocast inference; cached residue embeddings are float16.

### Ligand

Use the existing audited `src.data.dualcold.ligand_features` contract:

- Morgan radius 2, 1,024 bits;
- ten physicochemical descriptors;
- RDKit `2023.09.6`.

Morgan bits remain binary. Descriptor mean and standard deviation are fit on
fit-panel ligand cells only and then frozen for validation/audit.

## Exact additive null and interaction arm

Both arms have identical trainable parameters and initialization. Let frozen
residue input be \(e_{ti}\), ligand input be \(x_l\), and trainable towers
produce \(u_{ti},v_l\in\mathbb R^{64}\). Both arms contain the same tower
layers, LayerNorm, dropout, and scalar heads:

\[
A(t,l,i)=a(u_{ti})+b(v_l).
\]

The null prediction is exactly:

\[
P_{\mathrm{null}}(t,l,i)=A(t,l,i).
\]

The interaction prediction adds a parameter-free multiplicative term:

\[
P_{\mathrm{cross}}(t,l,i)
=A(t,l,i)+\frac{u_{ti}^{\mathsf T}v_l}{\sqrt{64}}.
\]

Thus the cross arm has no extra trainable parameter. Setting the
multiplicative term to zero recovers the exact null. A concatenated
target/ligand retrieval embedding is not an interaction certificate.

## Frozen optimization

- seeds: 1729, 1730, 1731;
- identical initial state for null/cross within each seed;
- AdamW, learning rate `3e-4`, weight decay `1e-4`;
- batch size 8 panels, 30 fixed epochs, no early stopping;
- dropout `0.10`, gradient clip norm `1.0`;
- CUDA float32 student training;
- identical panel order and update count for both arms.

For valid cells, use:

\[
\mathcal L
=\mathcal L_{\text{balanced BCE}}(Y,P)
+4\,\operatorname{MSE}
\left(
\widetilde Y,\,
\operatorname{center}_{l\in p}\sigma(P)
\right),
\]

where the centered MSE is evaluated only at contrast residues. Balanced BCE
is the mean positive loss plus mean negative loss divided by two. No audit
metric selects an epoch or seed.

## Frozen controls

Evaluate each trained arm on the same audit panels:

1. additive null with correct inputs;
2. cross arm with correct inputs;
3. cross arm with ligand features cyclically deranged within each panel;
4. cross arm with every residue embedding replaced by that target's
   sequence-mean embedding (`protein-free-position`).

The two destruction controls do not retrain the model.

## Frozen metrics

For every panel, center predicted probabilities over ligands.

1. **Directional accuracy:** for every unordered ligand pair and every residue
   where their labels differ, score whether the sign of the predicted
   difference matches the observed difference; an exact tie receives 0.5.
2. **Centered cosine:** cosine between flattened centered prediction and
   centered label over contrast residues.
3. **Assignment Recall@1:** cosine-match every predicted ligand residual to
   all observed ligand residuals in the same panel; award `1/k` for a correct
   label tied among `k` maxima.
4. **Centered MSE:** panel-macro squared error on contrast cells.

Report medians across the three seeds. For directional-accuracy deltas,
average each panel over seeds and use 2,000 seed-1729 panel-bootstrap
replicates. The audit panels are already conflict-free independent units.

## Frozen gates

All must pass:

1. **S1 substrate:** at least 800 fit, 50 validation, and 70 audit panels with
   contrast labels; all 88 audit panel identities remain unchanged.
2. **S2 nontrivial cross signal:** median-seed audit panel-macro directional
   accuracy of the cross arm is at least `0.60`, and centered cosine is at
   least `0.10`.
3. **S3 exact-null increment:** cross minus additive-null directional
   accuracy is at least `0.05`, its panel-bootstrap lower 95% bound is greater
   than zero, centered-cosine improvement is at least `0.05`, and at least
   two of three individual seeds have positive directional improvement.
4. **S4 ligand destruction:** correct minus cyclically deranged ligand
   directional accuracy is at least `0.05` with panel-bootstrap lower 95%
   bound greater than zero, and at least two of three seeds have positive
   directional destruction.
5. **S5 protein-position destruction:** correct minus
   protein-free-position directional accuracy is at least `0.05` with
   panel-bootstrap lower 95% bound greater than zero, and at least two of
   three seeds have positive directional destruction.
6. **S6 pair assignment:** cross assignment Recall@1 is at least `0.60` and
   exceeds the strongest of null, ligand-deranged, and
   protein-free-position by at least `0.05`.
7. **S7 execution/firewall:** all three seeds finish without nonfinite
   values; parameter counts are identical; CUDA is used for ESM inference and
   student training; no forbidden field or outcome is loaded.

Pass:
`REQUEST_UBSE_G2_CHEMBL_TRAIN_SEMANTIC_SUPPORT_PREREGISTRATION`.

Failure:
`STOP_UBSE_G1_NO_DEPLOYABLE_INTERACTION_RESIDUAL`.

Passing G1 only establishes a coarse residue-contact residual on held BioLiP
domains. Before any affinity fitting, G2 must show target-domain coverage and
nondegenerate predictions on ChEMBL TRAIN covariates under the same additive
null and destruction controls. G1 cannot authorize a 20k-50k coordinate
teacher, Stage-2 affinity readout, confirmation access, or the phrase
"ligand-induced causal transition".

## Expected compute

On the available RTX 4060 Laptop GPU:

- ESM2-8M residue cache: approximately 5-20 minutes;
- six fixed student fits (2 arms x 3 seeds): approximately 15-45 minutes;
- CPU parsing/evaluation/bootstrap: approximately 5-15 minutes.

Total expected wall time: 25-80 minutes, excluding code verification.
