# QPSMP Model Implementation Contract

## Primary Module

`model.qpsmp_meta.QPSMPBioModel` is the candidate learnable innovation. It consumes cached protein
language-model pooled/residue states and ligand molecular graphs, then applies shared trainable
protein and ligand encoders, ligand-conditioned residue localization, a crossed scalar potential,
and a centered neural support adapter.

The present implementation is the fixed standardized-Ki-context specialization. It does not encode
assay or panel context and therefore makes no cross-context generalization claim.

The primary output decomposes as

```text
zero_support = additive + cross_zero_shot
few_shot_prediction = calibrated_level + scaled_zero_shot_shape + scaled_SAR_adaptation
```

G2 must evaluate the crossed channel against additive, ligand-only, matched wrong-protein, shuffled
protein, and design-nuisance controls. Full prediction utility is a separate estimand.

## Structural Invariants

- Target identifiers are lookup keys only and never enter model tensors.
- Every support and query ligand uses the same ligand encoder.
- The interaction heads receive only a crossed protein-ligand representation.
- The support adapter is permutation invariant and receives centered residual evidence.
- Zero centered evidence implies an exactly zero neural SAR state.
- Query shape and SAR have separate source-learned reliabilities; neither depends on support evidence.
- At `k=1`, zero-shot query shape is retained and only target-level calibration can adapt.
- Delta and rectangle predictions are differences of retained scalar endpoint predictions.
- A foreign-support control may replace only the transient neural task state.

## Evidence Boundary

The implementation is trainable and its declared tensor interface is covered by focused tests. In a
consumed three-seed k=5 development diagnostic, the repaired complete predictor beat its matched
level baseline in every seed, but the stricter shared-checkpoint nested-k component test did not pass
its lower-confidence-bound criteria. SAR-only gain and protein-specific controls were not stable. No
G2, G3, biological, confirmation, safety, or integration claim is authorized. The compact 128-slot ESM2 bank is a
position-preserving representation, not a validated binding-pocket localizer.

Stable data, training, and governed nested-k evaluation entry points are in `scripts/qpsmp_data.py`,
`scripts/train_qpsmp.py`, and `scripts/evaluate_qpsmp.py`. Historical exploratory runners remain
under `research/` but are not part of the primary QPSMP path.
