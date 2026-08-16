# MEDIP S0 decision

**Verdict:** `MEDIP_S0_ENGINEERING_CALIBRATION_STOP`

## Gate results

- `S0_RECOVERY`: `True`
- `S0_SEPARATE_ENDPOINTS`: `True`
- `S0_METADATA_DESTRUCTION`: `False`
- `S0_SELECTIVITY_DESTRUCTION`: `False`
- `S0_INTERACTION_NULL`: `True`
- `S0_ARCHITECTURE`: `True`
- `S0_NUMERICS`: `True`
- `S0_REPRODUCIBLE`: `True`

## Primary results

- Correct mixed-difference correlation: `0.998710`.
- Correct same-ligand cross-target ordering accuracy:
  `0.984848`.
- Correct exact-mean RMSE: `0.133927`.
- Merge-minus-correct exact-mean RMSE:
  `1.707017`.
- Correct-minus-metadata-shuffle correlation:
  `0.005862`.
- Metadata-shuffle-minus-correct exact-mean RMSE:
  `1.677035`.
- Correct-minus-selectivity-shuffle ordering:
  `-0.000947`.
- Correct-minus-separable ordering:
  `0.397727`.

## Architecture result

The interaction encoder signature is
`forward(self, target_features, ligand_features)`. All named endpoint and
source parameters are below `observation.*`; changing metadata changes the
observation mean but leaves the latent score bit-identical.

## Ordering amendment and claim boundary

Amendment A1 was frozen before any run. S0 is allowed before `OMUT-I0` only by
the same independent synthetic estimator-falsification logic as R-MAON G0,
correcting the conservative ordering in
`open_evidence_pretraining_crossreview` step 6.

This result used generated coordinates and generated outcomes only. It read no
real affinity value or external dataset file. A pass is engineering evidence
for this module implementation only. It cannot unlock real outcomes, grant
protein or ligand representation credit, establish biological validity,
authorize open-data pretraining, or support a strict dual-cold predictive
claim. All real-data information, provenance, topology, firewall, and power
gates remain binding.
