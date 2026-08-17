# Stage Q report — decoupled frozen-feature level head: REJECTED (G3)

Development evidence, single seed, meta_val read once after freezing;
meta_test sealed. Authorities: Q_vs_T2.contrast.json, per-arm row summaries,
Q0_JOINT_FROZEN_IDENTIFIABILITY.json, PREREGISTRATION.md, per-arm RESULT.json.

## Verdict

**Rejected by G3; nothing promoted.** The Q0 probe met its preregistered
threshold (joint frozen-feature level MSE 1.3416, the best frozen predictor
on record), but the trained composition fails: k=0 MSE -0.1479 (unresolved)
while Spearman/CI degrade with RESOLVED intervals at k=0 (Spearman -0.0588
[-0.1174, -0.0008]), k=2 (-0.0501) and k=3 (-0.0496). Stop rule S2 fires.

## What this closes

The decoupling hypothesis is falsified: the head consumed only FROZEN
features (ESM bank vector, handcrafted panel statistics, journal table), so
it could not reshape the trunk through its own gradient path — yet the
trunk's ordering at k>=1 still degrades with resolved intervals. The cause
is the k=0 training signal itself: adding the head's scalar to the zero-shot
endpoint changes the k=0 objective for the trunk's own parameters
(protein_head and the interaction path), and that reshapes the shared
representation. Four compositions of a k=0-specialized level mechanism have
now failed the ranking gate: E (ungated, coupled features), J (ungated,
assay covariates), L (gated, coupled features), Q (gated, decoupled
features). The conflict is between the zero-shot level objective and
within-target ordering on one shared trunk, independent of gating and of
feature coupling.

Q-UNGATED confirms the gate is necessary but insufficient: always-on, the
head wrecks k=1 MSE (2.3095) through double level fitting, while its k=5
ranking (Spearman 0.3445, CI 0.635) is the best in the record — the
level/ranking trade-off of this trunk is fundamental.

## Final composition map

The only k=0-improving compositions that preserve ranking would need a
level mechanism whose training never touches the shared trunk — i.e. a
separately trained calibrator composed at inference, which the governing
contract excludes as a multi-stage regime. Within single-stage end-to-end
training, the k=0 frontier of this trunk family is T2 (2.5961) with the
Q0 probe as an unrealisable bound. The bounded conclusion
(report/BOUNDARY_20260817_NIGHT.md) stands.
