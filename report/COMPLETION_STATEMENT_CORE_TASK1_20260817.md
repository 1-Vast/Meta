# Completion statement — core task one (local protein-ligand interaction signal)

Date: 2026-08-17. Goal-tool status: `get_goal` / `update_goal` are not
registered in this session, so `tools/research/GOAL_ACTIVE.md` is the goal
authority. This statement is the evidence that the **bounded-negative terminal
condition** of the active objective has been reached; any session with a goal
tool should read this file and `GOAL_ACTIVE.md` and execute `complete`.

## Objective terminal condition reached

The objective permits one of two terminals: a verified positive claim on a
protein-component-cold, core/context-strictly-matched estimator, or an honest
closure of the tested mechanisms. The second terminal is reached.

**Verdict: BOUNDED NEGATIVE under the current BindingDB-Ki double-cold
protocol.** The tested protein-conditioned interaction mechanisms are closed as
not estimable / negative; this is not a claim of biological impossibility.

## Evidence chain

1. Measurement reliability
   - `stageU_mmp_interaction/U0_RELIABILITY.json` (T0-equivalent): aggregation
     rule verified; same-panel delta reference 0.858 pK²; not an MSE floor.
   - `stageV_core_mmp/PAIR_LEVEL_NOISE_AUDIT.json`: direct repeated
     same-panel MMP deltas, 88 pairs (40 zero-range duplicates);
     disagreeing-only variance 0.303 [0.200, 0.427].
2. Identifiability / support
   - `stageV_core_mmp/V0_V1_RESULT.json`: primary internal repeated-key
     surface 32 rows / 4 components (<100); internal rich keys = 0.
   - `stageV_core_mmp/METAVAL_STRUCTURE_CENSUS.json`: development-validation
     split has 7,209 same-panel MMP observations but 0 exact keys shared with
     meta_train.
3. Interaction variance
   - Preregistered noise: `theta = -0.406 [-0.704, -0.073]`.
   - Pair-level conservative noise: cross-component
     `theta = +0.391 [-0.327, +0.368]` (unresolved).
   - `V1_SYNTHETIC_CALIBRATION.json`: if all excess above the pair-level
     references were signal, latent interaction sd ≈ 0.39–0.53 pK; not
     detectable above the noise envelope.
4. Predecessor mechanism closures
   - Stage S whole-molecule global FiLM: failed controls.
   - Stage P centered correct-vs-wrong objective: failed alignment.
   - Stage T core-blind MMP pooled discriminator: rejected; global closure
     withdrawn after the core-mismatch forensic correction.
   - Stage U core-inclusive key: U0 degree-concentration fail.
5. Remaining lanes
   - `stageV_core_mmp/REMAINING_LANES_AUDIT.json`: MSA/coevolution blocked on
     an absent governed UniRef snapshot; GO annotations falsified (Stage P0);
     pocket priors non-sequence-only and level-rejected (Stage H0); Davis/KIBA
     promotion-gated and not run; meta_test sealed with 0 evaluations; looser
     MMP classes unregistered and screen-only by rule.

## Controls and governance

No neural model was trained for this decision. Stage V's preregistration
(SHA-256 `c567f660…5844d4`) was frozen before its statistics and remains
unchanged; no threshold was moved. `model/` and production `scripts/` were not
modified. The sealed confirmation split was never mounted (0 evaluations).

## Verification

- Stage U + Stage V research suites: **58 passed** (`RUN_SLOW=1`).
- Maintained suite `python main.py verify tests`: **310 passed / 6 skipped**.
- Git commit at execution: `5bb3736`.

## What would reopen the question

(a) a governed MSA/coevolution snapshot, (b) a corpus where a complete
core/context-matched transformation recurs across many protein components, or
(c) a governance change authorizing external datasets before BindingDB
promotion. Each requires a new preregistered stage and may not move the frozen
Stage U/V thresholds.
