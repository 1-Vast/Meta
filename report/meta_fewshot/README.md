# Cold-Target Few-Shot Experiment Lineage

This is the single location for zero-shot and few-shot DTA experiment outputs.
Each leaf directory is immutable provenance; `RESULT.json` is numerical
authority.

## Current retained incumbent (comparator for all R5+ gates)

```text
stageR3R4_level_shape_20260815/A0_incumbent_seed202608{15,16,17}/
```

`similarity_only` grammar checkpoints, 1200 steps, three seeds. Double-cold
`meta_val`: k=0 MSE 2.149 / CI 0.580 / Spearman 0.223 / calibration 1.236 /
shape 0.913. This is the development incumbent, not an admitted model.

## The R5-R8 cycle (2026-08-16): unified summary

Population: governed double-cold `bindingdb_ki_double_cold_v1`
(meta_train 5,643 cells / 346 targets; meta_val 41 targets / 19 components;
meta_test 22 targets / 10 components, physically sealed, never opened).

- **R5** (`stageR5_reltransport_20260816/`): repaired the experimental
  contract (same-split wrong-protein donors with meta_train-only whitening,
  aggregated gradient cosines, physical meta_test seal, full artifact
  records) and passed 23/23 structural + synthetic falsification gates for
  the relative-transport candidate.
- **R6** (`stageR6_reltransport_screening_20260816/`): three screening
  rounds, each eliminating a transport mechanism under its preregistered
  gates: R6a multiplicative saturating gate (deployment-inert, nogate gap
  0.000), R6b additive residual-identity correction (self-cancels by
  construction), R6c attention level readout (calibration needs the full
  budget). All artifacts retained in `R6a/R6b/R6c_archive/`.
- **R7** (`stageR7_reltransport_3seed_20260816/`): three-seed formal run,
  1200 steps. **Admission refused** — A2 k=0 2.420 vs A0 2.149 (-12.6%),
  CI 0.542 vs 0.580; the linear rho gate was again eval-inert while its
  training disturbed calibration (the seventh query-specific channel in the
  project with that signature). Positive: the shape-first training produced
  the project's first real shape gain (0.943 -> 0.895; k=0 cliff sign 0.536
  vs 0.512).
- **R8** (`stageR8_stronger_shape_20260816/`): stronger shape signal (A3
  configuration, shape_variance 1.5, relative 1.0, no gate). Best-ever shape
  0.896 and best-ever k=5 activity-cliff sign 0.768 (A0 0.675), k=0 2.167
  (tie) — but CI 0.535 vs 0.580. Both preregistered advance gates failed:
  **the relative-transport/gate model family is closed for the double-cold
  zero-shot target as a claimed core innovation.**

### Why the three gate designs failed (mechanism, measured)

1. A saturating gate cannot express the optimal per-query residual scaling,
   and it measured deployment-inert (nogate gap 0.000).
2. An additive `r_k + delta_hat - delta_f0` correction self-cancels: the
   trained relative potential converges to the endpoint's own implied
   relative, so the correction degenerates to the plain residual.
3. A linear zero-mean gate stays eval-inert (trained rho ~ 1) while its
   training gradients disturb the routed level's calibration.

### Why shape-first partially worked, and what did not transfer

The within-target pairwise ranking (ActFound/PBCNet-style, cliff-weighted)
plus direct relative supervision on an antisymmetric potential is the first
training method in this project to measurably move the shape term and the
activity-cliff ordering. What did not transfer: the shape gain concentrated
on cliff/high-gap pairs while the global concordance index regressed — the
open question the R9 pair audit answers stratum by stratum.

### What remains transferable to the next model

- the shape-first objective family (ranking-primary, routed, counterfactual)
  as Innovation B material;
- the attention-pooled routed level readout (converges to incumbent
  calibration at the full budget);
- the pair-level audit machinery and the closure discipline: any new
  query-specific channel must first beat `nogate` on the development split
  at 300 steps before any formal run.

## Rejected recent stages

- `stageC_*`: slot-gated endpoint readout.
- `stageD_*`: atom-aware localization.
- `stageF_*`, `stageG_*`: primitive kernel.
- `stageH_*` through `stageK_*`: context-kernel variants.
- `stageL_*`: label-noise smoke.
- R6a/R6b/R6c/R7/R8: relative-transport gate family (closed, see above).

Earlier `dmemt_*`, `cipf_term_*`, `lcipf_elmt_*`, `hypersar_*`, `lirms_*`,
and `main_v1_*` families are historical comparators. `bpsf_v2_research/`
contains outputs moved from the research source tree.

Read `../CURRENT_MODEL_EVIDENCE.md` before starting a new stage.
