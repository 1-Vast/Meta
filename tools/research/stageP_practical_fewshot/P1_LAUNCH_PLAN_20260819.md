# Stage P1 launch plan (budget discipline, 2026-08-19)

All launches wait until the Q2d-1e ladder and (if needed) the frozen
span-param diagnostic have finished using the GPU. Promotion ladder
per CORE1_FUNNEL_PLAN_20260819.md and P1_BAKEOFF_PREREGISTRATION.md
(single-seed screening on p_val first; promotion -> 3 seeds ->
independent p_test surface). p_test labels are never used for
selection; the checkpoint rule is the frozen p_val monitor.

## Order and estimated cost (RTX 4060 8GB; batch 256 cells)

1. arm 3 ordinary FT: incumbent learned baseline, ALWAYS runs 3 seeds.
   Seed ~10-15 min (train 6000 steps + eval 6376 records). Total ~0.6
   GPU-h. Estimated max memory: < 1.0 GB.
2. arm 4 FOMAML: single-seed screen ~40-55 min (5 inner steps x ~14
   tasks per outer step); promote only if p_val delta vs baselines
   clearly positive. Total screen ~0.9 GPU-h.
3. arm 5 Deep-Sets CNP: single-seed screen ~10-15 min. ~0.3 GPU-h.
4. arm 6 FS-CAP-style: single-seed screen ~8-12 min. ~0.2 GPU-h.
5. arm 7 ActFound-style: single-seed screen ~15-25 min. ~0.4 GPU-h.
6. P2 metrics: p2_arms.py on every promoted arm artifact (CPU, minutes).

Promoted arms rerun with seeds {1,2,3}; p_report.py v2 then computes
paired deltas and target-level bootstrap vs baselines and arm 3.
Expected stopping points: a non-promoted arm stops after its single
seed (recorded in commands.jsonl).

## Screens must be recorded with

- command, seeds, device, estimated and actual GPU time, max memory;
- p_val screen outcome vs the frozen promotion rule (no p_test access);
- promotion decision before any 3-seed or p_test run.
