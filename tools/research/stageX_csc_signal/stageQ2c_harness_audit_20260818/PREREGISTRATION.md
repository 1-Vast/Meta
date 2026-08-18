# Stage Q2c preregistration — corrected successor audit of the planted-signal harness (2026-08-18)

Frozen BEFORE any Q2c computation. SHA-256 recorded in PREREGISTRATION_SHA256.txt.
Baseline verdict carried at freeze time (from stageX0c artifacts, unchanged here):
Pipeline qualification FAILED at Q2; biological conclusion UNRESOLVED;
B1/B2/C/D NOT AUTHORIZED. This stage is synthetic-only; it reads no external
biological labels.

## Q2c-0 harness self-audit (runs first; precedes any representation change)

1. ORACLE ALIGNMENT: table comparing the diag21 P@U probe protocol against the
   Q2_PLANTED.json oracle_protein arm across: model class, input representation,
   train/eval graph, tau/rank/noise, endpoint transform, projection, seeds,
   selection rule. Until aligned, no claim that any oracle arm passes.
2. ANOVA projection is DEMOTED from primary gate to diagnostic. Primary
   metrics are the explicit interaction-head Spearman / dead-zone sign
   accuracy / gap vs ligand_only with unchanged thresholds 0.30 / 0.70 / 0.05
   (median over 3 seeds).
3. Projection negative controls N1-N7 (each must behave as stated, or the
   projection is removed from diagnostics entirely):
   N1 tau*=0: projected interaction correlation near-random
      (|pearson| < 0.15 AND permutation p > 0.05, 200 permutations, stable seed).
   N2 no_interaction_head: projected interaction correlation must fail
      (same bound).
   N3 ligand_only: projected interaction correlation must fail (same bound).
   N4 independent random bipartite graph with planted main effects only:
      projection null must satisfy the same bound.
   N5 the projection operator is fit on TRAIN cells only; eval labels never
      enter (asserted by construction and by code inspection, recorded).
   N6 the explicit interaction-head metric is computed without any
      projection (dead-zone sign accuracy on raw head output).
   N7 pre-projection and post-projection results both reported.
4. SIGMOID/ENDPOINT CONSISTENCY: rank-correlation between the latent
   interaction I and the interaction implied by the observable % activity
   (logit transform) on determinate eval cells, reported as an endpoint
   distortion budget.
5. SCALE SEPARATION: latent z-scale and observed % activity scale kept
   distinct in every artifact; truth saved separately for mu / pmain / lmain
   / I with SHA-256.
6. MINIMAL LINEAR NO-CENSORING Q2: linear interaction head, identity link,
   no censoring, same graph/tau/rank/seeds; recovery reported.
7. SIGMOID + CENSORING sensitivity branch: the full Q2 configuration
   (values already in stageX0c artifacts; recomputed here for alignment).

## Q2c-1 representation x learner matrix (after Q2c-0)

Grid: representations {one_hot_pocket, pair_centered_local_esm, klifs_pocket,
oracle_PU, random_protein, shuffled_protein} x learners {linear_interaction_probe,
mlp_head}; fixed (tau*=1.0, rank 4, dense, 3 seeds, same splits/graph/noise);
report dz/sp/gap per cell. Interpretation rules (frozen):
- oracle passes with a linear probe but the end-to-end model fails ->
  optimization / loss routing problem;
- oracle passes but local ESM fails -> representation capability gap;
- all probes fail -> harness definition / graph power / truth generation;
- all pass -> proceed to Q2c-2 with the best admissible representation.

## Q2c-2 frozen gate rerun (pair_centered_local_esm)

Same training graph, tau*=1.0, rank 4, same seeds, same budget/restart/
selection protocol, same negative controls, same gate 0.30 / 0.70 / 0.05.
No post-hoc threshold change. Q2c passes only if the gate passes AND every
negative control fails.

## Q3b pairability audit (parallel, no biological training)

exact WT-mutant matching, construct background, substrate compatibility,
ATP protocol, single-mutant vs fusion, saturated cells, responsive-window,
duplicate agreement, parent count, scaffold count, effective sample size.
License CC BY-NC-ND 4.0: code + value-free summaries only; no derivative
value matrices in Git.

## B1 authorization rule

B1 is authorized only after Q2c-2 passes. B1 compares arms: ligand-only,
target-ID, family-ID, nearest-pocket, global protein, local ESM, correct
protein, shuffled protein, family-preserving shuffle, random protein,
matched-parent wrong-variant. The correct arm must improve its own absolute
prediction, not merely make wrong arms worse.

## Governance

- SHA-256-seeded random everywhere; no Python hash().
- Artifacts carry schema / prereg SHA / input SHA / code commit.
- commands.jsonl appended for every run.
- Raw restricted data never committed; downloaded panels stay gitignored.
- stageX0c frozen artifacts and stageX0c/q2.py are read-only inputs to this
  stage; nothing in stageX0c is modified (doc-level corrections in
  REPORT.md per the independent review are allowed and logged).
- Separate commits for implementation and documentation.
