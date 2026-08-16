# Stage R11 preregistration: shape-first training on the incumbent grammar trunk

Frozen before any R11 run. Population, bank, seeds, budget, selection and the
wrong-protein contract are identical to R7-R10. `meta_test` remains sealed.
This experiment tests the **shape parameterization** question on the
architecture with the best-known calibration, with zero architecture change.

## Hypothesis

R9/R10 established that the shape-first training method is the project's
first real within-target shape source, and that the factorized
relative-transport trunk's anchor-mean-of-delta readout caps the global CI
below the incumbent's (best 0.562 vs 0.580) while its routed level
calibration matches the incumbent only at the full budget. The incumbent
trunk itself (`similarity_only` grammar, per-atom interaction head, the
project's best k=0 calibration 1.236 and CI 0.580) has never been trained
with the shape-first method. Hypothesis: applying the shape-first objectives
(ranking-primary, cliff weight 1.0 per R9, relative-free, routed,
counterfactual) to the incumbent trunk lifts its within-target ordering —
CI and cliff ordering — without losing its calibration, because the
per-atom interaction head is a strictly more expressive shape readout than
the anchor-mean-of-delta form.

## The experiment

`scripts/train_grammar_shape.py` — the exact `SimilarityGrammarModel` (A0's
architecture, `use_learned_key=False`, same Tanimoto transport), trained
with: level term (mean squared error) routed to ligand/protein heads;
shape objectives = pairwise RankNet (cliff weight 1.0) + shape variance
(1.0) on the interaction-branch path; label-free identifiability pin on the
interaction mean (0.3); wrong-protein shape/level contrasts and
permuted-label (k>=2) / wrong-ligand (k=1) binding contrasts (0.25). One
backward pass. 1200 steps, 3 episodes/step, lr 6e-4 cosine, seeds
20260815/16/17. Executed by `scripts/run_stage.py` (smoke first).

## Arms

- **A0** frozen R3R4 incumbent checkpoints (identical architecture, ordinary
  training) — the same-architecture control;
- **G1** the R11 arm (grammar trunk + shape-first training), 3 seeds.

## Gates

- **H1** G1 k=0 MSE does not regress beyond A0's 2.149 (the training method
  must not cost calibration);
- **H2** G1 k=0 CI improves over A0's 0.580;
- **H3** G1 k=0 cliff sign accuracy improves over A0's 0.512;
- **H4** the D1-vs-A0 contrasts hold in all three seeds (direction) with a
  component bootstrap lower bound above zero for H2's CI change;
- **H5** interaction-cut (`ligand_only`) is clearly worse than full.

## Decision rule

H1-H3 must all pass for the incumbent-trunk + shape-first combination to
advance to the standing Z-gate re-evaluation (three seeds, k=0/1/2/3/5,
full control battery). A failure is recorded per gate; no gate moves. The
next single variable after a failure is selected from the pair audit of the
G1 checkpoints, not from a sweep.
