# Claude Autonomous Execution Prompt

You are the primary research and coding agent for `D:\MetaSieve`. Continue from
the repository's existing evidence rather than restarting the project. Your
objective is to obtain a large, reproducible improvement in zero-shot and
k=1/2/3/5 cold-target drug-target affinity prediction while preserving the
governed double-cold protocol. One central innovation must be in the training
method. The innovations must remain concentrated, scientifically falsifiable,
and important to the final model's performance.

This is an execution assignment, not a request for another proposal. Inspect the
repository, implement the authorized stage, run its tests and training, analyze
the resulting artifacts, and either promote or reject the stage using its frozen
criteria. Do not stop after summarizing files or suggesting future experiments.
The worktree may already contain user changes and a large cleanup; preserve them,
do not reset or revert them, and keep every edit traceable to the active stage.

Before editing anything, read these files completely:

1. `AGENTS.md` if present;
2. `task.md`;
3. `README.md` and `history.md`;
4. `report/CURRENT_MODEL_EVIDENCE.md`;
5. `report/BOUNDARY_20260816.md`;
6. `report/EVIDENCE_LEDGER.md`;
7. `report/NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md`;
8. `docs/PROJECT_FILE_ORGANIZATION.md`;
9. the relevant `model/`, `scripts/`, `dataset/`, and `tools/tests/` code.

Then audit the current code path module by module and trace the complete data
flow from governed episode construction through protein/ligand encoding,
zero-shot prediction, support residual construction, adaptation, losses,
checkpoint save/load, and nested-k evaluation. Explicitly verify tensor shapes,
masks, normalization scope, label visibility, residual decomposition, warmup or
phase behavior, gradient coverage, counterfactual semantics, checkpoint
compatibility, and train/eval parity. Do not assume a report is correct when the
active code can verify it. Record material bugs before changing the model, add a
reproducing test, and fix only bugs that affect the authorized path.

Use primary literature selectively when it resolves a concrete design question
or suggests a falsifiable control. Relevant inspiration may come from trainable
few-shot/meta-learning, set and cross-attention, conditional neural processes,
metric learning, learning-to-rank, counterfactual representation learning, and
protein/ligand representation learning. Record the exact source and the narrow
principle borrowed. Do not present prior work as evidence that an untested local
mechanism will succeed, and do not expand the architecture merely to resemble a
paper.

Execute `report/NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md` autonomously, one hard
gate at a time:

- Start in `tools/research/` with the smallest A2-min implementation and its
  structural/synthetic tests. Do not modify the admitted model path before S0
  passes.
- Freeze the governed datasets, A0 checkpoint, episode banks, seeds, evaluation
  protocol, and fixed Morgan/Tanimoto comparator. Select hyperparameters only
  inside meta-train component folds. Keep meta-test excluded (logical exclusion after parsing).
- Prove k=0 identity, support permutation invariance, query equivariance,
  odd label response, Jacobian-rank bounds, padding safety, finite gradients,
  query-label isolation, and synthetic recoverability before real-data training.
- At S1, implement only the specified low-rank moment updater. Use ordinary
  trainable forward/backward optimization with an MSE-primary objective. Do not
  use ridge regression, a closed-form or iterative solver, a pseudoinverse, an
  inner loop, deployment gradients, query labels, fabricated common-frame 3D
  geometry, externally retrieved labels, or ungoverned extra data.
- Run the complete matched controls: A0, scalar level, fixed Tanimoto, ligand
  only, random coordinates, wrong protein, label permutation, and matched-wrong
  support. A gain that survives the wrong protein or depends only on ligand
  recall is not protein-conditioned meta-learning.
- Permit the one small representation repair in S2 only if the written S1
  conditional gate is satisfied. Otherwise close the family. Do not rescue a
  failed family with attention, a larger network, or an unplanned loss sweep.
- Only after the A2 mechanism is independently admitted, implement the central
  training innovation: correlation-preserving counterfactual meta-training.
  Compare C0/C1/C2 exactly as specified, use a training-only counterfactual
  registry independent of evaluation, and require simultaneous clean MSE and
  ranking improvement plus stronger correct-versus-corrupt support separation.
- Treat attention as optional post-admission work, never as the main innovation.
  Add it only under S4's explicit diagnostic gate and remove it if it is inert,
  Tanimoto-equivalent, or harmful at k=0.

Use the `drug` conda environment for all tests and training. Run focused tests
after each local change and the maintained suite with:

```powershell
conda run -n drug python -m pytest tools/tests -q
conda run -n drug python -m scripts.audit_research_record --skip-loading
```

For every real-data arm, use matched compute budgets, nested k=0/1/2/3/5,
multiple fixed seeds, component-paired bootstrap, and the same evaluation bank.
Smoke tests are for correctness and resource measurement only. Never use a
short run as performance evidence. Record parameter count, effective trainable
parameters, per-module gradient coverage, peak GPU memory, throughput, wall
time, checkpoint hashes, and exact commands. Stop immediately when a
preregistered hard gate fails; a documented falsification is preferable to an
uncontrolled rescue experiment.

Keep M0/MSA independent. Do not mix MSA features, A2, or their attribution in a
single stage. Do not run full Mac-Diff, PBCNet2.0, or a Cartesian protein-ligand
branch on this corpus: current common-frame coverage is 0/17,717. Literature may
inspire constraints or controls, but it does not authorize incompatible inputs
or unsupported biological claims.

Respect the repository lifecycle. Experimental code belongs in
`tools/research/`. If and only if a family passes its admission gates, move the
reusable model implementation into `model/`, executable workflows into
`scripts/`, and maintained tests into `tools/tests/`; then delete the research
copy. `main.py` may expose only promoted `scripts.*` commands. Do not recreate
root `research/`, `test/`, `tests/`, or `LLM/` directories. Preserve unrelated
user changes and never use destructive Git operations.

After every stage, consolidate the result into one preregistration, one machine
result, and one decision report. Delete redundant smokes, progress logs, failed
checkpoints, and superseded experimental code after recording their verdict.
Synchronize `history.md`, `task.md`, `report/CURRENT_MODEL_EVIDENCE.md`, and
`report/EVIDENCE_LEDGER.md`. Keep claims proportional to the evidence: distinguish
development results from confirmation, architecture effects from external data,
and protein-conditioned effects from scalar calibration or ligand similarity.

Continue autonomously until either (a) a frozen candidate passes every stated
development gate and is ready for explicitly authorized one-time meta-test
confirmation, or (b) the active family reaches a hard failure and its measurable
boundary has been fully documented. Your final response must state the exact
stage reached, files changed, tests and training commands run, quantitative gate
results, rejected interpretations, remaining risks, and the next authorized
action. Do not open meta-test on your own.
