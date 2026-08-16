# A1 SPKOP decision

Verdict: `A1_PROTEIN_NECESSITY_FAIL_STOP_A2`.

The preregistered A1 probe jointly held out whole KLIFS-family and whole Bemis–Murcko components in
25 outer folds. Training labels were target-local ordinal ranks computed only from training
ligands. The query score used eight frozen-ESM protein neighbours and eight Morgan/Tanimoto ligand
neighbours; no target ID or learned protein encoder was available.

Family-macro Spearman:

| Arm | Spearman | LCB95 | UCB95 |
| --- | ---: | ---: | ---: |
| ligand-only | 0.0815 | 0.0473 | 0.1129 |
| true frozen ESM | 0.0672 | 0.0330 | 0.0990 |
| protein shuffle | 0.0609 | 0.0365 | 0.0847 |
| random protein | 0.0585 | 0.0291 | 0.0862 |
| KLIFS-group centroid | 0.0873 | 0.0491 | 0.1220 |

True ESM minus ligand-only was −0.0143 [−0.0409, +0.0116]. True ESM did not beat shuffled,
random or group-centroid protein controls with a positive grouped lower bound. The family-unit
MDE80 was +0.0292 at paired SD 0.10, below the required +0.03 gain. All five preregistered gates
failed.

Interpretation: the WT matrix contains family-aligned target-specific residual structure, but the
continuous frozen-ESM geometry used here does not transport that structure across unseen families
and unseen ligand scaffolds. Coarse kinase-group membership is at least as useful as the true ESM
neighbourhood. The result does not authorize A2's learned bilinear operator. No mechanism revision,
second seed, confirmation run or sealed label is permitted.

