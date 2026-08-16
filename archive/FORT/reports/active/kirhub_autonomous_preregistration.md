# KirHub autonomous candidate round — preregistration

Date: 2026-07-26. The KirHub/OpenBind directions supplied in user attachments are external
instructions (U1–U3) and consume zero autonomous-candidate slots. The completed KirHub WT H0 is a
user-directed prerequisite audit, not an autonomous candidate.

## Frozen evidence entering this round

- KirHub H0: 358 eligible genes, 92 KLIFS family components, 87 Bemis–Murcko components, 46.22%
  5–95% non-saturated cells.
- Every one of the 25 balanced family-fold × scaffold-fold cells has at least 59 target profiles
  with five query ligands.
- Leave-family-out global-ligand Spearman is 0.3924, but median rho-squared is only 0.1858.
- Protein-family-aligned residual delta is +0.2860 [0.2370, 0.3364], with family-label permutation
  p=0.0005.
- Prospective family-macro MDE80 at paired SD 0.10 is +0.0292.

No neural-network training is authorized by H0. Frozen ESM-2 may be used only as an input to a
low-capacity necessity probe.

## Autonomous candidates (3/3 proposed before candidate results)

### A1 — Separable protein–ligand kernel ordinal probe (SPKOP)

Hypothesis: continuous frozen-protein similarity contains transferable information about the
target-specific ligand ranking residual after strict family/scaffold removal. For every outer
family-fold × scaffold-fold, training excludes all rows from the held target families and all
columns from the held ligand scaffolds. Within each training target, only training-ligand
5–95%-non-saturated observations are rank-normalized. A query score is a nonparametric weighted
average over the Cartesian product of the eight nearest frozen-ESM training genes and eight nearest
Morgan/Tanimoto training ligands. No target ID, learned encoder or affinity-pretrained representation
is allowed.

Required arms: ligand-only kernel, true frozen ESM, fixed target-feature shuffle, matched random
protein features and KLIFS-group centroid. Primary metric: target Spearman within the held scaffold
fold, collapsed first to target family. Gate: true-protein gain over ligand-only >=
max(0.03, MDE80), grouped LCB95 >0, and true protein must beat shuffle, random and group-centroid
with grouped LCB95 >0. One fixed seed (1729); no mechanism revision.

### A2 — Cross-fitted low-rank bilinear ordinal ridge (CLBOR)

Hypothesis: if A1 proves that continuous protein geometry is necessary, a rank-8 bilinear map
between train-only PCA projections of frozen ESM and Morgan features can extract a smoother
dual-cold interaction than the local product kernel. The same 25 outer folds, destructive controls
and family-macro gate apply. Rank, ridge strength and PCA dimensions must be selected only inside
training components. This candidate may enter one single-seed run only if A1 passes every
protein-necessity gate. It is substantively different from A1 because it learns a global
protein–ligand interaction operator rather than transporting labels locally.

### A3 — Mutation-direction transport calibration (MDTC)

Hypothesis: a predictor that passes WT family/scaffold dual cold should also preserve the sign of
matched WT→mutant double differences after replacing the WT protein vector by the mutant vector.
This is a perturbation calibration/falsification layer, not another WT predictor. It requires an
independent mutation panel or raw KirHub replicates because the current 34-gene aggregate subset
has MDE80 +0.048 at paired SD 0.10. A3 is registered but cannot enter experiment on current public
labels.

## Candidate and experiment budget

- Autonomous candidates proposed: 3/3.
- Candidates allowed to enter experiment: at most 2; A1 is first.
- A2 is contingent on a complete A1 pass. A3 is blocked before experiment by its independent-unit
  and replicate gate.
- Multi-seed, confirmation and sealed-test runs remain forbidden.

