# Section-Conditioned Mechanism Operator: exploratory preregistration

Status: frozen before any new KCGS outcome value is read. PKIS1 is source;
PKIS2 and Anastassiadis2011 have already been consumed by the preceding
`pkis_mechanism_pilot` and are development-only. This is not E-AFF-X1 and does
not authorize admission to `model/config.py`.

## Question

Can a new target's at-most-four-dimensional affinity section be identified
from `k <= 5` support measurements when the query ligands and target are both
cold, and does correct kinase-pocket side information improve that section?

The decomposition is

```text
y(s,t,l) = context(s) + ligand(l) + protein(t)
           + tau(t) + sum_q u_q(t) phi_q(t,l) + error.
```

`tau` is a location coordinate. There are at most three interaction
coordinates. Hence `d_adapt = 1 + r <= 4 <= k` for the primary `k=5` setting.
The pair basis must be deterministic, bounded after source-fitted scaling,
query-label-free and dependent on ligand chemistry and the aligned KLIFS
pocket. A learned support operator must be permutation invariant.

## Stages

### F0: attainable ceiling on consumed development panels

Fit a source-only low-rank interaction basis on PKIS1, predict its ligand and
protein factors from Morgan/scalar ligand features and aligned KLIFS pocket
features, and estimate `(tau,u)` from support by a positive-ridge solve. This is
only a ceiling/control, not the proposed learnable adapter. Select rank in
`{1,2,3,4}` and ridge values using source-only group/scaffold folds.

Run fixed support sizes `k in {5,20}` and seeds `20260808..20260827` on PKIS2
and Anastassiadis2011. Query rows exclude support rows. No target or ligand ID
is a feature.

F0 is viable only if, at `k=5`, the correct section improves target-macro query
MSE over both support-free and location-only adaptation with a target-bootstrap
95% lower bound above zero on PKIS2, and the point estimate is positive on
Anastassiadis2011. Correct support must beat wrong-target and label-permuted
support. If the source basis cannot meet this, the trainable version is not run.

### F1: learnable section operator

Train episodically on PKIS1 only. A hard-masked three-channel pair encoder uses
the KLIFS hinge, DFG/back-pocket and front/solvent-pocket regions together with
corresponding ligand pharmacophore-centred fingerprints. A DeepSets adapter
maps the unordered support set to `(tau,u_1,u_2,u_3,c)`, where `c` is a coverage
coordinate. The scalar diagnostic decoder is used only to test whether the
biological statistic exists; deployment must still pass through the frozen
`F -> B -> K` operator.

F1 must beat capacity-matched controls: no support, location only, no protein,
deranged protein, wrong-target support, support-label permutation, and a
closed-form ridge adapter. Success requires the PKIS2 F0 margins plus a positive
correct-minus-deranged lower bound and no degradation relative to the strongest
baseline in target-macro MSE. Five fixed training seeds are required.

### F2: untouched external transfer

KCGS is reserved until the complete F1 architecture, hyperparameters, parser,
identity mapping, exclusions, hashes and verdict code are frozen. Before F2,
only workbook names, dimensions, string metadata, compound identities and
structures may be inspected. No numeric assay outcome may be printed,
summarized or used for selection.

F2 requires strict source-scaffold cold ligands, exact-target cold targets,
`>= 20` targets, `>= 100` ligands and `>= 80%` finite cells. Main success:
correct beats support-free, location-only, wrong support and deranged protein
with positive target-cluster bootstrap lower bounds. Family-cold is reported
separately and cannot be substituted for the main stratum.

## Statistics and failure rules

- Primary: target-macro query MSE and its paired target-cluster bootstrap.
- Secondary: MAE, Pearson/Spearman, interaction-residual MSE, IC-index when its
  implementation is independently verified, coverage-risk curve and interval
  coverage after calibration.
- Bootstrap targets, never individual matrix cells. Fixed 10,000 draws for a
  final result; 2,000 may be used during development and must be labelled so.
- All support/query splits are nested within target and disjoint by molecule.
- Source/development/test generic Murcko scaffold overlaps are zero for the
  strict stratum.
- A positive PKIS2 result without Anastassiadis/KCGS transfer is development
  evidence only. A negative external result cannot be hidden by reporting an
  easier overlap stratum.
- No change to `model/`, the theory archive, CSMO, Band, anchors, positive
  ridge, recipient labels or DAVIS labels is authorized by this pilot.

