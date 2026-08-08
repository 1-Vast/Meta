# E0 Identifiability Research

## Two Separate Estimands

**Claim A** — does the correct protein provide affinity-*location* information
beyond population, ligand-only and protein-sequence-only baselines? Addressed by
`E-AFF-L0`.

**Claim B** — does a non-additive protein-by-ligand affinity *interaction*
exist? Addressed by `E-AFF-X0-B` and then `X1`.

**L0 and X1 are different estimands.** Neither replaces, implies nor authorizes
the other.

The P1C/P1R\*/E-AFF-P0/H0A/H0C evaluations used within-task concordance or
another rank-based metric. Those experiments remain valid and no verdict is
reclassified; their negative conclusions are limited to within-task ranking
information, because within-task concordance is invariant to task-wise
prediction shifts and positive rescaling.

This is the canonical terminal research package for the E0 identifiability
line. It preserves synthetic numerical closure, structural basis evidence and
the source-affinity falsification chain. There is no active model-training
stage and no unresolved synthetic numerical question.

`run_proposal_numerical_closure.py` is the separately registered E0R2
synthetic-only closure of the corrected residual/difference objective. It also
records which claims in the directional-potential proposal are established,
untested, or semantically overstated. A PASS remains research evidence and does
not authorize production integration.

It is not production code. Synthetic and structural stages must not read real
affinity labels; the separately registered E-AFF stages may read governed
ChEMBL37 source Ki/Kd labels but must never read DAVIS or recipient labels.
E0R2 closed the synthetic numerical boundary with the separately registered
verdict `SYNTHETIC_OBJECTIVE_DESIGN_SOLVER_IDENTIFIED`. This does not identify
real affinity energetics or authorize production integration.

`TDIR_P0_PREREGISTRATION.md`, `run_tdir_pilot.py`, `TDIR_P0_RESULT.md` and
`artifacts/tdir_p0_v1/` contain a small structure-only feasibility pilot. Its
verdict is `PILOT_LEARNABILITY_SIGNAL_NOT_OBSERVED`. The executed group-channel
mapper defect and omitted controls are recorded in `postrun_audit.json`; the
used test panel must not be rerun as untouched validation. No T-DIR code is
authorized for promotion.

`T_BASIS_R0_PREREGISTRATION.md`, `run_tbasis_radial.py`,
`T_BASIS_R0_RESULT.md`, and `artifacts/tbasis_r0_v1/` contain the fresh-panel
fixed radial basis study. Its verdict is
`RADIAL_BASIS_PARTNER_RECOVERABILITY_IDENTIFIED`. It authorizes only a new
structure-only angular/many-body basis registration; it does not authorize
affinity or production code.

`EAFF_P0_PREREGISTRATION.md`, `run_eaff_pilot.py`, `audit_eaff_pilot.py`,
`EAFF_P0_RESULT.md`, and `artifacts/eaff_p0_v1/` test one population-shared
288D residual-affinity direction. The result is
`SHARED_DIRECTION_NOT_OBSERVED_H0_DATA_SUPPORTED`.

`EAFF_H0A_PREREGISTRATION.md`, `run_eaff_h0a.py`, `audit_eaff_h0a.py`,
`EAFF_H0A_RESULT.md`, and `artifacts/eaff_h0a_v1/` test task-local radial
headroom on unseen ligands. Strong headroom was observed, but the correct-minus-
deranged contrast was only `+0.00864`, below the frozen `+0.03` requirement.
The verdict is `TASK_LOCAL_RADIAL_HEADROOM_WITHOUT_PARTNER_SPECIFICITY`; H0-B
and RFSA are not authorized.

`EAFF_H0C_PREREGISTRATION.md`, `run_eaff_h0c.py`, `audit_eaff_h0c.py`,
`EAFF_H0C_RESULT.md`, and `artifacts/eaff_h0c_v1_run2/` contain the fixed radial
interaction-residual follow-up. On 54 new scaffold-disjoint tasks, the
support-matched ligand nuisance improved ranking, but the centered radial
interaction residual did not. Verdict:
`FIXED_RADIAL_INTERACTION_RESIDUAL_NOT_OBSERVED`. This closes H0C without
authorizing RFSA or orientation.

`EAFF_X0_PREREGISTRATION.md`, `x0_metadata.sql`, `run_eaff_x0.py`,
`audit_eaff_x0.py`, `EAFF_X0_RESULT.md`, and `artifacts/eaff_x0_v1/` contain the
strictly label-blind crossed source census. Millions of nominal ChEMBL
rectangles reduce to only 36 Ki and 12 Kd dependency components, below the
frozen requirement of 245. Verdict: `STOP_SOURCE_INTERACTION_UNDERDETERMINED`.
X1/X2 are not authorized.

`EAFF_X0_FEAS_REGISTRATION.md`, `run_eaff_x0_feasibility.py` and
`EAFF_X0_FEAS_RESULT.md` audit the X0 estimand rather than the source. The
closure-component universe is `245`, while the frozen requirement is `245` *per
endpoint*; only `202` components carry Ki rows and `72` carry Kd rows. Because a
rectangle needs two proteins inside one document and D1 closure unions all
targets sharing a document, both proteins of every rectangle always lie in one
closure component, so crossing can never create a unit. Recomputing the closure
over the full governed D0 corpus and over shallower task populations raises the
best ceiling only to Ki `97` and Kd `56`, still far under `245`, and the Ki
ceiling falls as the corpus grows because added documents merge components.
Verdict:
`X0_UNIT_REQUIREMENT_UNATTAINABLE_BY_CONSTRUCTION`. This re-registers nothing
and authorizes nothing; it records that the X0 stop is a specification
consequence, not a measurement of ChEMBL crossing.

`EAFF_X0B_PREREGISTRATION.md`, `run_eaff_x0b.py`, `EAFF_X0B_RESULT.md` and
`artifacts/eaff_x0b_v1/` re-register the crossed-census unit as the cell-disjoint
rectangle with a design-effect effective sample size, leaving the frozen effect
size, alpha, power, `245` requirement and `+0.03` affinity margins unchanged. At
`rho = 1` the model reproduces X0's `36` and `12` exactly, so X0 was its
total-correlation corner. Ki packs `11,168` units in `36` clusters and Kd `1,041`
in `12`. Verdict: `X0B_CONDITIONAL_DESIGN_SUPPORTED_KI|KD`, conditional on
intra-cluster correlation at most `0.0915` (Ki) and `0.0164` (Kd). A separately
registered X1 must estimate `rho`, compare its upper confidence bound to that
threshold, and abstain if exceeded.

`EAFF_R0_REGISTRATION.md`, `run_eaff_r0_readout.py`, `EAFF_R0_RESULT.md` and
`artifacts/eaff_r0_v1/` diagnose the readout behind every affinity result rather
than the biology. Within-task concordance is **exactly** invariant to per-task
affinity shift and rescale (maximum deviation `0.0` on all four transforms), and
a simulated predictor that knows a task's affinity level perfectly scores exactly
`0.5000` even when the level holds `98.5%` of variance. H0C additionally removed
that channel upstream and shared the correct protein's support-fitted nuisance
with the deranged arm. Verdict:
`READOUT_BLIND_TO_TASK_LEVEL_AFFINITY_LOCATION|PERFECT_LEVEL_PREDICTOR_SCORES_CHANCE_AT_EVERY_VARIANCE_SHARE`.
No affinity label was read and no past result is overturned; their scope narrows
to within-task ranking.

`THEORY_BIOLOGY_INTEGRATION.md` is an unvalidated design proposal joining the
frozen operator to the validated geometry: stochastically ordered anchor bands so
the affinity sign is a deployment property instead of an estimated direction,
`kappa` carrying assay nuisance into `beta_0`, a small bounded pair-local
`z_bio`, a small coefficient simplex, and abstention as mass on `beta_0`. It
admits nothing and changes no Gate. Pre-execution amendment A1 recorded that L0
uses the deployment's own frozen `m = 7` rather than a truncated ladder, because
L0 performs no support-based adaptation and the `d_adapt <= k` limit is not
engaged; that limit still binds any future few-shot adapter, which stays frozen.

`EAFF_L0_PREREGISTRATION.md`, `EAFF_L0_DATA_CONTRACT.md`, `l0_contract.py`,
`audit_l0_operator_contract.py`, `run_l0_estimand_check.py`, `run_eaff_l0.py`,
`audit_eaff_l0.py`, `EAFF_L0_RESULT.md` and `artifacts/eaff_l0_*` are the Claim A
Gate. Preconditions passed: the operator/anchor contract froze
(`L0_OPERATOR_AND_ANCHOR_CONTRACT_FROZEN`, all `258` theory files hash-matched)
and the estimand check admitted Ki only (`C1=79`, `C2=218`, `C3=0.630`) while
excluding Kd (`C1=10 < 30`). `sigma_assay = 0.47971 [0.47034, 0.48946]` log
units from `4,261` replicate cells, giving `margin_L0 = 0.23985`.

The run executed on `115` tasks in `115` closure components, then failed closed:
`L0_NOT_RUN_NUMERICAL_PRECONDITION_FAILED`. Gate condition 3 used a
step-containment coverage statistic that is identically `0.0` for every arm on
the fixed `33`-point mesh, so one of three registered conditions did not execute
as specified; ligand-only additionally failed to beat population-only, so no
positive control was established. Claim A remains untested, the panel is
consumed, and a corrected Gate needs a new registration and a fresh panel.

`EAFF_L0R_PREREGISTRATION.md`, `run_eaff_l0r.py`, `EAFF_L0R_RESULT.md` and
`artifacts/eaff_l0r_v1/` are that corrected Gate, on `195` fresh tasks in `78`
closure components (`3,900` rows), excluding all `470` tasks consumed by P0, H0A,
H0C and L0. Both L0 defects were fixed and both fixes worked: mean-interval
coverage is informative (`0.0674`-`0.1664` across arms) and a registered
positive-control precondition was enforced. The control then failed - ligand-only
beat population-only by only `+0.03421` log units, CI `[-0.03304, 0.10793]`,
against a required `0.1 * sigma_assay = 0.04797` - so the Gate stopped with
`L0R_NOT_RUN_POSITIVE_CONTROL_ABSENT` and **no protein contrast is reported as
evidence**. Claim A remains untested. The reusable measurement is that the best
cross-component location signal in governed ChEMBL37 Ki is worth about `7%` of
assay noise and is not separable from zero at `78` components, which bounds any
L0-style Gate on this corpus.

## Evidence Disposition

| Class | Retained here | Meaning |
|---|---|---|
| Synthetic closure | E0S, E0R0, E0R1, E0R2 | objective/design/solver is numerically identified only |
| Structural negative | T-DIR-P0 | sparse PLIP event learnability was not observed on its consumed pilot panel |
| Structural positive | T-BASIS-R0 | fixed 288D radial basis is recoverable and partner-dependent |
| Source-affinity negative/mixed | E-AFF-P0, H0A, H0C | no admitted correct-partner affinity statistic |
| Source-data stop | E-AFF-X0 | ChEMBL crossed interaction is underdetermined at the frozen sensitivity |
| Estimand audit | E-AFF-X0-FEAS | the X0 requirement exceeds its own unit ceiling, so that stop is not a source measurement |
| Design re-registration | E-AFF-X0-B | cell-disjoint units reach 245 conditional on low intra-cluster correlation; X1 gated on rho |
| Readout diagnosis | E-AFF-R0 | the affinity readout is algebraically blind to task-level affinity location |
| Integration design | THEORY_BIOLOGY_INTEGRATION | ordered anchors make the affinity sign a deployment property; unvalidated |
| Claim A Gate | E-AFF-L0 | not run: a registered gate condition used a degenerate statistic; Claim A untested |

Preregistrations, result reports, final artifacts and independent audits are
retained because together they establish what was tested. Python/pytest caches,
failed-import selections and timeout-orphaned partial selections are disposable
runtime residue and are not scientific evidence.

`E0R1_PROPOSAL.md` and the original artifact manifests retain their pre-move
paths verbatim so their registered SHA-256 values remain auditable. Executable
defaults were updated to the paths under this package.

`RELOCATION_MANIFEST.json` binds the post-consolidation research code. It is
separate from the immutable manifests generated by the original experiments.

Focused verification:

```powershell
D:\anaconda\envs\drug\python.exe -m pytest -q research/e0_identifiability/tests
```
