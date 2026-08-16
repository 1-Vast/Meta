# UBSE-A1-v2 task-aligned feasible actions

Date: 2026-07-30

Status: binding queue clarification; no coordinate, CCD, parse, event,
student, affinity, confirmation, or sealed-data authorization

Authority: `task.md`

Current decision:
`NO_DECISION_A1V2_STAGE1_EXTRACTOR_TOPOLOGY_PENDING`

## 1. Binding recommendation

The current obstacle must be separated into three claims:

1. deposited complexes may add real ligand-conditioned 3D event information;
2. fixed-margin tetrads may make that information identifiable beyond
   residue and FG marginals;
3. a separately gated pose ensemble may add pair-conditioned geometry that
   is available to a strict dual-cold deployment student.

These claims must be tested in order. More model capacity, a flexible kernel,
or a different fusion block cannot replace any failed information gate.

The parent Stage-1 contract remains immutable. A binding mask/statistic
amendment closes the deployment-mask and primary-likelihood ambiguities. The
accepted innovation candidate remains **Fixed-Margin Tetrad Interaction
Distillation**:

```text
PLIP-only typed event measurement
-> event-blind candidate tetrads
-> no-P0A fixed-margin coupling
-> deployment-side coupling student
-> B0 + theta^T z_int
```

The tetrad statistic does not replace `theta`. It is Stage-2 evidence that a
teacher contains pair-specific coupling. `theta` is the later 224-dimensional
Stage-4A affinity readout.

## 2. Ordered action queue

| Gate | Action | Verification | Authorization state |
| --- | --- | --- | --- |
| P0 | Prepare the hash-first coordinate/CCD acquisition packet | Bind exact manifests, URLs, A1-R-only initial scope, hosts, byte limits, no-redirect streaming, write-once paths, recovery rules, and an explicit `PREPARED_NOT_AUTHORIZED` state | Allowed to prepare only |
| G1 | Acquire 421 A1-R coordinate bodies and the CCD snapshot | Hash-bound complete transfer, no overwrite, no unexpected body, and no parsing | Requires explicit coordinate/CCD GET authorization |
| G2 | Parse the acquired bodies and execute PLIP/ProLIF | Complete nullable legal-cell ledger; PLIP remains the sole teacher and ProLIF remains a non-voting audit | Requires a later parse/event authorization |
| G3 | Decide A1-R reliability | Frozen cross-deposition and cross-implementation thresholds, component bootstrap, and two-reviewer water/heuristic-charge audit | Requires A1-R event-read authorization |
| G4 | Close A1-S topology and power | Frozen members, event-blind legal-cycle rank, role-local components, five component-disjoint folds, and A1-S-TRAIN-only MDE calibration before development or A1-C | Requires separate A1-S authorization |
| G5 | Test Stage-2 coupling | No-P0A model beats exact null `Omega=0` and all registered destructions on held-component conditional tetrad deviance | Requires Stage-2 execution authorization |
| G6 | Test the deployable student | Coupling survives full-minus-null, cycle-direction, wrong-target, wrong-ligand, FG/position, and direct-baseline controls | Requires student authorization |
| G6B | Test a frozen pose-ensemble incremental arm | License, training membership, pocket/template closure, deployment availability, fixed affinity-blind weighting, interaction recovery, and pose-destruction gates all pass | Deferred separate Stage-3 amendment |
| G7 | Run minimal strict dual-cold DTA | `B0 + theta^T z_int` materially beats exact null `theta=0` without a raw target/ligand bypass | Requires Stage-4A and affinity authorization |

Each gate stops on failure. A downstream model may not rescue a failed
upstream information, reliability, topology, power, or deployment gate.

## 3. Accepted recommendations

- Keep no-P0A fixed-margin tetrads as the primary coupling hypothesis.
- Treat P0A only as a later incremental proposal after `M1 > M0`.
- Keep candidate tetrads event-blind and informative `1,1,0,0` or
  `0,0,1,1` patterns post-extraction only.
- Separate placement diagnostics from typed event-chemistry diagnostics.
  Actual pair distance cannot be a covariate in the placement null.
- Keep Level-0/Level-1 aggregates as diagnostics; they cannot substitute for
  the seven typed-channel decision.
- Use OT only as a fixed-margin solver. Its objective is not coupling
  evidence.
- Admit a pose ensemble only as a separately frozen Stage-3 incremental arm.
  Its value must be interaction recovery, not docking score or ligand RMSD.
- If public data fail realized rank or power, stop. The only route that
  clearly adds independent biological outcomes is the separately authorized
  prospective complete-block acquisition design.

## 4. Rejected or deferred recommendations

- PLIP/ProLIF consensus, union, voting, or backfill labels are rejected.
- Frozen A1-R/A1-S/A1-C members cannot be reselected using events,
  informative tetrads, log-determinants, or model outputs.
- The route is not called counterfactual because no intervention is observed.
- A single docking pose is not treated as true pair geometry.
- Increasing rank, attention, backbone size, epochs, ensembles, VIB, MMD,
  global OT, or kernel capacity is not an information remedy.
- `K in {1,2,4}` and uncertainty models require separate preregistration and
  cannot be selected on A1-C.
- The flexible kernel remains conditional Stage-4B and cannot execute unless
  Stages 1-3 and Stage-4A pass.

## 5. Current success criteria

Before any biological execution:

- the frozen contract SHA-256 remains
  `e2e66f1143c93ce886af86a976abe4dcb28f5677ab332451998b9d301e638c92`;
- `contract_valid=true` and `runtime_ready=true`;
- all 15 execution authorizations remain `false`;
- a standard `__pycache__` entry is excepted only when its same-name source is
  `RECORD`-listed and hash-valid; every other hashed payload remains checked;
- synthetic tetrad tests preserve exact row/column cancellation and reject
  non-evaluable or unknown cells.

After the required authorizations, A1-R must satisfy the frozen reliability
thresholds rather than merely produce many event rows. Coverage counts are
not power. The primary future evidence is now the unique component-held-out
conditional tetrad deviance gain against `Omega=0`, with A1-S-TRAIN-only MDE
calibration.

## 6. Estimate

| Work package | Expected elapsed time |
| --- | ---: |
| Acquisition authorization packet and independent review | 0.5-1 day |
| Acquisition runner and recovery tests after authorization | 0.5-1 day |
| A1-R acquisition, parse, and extraction | 0.5-1.5 days |
| Reliability, manual audit, and power decision | 0.5-1.5 days |
| A1-S extraction if A1-R passes | 1-3 days |
| Stage-2 coupling prototype if Stage 1 passes | 1-2 days |
| One-seed Stage-3 CUDA pilot if Stage 2 passes | 1-2 days |

These are engineering ranges, not a forecast that a scientific gate will
pass.

## 7. Binding next action

The acquisition packet is prepared at
`reports/active/ubse_a1v2_coordinate_acquisition_authorization.v1.json`,
SHA-256
`0b650526ba64160355e1da92d18219239bffbc7f44800b804f7744c26d37e7a2`.
It remains `PREPARED_NOT_AUTHORIZED`; no runner exists. The next action
requires separate user authorization for A1-R coordinate GET and CCD GET.
Parsing, event extraction, model training, affinity, confirmation, and sealed
access remain separately locked.
