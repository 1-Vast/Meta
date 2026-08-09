# Research boundary

## Implemented current path

`s7_l2b_r0r/` contains the completed S7/L2B Phase 0 integrity workflow, the
Phase 1 B4/B5 development implementation, and the Phase 2A audit-only
attribution stage (`pa0`–`pa5`).

Phase 2A trained nothing. It verified every artifact by SHA-256, censused the
corpus for within-construct multi-ligand structure, measured teacher ligand
conditionality against a replicate noise floor, decomposed the sealed logits
into weighted additive marginals plus a marginal-orthogonal coupling residual,
scored a degree-preserving rewiring null, and audited label semantics against a
locally built dense-distance teacher.

```text
TERMINAL VERDICT   LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING
NEXT ACTION        preregister one ligand-conditioned residue residual head
```

The labels are ligand-conditioned at the residue level (ΔJ `+0.258`
[LCB `+0.234`] over a measured replicate floor, chemistry association
ρ `+0.322`). B5's residue marginal is not — a wrong ligand retains 89% of it.
B5's ligand dependence sits entirely in its coupling term, which is real but
below the preregistered margin over both the rewiring null and the wrong-ligand
arm. The teacher's own edge coupling is likewise not identifiable.

## Phase 2B — implemented, run to its precondition, stopped

`p2b_residue_residual.py` and `p2b_run.py` implement
`PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md` (`5e6688f6…`, committed `b9753db`
**before** the code was written). That registration supersedes `ae6d1a01…`,
which was never executed and is kept byte-identical; its eleven defects are
itemised in `PHASE2B_DESIGN_AUDIT.md`.

```text
TERMINAL VERDICT  PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED
```

Every contract and numerical precondition passed. The stage stopped at its own
synthetic trainability control (`AP_bidir 0.3577` against a registered `>=
0.50`), so the real-label run was not executed and gates `R1`–`R6` were not
scored. No biological conclusion follows. Diagnostics put the shortfall in the
optimization budget, not the hypothesis class, the evaluation code or
generalisation.

The next stage may repair the **optimization contract only**: budget and
sampler caps from a measured synthetic scaling curve, an acceptance threshold
derived from that curve, a fresh synthetic teacher seed, and architecture,
projection, metric, controls and gates unchanged. It is not registered.

## Superseded design

`PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL.md` (SHA-256 `ae6d1a01…`), never
executed. It registered one head and nothing else:

```text
logit p_r(P, L) = b_r(P) + delta_r(P, L)
```

- `b_r(P)` is the frozen B5 residue-marginal prior; it is not retrained.
- `delta_r` is one low-rank bilinear residual, `K <= 8`, over existing frozen
  ESM2 residue states and existing atom features.
- `delta_r` is projected away from the constant, the pocket prior, and
  ligand-only/global directions, tolerance `1e-8`.
- Supervision is the same-protein ligand differential on symmetric-difference
  residues only, which removes the generic pocket marginal by construction.
- Gates `D1`–`D5`, a replicate-oracle ceiling, and a fail-closed
  module-participation audit are frozen in the document.

Explicitly excluded by that registration: any additional PLM, cross-attention
stack, typed-interaction head, geometry branch, knowledge graph, parallel SSL
module, affinity head, or PU loss. A pair-coupling head is excluded on evidence:
it would optimise a term worth `0.011` while leaving `0.32` of additive AP
unclaimed.

## Promotion rule

Research code may enter `model/` or `scripts/` only after its own preregistered
Gate, immutable inputs and hashes, closure-safe evaluation, shortcut controls,
regression tests, source affinity increment, and sealed transfer. No current
research statistic is admitted to production `z`.

Historical failed implementations are indexed in `history.md` and recoverable
from Git; they are not active research instructions.
