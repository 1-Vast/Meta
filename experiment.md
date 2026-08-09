# Active experiment protocol

## R-VERIFY external SSL evidence recovery

Before S5 execution, recover the exact source commits, immutable manifests,
checkpoints, per-unit predictions, control maps, bootstrap units and label-read
audit behind the supplied S5-S9 report. Missing evidence keeps those numbers at
`EXTERNAL_CLAIM_NOT_REPRODUCED`; it does not authorize a substitute run or a
new ensemble/gauge model. The existing S5 protocol below remains the only
registered experiment if recovery fails.

## Conditional G0 metadata census

`research/ssl_gauge_fixed/PREREG_G0_DATA_FEASIBILITY.md` may audit versioned
source metadata, licences, U1/U2/U3 coverage, apo/holo mappings, closure graphs
and compute feasibility. Numeric affinity values and A0/A1/A2 training remain
forbidden. G0 does not replace the active S5 experiment.

## S5 structural-mechanism observability

S5 is label-free with respect to affinity.  The 14,906 historical holo records
are development/training material.  The previously scored 1,118 S4 complexes
are development-exposed.  A new score-blind RCSB block, disjoint from both
registries, is required for final confirmation.

### Inputs

- frozen ligand GINE atom states;
- frozen protein ESM residue-slot states;
- explicit ligand atom chemistry and slot residue composition;
- frozen P1B contact probability;
- frozen P1B five-bin distance distribution.

### Target

A deterministic, mapped structural pseudo-teacher with pair-local contributions
for prespecified channels.  It is not called binding free energy and does not
observe waters, metals, complete protonation, solvation or entropy.

### Controls

- train mean;
- ligand-only local chemistry;
- historical pooled ESM+ECFP baseline;
- local chemistry without geometry;
- frozen P1B geometry;
- full frozen P1B local states;
- exact-coordinate oracle;
- capacity-matched random features;
- score-blind `<40%` deranged protein;
- within-complex pair-geometry shuffle.

### Statistics

Protein-homology/scaffold closure components are the inference units.  Pair rows
are never treated as IID.  At least one prespecified core channel must achieve
all four component-bootstrap criteria:

```text
R2(correct vs mean) >= 0.02, LCB95 > 0
R2(correct)-R2(ligand-only) >= 0.02, LCB95 > 0
R2(correct)-R2(deranged) >= 0.02, LCB95 > 0
R2(correct)-R2(pair-shuffle) >= 0.02, LCB95 > 0
```

### Conditional training

If the mapping, teacher, slot ceiling, observability and synthetic controls
pass, train only a 1–5M parameter pair-local head in `drug` on the RTX 4060.
ESM, GINE, ProteinEncoder and P1B remain frozen.  Use AMP, batch 2–4, gradient
accumulation and chunked atom-slot evaluation.  No affinity, identity or
wrong-protein training objective is allowed.

### Terminal outcomes

```text
S5_DATA_OR_MAPPING_CONTRACT_FAIL_CLOSED
S5_TEACHER_OR_SLOT_CONTRACT_DEFECT
S5_HEAD_OR_OPTIMIZATION_NOT_IDENTIFIED
P1B_PAIR_LOCAL_MECHANISM_NOT_OBSERVABLE
P1B_PAIR_LOCAL_STRUCTURAL_MECHANISM_OBSERVED
```
