# MetaSieve Project Summary

This file is the current status entry point. Numerical claims are controlled by the corresponding
`RESULT.json` and manifest files; mathematical and model claims are controlled by
`theory/CURRENT_THEORY/`.

## Objective

Learn a shared deep meta-learning model for few-shot affinity prediction on proteins absent from
training. Deployment may use only the protein sequence or legal structure representation, ligand
molecular information, context, and a disjoint support set. Target identifiers, query labels, and
target-specific parameter memory are prohibited.

## Core Innovation

The candidate core is the **Quotient-Preserving Section Meta-Potential (QPSMP)** in
`model/qpsmp_meta.py`. It is a trainable neural scalar potential with ligand-conditioned protein
localization and a permutation-invariant, centered, zero-preserving support adapter. Query loss
trains the localizer, crossed interaction map, scalar head, section basis, and adapter. The analytic
centered ridge in `model/qpsmp.py` is a comparator and section diagnostic, not the innovation.

## Current Status

`QPSMP_TRAINABILITY_INTERFACE_DEFINED=true`

The module interface and focused invariants are implemented and pass unit tests. Protein-specific
Cold Target admission is **not** established. The current frozen PLM-slot/T-BASIS diagnostics do
not authorize G2, G3, biological interpretation, or V1 integration.

The repaired endpoint contract improved over its matched level baseline in an already-consumed k=5
development diagnostic. A stricter shared-checkpoint nested-k evaluation over k={1,2,3,5} produced
positive point estimates against level at every k, but no component-bootstrap lower bound was above
zero, and SAR/protein-specific controls were not stable. The recipe has therefore not passed
development promotion. No G2, G3, biological, confirmation, performance-guarantee, or integration
claim is authorized. The governed manifest records only six independent homology components, and
the repository currently has no authorized untouched confirmation cohort.

## Authority Order

1. `theory/CURRENT_THEORY/` for current mathematics and model contract.
2. `task.md` for the active falsifiable work contract.
3. `report/**/RESULT.json` plus manifests for numerical evidence.
4. `report/EVIDENCE_LEDGER.md` for navigational summaries.
5. `history.md` and `archive/` for provenance only.

`theory/FINAL_FROZEN_THEORY/`, `PROJECT_HANDOFFS/`, and `archive/FORT/` are preserved historical
material and cannot override the current theory.
