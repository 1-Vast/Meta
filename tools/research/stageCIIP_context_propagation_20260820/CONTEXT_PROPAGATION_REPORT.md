# CIIP Contextual-Propagation Audit Report

Date: 2026-08-20  
Preregistration SHA-256: `cdd6e0a8fa1612efd6a7a9c23d93f6a38b699acf9eeb968a2acdfc700e9b5b19`

## Result

The read-only residue-level audit completed on all 49 ESM-covered verified
single-mutation pairs. Mutation erasure passed exactly: WT and variant inputs
became identical for 49/49 pairs after replacement of the verified residue
with `X`, and the largest absolute residue-embedding difference was `0.0`
(required maximum `<= 1e-5`). The original cache therefore reflects the
sequence edit rather than an input-alignment artifact.

| measurement | mean L2 delta norm |
|---|---:|
| mutation site | 4.0111 |
| radius-6 local window | 1.1605 |
| non-site context | 0.05749 |
| full sequence | 0.07392 |

The curve is highest at the mutation site and decays with sequence distance;
it remains nonzero outside the radius-6 window. Consequently, a distant ESM
window can carry distributed contextual information about the mutation. This
makes the prior random-window arm an explicit-coordinate control, not a proof
of complete mutation-information independence.

## What This Does Not Establish

No predictor was fitted in this stage, and no labels were used. Therefore
context-only centered-response prediction and the predictive increment of a
site-specific residual are **not evaluated**. Either would require a separate
preregistration with train-parent-only feature fitting and a new control matrix.
The audit neither reverses the CIIP-1A formal verdict
`ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED` nor authorizes CIIP-1B, BindingDB bridge
work, a deployable mutation-coordinate-free potential, or production changes.

## Artifacts

- `CONTEXT_PROPAGATION_RESULT.json` and `RESULT.json`: machine-readable result.
- `context_distance_curves.npz`: per-pair mean delta norm by distance.
- `commands.jsonl`: exact environment commands.
- `FAILURES.md`: no audit failure occurred.
