# CIIP-1A Control-Arm Formal Report

Date: 2026-08-20  
Frozen preregistration SHA-256: `39d02166f69acf235a34d351b649a4cdbf3b828491a0994901bf2378777463f7`  
Result SHA-256: `55b0de124dda9e092c5a6692c16211bdbffb575235aa585e844652c55bb38ffd`

## Formal Verdict

`ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED`

The adjudicator in `adjudicate_controls.py` reads the immutable seven-arm
result and applies the frozen verdict rules without retraining or changing any
threshold. Its machine-readable output is `CONTROL_ADJUDICATION.json`.

## Scope and Limits

This result applies only to the oracle mutation-coordinate, ESM-2 local-window
representation; the Duong-Ly centered percent-inhibition endpoint; the
49-pair ESM-covered subset; the current rank-8 low-capacity potential; and the
existing pair-level (not parent-disjoint) split. It is not a conclusion about
general cold-target DTA, deployable mutation-coordinate-free protein features,
or direct binding affinity (Ki, Kd, pK, or DeltaDeltaG).

## Frozen Rule Evidence

| comparison / requirement | observed result | adjudication |
|---|---:|---|
| Correct nonconstant coverage exceeds random-window coverage | 9/9 vs 9/9 | fail |
| Correct R2 minus random-window R2 | -0.1217; parent bootstrap [-0.4569, 0.0327] | fail |
| Correct R2 minus family-shuffle R2 | +0.0111; parent bootstrap [-0.3508, 0.2421] | fail |
| Correct R2 minus random-protein R2 | +0.0070; parent bootstrap [-0.3387, 0.2322] | fail |
| Correct R2 minus ligand-only R2 | +0.0075; parent bootstrap [-0.3487, 0.2350] | fail |
| Correct sign accuracy | 0.7025 | descriptive only; not sufficient |

The feature audit is positive but does not change the verdict: the correct-site
WT/variant delta norm is 0.5310 on average versus 0.0267 at matched random
windows, and 49/49 pairs are larger at the correct site. Thus ESM detects the
sequence edit at its annotated position, but this control does not show that
the detected information yields transferable ligand-conditioned prediction.

## Authorization

The deployable protein representation remains **not validated**. CIIP-1A PASS,
CIIP-1B, the BindingDB Potential Bridge, production `model/` changes, and
production `scripts/` changes remain **not authorized**. The only next action
authorized by the supplied task is a read-only contextual-propagation audit.
