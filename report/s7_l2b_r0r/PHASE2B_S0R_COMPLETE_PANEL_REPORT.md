# Phase 2B S0R complete-panel replay

Terminal verdict: `SURROGATE_AP_MISALIGNMENT_FULL_PANEL`

This was a synthetic-only replay over a pair universe constructed from
metadata fields only. No MONN residue-edge or affinity value was read.

## Contract

- Train hash panel: 14,333 pairs / 298 components
- Complete held-out A: 44,746 pairs / 81 components
- Candidate-path teacher AP: 0.997935
- Balanced train-only ray scale: 22.27278499

## Decision

| trajectory | AP(0) | AP(100) | drop | UCB95(delta) | misaligned |
|---|---:|---:|---:|---:|---|
| original_gauge | 0.997935 | 0.390110 | 0.607826 | -0.588935 | true |
| balanced_gauge | 0.997935 | 0.884815 | 0.113120 | -0.106646 | true |

The result localizes only the synthetic control. It does not identify
biology, affinity, few-shot adaptation or a biological z coordinate.
Real Phase 2B and the frozen law operator remain untouched.
