# Stage CIIP-1A 2x2 root-cause diagnostic — archived report (2026-08-19)

Artifact: RESULT_2X2_DIAG.json (prereg ee844b2b..., amendments A1-A4,
launch 4, PYEXIT=0). This stage is a root-cause diagnostic on the
oracle-covered subset; it is NOT CIIP-1A and cannot produce a PASS.

## Mandated status block

```text
KLIFS structural collapse: supported
oracle ESM restores nonconstant outputs: supported on covered subset
representation R2 gain: ambiguous
centered-only objective superiority: not supported
representation x objective interaction: ambiguous
biological protein-conditioned signal: UNRESOLVED
deployable representation: not validated
CIIP-1A PASS: not authorized
```

## Cell table (single seed, 9 covered test pairs, 6 parents)

| cell | R2 (pair-mean) | Spearman | sign-acc | MSE | nonconst | parents | collapsed | val MSE |
|---|---|---|---|---|---|---|---|---|
| KLIFS joint | -0.0745 | -0.038 | 0.479 | 247.9 | 3/9 | 2/6 | yes | 263.7 |
| KLIFS centered | +0.0136 | +0.115 | 0.515 | 229.5 | 3/9 | 2/6 | yes | 251.8 |
| ESM joint | +0.0436 | +0.257 | 0.628 | 190.0 | 9/9 | 6/6 | no | 211.1 |
| ESM centered | +0.0263 | +0.278 | 0.674 | 194.3 | 9/9 | 6/6 | no | 200.7 |

Structural cap: KLIFS has nonzero input on only 3/9 test pairs, so the
KLIFS cells are collapsed BY CONSTRUCTION (their collapse is a
representation statement, not an optimization statement). The ESM cells
restore nonconstant outputs on 9/9 test pairs across all 6 parents —
but their R2 is near zero: restored VARIATION is not restored SIGNAL.

## Effects (observed pair-mean R2; parent-cluster bootstrap, 2000 draws)

| effect | observed pair | observed parent | CI lo2.5..hi97.5 | status | LOPO stable |
|---|---|---|---|---|---|
| representation (joint) | +0.1181 | +0.0622 | -0.117..+0.300 | ambiguous | yes |
| representation (centered) | +0.0127 | -0.0676 | -0.281..+0.232 | absent | no |
| objective (KLIFS) | +0.0881 | +0.0971 | 0.000..+0.198 | ambiguous | yes |
| objective (ESM) | -0.0174 | -0.0328 | -0.165..+0.120 | absent | no |
| interaction | -0.1054 | -0.1299 | -0.337..+0.091 | ambiguous | yes |

## Interpretation (bounded per the frozen interpretation_bounds)

- Objective dominance persists at step 1 (R_g 616.8 KLIFS / 938.1 ESM,
  C_g ~ 0) but the centered-only objective does NOT improve ESM results
  (objective effect on ESM absent): dominance is confirmed, benefit of
  removing it is not established on this surface.
- The objective effect on KLIFS is ambiguous AND operates on
  structurally collapsed cells — not a fair estimate (frozen bound).
- Interaction ambiguous; every effect crosses or touches zero.
- oracle local ESM lifted the structural zero-space, but the only
  signal-shaped quantity (sign accuracy 0.67, Spearman 0.28) needs
  control-arm attribution before ANY signal claim; R2 <= 0.05 means
  even the correct arm explains little centered variance.

## Objective-sampling report (joint cells)

KLIFS joint L_abs pool 8,052 valid cells (WT+variant rows); ESM joint
5,837 (variant rows only) — the frozen objective-sampling confound is
confirmed and recorded; the stage is a matched representation
diagnostic, not a pure factorial estimate.

## Consequence

The attribution of the restored nonconstant response requires the
control-arm stage (annotation-shortcut audit): correct mutation window
vs random matched window vs family shuffle vs random protein vs zero
floors vs free pairwise. Until then the oracle ESM result remains a
representation-capability observation, NOT protein-conditioned signal
evidence.
