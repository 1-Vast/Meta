# Stage Q2d-1d final report — span-restricted identifiable interaction (2026-08-19)

## Verdict

**GATE FAIL.** The Q2d-1d ladder (never-moved gate) fails at its primary
surface M1 level A, double-cold:

| Quantity | Observed (median over 3 truth seeds) | Threshold | Result |
|---|---|---|---|
| correct-arm dead-zone dz (double-cold) | 0.5616 | ≥ 0.70 | FAIL |
| correct-arm Spearman (double-cold) | 0.1278 | ≥ 0.30 | FAIL |
| gap vs ligand_only (dz) | 0.0374 | ≥ 0.05 | FAIL |
| correct dz vs best negative (family_preserving 0.6121) | −0.0505 | ≥ +0.03 | FAIL |

Adjudicator: `stageQ2d1d.../adjudicate_d.py` → `Q2D1D_GATE.json`
(GATE_PASS=false, schema GATE.v1, prereg baf4bb72...). Every negative arm
itself fails the gate (every_negative_fails=true), but the correct arm
fails it too, and the family-preserving shuffle beats the correct arm on
double-cold dz — the learned map carries family-identity structure rather
than identifiable interaction structure.

## Execution record

- Ladder `runner_d.py` (PID 20348) ran M1 levels A-E x 3 seeds x 8 arms
  to completion and then **crashed at M2 truth generation**:
  `NameError: name 'PCA_VT' is not defined` — the frozen M2 definition
  requires the pre-compression PCA basis, which exists in the frozen
  label-free feature artifact `q2d1d_features.npz` (key PCA_VT) but was
  never loaded by truth_d.py. M1 results are preserved in runner_d.log.
- Because the crash lost the ladder JSON (which alone carries Spearman
  and sign-accuracy per arm), a **M1-A-only recovery run**
  (`recover_m1a.py`) re-executed the exact frozen runner_d code path for
  M1 level A and wrote `Q2D1D_LADDER.json` (schema LADDER.v1,
  provenance flags, `recovered: true`).
- Recovery cross-check vs the original runner_d.log: 7 of 8 arms
  reproduced bitwise; the family_preserving_shuffle arm differs (e.g.
  seed 0 dc dz 0.683 recovered vs 0.706 in the original log) because the
  frozen runner iterates `set(fams.tolist())`, whose order is
  PYTHONHASHSEED-dependent. The gate verdict is insensitive to this
  (correct dz 0.562 < family_preserving dz in both instantiations).
- Stage GPU regression tests: 17 passed (censored-loss tensor path,
  minibatch-order pin, 1e runner tests, AD1 tests).

## M1 ladder (double-cold dz, median over seeds; from the original log)

| Level | correct | ligand_only | shuffled | family_pres. | random | oracle |
|---|---|---|---|---|---|---|
| A | 0.562 (0.562/0.570/0.499) | 0.524 | 0.528 | 0.605 (0.706/0.605/0.536) | 0.545 | 0.992 |
| B | 0.519 (0.454/0.552/0.519) | 0.492 | 0.535 | 0.577 | 0.545 | 0.990 |
| C | 0.564 (0.734/0.564/0.491) | 0.472 | 0.518 | 0.609 | 0.555 | 0.961 |
| D | 0.432 (0.486/0.432/0.414) | 0.472 | 0.515 | 0.619 | 0.563 | 0.968 |
| E | 0.618 (0.655/0.618/0.501) | 0.454 | 0.500 | 0.459 | 0.447 | 0.761 |

At no level does the correct arm reach dz 0.70; at A-D a negative arm
beats it; at D the correct arm is below ligand_only. Oracle ceilings
0.76-1.0 confirm the truth is recoverable in principle — the trained
low-rank bilinear learner does not recover it under double-cold.

## Adjudicator/preregistration wording note

The Q2d-1d prereg gate text (baf4bb72...) lists for the correct arm:
Spearman ≥ 0.30, dead-zone sign accuracy ≥ 0.70, gap vs ligand_only ≥
0.05, margin over best negative ≥ 0.03. The committed adjudicator
additionally requires raw double-cold dz ≥ 0.70 for every arm (a rule
written in the Q2d-1e prereg and inherited here). The verdict does not
depend on this wording difference: correct dz 0.5616 < 0.70, sp 0.1278 <
0.30, and the margin is negative under either reading. Flagged for the
record; no gate was moved.

## Downstream

Per the bounded chain: Q2d-1d FAIL → the frozen successor **Q2d-1e**
(span-initialized A in the train-row feature span + L2 1e-3 on factor
maps) is the last authorized stage. It is running with AD1 repairs
(truth_e.py; M1/M2/M3 streams bit-identical to truth_d, tests green).
If Q2d-1e fails, exactly one finite diagnostic (A = V_train·G
reparameterization, frozen) runs, then a terminal PASS/FAIL on the
low-rank bilinear learner family is recorded. Q2d-2/Q2d-3/B1 remain
unauthorized until then.
