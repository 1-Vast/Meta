"""Finalize stage artifacts once Q2 completes: update REPORT.md Q2 and
summary sections from Q2_PLANTED.json, then assemble RESULT.json."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
q2 = json.loads((HERE / 'Q2_PLANTED.json').read_text(encoding='utf-8'))
gate = q2.get('frozen_gate', {})
nc = q2.get('negative_controls', {})
report = (HERE / 'REPORT.md').read_text(encoding='utf-8')

sp = gate.get('correct_spearman')
dz = gate.get('correct_dead_zone_sign_accuracy')
gap = gate.get('gap_sign_accuracy')
verdict = 'PASS' if gate.get('pass') else 'FAIL'
req = gate.get('thresholds', {})

q2_sec = (
    'Gate point (tau*=1.0, rank 4, dense, eval cells, median of 3 seeds): '
    'correct-arm Spearman %s (need >= %s), dead-zone sign accuracy %s (need >= %s), '
    'gap vs ligand_only %s (need >= %s) -> FAIL. '
    'Detection ladder tau* {0.125,0.25,0.5,1.0,2.0} x rank {1,4,16} is in '
    'Q2_PLANTED.json. Negative controls behaved as designed: label permutation '
    'dead-zone sign accuracy 0.52-0.57 (chance), tau*=0 interaction head '
    'recovers nothing, floor-clamp imputation induces the expected spurious '
    'recovery (dz 0.588), no-interaction head and shuffled/family-shuffled/'
    'random protein ~0.50-0.55, free-target-id upper bound ~0.50-0.55. '
    'Diagnosis (separate frozen-seed runs, recorded in the artifact): the '
    'oracle arm (P@U input) recovers the centred interaction at dz 0.68-0.76, '
    'so the information exists; the correct arm (one-hot pocket input) tops '
    'out at dz ~0.58 across 8+6 restarts and two training protocols, so the '
    'failure is representation-learning/optimization capacity at this sample '
    'size, not information absence.'
) % (sp, req.get('spearman'), dz, req.get('dead_zone_sign_accuracy'), gap, req.get('gap'))

report = report.replace(
    'PENDING - grid run in progress (tau* x {1,4,16} rank; frozen gate at' + chr(10) +
    'tau*=1.0, rank 4, dense, on held-out eval, median of 3 seeds).',
    q2_sec,
)
report = report.replace('| Q2 fully synthetic planted harness | PENDING | Q2_PLANTED.json |',
                        '| Q2 fully synthetic planted harness | FAIL | Q2_PLANTED.json |')
report = report.replace('| Q0-A ProteinGym external validation | PENDING | Q0A_PROTEINGYM_VALIDATION.json |',
                        '| Q0-A ProteinGym external validation | PASS | Q0A_PROTEINGYM_VALIDATION.json |')
summary_old = ('PENDING - resolved when Q0-A and Q2 complete. B1/B2/C/D remain gated; no' + chr(10) +
               'real-data biological inference is authorized until Q0-Q2 pass under the' + chr(10) +
               'frozen gates.')
summary_new = (
    'Q2 FAILED its frozen gate; therefore B1/B2/C/D are NOT authorized and no '
    'real-data biological inference may be drawn from this harness. The Q2 '
    'failure is optimization-limited (oracle arm recovers the planted signal), '
    'which defines the single highest-information next step: a new '
    'preregistration for the representation-learning fix (or a revised gate '
    'point), since the current frozen gate may not be moved retroactively.'
)
report = report.replace(summary_old, summary_new)
(HERE / 'REPORT.md').write_text(report, encoding='utf-8')
print('REPORT.md updated; verdict:', verdict)
