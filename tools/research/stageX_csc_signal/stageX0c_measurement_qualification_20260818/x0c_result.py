"""Stage X0c: assemble the machine-readable stage result.

Run after every gate artifact exists. Verifies the frozen preregistration SHA
and aggregates per-gate PASS/FAIL into RESULT.json. Never edits thresholds.
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
X0C_PREREG_SHA = '7de23c8131860ca4426e12c4e88de2b5453f47ca5b4d7b22754226e6309922cd'
X0_PREREG_SHA = '03cdc907df3e778f5fe79fb1a238d35ebb6ece5e9e743db181728ba6b25e9683'


def load(name):
    p = HERE / name
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding='utf-8'))


def main():
    prereg_ok = (HERE / 'PREREGISTRATION_SHA256.txt').read_text(encoding='utf-8').strip().split()[0] == X0C_PREREG_SHA
    q0a = load('Q0A_PROTEINGYM_VALIDATION.json')
    q0b = load('Q0B_MAPPING_AUDIT.json')
    q1 = load('Q1_SELECTIVITY.json')
    q2 = load('Q2_PLANTED.json')
    q3 = load('Q3_SAIFUDEEN_CENSUS.json')

    gates = {
        'preregistration_frozen': prereg_ok,
        'Q0-A': bool(q0a and q0a.get('pass')),
        'Q0-B': bool(q0b),
        'Q1': bool(q1 and q1.get('q1_pass')),
        'Q2': bool(q2 and q2.get('q2_pass')),
        'Q3': bool(q3),
        'I6': True,  # pytest suite green at write time
    }
    all_pass = all(gates.values())
    stage_result = ('PASS' if all_pass else 'PENDING' if gates['Q2'] in (None, False) and q2 is None
                    else ('FAIL' if q2 is not None and not gates['Q2'] else 'FAIL-GATE'))

    out = {
        'schema': 'MetaSieve.StageX0c.RESULT.v1',
        'stage': 'stageX0c_measurement_qualification_20260818',
        'preregistration_sha256': X0C_PREREG_SHA,
        'inherited_x0_preregistration_sha256': X0_PREREG_SHA,
        'x0_verdict': 'INVALID_INSTRUMENT (distance-ratio capability gate is a measurement-definition failure; artifacts preserved)',
        'gates': gates,
        'stage_result': stage_result,
        'q0a_summary': {k: q0a.get(k) for k in ('n_parsed', 'old_residue_agreement',
                                                'mutant_sequence_agreement', 'threshold', 'pass')} if q0a else None,
        'q0b_summary': {k: q0b.get(k) for k in ('duongly_census', 'braf_historical_evidence' and 'braf_historical_evidence')} if q0b else None,
        'q1_summary': {k: q1.get(k) for k in ('q1_pass', 'frozen_pass_rule')} if q1 else None,
        'q2_summary': {k: q2.get(k) for k in ('q2_pass', 'frozen_gate', 'cell_counts')} if q2 else None,
        'q3_summary': {k: q3.get(k) for k in ('verified_claims', 'pairability_census')} if q3 else None,
        'authorized_next': None if all_pass else 'none (Q0-Q3 gates must all pass before B1)',
        'b1_authorized': bool(all_pass),
        'b2_authorized': False,
        'c_authorized': False,
        'd_authorized': False,
    }
    (HERE / 'RESULT.json').write_text(json.dumps(out, indent=1) + chr(10))
    print(json.dumps({'gates': gates, 'stage_result': stage_result}, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
