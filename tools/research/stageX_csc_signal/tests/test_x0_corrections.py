"""Stage X0 round-2 correction tests: regression + contract tests for the
corrected pair table and I2 instrument.

Historical-defect tests (test_old_*) document that the round-1 artifacts
encode the defects that were reproduced; they assert on the preserved
negative evidence. Contract tests (test_i2_*) assert the corrected behavior.
"""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from x0_common import (PREREG_SHA, sha256_seed, stable_rng, parse_mutation_list,
                       parse_construct_range, map_canonical_to_construct)
from x0_i2 import window_onehot, ESM_WINDOW_RADIUS, ESM_MAX_LEN

PAIR_TABLE = json.loads((HERE / 'X0_PAIR_TABLE.json').read_text(encoding='utf-8'))
I2 = json.loads((HERE / 'X0_I2.json').read_text(encoding='utf-8'))
OLD = json.loads((HERE / 'X0_INSTRUMENTS.json').read_text(encoding='utf-8'))


def pair_rows(status=None):
    rows = PAIR_TABLE['pairs']
    return [r for r in rows if status is None or r['admission_status'] == status]


# ------------------------------------------------- historical defect evidence
def test_old_artifact_documents_mutation_token_zero_denominator():
    assert OLD['representation_capability']['mutation_token']['ratio'] > 1e11


def test_old_artifact_documents_midpoint_window_artifact():
    assert abs(OLD['representation_capability']['local_window']['ratio'] - 1.0) < 1e-6


def test_old_artifact_documents_hash_based_random_ratio():
    assert abs(OLD['representation_capability']['random']['ratio'] - 1.0269) < 1e-3


def test_prereg_still_frozen():
    import hashlib
    digest = hashlib.sha256((HERE / 'STAGE_X0_PREREGISTRATION.md').read_bytes()).hexdigest()
    assert digest == PREREG_SHA


# ---------------------------------------------------------------- pair table
def test_braf_historical_renumbering_admitted_with_canonical_coordinate():
    row = next(r for r in pair_rows() if r['reported_construct'] == 'BRAF (V599E)')
    assert row['admission_status'] == 'admitted_point_pair'
    m = row['mutations'][0]
    assert m['reported_pos'] == 599 and m['canonical_coordinate'] == 600
    assert m['expected_old_residue'] == 'V' and m['observed_reference_residue'] == 'V'
    assert m['new'] == 'E'
    assert 'nature00766' in ' '.join(m['mapping_evidence'])


def test_pdgfra_uses_human_p16234_not_zebrafish():
    rows = [r for r in pair_rows() if r['parent_kinase'] == 'PDGFRA']
    assert all(r['canonical_accession'] == 'P16234' for r in rows)
    assert all(r['species'] == 'Homo sapiens' for r in rows)
    t674 = next(r for r in rows if 'T674I' in r['reported_construct'])
    v561 = next(r for r in rows if 'V561D' in r['reported_construct'])
    d842 = next(r for r in rows if 'D842V' in r['reported_construct'])
    assert t674['admission_status'] == 'admitted_point_pair'
    assert v561['admission_status'] == 'admitted_point_pair'
    assert d842['admission_status'] == 'excluded_construct_unresolved'
    assert '1210' in d842['exclusion_reason']  # reported construct exceeds canonical length
    assert t674['notation_fix_note'] is not None  # T6741I -> T674I recorded


def test_no_admitted_pair_has_residue_mismatch():
    for r in pair_rows('admitted_point_pair'):
        for m in r['mutations']:
            assert m['residue_verified'] is True, r['reported_construct']
            assert m['construct_residue_verified'] is True, r['reported_construct']


def test_excluded_classes_carry_reasons():
    census = PAIR_TABLE['census']['admission_counts']
    assert census['excluded_multi_point'] == 3
    assert census['excluded_deletion'] == 2
    assert census['excluded_deletion_plus_point'] == 3
    assert census['excluded_insertion'] == 1
    assert census['excluded_construct_unresolved'] == 2
    assert census['admitted_point_pair'] == 65
    for r in pair_rows():
        if r['admission_status'] != 'admitted_point_pair':
            assert r['exclusion_reason'], r['reported_construct']


def test_construct_coordinate_mapping():
    f317 = next(r for r in pair_rows() if r['reported_construct'] == 'ABL1 (F317I)')
    assert f317['mutations'][0]['construct_coordinate'] == 200  # 118-535 -> 317-117
    alk = next(r for r in pair_rows() if r['reported_construct'] == 'ALK (C1156Y)')
    assert alk['mutations'][0]['construct_coordinate'] == 99  # 1058-1620 -> 1156-1057
    kit = next(r for r in pair_rows() if r['reported_construct'] == 'C-KIT (V559D)')
    assert kit['mutations'][0]['construct_coordinate'] == 16  # 544-976 -> 559-543


def test_egfr_deletion_segments_recorded():
    row = next(r for r in pair_rows() if r['reported_construct'] == 'EGFR (d746-750)')
    del_ = [m for m in row['mutations'] if m['kind'] == 'deletion'][0]
    assert del_['deleted_segment_observed'] == 'ELREA'
    assert row['admission_status'] == 'excluded_deletion'


def test_pair_table_required_columns_present():
    required = {'dataset', 'assay_row_s2', 'parent_kinase', 'reported_construct',
                'canonical_accession', 'isoform', 'species', 'reported_mutation_notation',
                'construct_coordinate', 'admission_status', 'exclusion_reason',
                'confidence', 'provenance'}
    for r in pair_rows():
        assert required <= set(r.keys()), r['reported_construct']
        pts = [m for m in r['mutations'] if m['kind'] == 'point']
        if pts:
            assert any(m.get('mapping_evidence') for m in pts), r['reported_construct']


def test_all_parents_have_species_and_accession():
    assert len(PAIR_TABLE['parents']) == 21
    for p in PAIR_TABLE['parents']:
        assert p['species'] == 'Homo sapiens'
        assert p['canonical_length'] > 0


# ------------------------------------------------------------------------ I2
def test_pair_centered_window_distance_is_exactly_sqrt2():
    """A single substitution changes exactly two one-hot positions."""
    rep = I2['representations']['pair_centered_local_window']
    for p in rep['pairs']:
        assert abs(p['d'] - np.sqrt(2.0)) < 1e-5


def test_window_onehot_centered_at_mutation_not_midpoint():
    # ABL1 E255K on P00519: window at 255 must differ from window at midpoint
    seq = (HERE / 'uniprot' / 'P00519.fasta').read_text(encoding='utf-8').splitlines()
    seq = ''.join(l.strip() for l in seq if not l.startswith('>'))
    w_mut = window_onehot(seq, 255, ESM_WINDOW_RADIUS)
    w_mid = window_onehot(seq, len(seq) // 2, ESM_WINDOW_RADIUS)
    assert not np.array_equal(w_mut, w_mid)
    mt = seq[:254] + 'K' + seq[255:]
    assert abs(np.linalg.norm(window_onehot(seq, 255, ESM_WINDOW_RADIUS)
                              - window_onehot(mt, 255, ESM_WINDOW_RADIUS)) - np.sqrt(2)) < 1e-5


def test_no_degenerate_denominators():
    for name, rep in I2['representations'].items():
        if name == 'random' or rep.get('denominator') is None:
            continue
        assert rep['denominator'] > 0 and np.isfinite(rep['denominator'])
        assert rep['ratio'] < 10, (name, rep['ratio'])


def test_mutation_token_excluded_from_gate():
    gate = I2['gate']
    assert gate['excluded_from_gate']['mutation_token'] == 'excluded_from_representation_gate'
    assert 'random' in gate['excluded_from_gate']
    assert I2['edit_descriptor']['denominator'] is None
    assert 'mutation_token' not in gate['eligible_representations']


def test_random_control_labelled_sensitivity_only():
    rep = I2['representations']['random']
    assert 'sensitivity' in rep['interpretation']
    assert 'NOT' in rep['interpretation']
    assert rep['type'] == 'control'


def test_esm_window_admission_rules():
    excluded = {e['construct'] for e in I2['census']['esm_excluded_pairs']}
    assert 'ALK (C1156Y)' in excluded and 'LRRK2 (G2019S)' in excluded
    assert 'C-MET (F1200I)' in excluded and 'TIE2 (Y1108F)' in excluded
    rep = I2['representations']['global_esm']
    for p in rep['pairs']:
        row = next(r for r in pair_rows() if r['reported_construct'] == p['construct'])
        assert row['mutations'][0]['canonical_coordinate'] <= ESM_MAX_LEN


def test_klifs_census_fields_and_reasons():
    census = I2['census']['klifs_parent_coverage']
    assert len(census) == 21
    for c in census:
        assert 'parent_kinase' in c and 'pocket_available' in c and 'reason' in c
    lrrk = next(c for c in census if c['parent_kinase'] == 'LRRK2')
    assert lrrk['pocket_available'] and lrrk['alignment'] is None
    rep = I2['representations']['klifs_pocket']
    for p in rep['pairs']:
        if p.get('excluded'):
            assert p.get('reason'), p
        else:
            assert p.get('pocket_index') is not None
            assert abs(p['d'] - np.sqrt(2)) < 1e-5


def test_klifs_gatekeeper_maps_to_pocket_index_45():
    rep = I2['representations']['klifs_pocket']
    gatekeepers = {'EGFR (T790M)', 'C-KIT (T670I)', 'RET (V804L)', 'FGFR1 (V561M)',
                   'C-SRC (T341M)', 'PDGFRα (T674I)'}
    found = {p['construct']: p['pocket_index'] for p in rep['pairs'] if not p.get('excluded')}
    for g in gatekeepers:
        assert found.get(g) == 45, (g, found.get(g))


def test_i2_gate_pass_with_three_admissible():
    gate = I2['gate']
    assert gate['n_passing'] == 3 and gate['n_local_passing'] == 3
    assert gate['pass'] is True


def test_i2_bootstrap_cluster_based():
    for name in ('global_esm', 'esm_local_window', 'klifs_pocket', 'composition',
                 'pair_centered_local_window'):
        b = I2['representations'][name]['bootstrap']
        assert b['draws'] == 2000 and b['n_clusters'] > 0
        assert b['ci_lo'] <= b['estimate'] <= b['ci_hi']


# ------------------------------------------------------------------ seeds
def test_sha256_seed_cross_process_stable():
    code = ("import sys; sys.path.insert(0, '.'); "
            "from x0_common import stable_rng; "
            "print(','.join(f'{v:.6f}' for v in stable_rng('stageX0','random_representation','KIT','dim',128).normal(0,1,size=4)))")
    outs = set()
    for _ in range(2):
        p = subprocess.run([sys.executable, '-c', code], cwd=str(HERE),
                           capture_output=True, text=True, timeout=120)
        assert p.returncode == 0, p.stderr
        outs.add(p.stdout.strip())
    assert len(outs) == 1


def test_sha256_seed_deterministic_values():
    a = sha256_seed('stageX0', 'random_representation', 'KIT')
    b = sha256_seed('stageX0', 'random_representation', 'KIT')
    assert a == b
    rng1 = stable_rng('x', 1); rng2 = stable_rng('x', 1)
    assert np.array_equal(rng1.normal(size=10), rng2.normal(size=10))


# ----------------------------------------------------------------- parsing
def test_parse_mutation_notation_quirks():
    assert parse_mutation_list('T6741I')[0]['kind'] == 'point'
    assert parse_mutation_list('V559D, T670I')[1]['old'] == 'T'
    itd = parse_mutation_list('internal tandem duplication aa591-601')
    assert itd[0]['kind'] == 'insertion' and itd[0]['start'] == 591
    assert parse_mutation_list('-') == []


def test_parse_construct_ranges():
    assert parse_construct_range('aa 27-end', 1130)['kind'] == 'range'
    assert parse_construct_range('aa 27-end', 1130)['end'] == 1130
    assert parse_construct_range('aa 544-976 (end)', 976)['start'] == 544
    r = parse_construct_range('cytoplasmic domain [669-745, 751-1210(end)', 1210)
    assert r['kind'] == 'parts' and r['parts'] == [[669, 745], [751, 1210]]
    assert parse_construct_range('full-length', 500)['kind'] == 'full'
    assert map_canonical_to_construct(600, parse_construct_range('aa 416-766', 766)) == 185
