"""Stage X0c Q0 unit tests: VariantRecord mutation application,
admission rules and deterministic serialization."""
import json, sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from variant_record import (Substitution, VariantRecord, apply_mutations,
                            apply_substitution, parse_proteingym_mutant)

SEQ = 'ACDEFGHIKLMNPQRSTVWYACDEFGHIK'


def test_old_residue_hard_fail():
    with pytest.raises(ValueError):
        apply_substitution(SEQ, Substitution('W', 3, 'G', 'canonical'), verify_old=True)


def test_old_residue_ok():
    out, obs = apply_substitution(SEQ, Substitution('D', 3, 'G', 'canonical'), verify_old=True)
    assert obs == 'D' and out == 'AC' + 'G' + SEQ[3:]


def test_multi_mutation_applied_atomic():
    subs = [Substitution('A', 1, 'W', 'canonical'), Substitution('E', 4, 'K', 'canonical')]
    mt, obs = apply_mutations(SEQ, subs)
    assert mt == 'WCD' + 'K' + SEQ[4:]
    assert obs == ['A', 'E']


def test_multi_mutation_fails_atomically_on_bad_residue():
    subs = [Substitution('C', 1, 'W', 'canonical'), Substitution('X', 4, 'K', 'canonical')]
    with pytest.raises(ValueError):
        apply_mutations(SEQ, subs)


def test_proteingym_parser():
    subs = parse_proteingym_mutant('A1P:D2N')
    assert [s.old for s in subs] == ['A', 'D'] and [s.pos for s in subs] == [1, 2]
    assert parse_proteingym_mutant('A1P:garbage') is None
    assert parse_proteingym_mutant('*5A') is not None  # stop-gain notation


def test_canonical_serialization_deterministic():
    vr1 = VariantRecord(dataset='duongly_2016', parent_gene='ABL1',
                        substitutions=(Substitution('E', 255, 'K', 'canonical'),))
    vr2 = VariantRecord(dataset='duongly_2016', parent_gene='ABL1',
                        substitutions=(Substitution('E', 255, 'K', 'canonical'),))
    j1 = json.dumps({k: vr1.canonical_json() for k in ['a']}, sort_keys=True)
    j2 = json.dumps({k: vr2.canonical_json() for k in ['a']}, sort_keys=True)
    assert vr1.canonical_json() == vr2.canonical_json()


def test_record_hash_changes_with_content():
    vr = VariantRecord(dataset='x', parent_gene='ABL1').with_hash()
    vr2 = VariantRecord(dataset='x', parent_gene='ABL2').with_hash()
    assert vr.record_hash != vr2.record_hash


def test_mutation_classes_vocabulary():
    from variant_record import MUTATION_CLASSES, ADMISSION_STATUSES
    assert 'historical_numbering_alias' in MUTATION_CLASSES
    assert 'construct_domain_offset' in MUTATION_CLASSES
    assert 'sequence_exceeds_plm_limit' in MUTATION_CLASSES
    assert 'quarantined' in ADMISSION_STATUSES
