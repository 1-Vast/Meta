"""Stage X0c Q0: typed variant-coordinate layer (VariantRecord) and
Q0-A external validation against ProteinGym DMS_substitutions.

Hard rules implemented:
  - old residue verified BEFORE applying any mutation; mismatch => hard fail
    or quarantine, never a silent fix.
  - coordinate transforms are explicit typed objects with per-record
    evidence; BRAF V599E->V600E is one cited alias, never generalized.
  - mutation classes are an enumerated vocabulary.
  - deterministic canonical serialization (sorted keys); record hash =
    SHA-256 of the canonical JSON.
  - all seeds SHA-256-derived; Python hash() never used.
"""
from __future__ import annotations
import hashlib, json, re, sys, zipfile, csv, io
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
sys.path.insert(0, str(PARENT))
from x0_common import sha256_seed, stable_rng, sha256_file, sha256_text, write_artifact

PREREG_SHA = '7de23c8131860ca4426e12c4e88de2b5453f47ca5b4d7b22754226e6309922cd'

MUTATION_CLASSES = [
    'canonical_single_substitution',
    'historical_numbering_alias',
    'mature_protein_offset',
    'signal_peptide_offset',
    'construct_domain_offset',
    'isoform_offset',
    'multi_substitution',
    'deletion',
    'insertion',
    'truncation',
    'fusion',
    'phosphorylation_state_construct',
    'unknown_notation',
    'sequence_exceeds_plm_limit',
    'mutation_position_outside_default_truncation',
]

ADMISSION_STATUSES = ['admitted', 'quarantined', 'excluded', 'unresolved']


@dataclass(frozen=True)
class Substitution:
    old: str
    pos: int
    new: str
    coordinate_kind: str  # canonical | construct


@dataclass(frozen=True)
class CoordinateTransform:
    kind: str
    reported_position: int
    canonical_position: int
    basis: str
    evidence: tuple


@dataclass(frozen=True)
class VariantRecord:
    schema: str = 'MetaSieve.StageX0c.VariantRecord.v1'
    dataset: str = ''
    source_row: str = ''
    parent_gene: str = ''
    reported_construct: str = ''
    species: str = ''
    canonical_accession: str = ''
    canonical_version: str = ''
    isoform: str = ''
    reported_mutation_notation: str = ''
    mutation_class: str = 'canonical_single_substitution'
    substitutions: tuple = ()
    coordinate_transforms: tuple = ()
    construct_start: Optional[int] = None
    construct_end: Optional[int] = None
    wt_residue: str = ''
    mutant_residue: str = ''
    parent_sequence: str = ''
    mutant_sequence: str = ''
    reference_release: str = ''
    provenance_url: str = ''
    evidence_grade: str = ''
    admission_status: str = 'unresolved'
    exclusion_reason: str = ''
    sequence_hash: str = ''
    record_hash: str = ''

    def canonical_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True, ensure_ascii=False)

    def with_hash(self):
        d = dict(self.__dict__)
        d['record_hash'] = sha256_text(self.canonical_json())
        return VariantRecord(**d)


def apply_substitution(seq, sub, verify_old=True, ref_pos=None):
    pos = ref_pos if ref_pos is not None else sub.pos
    if not (1 <= pos <= len(seq)):
        raise ValueError(f'position {pos} outside sequence length {len(seq)}')
    observed = seq[pos - 1]
    if verify_old and observed != sub.old:
        raise ValueError(
            f'old-residue mismatch at {pos}: expected {sub.old}, observed {observed}')
    return seq[:pos - 1] + sub.new + seq[pos:], observed


def apply_mutations(seq, subs, verify_old=True):
    observed = []
    for s in subs:
        _, obs = apply_substitution(seq, s, verify_old=verify_old)
        observed.append(obs)
    out = seq
    for s in subs:
        out = out[:s.pos - 1] + s.new + out[s.pos:]
    return out, observed


def parse_proteingym_mutant(mutant):
    subs = []
    for part in mutant.split(':'):
        m = re.match(r'^([A-Z*])([0-9]+)([A-Z*])$', part)
        if not m:
            return None
        subs.append(Substitution(m.group(1), int(m.group(2)), m.group(3), 'canonical'))
    return subs


def run_q0a(proteingym_zip, out_path, max_records=50000):
    rng = stable_rng('stageX0c/q0a/sample')
    # reference sequences come from the official ProteinGym reference file
    ref_csv = proteingym_zip.parent / 'ProteinGym_DMS_substitutions_reference.csv'
    ref_seqs = {}
    with open(ref_csv, encoding='utf-8') as f:
        for r in csv.DictReader(f):
            ref_seqs[r.get('DMS_id', '')] = (r.get('target_seq') or '').strip()
    result_rows = []
    n_total = 0
    n_parsed = 0
    agree_old = 0
    agree_seq = 0
    failures = []
    datasets_seen = []
    with zipfile.ZipFile(proteingym_zip) as zf:
        names = sorted(n for n in zf.namelist() if n.endswith('.csv'))
        for name in names:
            with zf.open(name) as f:
                reader = csv.DictReader(io.TextIOWrapper(f, encoding='utf-8', errors='replace'))
                rows = list(reader)
            if not rows:
                continue
            datasets_seen.append(name)
            for row in rows:
                n_total += 1
                if rng.random() > (max_records / 2700000.0):
                    continue
                result_rows.append((name, row))
                if len(result_rows) >= max_records:
                    break
            if len(result_rows) >= max_records:
                break
    for name, row in result_rows:
        dms_id = Path(name).stem
        seq = ref_seqs.get(dms_id, '')
        mutant_ann = (row.get('mutant') or '').strip()
        ref_mutseq = (row.get('mutated_sequence') or '').strip()
        subs = parse_proteingym_mutant(mutant_ann)
        if not subs or not seq:
            continue
        n_parsed += 1
        try:
            mt, obs = apply_mutations(seq, subs, verify_old=True)
        except ValueError as e:
            failures.append({'dataset': name, 'mutant': mutant_ann, 'error': str(e)})
            continue
        if all(o == s.old for o, s in zip(obs, subs)):
            agree_old += 1
        if mt == ref_mutseq:
            agree_seq += 1
        else:
            failures.append({'dataset': name, 'mutant': mutant_ann,
                             'error': 'mutated_sequence mismatch'})
    result = {
        'schema': 'MetaSieve.StageX0c.Q0A.v1',
        'preregistration_sha256': PREREG_SHA,
        'source': 'ProteinGym v1.3 DMS_ProteinGym_substitutions.zip (official benchmark)',
        'source_url': 'https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip',
        'license_note': 'ProteinGym benchmark data (MIT license); local analysis only, raw files not committed',
        'sample_seed': 'sha256(stageX0c/q0a/sample)',
        'n_datasets_scanned': len(datasets_seen),
        'n_total_rows_seen': n_total,
        'n_sampled': len(result_rows),
        'n_parsed': n_parsed,
        'old_residue_agreement': agree_old / n_parsed if n_parsed else None,
        'mutant_sequence_agreement': agree_seq / n_parsed if n_parsed else None,
        'threshold': 0.995,
        'pass_old_residue': (agree_old / n_parsed >= 0.995) if n_parsed else False,
        'pass_sequence': (agree_seq / n_parsed >= 0.995) if n_parsed else False,
        'pass': bool(n_parsed and agree_old / n_parsed >= 0.995 and agree_seq / n_parsed >= 0.995),
        'failures': failures[:20],
        'n_failures': len(failures),
    }
    write_artifact(out_path, result, [proteingym_zip])
    return result


if __name__ == '__main__':
    zip_path = PARENT / 'downloads' / 'DMS_ProteinGym_substitutions.zip'
    if not zip_path.exists():
        print('zip missing; wait for the detached download')
        raise SystemExit(1)
    res = run_q0a(zip_path, HERE / 'Q0A_PROTEINGYM_VALIDATION.json')
    print(json.dumps(res, indent=1))
