"""Stage X0 round-2: corrected WT-mutant pair table (mutation/construct mapping audit).

Builds X0_PAIR_TABLE.json from the governed local caches. Rules (see
X0_CORRECTION_PLAN.md):
- old residue at the canonical coordinate must match the reference sequence
  BEFORE any substitution; mismatches are excluded, never silently fixed.
- historical numbering is an explicit, cited mapping (BRAF V599E -> V600E).
- wrong-species accessions are corrected with cited mapping basis
  (PDGFRalpha: S1 lists Q9DE49 = Danio rerio; S1 GenBank NP_006197 and KLIFS
  resolve human PDGFRA = P16234).
- multi-mutation / deletion / insertion / unknown notations get explicit
  admission statuses and exclusion reasons; every excluded row stays in the
  coverage census.
"""
from __future__ import annotations
import json, re
from pathlib import Path
import numpy as np
import pandas as pd

from x0_common import (
    PREREG_SHA, AAS, HERE, load_duongly, normalize_construct_name,
    normalize_parent_name, parse_mutation_list, parse_point_notation,
    parse_construct_range, construct_sequence, map_canonical_to_construct,
    write_artifact, sha256_file)

# --- cited mapping overrides (historical numbering / wrong accessions) ---
BRAF_EVIDENCE = [
    'https://doi.org/10.1038/nature00766 (Davies et al. 2002, Nature 417:949, reported BRAF V599E using the shorter reference)',
    'https://www.uniprot.org/uniprotkb/P15056 (canonical P15056: V600; UniProt variant VAR_018629 p.Val600Glu; residue 599 is T)',
    'https://klifs.net/api/kinase_ID?kinase_name=BRAF&species=Human (KLIFS BRAF Human -> P15056)',
]
PDGFRA_EVIDENCE = [
    'Duong-Ly Table S1 lists GenBank NP_006197 (human PDGFRA isoform) and Protein Accession Q9DE49 (UniProt: Danio rerio pdgfra) for the same rows',
    'https://rest.uniprot.org/uniprotkb/P16234.fasta (human PDGFRA canonical)',
    'https://klifs.net/api/kinase_ID?kinase_name=PDGFRa&species=Human (KLIFS PDGFRa Human -> P16234)',
]

# accession override: parent -> canonical accession
ACCESSION_OVERRIDES = {'PDGFRA': 'P16234'}
# reported -> canonical coordinate remap (historical numbering), cited.
HISTORICAL_RENUMBERING = {
    'BRAF': {599: {'canonical_pos': 600,
                   'basis': 'historical numbering V599E renumbered to canonical V600E '
                            '(shorter reference used in Davies 2002); P15056 residue 599 is T, 600 is V',
                   'evidence': BRAF_EVIDENCE}},
}
# reported notation quirks (S1 Mutation column) -> canonical notation
NOTATION_FIXES = {
    'T6741I': {'canonical': 'T674I',
               'basis': 'S1 Mutation column typo "T6741I"; construct name column reads T674I'},
}


def parse_name_mutations(name: str):
    """Parse mutation notation embedded in a construct name like 'EGFR(d747-749/A750P)'."""
    body = name.split('(', 1)[1].rsplit(')', 1)[0] if '(' in name else ''
    if not body:
        return []
    return parse_mutation_list(body.replace('/', ','))


def build_parent_table(info, seqs):
    """One row per parent kinase (21)."""
    parents = {}
    for _, row in info.iterrows():
        name = str(row['Kinase (Mutation)']).strip()
        parent = normalize_parent_name(name)
        acc = str(row['Protein Accession #']).strip()
        if parent not in parents:
            parents[parent] = {'accessions_seen': [], 'genbanks_seen': []}
        if acc not in parents[parent]['accessions_seen']:
            parents[parent]['accessions_seen'].append(acc)
        gb = str(row['Genbank Accession #']).strip()
        if gb not in parents[parent]['genbanks_seen']:
            parents[parent]['genbanks_seen'].append(gb)
    rows = []
    for parent in sorted(parents):
        acc = ACCESSION_OVERRIDES.get(parent, parents[parent]['accessions_seen'][0])
        seq = seqs.get(acc, {}).get('sequence', '')
        header = seqs.get(acc, {}).get('header', '')
        m = re.search(r'OS=([^ ]+ [^ ]+)', header)
        species = m.group(1) if m else 'unresolved'
        iso = re.search(r'SV=([0-9]+)', header)
        rows.append({
            'parent_kinase': parent,
            'canonical_accession': acc,
            'accessions_reported_in_s1': parents[parent]['accessions_seen'],
            'genbank_reported_in_s1': parents[parent]['genbanks_seen'],
            'species': species,
            'isoform': f"SV={iso.group(1)}" if iso else 'unresolved',
            'canonical_length': len(seq),
            'reference_fasta_sha256': sha256_file(HERE / 'uniprot' / f'{acc}.fasta'),
            'provenance': 'UniProt canonical reference fetched locally; S1-reported accessions listed separately',
        })
    return rows


def build_pair_rows(info, matrix, seqs):
    s2_labels = [str(x).strip() for x in matrix.iloc[:, 0].tolist()]
    s2_norm = [normalize_construct_name(x) for x in s2_labels]

    rows = []
    for s1_idx, (_, s1) in enumerate(info.iterrows()):
        name = str(s1['Kinase (Mutation)']).strip()
        norm = normalize_construct_name(name)
        parent = normalize_parent_name(name)
        acc_reported = str(s1['Protein Accession #']).strip()
        acc = ACCESSION_OVERRIDES.get(parent, acc_reported)
        canon = seqs.get(acc, {})
        seq = canon.get('sequence', '')
        header = canon.get('header', '')
        species = (re.search(r'OS=([^ ]+ [^ ]+)', header).group(1)
                   if re.search(r'OS=([^ ]+ [^ ]+)', header) else 'unresolved')
        iso = re.search(r'SV=([0-9]+)', header)
        clone_raw = str(s1['Clone']).strip()
        crange = parse_construct_range(clone_raw, len(seq))
        cseq = construct_sequence(seq, crange) if seq else None

        if norm in s2_norm:
            s2_row_idx = s2_norm.index(norm)
        else:
            s2_row_idx = None

        reported_mut_raw = str(s1['Mutation']).strip()
        fix = NOTATION_FIXES.get(reported_mut_raw)
        if fix:
            reported_mut_raw = fix['canonical']
            notation_fix_note = fix['basis']
        else:
            notation_fix_note = None
        muts = parse_mutation_list(reported_mut_raw)
        name_muts = parse_name_mutations(name)
        # the construct name is authoritative for construct identity: merge
        # deletions/insertions that only appear in the name (e.g. d746-750/T790M
        # whose S1 Mutation column lists only T790M).
        have_kinds = {m['kind'] for m in muts}
        for nm in name_muts:
            if nm['kind'] in ('deletion', 'insertion') and nm['kind'] not in have_kinds:
                muts.append(nm)
                have_kinds.add(nm['kind'])
        if not muts:
            muts = [{'kind': 'unknown', 'reported': reported_mut_raw}]

        for m in muts:
            if m['kind'] == 'point':
                if parent in HISTORICAL_RENUMBERING and m['pos'] in HISTORICAL_RENUMBERING[parent]:
                    h = HISTORICAL_RENUMBERING[parent][m['pos']]
                    m['reported_pos'] = m['pos']
                    m['pos'] = h['canonical_pos']
                    m['mapping_basis'] = h['basis']
                    m['mapping_evidence'] = h['evidence']
                else:
                    m['mapping_basis'] = 'direct: reported coordinate is canonical UniProt numbering'
                    m['mapping_evidence'] = [f'https://www.uniprot.org/uniprotkb/{acc}']
                expected = m['old']
                observed = seq[m['pos'] - 1] if seq and 1 <= m['pos'] <= len(seq) else None
                m['expected_old_residue'] = expected
                m['observed_reference_residue'] = observed
                m['canonical_coordinate'] = m['pos']
                m['construct_coordinate'] = map_canonical_to_construct(m['pos'], crange)
                m['residue_verified'] = (observed == expected)
                if m['construct_coordinate'] is not None and cseq:
                    c_obs = cseq[m['construct_coordinate'] - 1]
                    m['construct_residue'] = c_obs
                    m['construct_residue_verified'] = (c_obs == expected)
                else:
                    m['construct_residue'] = None
                    m['construct_residue_verified'] = None
            elif m['kind'] == 'deletion':
                seg = seq[m['start'] - 1:m['end']] if seq and m['end'] <= len(seq) else None
                m['deleted_segment_observed'] = seg

        kinds = [m['kind'] for m in muts]
        status, reason = decide_admission(muts, kinds, crange, seq, cseq, acc, acc_reported, parent, notation_fix_note)

        if acc_reported != acc:
            accession_note = (f"S1 Protein Accession {acc_reported} superseded by {acc}: "
                              f"{acc_reported} is not the human canonical for {parent}")
        else:
            accession_note = f"S1 accession {acc_reported} used as canonical"

        first_point = next((m for m in muts if m['kind'] == 'point'), None)
        rows.append({
            'dataset': 'duongly_2016',
            'assay_row_s2': s2_labels[s2_row_idx] if s2_row_idx is not None else None,
            'assay_row_s2_index': s2_row_idx,
            's1_row_index': int(s1_idx),
            'construct_coordinate': (first_point.get('construct_coordinate')
                                     if first_point else None),
            'parent_kinase': parent,
            'reported_construct': name,
            'canonical_accession': acc,
            's1_reported_protein_accession': acc_reported,
            's1_reported_genbank': str(s1['Genbank Accession #']).strip(),
            'isoform': f"SV={iso.group(1)}" if iso else 'unresolved',
            'species': species,
            'reported_mutation_notation': str(s1['Mutation']).strip(),
            'construct_name_notation': name,
            'clone_field': clone_raw,
            'construct_range': crange,
            'construct_sequence_length': len(cseq) if cseq is not None else None,
            'construct_sequence_resolvable': cseq is not None,
            'mutations': muts,
            'mutation_class': classify(kinds),
            'admission_status': status,
            'exclusion_reason': reason,
            'accession_mapping_note': accession_note,
            'notation_fix_note': notation_fix_note,
            'confidence': 'high' if status == 'admitted_point_pair' and
                          all(m.get('construct_residue_verified', True) is not False for m in muts)
                          else ('medium' if status == 'admitted_point_pair' else 'excluded'),
            'provenance': 'Duong-Ly 2016 Cell Reports Table S1/S2 + UniProt canonical reference '
                          f'{acc} (local SHA-256 in input_sha256)',
        })
    return rows


def classify(kinds):
    if kinds == ['point']:
        return 'point'
    if all(k == 'point' for k in kinds):
        return 'multi_point'
    if kinds == ['deletion']:
        return 'deletion'
    if 'deletion' in kinds and 'point' in kinds:
        return 'deletion_plus_point'
    if kinds == ['insertion']:
        return 'insertion'
    if 'unknown' in kinds:
        return 'unknown_notation'
    return 'mixed_other'


def decide_admission(muts, kinds, crange, seq, cseq, acc, acc_reported, parent,
                      notation_fix_note=None):
    cls = classify(kinds)
    if cls == 'point':
        m = muts[0]
        if not m.get('residue_verified'):
            return ('excluded_old_residue_mismatch',
                    f"expected {m['expected_old_residue']} at canonical position {m['pos']} "
                    f"of {acc}, observed {m.get('observed_reference_residue')}; reported "
                    f"notation {m.get('reported', '')}")
        if m.get('construct_coordinate') is None:
            if crange.get('kind') == 'unresolved':
                return ('excluded_construct_unresolved',
                        f"construct range could not be parsed from clone field {crange.get('note')}; "
                        f"canonical position {m['pos']} cannot be mapped to a construct coordinate")
            return ('excluded_outside_construct',
                    f"canonical position {m['pos']} not inside parsed construct {crange}")
        if not cseq:
            return ('excluded_construct_unresolved',
                    f"construct range not resolvable on {acc}: reported range {crange} extends "
                    f"beyond canonical length {len(seq)}")
        if m.get('construct_residue_verified') is False:
            return ('excluded_construct_offset_error',
                    f"construct residue at construct coordinate {m['construct_coordinate']} "
                    f"is {m.get('construct_residue')}, expected {m['expected_old_residue']}")
        return ('admitted_point_pair', None)
    if cls == 'multi_point':
        return ('excluded_multi_point', 'two or more point mutations in one construct')
    if cls == 'deletion':
        return ('excluded_deletion', 'deletion-only construct')
    if cls == 'deletion_plus_point':
        return ('excluded_deletion_plus_point', 'deletion combined with point mutation')
    if cls == 'insertion':
        return ('excluded_insertion', 'internal tandem duplication / insertion construct')
    if cls == 'unknown_notation':
        return ('excluded_unknown_notation', 'mutation notation could not be parsed')
    return ('excluded_other', 'mixed or unclassified mutation notation')


def census(rows):
    out = {'total_assay_rows_s2': 97, 'wt_parent_rows': 21, 'mutant_construct_rows': len(rows)}
    counts = {}
    for r in rows:
        counts[r['admission_status']] = counts.get(r['admission_status'], 0) + 1
    out['admission_counts'] = counts
    out['admitted_point_pairs'] = counts.get('admitted_point_pair', 0)
    return out


def main():
    info, matrix, seqs = load_duongly()
    parents = build_parent_table(info, seqs)
    rows = build_pair_rows(info, matrix, seqs)
    out = {
        'schema': 'MetaSieve.StageX.X0PairTable.v2',
        'stage': 'stageX_csc_signal',
        'preregistration_sha256': PREREG_SHA,
        'parents': parents,
        'pairs': rows,
        'census': census(rows),
        'historical_numbering_map': HISTORICAL_RENUMBERING,
        'notation_fixes': NOTATION_FIXES,
        'accession_overrides': ACCESSION_OVERRIDES,
        'admission_policy': {
            'admitted_point_pair': 'single point mutation, old residue verified at canonical coordinate, construct coordinate resolvable and residue-verified',
            'excluded_multi_point': 'two or more point mutations in one construct (no single-site window)',
            'excluded_deletion': 'deletion-only construct',
            'excluded_deletion_plus_point': 'deletion combined with point mutation',
            'excluded_insertion': 'internal tandem duplication / insertion construct',
            'excluded_unknown_notation': 'mutation notation could not be parsed',
            'excluded_old_residue_mismatch': 'reference residue at canonical coordinate does not equal expected old residue',
            'excluded_outside_construct': 'mutation position not inside the reported construct',
            'excluded_construct_unresolved': 'construct range could not be mapped onto the canonical sequence',
            'excluded_construct_offset_error': 'construct residue at mapped construct coordinate differs from expected old residue',
        },
    }
    inputs = [HERE / 'downloads/duongly_mmc2.xlsx', HERE / 'downloads/duongly_mmc3.xlsx']
    inputs += sorted((HERE / 'uniprot').glob('*.fasta'))
    write_artifact(HERE / 'X0_PAIR_TABLE.json', out, inputs)
    print(json.dumps({'schema': out['schema'],
                      'census': out['census']}, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
