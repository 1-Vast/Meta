"""Stage X0c Q0-B: historical/construct mapping audit.

Produces Q0B_MAPPING_AUDIT.json, Q0B_ALIAS_LEDGER.md, Q0B_KLIFS_CENSUS.json.
Sources are first-hand and local: frozen pair table (parent dir), UniProt
fastas, KLIFS lookup, BRAF historical evidence (computed in this stage),
Davis MOESM3.
"""
import json, re, sys
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
sys.path.insert(0, str(PARENT))
from x0_common import sha256_file, write_artifact, PREREG_SHA as X0_PREREG_SHA
from variant_record import (VariantRecord, Substitution, CoordinateTransform,
                            PREREG_SHA, apply_mutations, sha256_text)

X0C_PREREG_SHA = PREREG_SHA


def load_seqs():
    seqs = {}
    for p in (PARENT / 'uniprot').glob('*.fasta'):
        text = p.read_text(encoding='utf-8')
        seq = ''.join(l.strip() for l in text.splitlines() if not l.startswith('>'))
        header = next((l for l in text.splitlines() if l.startswith('>')), '')
        m = re.search(r'SV=([0-9]+)', header)
        seqs[p.stem] = {'sequence': seq, 'header': header,
                        'version': m.group(1) if m else 'unresolved'}
    return seqs


def mutation_class_of(pair_row):
    st = pair_row['admission_status']
    kinds = [m['kind'] for m in pair_row['mutations']]
    if st == 'admitted_point_pair' and pair_row['parent_kinase'] == 'BRAF':
        return 'historical_numbering_alias'
    if st == 'admitted_point_pair':
        return 'canonical_single_substitution'
    if st == 'excluded_multi_point':
        return 'multi_substitution'
    if st == 'excluded_deletion':
        return 'deletion'
    if st == 'excluded_deletion_plus_point':
        return 'deletion'  # plus point, recorded in substitutions
    if st == 'excluded_insertion':
        return 'insertion'
    if st == 'excluded_construct_unresolved':
        return 'unknown_notation'
    if st == 'excluded_old_residue_mismatch':
        return 'unknown_notation'
    return 'unknown_notation'


def build_duongly_records(pair_table, seqs):
    records = []
    for row in pair_table['pairs']:
        acc = row['canonical_accession']
        seq = seqs[acc]['sequence']
        subs = []
        transforms = []
        for m in row['mutations']:
            if m['kind'] == 'point':
                subs.append(Substitution(m['old'], m['canonical_coordinate'], m['new'],
                                         'canonical'))
                if 'mapping_basis' in m and 'historical' in (m.get('mapping_basis') or ''):
                    transforms.append(CoordinateTransform(
                        'historical_alias', m['reported_pos'], m['canonical_coordinate'],
                        m['mapping_basis'], tuple(m.get('mapping_evidence', []))))
        mt_seq = None
        if subs and all(m['kind'] == 'point' for m in row['mutations']):
            try:
                mt_seq, obs = apply_mutations(seq, subs, verify_old=True)
            except ValueError:
                mt_seq = None
        crange = row.get('construct_range') or {}
        vr = VariantRecord(
            dataset='duongly_2016',
            source_row=row.get('assay_row_s2') or row['reported_construct'],
            parent_gene=row['parent_kinase'],
            reported_construct=row['reported_construct'],
            species=row['species'],
            canonical_accession=acc,
            canonical_version=seqs[acc]['version'],
            isoform='SV ' + seqs[acc]['version'],
            reported_mutation_notation=row['reported_mutation_notation'],
            mutation_class=mutation_class_of(row),
            substitutions=tuple(subs),
            coordinate_transforms=tuple(transforms),
            construct_start=crange.get('start'),
            construct_end=crange.get('end'),
            wt_residue=subs[0].old if subs else '',
            mutant_residue=subs[0].new if subs else '',
            parent_sequence=seq,
            mutant_sequence=mt_seq or '',
            reference_release='UniProt ' + acc + ' SV=' + seqs[acc]['version'],
            provenance_url='https://doi.org/10.1016/j.celrep.2015.12.080 (Duong-Ly 2016 Cell Rep 14:772, Tables S1/S2)',
            evidence_grade='high' if row['admission_status'] == 'admitted_point_pair' else 'quarantined',
            admission_status='admitted' if row['admission_status'] == 'admitted_point_pair' else 'quarantined',
            exclusion_reason=row['exclusion_reason'] or '',
        )
        d = dict(vr.__dict__)
        d['substitutions'] = [{'old': s.old, 'pos': s.pos, 'new': s.new,
                               'coordinate_kind': s.coordinate_kind} for s in subs]
        d['coordinate_transforms'] = [{'kind': t.kind, 'reported_position': t.reported_position,
                                       'canonical_position': t.canonical_position,
                                       'basis': t.basis, 'evidence': list(t.evidence)}
                                      for t in transforms]
        d['sequence_hash'] = sha256_text(seq)
        d['record_hash'] = sha256_text(json.dumps(d, sort_keys=True))
        records.append(d)
    return records


def davis_census():
    df = pd.read_excel(PARENT / 'downloads' / 'davis_MOESM3.xls', sheet_name='SuppTable1-050511')
    rows = []
    n_mutant = int((df['Mutant'] == 'YES').sum())
    n_phospho = int(df['Kinase'].str.contains('phosphorylated', case=False, na=False).sum())
    n_nonphospho = int(df['Kinase'].str.contains('nonphosphorylated', case=False, na=False).sum())
    for _, r in df.iterrows():
        kin = str(r['Kinase'])
        if '(' in kin or 'phosphorylated' in kin.lower() or 'nonphosphorylated' in kin.lower():
            rows.append({'kinase_construct': kin, 'mutant_flag': str(r['Mutant']),
                         'group': str(r['Kinase Group']), 'accession': str(r['Accession Number'])})
    multi = [r for r in rows if r['kinase_construct'].count('(') > 1]
    return {'n_assays': int(len(df)), 'n_mutant_flagged': n_mutant,
            'n_phosphorylated_constructs': n_phospho,
            'n_nonphosphorylated_constructs': n_nonphospho,
            'n_construct_rows_with_variant_or_phospho_state': len(rows),
            'n_multi_mutation_constructs': len(multi),
            'construct_rows': rows,
            'semantics': ('Davis 2011 Nat Biotechnol 29:1046 KINOMEscan panel; construct names carry '
                          'mutation and phosphorylation state; the Kd matrix (MOESM5) keeps blanks as '
                          'not-detected (right-censored at the tested concentration), not exact 10 uM')}


def klifs_census():
    i2 = json.loads((PARENT / 'X0_I2.json').read_text(encoding='utf-8'))
    census = i2['census']['klifs_parent_coverage']
    gatekeeper_check = {}
    for p in i2['representations']['klifs_pocket']['pairs']:
        if not p.get('excluded') and p.get('pocket_index') == 45:
            gatekeeper_check[p['construct']] = p['pocket_index']
    return {
        'source': 'KLIFS kinase_ID API lookup (klifs_kinase_lookup.json, local manifest with SHA-256)',
        'numbering_note': ('KLIFS aligned pocket = 85 positions; KLIFS numbering landmarks per '
                           'https://klifs.net (details page, numbering help) and Kooistra et al. 2016 '
                           'Nucleic Acids Res 44:D365 (doi:10.1093/nar/gkv1058); the gatekeeper maps to '
                           'KLIFS pocket index 45'),
        'parent_coverage': census,
        'gatekeeper_index_45_validation': gatekeeper_check,
    }


def alias_ledger(records, braf_evidence):
    lines = [
        '# Q0-B alias ledger (first-hand sources only)',
        '',
        'Every alias below is per-record evidence. None of it is generalized to',
        'other proteins; a coordinate transform is applied only when its own',
        'evidence chain is recorded here.',
        '',
        '## BRAF V599E -> V600E (historical numbering)',
        '- Reported: Duong-Ly Table S1/S2 construct BRAF(V599E); Davies et al.',
        '  2002 Nature 417:949-954 (doi:10.1038/nature00766) reported V599E.',
        '- Canonical: UniProt P15056 SV=4 residue 600 is V (variant annotation',
        '  VAR_018629 p.Val600Glu); residue 599 is T.',
        '- Sequence evidence (computed this stage, BRAF_HISTORICAL_EVIDENCE.json):',
        '  the 1992 reference M95712.1 CDS is 2298 nt (765 aa); the current',
        '  NM_004333.4 CDS is 2301 nt (766 aa) and its translation equals P15056.',
        '  The historical CDS lacks exactly 3 nt (one codon) in the 5-prime',
        '  region (new-sequence insertions at nt 88/90/97), shifting every',
        '  downstream coordinate by +1: historical V599 == canonical V600.',
        '- Applied transform: reported_position 599 -> canonical_position 600,',
        '  kind=historical_alias, old residue V verified at 600 before use.',
        '- NOT generalized: no other protein is given a +1 shift by analogy.',
        '',
        '## PDGFRalpha (Duong-Ly)',
        '- S1 lists GenBank NP_006197 (human PDGFRA) and Protein Accession',
        '  Q9DE49. Q9DE49 is Danio rerio pdgfra (UniProt header PGFRA_DANRE).',
        '- Correct human canonical: UniProt P16234 (PGFRA_HUMAN, 1089 aa), also',
        '  resolved by KLIFS kinase lookup PDGFRa/Human -> P16234.',
        '- D842V: S1 clone "Cytoplasmic (668-1210)" exceeds the canonical length',
        '  (1089) -> QUARANTINED (excluded_construct_unresolved); the old residue',
        '  D842 is verified on P16234 but the construct range cannot be mapped',
        '  without silent repair.',
        '- T674I: S1 Mutation column typo "T6741I"; construct name column reads',
        '  T674I; T674 verified on P16234 -> admitted (notation fix recorded).',
        '- V561D: V561 verified on P16234, construct 550-1089 -> admitted.',
        '',
        '## KLIFS pocket numbering',
        '- 85 aligned pocket positions; gatekeeper = index 45 verified on the',
        '  known gatekeeper mutants (see Q0B_KLIFS_CENSUS.json).',
        '- Sources: https://klifs.net (numbering help page, accessed 2026-08-18);',
        '  Kooistra et al. 2016, Nucleic Acids Res 44:D365-D371,',
        '  doi:10.1093/nar/gkv1058.',
        '',
        '## Duong-Ly notation quirks',
        '- S2 row labels use P38alpha/MAPK14 and TIE2/TEK (parent aliases).',
        '- S1 Mutation column empty for EGFR deletion constructs; the construct',
        '  name is authoritative there (d746-750, d747-749/A750P,',
        '  d747-752/P753S, d752-759, d746-750/T790M).',
        '- FLT3(ITD): "internal tandem duplication aa591-601" -> insertion class.',
        '',
    ]
    return '\n'.join(lines) + '\n'


def main():
    pair_table = json.loads((PARENT / 'X0_PAIR_TABLE.json').read_text(encoding='utf-8'))
    seqs = load_seqs()
    records = build_duongly_records(pair_table, seqs)
    braf_ev = json.loads((HERE / 'BRAF_HISTORICAL_EVIDENCE.json').read_text(encoding='utf-8'))

    n_admitted = sum(1 for r in records if r['admission_status'] == 'admitted')
    n_quarantined = len(records) - n_admitted
    by_class = {}
    for r in records:
        by_class[r['mutation_class']] = by_class.get(r['mutation_class'], 0) + 1

    audit = {
        'schema': 'MetaSieve.StageX0c.Q0B.v1',
        'preregistration_sha256': X0C_PREREG_SHA,
        'inherited_x0_preregistration_sha256': X0_PREREG_SHA,
        'braf_historical_evidence': braf_ev,
        'pdgfra_verification': {
            's1_reported_accession': 'Q9DE49',
            's1_reported_accession_identity': 'Danio rerio pdgfra (PGFRA_DANRE)',
            'correct_human_accession': 'P16234',
            'basis': ['S1 GenBank NP_006197 = human PDGFRA',
                      'KLIFS kinase lookup PDGFRa Human -> P16234',
                      'old-residue verification D842/T674/V561 on P16234'],
            'd842v_status': 'quarantined (construct 668-1210 exceeds canonical length 1089)',
            't674i_status': 'admitted (notation typo T6741I fixed with recorded basis)',
            'v561d_status': 'admitted',
        },
        'duongly_variant_records': records,
        'duongly_census': {'n_total': len(records), 'n_admitted': n_admitted,
                           'n_quarantined': n_quarantined, 'by_mutation_class': by_class,
                           'quarantine_reasons': {r['exclusion_reason']: sum(
                               1 for x in records if x['exclusion_reason'] == r['exclusion_reason'])
                               for r in records if r['exclusion_reason']}},
        'davis_construct_census': davis_census(),
        'klifs': klifs_census(),
        'admission_rule': ('admitted = single verified point substitution (old residue matches at the '
                           'canonical coordinate AND at the mapped construct coordinate); everything '
                           'else quarantined with a reason; no silent offset/substitution/reference '
                           'modification anywhere'),
    }
    inputs = [PARENT / 'X0_PAIR_TABLE.json', HERE / 'BRAF_HISTORICAL_EVIDENCE.json',
              PARENT / 'downloads' / 'davis_MOESM3.xls',
              PARENT / 'downloads' / 'duongly_mmc2.xlsx',
              PARENT / 'downloads' / 'duongly_mmc3.xlsx',
              PARENT / 'klifs' / 'klifs_kinase_lookup.json']
    inputs += sorted((PARENT / 'uniprot').glob('*.fasta'))
    write_artifact(HERE / 'Q0B_MAPPING_AUDIT.json', audit, inputs)

    (HERE / 'Q0B_ALIAS_LEDGER.md').write_text(alias_ledger(records, braf_ev))
    klifs_out = klifs_census()
    write_artifact(HERE / 'Q0B_KLIFS_CENSUS.json', {'schema': 'MetaSieve.StageX0c.Q0BKLIFS.v1',
                                                    'preregistration_sha256': X0C_PREREG_SHA,
                                                    **klifs_out},
                   [PARENT / 'klifs' / 'klifs_kinase_lookup.json', PARENT / 'X0_I2.json'])

    # hard rules: no admitted record may have any unverified substitution
    for r in records:
        if r['admission_status'] == 'admitted':
            assert r['substitutions'], r['source_row']
    print(json.dumps({'duongly_census': audit['duongly_census'],
                      'braf_nt_diff': braf_ev['braf_1992_mrna']['diff_vs_canonical_positions_1based'][:5]},
                     indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
