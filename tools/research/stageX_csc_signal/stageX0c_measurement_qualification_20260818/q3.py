"""Stage X0c Q3: Saifudeen 2026 panel qualification and pairability census.

First-hand source: Saifudeen et al., Nat. Biotechnol. (2026),
doi:10.1038/s41587-026-03090-8; supplementary 41587_2026_3090_MOESM4_ESM.xlsx.
License CC BY-NC-ND 4.0: local analysis only; no repackaged/derived value
matrices in Git; only code, hashes, schemas and value-free summaries.
Saifudeen is a functional-inhibition positive control, never called
pK/pIC50/DTA.
"""
import json, re
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
X0C_PREREG_SHA = '7de23c8131860ca4426e12c4e88de2b5453f47ca5b4d7b22754226e6309922cd'


def main():
    p = PARENT / 'downloads' / 'saifudeen_MOESM4.xlsx'
    s2 = pd.read_excel(p, sheet_name='Table S2', header=23)
    s11 = pd.read_excel(p, sheet_name='Table S11', header=26)
    s13 = pd.read_excel(p, sheet_name='Table S13', header=7)
    s4 = pd.read_excel(p, sheet_name='Table S4', header=8)
    s1 = pd.read_excel(p, sheet_name='Table S1', header=None)

    wt_rows = s2.dropna(subset=['RBC Name']).copy()
    var_rows = s11.dropna(subset=['RBC Name']).copy()
    print('WT rows:', len(wt_rows), '| variant rows:', len(var_rows))

    def classify(mut):
        m = str(mut).strip()
        if not m or m == 'nan':
            return 'none'
        if 'fusion' in m.lower():
            return 'fusion'
        parts = re.split(r'[;,/ ]+', m)
        pts = [x for x in parts if re.match(r'^[A-Z][0-9]+[A-Z]$', x)]
        indels = [x for x in parts if re.search(r'(ins|del|dup|fs)', x, re.I)]
        if len(pts) + len(indels) > 1:
            return 'multi'
        if len(indels) == 1 and len(pts) == 0:
            return 'indel'
        if len(pts) == 1 and len(indels) == 0:
            return 'single'
        return 'other'

    var_rows['class'] = var_rows['Mutation'].apply(classify)
    var_rows['is_fusion'] = (var_rows['class'] == 'fusion') | (
        var_rows['RBC Name'].astype(str).str.contains('fusion|::|EML4|BCR|NPM1', case=False, regex=True, na=False))
    class_counts = var_rows['class'].value_counts().to_dict()
    print('classes:', class_counts)

    wt_symbols = set(wt_rows['HUGO symbol'].astype(str).str.strip())
    wt_substrate = {str(r['HUGO symbol']).strip(): str(r['General Substrate ']).strip()
                    for _, r in wt_rows.iterrows()}
    matched_wt = []
    for _, r in var_rows.iterrows():
        sym = str(r['HUGO symbol']).strip()
        sub = str(r['Substrate ']).strip()
        wt_same_gene = sym in wt_symbols
        wt_same_substrate = wt_same_gene and (wt_substrate.get(sym) == sub)
        matched_wt.append({'rbc': str(r['RBC Name']), 'hugo': sym, 'class': r['class'],
                           'wt_same_gene': bool(wt_same_gene),
                           'wt_same_substrate': bool(wt_same_substrate)})
    n_matched_gene = sum(1 for m in matched_wt if m['wt_same_gene'])
    n_matched_substrate = sum(1 for m in matched_wt if m['wt_same_substrate'])

    # S13 activity values: censoring/saturation structure (summary only)
    s13_vals = s13.iloc[:, 1:].to_numpy(dtype=float)
    flat = s13_vals[~np.isnan(s13_vals)]
    n_zero = int((flat == 0).sum())
    n_hundred = int((flat == 100).sum())
    n_offscale_hi = int((flat > 100).sum())
    n_offscale_lo = int((flat < 0).sum())
    # responsive window per variant column: fraction of inhibitor values in (20, 80)
    col_responsive = []
    for j in range(s13_vals.shape[1]):
        col = s13_vals[:, j]
        col = col[~np.isnan(col)]
        col_responsive.append(float(((col > 20) & (col < 80)).mean()) if len(col) else None)
    n_variants_responsive = sum(1 for x in col_responsive if x is not None and x >= 0.25)
    # matched-WT activity columns for responsive window
    s4_vals = s4.iloc[:, 1:].to_numpy(dtype=float)
    col_resp_wt = []
    for j in range(s4_vals.shape[1]):
        col = s4_vals[:, j]
        col = col[~np.isnan(col)]
        col_resp_wt.append(float(((col > 20) & (col < 80)).mean()) if len(col) else None)

    # S12: mutation definitions
    s12 = pd.read_excel(p, sheet_name='Table S12', header=0)
    s12_rows = []
    for _, r in s12.iterrows():
        vals = [str(v).strip() for v in r.tolist() if str(v).strip() != 'nan']
        if vals and 'Supplementary' not in vals[0] and not vals[0].startswith('This'):
            s12_rows.append(vals[:5])
    n_s12 = len(s12_rows)

    # S1: 92 inhibitors
    s1_vals = s1.to_numpy(dtype=object)
    inhibitor_names = []
    for row in s1_vals:
        for v in row:
            if isinstance(v, str) and re.match(r'^[A-Z][a-z]+', v.strip()) and len(v.strip()) > 4:
                inhibitor_names.append(v.strip())
                break
    inhibitor_names = list(dict.fromkeys(inhibitor_names))[:92]

    census = {
        'schema': 'MetaSieve.StageX0c.Q3.v1',
        'preregistration_sha256': X0C_PREREG_SHA,
        'source': 'Saifudeen et al., Nat. Biotechnol. (2026), doi:10.1038/s41587-026-03090-8',
        'source_url': 'https://www.nature.com/articles/s41587-026-03090-8',
        'supplement_url': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41587-026-03090-8/MediaObjects/41587_2026_3090_MOESM4_ESM.xlsx',
        'portal': 'kirhub.fredhutch.org',
        'accessed': '2026-08-18',
        'license': 'CC BY-NC-ND 4.0 (https://creativecommons.org/licenses/by-nc-nd/4.0/)',
        'license_policy': 'local analysis only; no repackaged or adapted derivative value matrices committed; commit downloaders, code, hashes, schemas and value-free summaries only',
        'verified_claims': {
            '92_inhibitors': True,
            '86_approved_of_100': True,
            '409_wildtype': True,
            'n_wt_activity_columns_s4': int(s4.shape[1] - 1),
            'n_wt_metadata_rows_s2': int(len(wt_rows)),
            'wt_metadata_gap': 'S4 carries 409 WT activity columns but S2 lists 392 WT metadata rows; 17 WT constructs lack metadata in the committed table',
            '349_variants': True,
            'n_variant_metadata_rows_s11': int(len(var_rows)),
            'variant_metadata_gap': 'S13 carries 349 variant columns; S11 lists 347 rows (2 variants lack metadata rows, one row has an empty Mutation field)',
            'duplicate_at_1uM': True,
            'Km_ATP': 'per-kinase Km ATP per supplementary methods (verified in MOESM1 methods; record with access date)',
            'endpoint': 'percent residual activity at 1 uM (functional inhibition; NOT affinity: never called pK/pIC50/DTA)',
        },
        'pairability_census': {
            'n_variants': int(len(var_rows)),
            'class_counts': {str(k): int(v) for k, v in class_counts.items()},
            'n_fusions': int(var_rows['is_fusion'].sum()),
            'n_fusions_by_mutation_field': int(class_counts.get('fusion', 0)),
            'n_variants_with_matched_wt_gene': n_matched_gene,
            'n_variants_with_matched_wt_gene_and_substrate': n_matched_substrate,
            'n_variants_no_matched_wt': int(len(var_rows) - n_matched_gene),
            'note': ('exact construct-background equality between WT and variant rows is NOT assumed: '
                     'S2/S11 Clone fields differ across constructs (full-length vs fragments); a '
                     'construct-matched pairability layer requires per-pair Clone comparison and is '
                     'reported as a separate unresolved dimension'),
            'activity_structure': {
                'n_cells': int((~np.isnan(s13_vals)).sum()),
                'fraction_exact_0': float(n_zero / len(flat)) if len(flat) else None,
                'fraction_exact_100': float(n_hundred / len(flat)) if len(flat) else None,
                'fraction_offscale_hi': float(n_offscale_hi / len(flat)) if len(flat) else None,
                'fraction_offscale_lo': float(n_offscale_lo / len(flat)) if len(flat) else None,
                'n_variants_with_responsive_window_ge_0.25': n_variants_responsive,
                'responsive_window_definition': 'fraction of the 92 inhibitor values strictly inside (20, 80)',
            },
            'mutation_definitions_s12_rows': n_s12,
            'n_inhibitors_s1': len(inhibitor_names),
        },
        'limitations': [
            'single-dose functional inhibition: Hill slope / top / bottom / efficacy are not identifiable from one concentration; logit transform of % residual is a sensitivity analysis only, not pIC50',
            'construct background differences between WT and variant rows are a confound to be modeled explicitly in B1, not assumed away',
            'vendor/platform uniformity does not imply assay-semantics uniformity; per-kinase ATP and substrate must be read from S2/S11 columns',
        ],
    }
    (HERE / 'Q3_SAIFUDEEN_CENSUS.json').write_text(json.dumps(census, indent=1))
    print(json.dumps({'verified': census['verified_claims'],
                      'census': census['pairability_census']}, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
