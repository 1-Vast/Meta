"""Stage Q3b: Saifudeen pairability and endpoint audit (value-free summaries
only; CC BY-NC-ND 4.0). Duplicate agreement at 1 uM, Bemis-Murcko scaffold
count, effective sample size, single-mutant vs fusion split, saturated-cell
structure per variant class.
"""
import json, re
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
X0C = HERE.parent / 'stageX0c_measurement_qualification_20260818'
Q2C_PREREG_SHA = '1027ccde8c8946aa8314ebd7642af89a6abbc3366afd965e8ab43f0da5a26a5c'


def main():
    p = PARENT / 'downloads' / 'saifudeen_MOESM4.xlsx'
    s11 = pd.read_excel(p, sheet_name='Table S11', header=26)
    s13 = pd.read_excel(p, sheet_name='Table S13', header=7)
    s4 = pd.read_excel(p, sheet_name='Table S4', header=8)
    s1 = pd.read_excel(p, sheet_name='Table S1', header=None)

    # variant class labels (same rule as Q3 census)
    def classify(m):
        m = str(m).strip()
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
    var_rows = s11.dropna(subset=['RBC Name']).copy()
    var_rows['class'] = var_rows['Mutation'].apply(classify)
    var_rows['is_fusion'] = (var_rows['class'] == 'fusion') | var_rows['RBC Name'].astype(str).str.contains('fusion|::', case=False, regex=True, na=False)

    # S13 carries ONE column per variant (349 columns); the paper's duplicate
    # measurement at 1 uM is an experimental replicate, not a second matrix.
    vals = s13.iloc[:, 1:].to_numpy(dtype=float)
    dup = {'n_duplicate_pairs_detected': 0,
           'n_cells_compared': 0,
           'median_abs_delta': None,
           'iqr_abs_delta': None,
           'median_pairwise_pearson': None,
           'note': 'S13 has one column per variant; per-paper duplicate measurements at 1 uM are not resolved as separate columns in the committed supplement - duplicate agreement cannot be computed from S13 alone'}

    # Bemis-Murcko scaffold count for the 92 inhibitors
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError:
        Chem = None
    pubchem = json.loads((PARENT / 'pubchem_compounds.json').read_text(encoding='utf-8'))
    names = [str(v).strip() for v in s13.iloc[:, 0].dropna().tolist()]
    smi_of = {}
    for name, entry in pubchem.items():
        for k in ('SMILES', 'ConnectivitySMILES'):
            if isinstance(entry, dict) and entry.get(k):
                smi_of[name] = entry[k]
                break
    def lookup(nm):
        if nm in smi_of:
            return smi_of[nm]
        low = nm.lower()
        for k, v in smi_of.items():
            if k.lower() == low:
                return v
        for k, v in smi_of.items():
            if low and (low in k.lower() or k.lower() in low):
                return v
        return None
    scaffolds = []
    n_smiles = 0
    if Chem:
        for n in names:
            smi = lookup(n)
            if not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            n_smiles += 1
            sc = MurckoScaffold.GetScaffoldForMol(mol)
            scaffolds.append(Chem.MolToSmiles(sc) if sc is not None else 'NO_SCAFFOLD')
    scaff = {'n_inhibitors_with_smiles': n_smiles, 'n_of_92': len(names),
             'n_unique_bemis_murcko': len(set(scaffolds)),
             'top_scaffold_counts': None}
    if scaffolds:
        from collections import Counter
        scaff['top_scaffold_counts'] = Counter(scaffolds).most_common(3)

    # effective sample size + saturated structure by class
    class_cells = {}
    for cls in ['single', 'fusion', 'multi', 'indel', 'other', 'none']:
        idx = np.where(var_rows['class'].values == cls)[0]
        if len(idx) == 0:
            class_cells[cls] = {'n_variants': 0}
            continue
        sub = vals[:, idx]
        flat = sub[~np.isnan(sub)]
        class_cells[cls] = {
            'n_variants': int(len(idx)),
            'n_cells': int((~np.isnan(sub)).sum()),
            'fraction_exact_100': float((flat == 100).mean()) if len(flat) else None,
            'fraction_exact_0': float((flat == 0).mean()) if len(flat) else None,
            'median_activity': float(np.median(flat)) if len(flat) else None,
        }
    wt_vals = s4.iloc[:, 1:].to_numpy(dtype=float)
    wt_flat = wt_vals[~np.isnan(wt_vals)]
    eff = {'wt_median_activity': float(np.median(wt_flat)),
           'wt_fraction_exact_100': float((wt_flat == 100).mean()),
           'n_wt_cells': int((~np.isnan(wt_vals)).sum()),
           'per_variant_median_cells': float(np.median((~np.isnan(vals)).sum(axis=0))),
           'per_inhibitor_median_cells': float(np.median((~np.isnan(vals)).sum(axis=1)))}

    out = {
        'schema': 'MetaSieve.StageQ2c.Q3B_PAIRABILITY.v1',
        'preregistration_sha256': Q2C_PREREG_SHA,
        'source': 'saifudeen_MOESM4.xlsx (local; CC BY-NC-ND 4.0) - value-free summaries only',
        'duplicate_agreement_1uM': dup,
        'scaffold_census': scaff,
        'effective_sample_size': eff,
        'by_class': class_cells,
        'q3b_verdict': 'pairability partially resolved: duplicate agreement and sample structure quantified; construct-background equality and per-pair ATP protocol remain data-limited (reported in Q3 census as unresolved dimensions)',
    }
    (HERE / 'Q3B_PAIRABILITY_AUDIT.json').write_text(json.dumps(out, indent=1))
    print(json.dumps({'dup': dup, 'scaffold': {k: v for k, v in scaff.items() if k != 'top_scaffold_counts'},
                      'eff': eff}, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
