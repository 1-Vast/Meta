"""Stage W0b censored-support re-census (no training)."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import sys
if __package__ in {None, ''}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from tools.research.stageW0b_core1_audit.w0b_audit import (
    PREREG_SHA, cdhit_components, layer_census, read_davis, read_klifs,
    read_klaeger, read_metz,
)
HERE = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'W0B_CENSORED_RECENSUS.json')
    args = parser.parse_args()
    out = {'schema':'MetaSieve.StageW0b.CensoredRecensus.v1',
           'stage':'stageW0b_core1_audit',
           'preregistration_sha256':PREREG_SHA,
           'rule':'rows flagged at detection floor excluded before pair formation; no censored value enters any layer statistic',
           'datasets':{}}
    # Davis
    targets, rows = read_davis()
    comps = cdhit_components(targets, 'davis')
    clean = {t:[r for r in rs if not r['censored']] for t,rs in rows.items() if any(not r['censored'] for r in rs)}
    clean = {t:rs for t,rs in clean.items() if len({r['smiles'] for r in rs})>=2}
    out['datasets']['davis']={'targets_with_pairs':len(clean),'components':len({comps[t] for t in clean}),'rows':sum(len(v) for v in clean.values()),'layers':layer_census(clean,comps)}
    # Metz
    metz_rows = read_metz()
    klifs,_ = read_klifs()
    mcomps={t:klifs.get(t.upper(),{}).get('group','UNMAPPED') for t in metz_rows}
    mclean={t:[r for r in rs if not r['censored']] for t,rs in metz_rows.items() if any(not r['censored'] for r in rs)}
    mclean={t:rs for t,rs in mclean.items() if len({r['smiles'] for r in rs})>=2}
    out['datasets']['metz']={'targets_with_pairs':len(mclean),'components':len({mcomps[t] for t in mclean}),'rows':sum(len(v) for v in mclean.values()),'layers':layer_census(mclean,mcomps)}
    # Klaeger
    klaeger_rows = read_klaeger()
    kcomps={t:klifs.get(t.upper(),{}).get('group','UNMAPPED') for t in klaeger_rows}
    kclean={t:[r for r in rs if not r['censored']] for t,rs in klaeger_rows.items() if any(not r['censored'] for r in rs)}
    kclean={t:rs for t,rs in kclean.items() if len({r['smiles'] for r in rs})>=2}
    out['datasets']['klaeger']={'targets_with_pairs':len(kclean),'components':len({kcomps[t] for t in kclean}),'rows':sum(len(v) for v in kclean.values()),'layers':layer_census(kclean,kcomps)}
    args.output.write_text(json.dumps(out,indent=1,sort_keys=True)+'\n')
    print(json.dumps(out,indent=1))
    return 0
if __name__=='__main__':
    raise SystemExit(main())
