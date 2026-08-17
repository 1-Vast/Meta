"""Stage X0-D: record acquired external data with provenance and semantics."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PREREG_SHA = '03cdc907df3e778f5fe79fb1a238d35ebb6ece5e9e743db181728ba6b25e9683'
D = HERE / 'downloads'

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def shape_info(path, sheet=None):
    try:
        if path.suffix=='.xlsx' or path.suffix=='.xls':
            xl=pd.ExcelFile(path)
            sheet=sheet or xl.sheet_names[0]
            df=xl.parse(sheet)
            return {'sheet':sheet,'shape':list(df.shape),'columns':[str(c) for c in df.columns[:10]],'non_null_fraction':float(df.notna().mean().mean())}
        if path.suffix=='.file':
            head=path.read_bytes()[:4]
            return {'magic':head.hex(),'bytes':path.stat().st_size}
    except Exception as e:
        return {'error':str(e)}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,default=HERE/'X0D_DATA_AUDIT.json'); args=ap.parse_args()
    files={
     'duongly_2016':[
       ('duongly_mmc2.xlsx','Table S1 construct/mutation metadata','Cell Reports CC BY-NC-ND 4.0; no derivative redistribution; local cache only','https://ars.els-cdn.com/content/image/1-s2.0-S2211124715015363-mmc2.xlsx','Table S1'),
       ('duongly_mmc3.xlsx','Table S2 183 inhibitors x 76 mutant kinases, % remaining activity','CC BY-NC-ND 4.0','https://ars.els-cdn.com/content/image/1-s2.0-S2211124715015363-mmc3.xlsx','Table S2'),
       ('duongly_mmc4.xlsx','Table S3 platform comparison','CC BY-NC-ND 4.0','https://ars.els-cdn.com/content/image/1-s2.0-S2211124715015363-mmc4.xlsx','Table S3'),
       ('duongly_mmc5.xlsx','Table S5 WT/mutant comparison','CC BY-NC-ND 4.0','https://ars.els-cdn.com/content/image/1-s2.0-S2211124715015363-mmc5.xlsx','Table S5'),
       ('duongly_mmc6.xlsx','Table S6 compound list','CC BY-NC-ND 4.0','https://ars.els-cdn.com/content/image/1-s2.0-S2211124715015363-mmc6.xlsx','Table S6'),
     ],
     'anastassiadis_2011':[
       ('anastassiadis_MOESM23.xls','Supplementary Table 3 complete pairwise kinase-compound % remaining activity, 178 compounds x 300 kinases','Nature Biotechnology supplementary; local analysis only','https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.2017/MediaObjects/41587_2011_BFnbt2017_MOESM23_ESM.xls','Sheet1'),
       ('anastassiadis_MOESM22.pdf','Supplementary text and figures + tables 1,2,4,5','Nature Biotechnology supplementary','https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.2017/MediaObjects/41587_2011_BFnbt2017_MOESM22_ESM.pdf',None),
     ],
     'davis_2011':[
       ('davis_MOESM3.xls','Table S1 kinase assay list with mutant and phospho-state annotations','Nature Biotechnology supplementary','https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.1990/MediaObjects/41587_2011_BFnbt1990_MOESM3_ESM.xls','SuppTable1-050511'),
       ('davis_MOESM5.xls','Full Kd matrix 72 compounds x 442 kinase assays, blanks preserved (tested, Kd>10 uM or not detected)','Nature Biotechnology supplementary','https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.1990/MediaObjects/41587_2011_BFnbt1990_MOESM5_ESM.xls','Sheet1'),
       ('davis_MOESM6.xls','Reference compound Kd list and binding-mode notes','Nature Biotechnology supplementary','https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.1990/MediaObjects/41587_2011_BFnbt1990_MOESM6_ESM.xls','RefCmpds73Kdsvs442Complete10uMN'),
     ],
     'pkis2_2017':[
       ('pkis2_s004.file','Supplementary data archive (likely xlsx/zip)','PLOS ONE CC BY 4.0 / public domain release','https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0181585.s004&type=supplementary',None),
     ],
     'local_reanalysis_sources':[
       ('metz_matrix.csv','Metz 2011 pKi matrix 704 compounds x 172 kinases (local)','local raw','dataset/raw/crossed_panels/kinase_panels/metz_matrix.csv',None),
       ('klaeger_matrix.csv','Klaeger 2017 Kinobeads apparent pKd matrix 222 drugs x 343 proteins (local)','local raw','dataset/raw/crossed_panels/kinase_panels/klaeger_matrix.csv',None),
       ('kiba.tab','KIBA score matrix (local)','local raw','dataset/raw/dta/kiba.tab',None),
     ],
    }
    out={'schema':'MetaSieve.StageX.X0DDataAudit.v1','stage':'stageX_csc_signal','preregistration_sha256':PREREG_SHA,'sources':{}}
    for source,entries in files.items():
        out['sources'][source]=[]
        for fname,semantics,license_note,url,sheet in entries:
            p=D/fname if source not in ('local_reanalysis_sources',) else ROOT/url
            exists=p.exists()
            out['sources'][source].append({'file':fname,'semantics':semantics,'license_note':license_note,'url':url,'exists':exists,'sha256':sha(p) if exists else None,'info':shape_info(p,sheet) if exists else None})
    args.output.write_text(json.dumps(out,indent=1,sort_keys=True)+'\n')
    print(json.dumps(out,indent=1))
    return 0
if __name__=='__main__':
    raise SystemExit(main())
