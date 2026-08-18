"""X0 round-2 data fetch: P16234 fasta, KLIFS pocket sequences, PubChem SMILES.
Records every fetched asset with URL, access date and SHA-256 into
x0_fetch_manifest.json. Network-available; idempotent.
"""
from __future__ import annotations
import hashlib, json, time, urllib.request, urllib.parse
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
PREREG_SHA = '03cdc907df3e778f5fe79fb1a238d35ebb6ece5e9e743db181728ba6b25e9683'
OUT = HERE / 'x0_fetch_manifest.json'
UNIPROT = HERE / 'uniprot'
KLIFS = HERE / 'klifs'
KLIFS.mkdir(exist_ok=True)

UA = {'User-Agent': 'MetaSieve-StageX0/2.0 (research audit; single-pass fetch)'}

def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()

def http_get(url, tries=3):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception as e:
            last = e
            time.sleep(2 * (i + 1))
    raise last

manifest = {'schema': 'MetaSieve.StageX.X0FetchManifest.v1',
            'preregistration_sha256': PREREG_SHA,
            'fetched_at': datetime.now(timezone.utc).isoformat(),
            'entries': []}
def record(kind, name, url, data, note=''):
    digest = sha256_bytes(data)
    manifest['entries'].append({'kind': kind, 'name': name, 'url': url,
                                'sha256': digest, 'bytes': len(data), 'note': note})
    return digest

# 1. P16234 human PDGFRA fasta
try:
    url = 'https://rest.uniprot.org/uniprotkb/P16234.fasta'
    data = http_get(url)
    p = UNIPROT / 'P16234.fasta'
    p.write_bytes(data)
    record('uniprot_fasta', 'P16234.fasta', url, data,
           'human PDGFRA; correct accession for Duong-Ly PDGFRalpha rows (S1 lists Q9DE49 = zebrafish)')
    print('P16234 fetched, header:', data.splitlines()[0].decode())
except Exception as e:
    print('P16234 fetch failed:', e)

# 2. KLIFS kinase_ID + pocket for the 21 parents
kinases = ['ABL1','ALK','BRAF','BTK','KIT','MET','SRC','CHK2','EGFR','FGFR1','FGFR2',
           'FGFR3','FGFR4','FLT3','JAK2','LRRK2','MAP2K1','MAPK14','PDGFRa','RET','TEK']
klifs_rows = {}
for k in kinases:
    try:
        url = 'https://klifs.net/api/kinase_ID?kinase_name=' + urllib.parse.quote(k) + '&species=Human'
        data = http_get(url)
        obj = json.loads(data)
        klifs_rows[k] = obj
        record('klifs_kinase_json', f'klifs_{k}.json', url, data,
               f'KLIFS kinase_ID lookup for {k}; {len(obj)} match(es)')
        print(k, '->', [(m.get('kinase_ID'), m.get('pocket','')[:12]) for m in obj])
        time.sleep(0.4)
    except Exception as e:
        print('KLIFS', k, 'failed:', e)
        klifs_rows[k] = None
        record('klifs_kinase_json', f'klifs_{k}.json', url, b'',
               f'FETCH FAILED: {e}')
with open(KLIFS / 'klifs_kinase_lookup.json', 'w', encoding='utf-8') as f:
    json.dump(klifs_rows, f, indent=1)

# 3. PubChem SMILES + InChIKey for the 183 Duong-Ly compounds
matrix = pd.read_excel(HERE / 'downloads/duongly_mmc3.xlsx', sheet_name='Table S2')
compounds = [str(c).strip() for c in matrix.columns[1:]]
print('n compounds:', len(compounds))
pubchem = {}
for i, name in enumerate(compounds):
    try:
        url = ('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
               + urllib.parse.quote(name) + '/property/CanonicalSMILES,InChIKey,IsomericSMILES/JSON')
        data = http_get(url, tries=2)
        obj = json.loads(data)
        props = obj['PropertyTable']['Properties'][0]
        pubchem[name] = props
        record('pubchem_json', f'pc_{i:03d}.json', url, data, f'lookup for compound {name!r}')
        print(i, repr(name), '->', (props.get('InChIKey') or '')[:27], (props.get('CanonicalSMILES') or '')[:24])
        time.sleep(0.35)
    except Exception as e:
        pubchem[name] = {'error': str(e)}
        print(i, repr(name), 'FAILED', e)
with open(HERE / 'pubchem_compounds.json', 'w', encoding='utf-8') as f:
    json.dump(pubchem, f, indent=1)

OUT.write_text(json.dumps(manifest, indent=1) + '\n')
print('manifest written:', OUT, len(manifest['entries']), 'entries')
