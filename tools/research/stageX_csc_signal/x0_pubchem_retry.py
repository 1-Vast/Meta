"""Round-2 pass: retry the 45 failed PubChem lookups with cleaned names.
Idempotent; updates pubchem_compounds.json and the fetch manifest.
"""
from __future__ import annotations
import json, re, time, urllib.request, urllib.parse, hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
UA = {'User-Agent': 'MetaSieve-StageX0/2.0 (research audit; retry pass)'}

def http_get(url, tries=2):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=45) as r:
                return r.read()
        except Exception as e:
            last = e
            time.sleep(1.5 * (i + 1))
    raise last

pc = json.loads((HERE / 'pubchem_compounds.json').read_text(encoding='utf-8'))
def clean(name):
    s = re.sub(r'\s+', ' ', str(name)).strip()
    return s

candidates = {}
for name, v in pc.items():
    if 'InChIKey' in v and v.get('InChIKey'):
        continue
    cands = []
    c = clean(name)
    cands.append(c)
    parts = [p.strip() for p in c.split(',') if p.strip()]
    if len(parts) > 1:
        cands.append(parts[-1])
        cands.append(parts[0])
    cands = list(dict.fromkeys(cands))
    candidates[name] = cands

resolved, still_fail = 0, []
manifest = json.loads((HERE / 'x0_fetch_manifest.json').read_text(encoding='utf-8'))
for name, cands in candidates.items():
    ok = False
    for cand in cands:
        try:
            url = ('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
                   + urllib.parse.quote(cand) + '/property/CanonicalSMILES,InChIKey,IsomericSMILES/JSON')
            data = http_get(url)
            props = json.loads(data)['PropertyTable']['Properties'][0]
            pc[name] = props
            pc[name]['_name_used'] = cand
            manifest['entries'].append({'kind': 'pubchem_json', 'name': f'retry:{name}',
                                        'url': url, 'sha256': hashlib.sha256(data).hexdigest(),
                                        'bytes': len(data),
                                        'note': f'retry lookup via cleaned name {cand!r}'})
            resolved += 1
            ok = True
            break
        except Exception:
            continue
    if not ok:
        still_fail.append((name, cands))
    time.sleep(0.3)

(HERE / 'pubchem_compounds.json').write_text(json.dumps(pc, indent=1) + '\n')
(HERE / 'x0_fetch_manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
print('resolved now:', resolved)
print('still failing:', len(still_fail))
for name, cands in still_fail[:20]:
    print('  ', repr(name), 'tried', cands)
