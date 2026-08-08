"""Resolve Klaeger 2017 drug names to structures via PubChem PUG-REST.

Label-blind: reads drug NAMES from the pinned matrix header only.  No affinity
value is read by this script.  Results are cached and hashed.
"""
from __future__ import annotations

import hashlib
import json
import os
import ssl
import time
import urllib.parse
import urllib.request

import pandas as pd

RAW = r"D:\MetaSieve\dataset\raw\crossed_panels\kinase_panels"
CACHE = r"D:\MetaSieve\dataset\processed\crossed_panels_xp2"
os.makedirs(CACHE, exist_ok=True)
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
UA = {"User-Agent": "Mozilla/5.0"}
BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/{}/TXT"


def fetch(name, prop="CanonicalSMILES"):
    url = BASE.format(urllib.parse.quote(name), prop)
    req = urllib.request.Request(url, headers=UA)
    txt = urllib.request.urlopen(req, timeout=45, context=CTX).read().decode().strip()
    return txt.splitlines()[0].strip() if txt else None


def main():
    hdr = pd.read_csv(os.path.join(RAW, "klaeger_matrix.csv"), nrows=0)
    names = pd.read_csv(os.path.join(RAW, "klaeger_matrix.csv"),
                        usecols=[hdr.columns[0]]).iloc[:, 0].astype(str).str.strip().tolist()
    print(f"drug names read from the pinned matrix header column: {len(names)}")
    path = os.path.join(CACHE, "klaeger_smiles.json")
    got = json.load(open(path)) if os.path.exists(path) else {}
    todo = [n for n in names if n not in got]
    print(f"to resolve: {len(todo)}")
    for i, n in enumerate(todo):
        for cand in (n, n.replace("-", ""), n.split("(")[0].strip()):
            try:
                s = fetch(cand)
                if s:
                    got[n] = {"smiles": s, "query": cand}
                    break
            except Exception:
                pass
            time.sleep(0.25)
        else:
            got[n] = None
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(todo)} resolved={sum(1 for v in got.values() if v)}")
            json.dump(got, open(path, "w"))
        time.sleep(0.2)
    json.dump(got, open(path, "w"))
    ok = sum(1 for v in got.values() if v)
    print(f"resolved {ok}/{len(names)}")
    blob = json.dumps({k: (v or {}).get("smiles") for k, v in sorted(got.items())})
    print("resolution_sha256:", hashlib.sha256(blob.encode()).hexdigest())
    json.dump({"n_names": len(names), "n_resolved": ok, "source": "PubChem PUG-REST",
               "resolution_sha256": hashlib.sha256(blob.encode()).hexdigest(),
               "affinity_values_read": 0},
              open(os.path.join(CACHE, "klaeger_structure_manifest.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
