"""Acquire the release-pinned BindingDB curated-articles subset.

License: BindingDB is distributed under CC BY 3.0.  Acquisition only; no
affinity value is read by this script.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import os
import ssl
import urllib.request

OUT = r"D:\MetaSieve\dataset\raw\crossed_panels\bindingdb"
os.makedirs(OUT, exist_ok=True)
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
UA = {"User-Agent": "Mozilla/5.0"}

RELEASE = "202608"
CANDIDATES = [
    f"https://www.bindingdb.org/rwd/bind/downloads/BindingDB_BindingDB_Articles_{RELEASE}_tsv.zip",
    f"https://www.bindingdb.org/bind/downloads/BindingDB_BindingDB_Articles_{RELEASE}_tsv.zip",
    f"https://www.bindingdb.org/rwd/bind/downloads/BindingDB_Assays_{RELEASE}_tsv.zip",
]


def get(url, dest):
    req = urllib.request.Request(url, headers=UA)
    r = urllib.request.urlopen(req, timeout=900, context=CTX)
    data = r.read()
    with open(dest, "wb") as f:
        f.write(data)
    return len(data), hashlib.sha256(data).hexdigest()


man = {"acquired_utc": datetime.datetime.utcnow().isoformat() + "Z",
       "release": RELEASE, "license": "CC BY 3.0 (BindingDB)",
       "affinity_values_read_by_this_script": 0, "files": {}}
for url in CANDIDATES:
    name = url.rsplit("/", 1)[-1]
    dest = os.path.join(OUT, name)
    if os.path.exists(dest):
        print("cached", name)
        continue
    try:
        n, h = get(url, dest)
        man["files"][name] = {"url": url, "bytes": n, "sha256": h}
        print(f"{name:56s} {n:10d} bytes sha256={h[:16]}...")
    except Exception as e:
        print(f"FAIL {name}: {type(e).__name__}: {str(e)[:80]}")
json.dump(man, open(os.path.join(OUT, "acquisition_manifest.json"), "w"), indent=2)
print("manifest written")
