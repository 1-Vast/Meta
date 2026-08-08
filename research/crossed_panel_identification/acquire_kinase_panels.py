"""Acquire independent crossed kinase panels.

Explicitly EXCLUDED by MetaSieve governance:
  - davis_affinity.csv / davis_proteins.csv  (DAVIS access prohibited)
  - anastassiadis_matrix.csv / anastassiadis.xls (consumed development panel)
"""
import urllib.request, ssl, hashlib, os, json, datetime

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

REPO = "polinavino/kinase-selectivity-definitions"
SHA = "8ab79cae31c18e49007dcce6dd11f93d2667ab14"
BASE = f"https://raw.githubusercontent.com/{REPO}/{SHA}/"

FILES = [
    "metz.xls",                 # Metz 2011 Nat Chem Biol original supplement
    "metz_matrix.csv",          # derived kinase x compound pKi matrix
    "aan4368_Table_S2.xlsx",    # Klaeger 2017 Science Table S2 original
    "klaeger_matrix.csv",       # derived kinase x drug apparent Kd matrix
]

OUT = r"D:\MetaSieve\dataset\raw\crossed_panels\kinase_panels"
os.makedirs(OUT, exist_ok=True)

manifest = {
    "acquired_utc": datetime.datetime.utcnow().isoformat() + "Z",
    "mirror_repo": REPO,
    "mirror_commit": SHA,
    "excluded_by_governance": [
        "davis_affinity.csv", "davis_proteins.csv",
        "anastassiadis_matrix.csv", "anastassiadis.xls",
    ],
    "files": {},
}

for name in FILES:
    url = BASE + name
    req = urllib.request.Request(url, headers=UA)
    r = urllib.request.urlopen(req, timeout=300, context=ctx)
    data = r.read()
    h = hashlib.sha256(data).hexdigest()
    p = os.path.join(OUT, name)
    with open(p, "wb") as f:
        f.write(data)
    manifest["files"][name] = {"url": url, "bytes": len(data), "sha256": h}
    print(f"{name:28s} {len(data):9d} bytes  sha256={h}")

with open(os.path.join(OUT, "acquisition_manifest.json"), "w") as f:
    json.dump(manifest, f, indent=2)
print("manifest written")
