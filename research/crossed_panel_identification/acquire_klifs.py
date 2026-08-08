import urllib.request, ssl, hashlib, os, json, datetime

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
UA = {"User-Agent": "Mozilla/5.0"}
OUT = r"D:\MetaSieve\dataset\raw\crossed_panels\protein_annotation"
os.makedirs(OUT, exist_ok=True)

man = {"acquired_utc": datetime.datetime.utcnow().isoformat() + "Z", "files": {}}
for name, url in [
    ("klifs_kinase_information_human.json", "https://klifs.net/api/kinase_information?species=HUMAN"),
    ("klifs_kinase_groups.json", "https://klifs.net/api/kinase_groups?species=HUMAN"),
]:
    req = urllib.request.Request(url, headers=UA)
    data = urllib.request.urlopen(req, timeout=180, context=ctx).read()
    h = hashlib.sha256(data).hexdigest()
    with open(os.path.join(OUT, name), "wb") as f:
        f.write(data)
    man["files"][name] = {"url": url, "bytes": len(data), "sha256": h}
    print(f"{name:42s} {len(data):8d} sha256={h[:16]}...")

rec = json.loads(open(os.path.join(OUT, "klifs_kinase_information_human.json"), encoding="utf-8").read())
print("human kinases:", len(rec))
print("with 85-mer pocket:", sum(1 for r in rec if r.get("pocket") and len(r["pocket"]) == 85))
plens = {}
for r in rec:
    plens[len(r.get("pocket") or "")] = plens.get(len(r.get("pocket") or ""), 0) + 1
print("pocket length histogram:", dict(sorted(plens.items())))
print("groups:", sorted({r["group"] for r in rec}))
print("families:", len({r["family"] for r in rec}))

with open(os.path.join(OUT, "acquisition_manifest.json"), "w") as f:
    json.dump(man, f, indent=2)
