"""S1b — acquire an independent RCSB CC0 structural set disjoint from P1B exposure."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import ssl
import time
import urllib.parse
import urllib.request

import numpy as np

ROOT = r"D:\MetaSieve"
OUT = os.path.join(ROOT, "dataset", "raw", "ssl_b2_independent")
MMCIF = os.path.join(OUT, "mmcif")
os.makedirs(MMCIF, exist_ok=True)
REPORT = os.path.join(ROOT, "report", "ssl_b2_structural_observability")
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
UA = {"User-Agent": "Mozilla/5.0"}
TARGET_N = 1500
SEED = 20260808

QUERY = {
    "query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.experimental_method",
            "operator": "exact_match", "value": "X-ray"}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.resolution_combined",
            "operator": "less_or_equal", "value": 2.5}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.nonpolymer_entity_count",
            "operator": "greater", "value": 0}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.polymer_entity_count_protein",
            "operator": "equals", "value": 1}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_accession_info.initial_release_date",
            "operator": "greater", "value": "2024-01-01T00:00:00Z"}},
    ]},
    "return_type": "entry",
    "request_options": {"paginate": {"start": 0, "rows": 10000},
                        "results_content_type": ["experimental"]},
}


def search_ids():
    ids = []
    start = 0
    while True:
        q = json.loads(json.dumps(QUERY))
        q["request_options"]["paginate"] = {"start": start, "rows": 5000}
        url = ("https://search.rcsb.org/rcsbsearch/v2/query?json="
               + urllib.parse.quote(json.dumps(q)))
        r = urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                   timeout=180, context=CTX)
        d = json.loads(r.read())
        got = [x["identifier"].lower() for x in d.get("result_set", [])]
        ids += got
        print(f"  fetched {len(ids)}/{d.get('total_count')}")
        if len(got) < 5000 or len(ids) >= d.get("total_count", 0):
            break
        start += 5000
        time.sleep(0.4)
    return sorted(set(ids))


def main():
    cand_p = os.path.join(OUT, "candidate_ids.json")
    if os.path.exists(cand_p):
        cand = json.load(open(cand_p))
    else:
        cand = search_ids()
        json.dump(cand, open(cand_p, "w"))
    exposed = set(json.load(open(os.path.join(
        ROOT, "dataset", "processed", "ssl_b2_exposed_pdb_ids.json"))))
    fresh = sorted(set(cand) - exposed)
    print(f"candidates {len(cand)}, exposed overlap {len(set(cand) & exposed)}, "
          f"fresh {len(fresh)}")
    rng = np.random.default_rng(SEED)
    pick = sorted(rng.choice(fresh, size=min(TARGET_N, len(fresh)), replace=False))

    man = {"schema": "MetaSieve.IndependentStructureAcquisition.v1",
           "coordinate_license": "CC0-1.0",
           "coordinate_service": "https://files.rcsb.org/download/",
           "query": QUERY["query"], "seed": SEED,
           "candidates": len(cand), "fresh_after_exposure_subtraction": len(fresh),
           "requested": len(pick), "files": {}, "failures": []}
    ok = 0
    for i, pid in enumerate(pick):
        dest = os.path.join(MMCIF, f"{pid}.cif.gz")
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            ok += 1
            continue
        try:
            u = f"https://files.rcsb.org/download/{pid}.cif.gz"
            data = urllib.request.urlopen(urllib.request.Request(u, headers=UA),
                                          timeout=120, context=CTX).read()
            with open(dest, "wb") as f:
                f.write(data)
            man["files"][pid] = {"bytes": len(data),
                                 "sha256": hashlib.sha256(data).hexdigest()}
            ok += 1
        except Exception as e:
            man["failures"].append({"pdb": pid, "error": f"{type(e).__name__}"})
        if (i + 1) % 100 == 0:
            print(f"  downloaded {ok}/{i+1}")
            json.dump(man, open(os.path.join(OUT, "acquisition_manifest.json"), "w"),
                      indent=2)
        time.sleep(0.05)
    man["downloaded"] = ok
    json.dump(man, open(os.path.join(OUT, "acquisition_manifest.json"), "w"), indent=2)
    print(f"downloaded {ok} of {len(pick)}; failures {len(man['failures'])}")


if __name__ == "__main__":
    main()
