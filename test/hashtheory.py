"""Verify the frozen theory tree is untouched, and record its hashes.

The theory directory is READ-ONLY for every implementation phase. This script
never writes into it.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
BASE = os.path.join(ROOT, "theory")
CORE = "FINAL_FROZEN_THEORY/00_CORE_THEORY/FINAL_THEORY_COMPLETE.md"
EXPECTED = "3d660448a585662083979c198d42258466cdcca7e0aab197095800cc2d42501e"


def main():
    recs = {}
    for root, _dirs, files in os.walk(BASE):
        for f in files:
            p = os.path.join(root, f)
            try:
                if os.path.getsize(p) > 5_000_000:
                    continue
                h = hashlib.sha256(open(p, "rb").read()).hexdigest()
            except OSError:
                continue
            rel = os.path.relpath(p, BASE).replace(os.sep, "/")
            recs[rel] = h
    ok = recs.get(CORE) == EXPECTED
    print(f"theory files hashed : {len(recs)}")
    print(f"core theory sha256  : {recs.get(CORE)}")
    print(f"matches frozen ref  : {ok}")
    out = os.path.join(ROOT, "theory", "FINAL_FROZEN_THEORY", "THEORY_HASHES.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump({"n_files": len(recs), "core_theory": CORE,
                   "core_sha256": recs.get(CORE), "matches_reference": ok,
                   "files": recs,
                   "note": "theory tree is READ-ONLY; hashed after implementation"},
                  fh, indent=2)
    print(f"wrote {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
