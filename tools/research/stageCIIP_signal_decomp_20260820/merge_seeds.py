"""Merge per-seed results into RESULT.json (multi-seed dict)."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
seeds = sys.argv[1:] if len(sys.argv) > 1 else ["1", "2", "3"]
out = {}
for s in seeds:
    out[s] = json.loads((HERE / f"SEED{s}_RESULT.json").read_text(encoding="utf-8"))
(HERE / "RESULT.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
print("merged seeds:", sorted(out.keys()))
