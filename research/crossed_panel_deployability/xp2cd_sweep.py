"""XP2-C and XP2-D sweeps: support size k <= 5, section rank d <= 5, both closures."""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from xp2cd_section import OUT, pretty, run  # noqa: E402

ALL = {}


def go(tag, **kw):
    r = run(**kw)
    pretty(r)
    ALL[tag] = r
    json.dump(ALL, open(os.path.join(OUT, "xp2cd_sweeps.json"), "w"),
              indent=2, default=float)
    return r


print("#" * 78)
print("# XP2-C  protein-group closure only, ligands reused (XP1-comparable regime)")
for k in (1, 2, 3, 4, 5):
    go(f"C_k{k}", closure="protein", rank=3, k=k, ligand_arm="L-ECFP")
for d in (1, 2, 5):
    go(f"C_d{d}", closure="protein", rank=d, k=5, ligand_arm="L-ECFP")

print("\n" + "#" * 78)
print("# XP2-D  DOUBLE held-out: protein group AND ligand scaffold component")
for k in (1, 2, 3, 4, 5):
    go(f"D_k{k}", closure="double", rank=3, k=k, ligand_arm="L-ECFP")
for d in (1, 2, 5):
    go(f"D_d{d}", closure="double", rank=d, k=5, ligand_arm="L-ECFP")
for feat in ("L-CHEMBERTA", "L-DESC", "L-RANDOM"):
    go(f"D_feat_{feat}", closure="double", rank=3, k=5, ligand_arm=feat)

print("\n" + "=" * 78)
print("SUMMARY  (R2_gamma vs ADD on held-out cells)")
hdr = f"{'run':16s} {'SEC':>26s} {'Delta_specific_foreign':>26s} {'ORACLE':>22s}"
print(hdr)
for tag, r in ALL.items():
    s = r["arms"].get("SEC", {}).get("r2_gamma_vs_ADD")
    o = r["arms"].get("ORACLE-TRSC", {}).get("r2_gamma_vs_ADD")
    f = r["contrasts"].get("Delta_specific_foreign__FOREIGN_minus_SEC")
    def fmt(v, k="point"):
        return f"{v[k]:+.4f}[{v['ci95'][0]:+.4f},{v['ci95'][1]:+.4f}]" if v else "--"
    print(f"{tag:16s} {fmt(s):>26s} {fmt(f):>26s} {fmt(o):>22s}")
print("\nwrote", os.path.join(OUT, "xp2cd_sweeps.json"))
