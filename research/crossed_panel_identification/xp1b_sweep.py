"""XP1-B sweeps: closure level, section rank, support size, censoring block."""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xp1b_transfer import REPORT, pretty, run  # noqa: E402

ALL = {}


def go(tag, **kw):
    res = run(verbose=True, **kw)
    pretty(res)
    ALL[tag] = res
    return res


print("#" * 78)
print("# S1  closure level (group = strict, family, pocket60 = permissive)")
for cl in ("group", "family", "pocket60"):
    go(f"closure_{cl}", closure=cl, rank=8, k_support=16)

print("\n" + "#" * 78)
print("# S2  section rank r  (= required dim of the biological coordinate)")
for r in (1, 2, 3, 5, 8, 12):
    go(f"rank_{r}", closure="group", rank=r, k_support=16)

print("\n" + "#" * 78)
print("# S3  support size k")
for k in (4, 8, 16, 32, 64):
    go(f"k_{k}", closure="group", rank=8, k_support=k)

print("\n" + "#" * 78)
print("# S4  censoring sensitivity: BLK-METZ-70")
go("density_70", closure="group", rank=8, k_support=16, density=0.70)

p = os.path.join(REPORT, "xp1b_sweeps.json")
with open(p, "w") as f:
    json.dump(ALL, f, indent=2, default=float)
print("\nwrote", p)

print("\n" + "=" * 78)
print("SUMMARY TABLE  (R2_gamma vs A2, held-out cells of held-out proteins)")
print(f"{'run':16s} {'A4 few-shot':>24s} {'A3 zero-shot kernel':>24s} "
      f"{'A3B knn':>22s} {'AO1 oracle':>22s}")
for tag, r in ALL.items():
    def g(a):
        v = r["arms"].get(a, {}).get("r2_gamma_vs_A2")
        return f"{v['point']:+.3f}[{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]" if v else "--"
    print(f"{tag:16s} {g('A4'):>24s} {g('A3::pocket_identity_kernel'):>24s} "
          f"{g('A3B::knn_pocket'):>22s} {g('AO1'):>22s}")
