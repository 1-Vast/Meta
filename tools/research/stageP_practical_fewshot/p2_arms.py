"""P2 screening metrics from learned-arm artifacts.

For each arm artifact: per (split, k) rank query cells within each
record by yhat, compute PR-AUC (active = pKi >= 6.0), EF1%/EF5%/EF10%,
top-1 hit rate, and active fraction (record-mean, matching
p_analysis.py's P2_SCREENING definitions). Writes P2_ARMS.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "artifacts"
ACTIVE_THRESH = 6.0
K_LIST = (0, 5)
N_BOOT = 2000


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def stable_rng(*parts):
    raw = "|".join(str(x) for x in parts).encode("utf-8")
    return np.random.default_rng(int(hashlib.sha256(raw).hexdigest()[:16], 16))


def pr_auc(yhat, y):
    order = np.argsort(-np.asarray(yhat))
    ys = np.asarray(y)[order]
    act = ys >= ACTIVE_THRESH
    n_act = int(act.sum())
    if n_act == 0 or n_act == len(ys):
        return float("nan")
    precs = []
    seen = 0
    hits = 0
    for i, a in enumerate(act, start=1):
        if a:
            hits += 1
            precs.append(hits / i)
    return float(np.mean(precs))


def ef(yhat, y, frac):
    n = len(yhat)
    k = max(1, int(round(n * frac)))
    order = np.argsort(-np.asarray(yhat))[:k]
    act = np.asarray(y) >= ACTIVE_THRESH
    n_act = int(act.sum())
    if n_act == 0:
        return float("nan")
    return float(act[order].sum() / k / (n_act / n))


def top1_hit(yhat, y):
    order = np.argsort(-np.asarray(yhat))[:1]
    return float((np.asarray(y) >= ACTIVE_THRESH)[order].sum())


def arm_metrics(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    per = {}
    for seed, sdata in data["seeds"].items():
        for r in sdata["records"]:
            key = (r["split"], r["k"])
            per.setdefault(key, []).append({
                "pr_auc": pr_auc(r["yhat"], r["y"]),
                "ef1": ef(r["yhat"], r["y"], 0.01),
                "ef5": ef(r["yhat"], r["y"], 0.05),
                "ef10": ef(r["yhat"], r["y"], 0.10),
                "top1_hit": top1_hit(r["yhat"], r["y"]),
            })
    out = {}
    for key, rows in per.items():
        out[f"{key[0]}:k{key[1]}"] = {
            "n_records": len(rows),
            "pr_auc": float(np.nanmean([x["pr_auc"] for x in rows])),
            "ef1": float(np.nanmean([x["ef1"] for x in rows])),
            "ef5": float(np.nanmean([x["ef5"] for x in rows])),
            "ef10": float(np.nanmean([x["ef10"] for x in rows])),
            "top1_hit": float(np.mean([x["top1_hit"] for x in rows])),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm3", default=str(OUT / "P1_ARM3_ORDINARYFT.json"))
    ap.add_argument("--arm4", default="")
    ap.add_argument("--arm5", default="")
    ap.add_argument("--arm6", default="")
    ap.add_argument("--arm7", default="")
    args = ap.parse_args()
    out = {"schema": "MetaSieve.StageP.P2Arms.v1",
           "active_threshold": ACTIVE_THRESH,
           "arms": {}}
    for arm_name, path in (("arm3", args.arm3), ("arm4", args.arm4),
                           ("arm5", args.arm5), ("arm6", args.arm6),
                           ("arm7", args.arm7)):
        if not path:
            continue
        out["arms"][arm_name] = {
            "artifact_sha256": sha256_file(Path(path)),
            "metrics": arm_metrics(Path(path)),
        }
    text = json.dumps(out, indent=1, sort_keys=True)
    path = OUT / "P2_ARMS.json"
    path.write_text(text, encoding="utf-8")
    art_sha = sha256_file(path)
    (OUT / "P2_ARMS.json.manifest.json").write_text(json.dumps({
        "schema": "MetaSieve.StageP.P2Arms.v1.Manifest",
        "file": "P2_ARMS.json", "sha256": art_sha}, indent=1), encoding="utf-8")
    print("wrote", path, "sha", art_sha)
    for arm_name, arm in out["arms"].items():
        for key, row in sorted(arm["metrics"].items()):
            print(f"{arm_name:<6} {key:<12} PR-AUC {row['pr_auc']:.3f} "
                  f"EF1 {row['ef1']:.2f} EF5 {row['ef5']:.2f} "
                  f"top1 {row['top1_hit']:.2f} n={row['n_records']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
