"""P1 paired comparison report: learned arms vs frozen baselines and
vs each other.

Pairs records by (split, k, draw, target_id) across arms; computes
record-mean deltas (arm MSE - baseline MSE) and paired target-level
bootstrap 95% intervals; also CI/Spearman deltas. Reads:
- P1_BASELINES_RECORDS.json (ligand_only, tanimoto per-record metrics)
- any learned-arm artifact (P1_ARM3_ORDINARYFT.json, P1_ARM4_MAML.json,
  P1_ARM5_CNP.json, P1_ARM6_FSCAP.json, P1_ARM7_ACTFOUND.json)
Writes P1_COMPARISON.json + console table.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "artifacts"
K_LIST = (0, 1, 2, 3, 5, 10, 20, 40)
N_BOOT = 2000
BOOT_SEED = 20260862


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def stable_rng(*parts):
    raw = "|".join(str(x) for x in parts).encode("utf-8")
    return np.random.default_rng(int(hashlib.sha256(raw).hexdigest()[:16], 16))


def load_baseline_records():
    data = json.loads((OUT / "P1_BASELINES_RECORDS.json").read_text(encoding="utf-8"))
    out = {}
    for arm, by_k in data["records"].items():
        for k, rows in by_k.items():
            for r in rows:
                # pseudo-seed 0: baselines are deterministic (no seed loop)
                out[(arm, 0, r["split"], int(k), r["target_id"], r["draw"])] = r
    return out


def load_arm_records(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    out = {}
    for seed, sdata in data["seeds"].items():
        for r in sdata["records"]:
            out[(int(seed), r["split"], r["k"], r["target_id"], r["draw"])] = r
    return out


def bootstrap_delta(deltas, targets, rng, n=N_BOOT):
    ut = np.unique(targets)
    per_t = {t: np.asarray([d for d, tt in zip(deltas, targets) if tt == t])
             for t in ut}
    means = []
    for _ in range(n):
        idx = rng.integers(len(ut), size=len(ut))
        sample = np.concatenate([per_t[ut[i]] for i in idx])
        means.append(float(sample.mean()))
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def paired_rows(arm_records, other, split_name, k, seeds=(1, 2, 3),
                 ref_is_baseline=False):
    """[(arm_record, other_record, target_id)] matched pairs."""
    pairs = []
    for (seed, sp, kk, t, d), r in arm_records.items():
        if sp != split_name or kk != k or seed not in seeds:
            continue
        ref_seed = 0 if ref_is_baseline else seed
        o = other.get((ref_seed, sp, k, t, d))
        if o is not None:
            pairs.append((r, o, t))
    return pairs


def compare(arm_name, arm_records, ref_name, ref_records, rng,
             splits=("p_val", "p_test"), ref_is_baseline=False):
    out = {}
    for split_name in splits:
        for k in K_LIST:
            pairs = paired_rows(arm_records, ref_records, split_name, k,
                                ref_is_baseline=ref_is_baseline)
            if not pairs:
                continue
            mse_a = np.mean([a["mse"] for a, _, _ in pairs])
            mse_b = np.mean([b["mse"] for _, b, _ in pairs])
            d_mse = np.asarray([a["mse"] - b["mse"] for a, b, _ in pairs])
            d_ci = np.asarray([a["ci"] - b["ci"] for a, b, _ in pairs])
            d_rho = np.asarray([a["spearman"] - b["spearman"] for a, b, _ in pairs])
            targets = [t for _, _, t in pairs]
            lo, hi = bootstrap_delta(d_mse, targets, rng)
            out[f"{split_name}:k{k}"] = {
                "n_pairs": len(pairs),
                "arm_mse": float(mse_a), "ref_mse": float(mse_b),
                "delta_mse_mean": float(d_mse.mean()),
                "delta_mse_ci95": [float(lo), float(hi)],
                "delta_ci_mean": float(np.nanmean(d_ci)),
                "delta_spearman_mean": float(np.nanmean(d_rho)),
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
    bl = load_baseline_records()
    rng = stable_rng("stageP", "pcompare", BOOT_SEED)
    comparisons = {}
    for arm_name, path, is_learned in (("arm3", args.arm3, True),
                                       ("arm4", args.arm4, True),
                                       ("arm5", args.arm5, True),
                                       ("arm6", args.arm6, True),
                                       ("arm7", args.arm7, True),
                                       ("ligand_only", "", False),
                                       ("tanimoto", "", False)):
        if not path:
            continue
        recs = load_arm_records(Path(path))
        bl_lig = {k[1:]: v for k, v in bl.items() if k[0] == "ligand_only"}
        bl_tan = {k[1:]: v for k, v in bl.items() if k[0] == "tanimoto"}
        comparisons[f"{arm_name}_vs_ligand_only"] = compare(
            arm_name, recs, "ligand_only", bl_lig, rng, ref_is_baseline=True)
        comparisons[f"{arm_name}_vs_tanimoto"] = compare(
            arm_name, recs, "tanimoto", bl_tan, rng, ref_is_baseline=True)
        if is_learned and args.arm3:
            a3 = load_arm_records(Path(args.arm3))
            comparisons[f"{arm_name}_vs_arm3"] = compare(
                arm_name, recs, "arm3", a3, rng)
    out = {"schema": "MetaSieve.StageP.P1Comparison.v2",
           "artifacts": {"baselines": sha256_file(OUT / "P1_BASELINES_RECORDS.json"),
                         "arm3": sha256_file(Path(args.arm3)) if args.arm3 else None,
                         "arm4": sha256_file(Path(args.arm4)) if args.arm4 else None,
                         "arm5": sha256_file(Path(args.arm5)) if args.arm5 else None,
                         "arm6": sha256_file(Path(args.arm6)) if args.arm6 else None,
                         "arm7": sha256_file(Path(args.arm7)) if args.arm7 else None},
           "comparisons": comparisons}
    text = json.dumps(out, indent=1, sort_keys=True)
    path = OUT / "P1_COMPARISON.json"
    path.write_text(text, encoding="utf-8")
    art_sha = sha256_file(path)
    (OUT / "P1_COMPARISON.json.manifest.json").write_text(json.dumps({
        "schema": "MetaSieve.StageP.P1Comparison.v2.Manifest",
        "file": "P1_COMPARISON.json", "sha256": art_sha}, indent=1), encoding="utf-8")
    print("wrote", path, "sha", art_sha)
    for key, cmp in sorted(comparisons.items()):
        for sk, row in cmp.items():
            if sk.endswith("k5") or sk.endswith("k10"):
                print(f"{key:<28} {sk:<10} dMSE {row['delta_mse_mean']:+.4f} "
                      f"[{row['delta_mse_ci95'][0]:+.3f},{row['delta_mse_ci95'][1]:+.3f}] "
                      f"n={row['n_pairs']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
