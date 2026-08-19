"""Stage CIIP-1A data contract builder (prereg 31d3eeaf...).

Frozen pair table = Q0B-admitted single-point Duong-Ly variant records
(65). Rows/ligands/features = x0_i1.load_features(). Labels = % inhibition
from duongly_mmc3 Table S2, raw values kept (endpoint never relabeled).
Targets: c_vl = d_vl - mean_l(d_vl), d_vl = y[v,L] - y[WT,L]. Pair-level
split 60/20/20 stratified per parent, SHA-256 keyed. Writes DATA1A.json
+ updates SHA256SUMS. Read-only with respect to every input.
"""
from __future__ import annotations

import hashlib
import json
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SIG = HERE.parent / "stageX_csc_signal"
sys.path.insert(0, str(SIG))
sys.path.insert(0, str(SIG / "stageX0c_measurement_qualification_20260818"))

import x0_i1  # noqa: E402
from x0_common import load_duongly, normalize_construct_name, stable_rng, sha256_file  # noqa: E402

PREREG_SHA = "31d3eeaf6a0d77c46b3bbbee0fe9d2ff667aadeeb7d9dcabd26ca59ec48d5196"
SCHEMA = "MetaSieve.StageCIIP1A.Data.v1"
SPLIT_SEED = 20260819
VAL_FRAC = 0.20
TEST_FRAC = 0.20
SPLITS = ("train", "val", "test")


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def main() -> int:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        rows, compounds, prot_feats, lig_feats, scaffolds, row_meta = x0_i1.load_features()
        _info, matrix, _seqs = load_duongly()
    Y = matrix.iloc[1:, 1:].to_numpy(dtype=np.float64)  # (97, 183) raw % inhibition
    s2_labels = [str(x).strip() for x in matrix.iloc[1:, 0].tolist()]
    assert len(s2_labels) == Y.shape[0] == len(rows) == prot_feats.shape[0]

    audit_path = (SIG / "stageX0c_measurement_qualification_20260818"
                  / "Q0B_MAPPING_AUDIT.json")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    admitted = [r for r in audit["duongly_variant_records"]
                if r["admission_status"] == "admitted"]
    assert len(admitted) == 65, len(admitted)

    def find_row(label: str):
        key = normalize_construct_name(label)
        hits = [i for i, s in enumerate(s2_labels)
                if normalize_construct_name(s) == key]
        assert len(hits) == 1, (label, key, hits)
        return hits[0]

    pairs = []
    for r in admitted:
        wt_i = find_row(r["parent_gene"])
        var_i = find_row(r["source_row"])
        assert wt_i != var_i
        pairs.append({
            "wt_row": wt_i, "var_row": var_i,
            "parent": r["parent_gene"],
            "mutation": r["reported_mutation_notation"],
            "pos": int(r["substitutions"][0]["pos"]),
            "wt_label": s2_labels[wt_i],
            "var_label": s2_labels[var_i],
        })
    # per-pair contrast targets (common finite ligands only)
    targets = []
    for p in pairs:
        mask = np.isfinite(Y[p["wt_row"]]) & np.isfinite(Y[p["var_row"]])
        d = Y[p["var_row"]][mask] - Y[p["wt_row"]][mask]
        c = d - d.mean()
        targets.append({
            "lig_idx": np.where(mask)[0].tolist(),
            "d": d.tolist(), "c": c.tolist(),
            "n_lig": int(mask.sum()),
        })
        assert abs(float(c.mean())) < 1e-9
    # per-parent stratified split 60/20/20
    rng = stable_rng("stageCIIP1A", "split", SPLIT_SEED)
    by_parent = {}
    for i, p in enumerate(pairs):
        by_parent.setdefault(p["parent"], []).append(i)
    split_of = np.zeros(len(pairs), dtype=np.int8)
    for parent, idxs in sorted(by_parent.items()):
        order = idxs.copy()
        rng.shuffle(order)
        n = len(order)
        n_test = max(0, int(round(n * TEST_FRAC)))
        n_val = max(0, int(round(n * VAL_FRAC)))
        n_train = n - n_test - n_val
        for j, i in enumerate(order):
            split_of[i] = 0 if j < n_train else (1 if j < n_train + n_val else 2)
    counts = {s: int((split_of == k).sum()) for k, s in enumerate(SPLITS)}
    assert sum(counts.values()) == len(pairs)

    out = {
        "schema": SCHEMA,
        "preregistration_sha256": PREREG_SHA,
        "inputs_sha256": {
            "q0b_audit": sha256_file(audit_path),
            "duongly_mmc3": sha256_file(SIG / "downloads" / "duongly_mmc3.xlsx"),
            "duongly_mmc2": sha256_file(SIG / "downloads" / "duongly_mmc2.xlsx"),
        },
        "endpoint": "percent inhibition (raw; never relabeled pK/Ki/Kd)",
        "shape": {"rows": int(Y.shape[0]), "ligands": int(Y.shape[1])},
        "label_stats": {
            "finite_fraction": float(np.isfinite(Y).mean()),
            "min": float(np.nanmin(Y)), "max": float(np.nanmax(Y)),
            "n_out_of_bounds_0_100": int(((Y < 0) | (Y > 100)).sum()),
        },
        "pairs": pairs,
        "targets": targets,
        "split": {"seed": SPLIT_SEED, "frac": [1 - VAL_FRAC - TEST_FRAC,
                                               VAL_FRAC, TEST_FRAC],
                  "pair_split": split_of.tolist(), "counts": counts,
                  "parents_per_split": {
                      s: int(len({pairs[i]["parent"] for i in np.where(
                          split_of == k)[0]}))
                      for k, s in enumerate(SPLITS)}},
        "feature_shapes": {"prot": list(prot_feats.shape),
                           "lig": list(lig_feats.shape)},
        "rows": rows,
        "ligands": compounds,
        "row_meta_admitted": [row_meta[p["var_row"]] for p in pairs],
    }
    art = HERE / "DATA1A.json"
    art.write_text(json.dumps(out, indent=1), encoding="utf-8")
    np.savez_compressed(HERE / "DATA1A.npz",
                        Y=Y, prot=prot_feats, lig=lig_feats,
                        pair_split=split_of)
    sums = {}
    for f in ("DATA1A.json", "DATA1A.npz"):
        sums[f] = sha256_file(HERE / f)
    (HERE / "SHA256SUMS").write_text(
        "".join(v + " *" + k + "\n" for k, v in sorted(sums.items())),
        encoding="utf-8")
    print(json.dumps({"pairs": len(pairs), "counts": counts,
                      "parents_per_split": out["split"]["parents_per_split"],
                      "finite_frac": out["label_stats"]["finite_fraction"],
                      "oob_0_100": out["label_stats"]["n_out_of_bounds_0_100"],
                      "sha": sums}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
