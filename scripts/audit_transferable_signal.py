"""Is there any transferable query-specific signal in support labels?

Stage 5 rejected two different query-specific transports. Before designing a
third, this measures — with labels and chemistry only, no training — how much a
*perfect* support-selection scheme could possibly gain over the support mean.

Estimators per episode (support labels are inputs, query labels are targets):

* `support_mean`      : the level baseline every learned transport collapsed to
* `nearest_support`   : the single most Tanimoto-similar support ligand
* `similarity_softmax`: softmax(beta * Tanimoto) weighted support labels
* `oracle_best_support`: the best single support *chosen with query labels* —
  not a model, an upper bound on any selection mechanism
* `oracle_target_mean`: the level ceiling

If `oracle_best_support` is close to `support_mean`, no selection mechanism can
help and the bottleneck is data, not architecture.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK


def morgan_fingerprints(smiles: dict[str, str], radius: int = 2, bits: int = 2048):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius, fpSize=bits)
    out: dict[str, np.ndarray] = {}
    for key, value in smiles.items():
        molecule = Chem.MolFromSmiles(value) if value else None
        if molecule is None:
            continue
        array = np.zeros(bits, dtype=np.uint8)
        for index in generator.GetFingerprint(molecule).GetOnBits():
            array[index] = 1
        out[key] = array
    return out


def tanimoto(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """`left` [n,b], `right` [m,b] -> [n,m]."""
    intersection = left.astype(np.float32) @ right.astype(np.float32).T
    left_sum = left.sum(-1, keepdims=True).astype(np.float32)
    right_sum = right.sum(-1, keepdims=True).astype(np.float32).T
    union = left_sum + right_sum - intersection
    return intersection / np.maximum(union, 1e-9)


def component_mean(rows: list[dict], field: str) -> float:
    by_target: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        value = row.get(field)
        if value is None or not np.isfinite(value):
            continue
        by_target.setdefault((row["component"], row["target"]), []).append(value)
    by_component: dict[str, list[float]] = {}
    for (component, _), values in by_target.items():
        by_component.setdefault(component, []).append(float(np.mean(values)))
    if not by_component:
        return float("nan")
    return float(np.mean([np.mean(v) for v in by_component.values()]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--split", default="meta_val",
                        choices=("meta_train", "meta_val", "meta_test"))
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--query-size", type=int, default=20)
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--beta", type=float, default=8.0)
    parser.add_argument("--bits", type=int, default=2048,
                        help="fingerprint width; the production episode "
                             "pipeline in scripts/qpsmp_data.py uses 1024")
    parser.add_argument("--radius", type=int, default=2)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    ligand_rows = [json.loads(line) for line in
                   (data.corpus / "ligands.jsonl").read_text(
                       encoding="utf-8").splitlines() if line.strip()]
    fingerprints = morgan_fingerprints(
        {row["drug_key"]: row.get("smiles") for row in ligand_rows},
        radius=args.radius, bits=args.bits)

    banks = data.fixed_nested_episode_banks(
        args.split, (0, 1, 2, 3, 5), args.query_size, args.draws,
        args.evaluation_seed, None)

    report: dict[str, dict] = {}
    similarity_pairs: list[tuple[float, float]] = []
    for k, specs in banks.items():
        if k == 0:
            continue
        rows = []
        for spec in specs:
            support_keys = [data.cells[i]["ligand_id"] for i in spec.support]
            query_keys = [data.cells[i]["ligand_id"] for i in spec.query]
            if not all(key in fingerprints for key in support_keys + query_keys):
                continue
            support_y = np.asarray([data.cells[i]["pK"] for i in spec.support])
            query_y = np.asarray([data.cells[i]["pK"] for i in spec.query])
            support_fp = np.stack([fingerprints[key] for key in support_keys])
            query_fp = np.stack([fingerprints[key] for key in query_keys])
            similarity = tanimoto(query_fp, support_fp)          # [Q,K]
            errors = (query_y[:, None] - support_y[None, :]) ** 2  # [Q,K]
            weights = np.exp(args.beta * (similarity - similarity.max(
                -1, keepdims=True)))
            weights = weights / weights.sum(-1, keepdims=True)
            nearest = similarity.argmax(-1)
            if k >= 2:
                pairs = tanimoto(query_fp, query_fp)
                gap = np.abs(query_y[:, None] - query_y[None, :])
                triu = np.triu_indices(len(query_y), 1)
                similarity_pairs.extend(zip(pairs[triu].tolist(), gap[triu].tolist()))
            rows.append({
                "component": spec.component, "target": spec.target,
                "support_mean": float(((query_y - support_y.mean()) ** 2).mean()),
                "nearest_support": float(errors[np.arange(len(query_y)), nearest].mean()),
                "similarity_softmax": float(
                    ((query_y - (weights * support_y[None, :]).sum(-1)) ** 2).mean()),
                "oracle_best_support": float(errors.min(-1).mean()),
                "oracle_target_mean": float(((query_y - query_y.mean()) ** 2).mean()),
                "mean_max_tanimoto": float(similarity.max(-1).mean()),
            })
        report[str(k)] = {
            field: component_mean(rows, field)
            for field in ("support_mean", "nearest_support", "similarity_softmax",
                          "oracle_best_support", "oracle_target_mean",
                          "mean_max_tanimoto")
        }
        report[str(k)]["episodes"] = len(rows)

    correlation = float("nan")
    if similarity_pairs:
        values = np.asarray(similarity_pairs)
        if values[:, 0].std() > 0 and values[:, 1].std() > 0:
            correlation = float(np.corrcoef(values[:, 0], values[:, 1])[0, 1])

    payload = {
        "schema": "MetaSieve.TransferableSignalAudit.v1",
        "split": args.split,
        "evaluation_seed": args.evaluation_seed,
        "draws_per_target": args.draws,
        "fingerprint": (f"Morgan r={args.radius} {args.bits} bits, Tanimoto, "
                        f"softmax beta={args.beta}"),
        "matches_production_pipeline": args.bits == 1024 and args.radius == 2,
        "estimator_mse_pk": report,
        "within_target_similarity_vs_abs_affinity_gap_pearson": correlation,
        "note": ("oracle_best_support selects the best support using query "
                 "labels; it is an upper bound on any selection mechanism, not "
                 "a deployable predictor"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
