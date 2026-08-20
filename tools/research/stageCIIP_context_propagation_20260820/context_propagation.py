"""Read-only residue-level ESM contextual-propagation audit."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BRIDGE = ROOT / "tools" / "research" / "stageCIIP_potential_bridge"
X0 = ROOT / "tools" / "research" / "stageX_csc_signal"
X0C = X0 / "stageX0c_measurement_qualification_20260818"
RADIUS = 6
ERASURE_TOLERANCE = 1e-5

sys.path.insert(0, str(X0))
sys.path.insert(0, str(X0C))
from x0_common import load_duongly, normalize_construct_name, normalize_parent_name  # noqa: E402
from x0_i2 import build_pair_records  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def erase_at(sequence: str, pos: int) -> str:
    if not 1 <= pos <= len(sequence):
        raise ValueError(f"position {pos} outside sequence length {len(sequence)}")
    return sequence[:pos - 1] + "X" + sequence[pos:]


def summarize_delta(wt: np.ndarray, variant: np.ndarray, pos: int) -> tuple[dict, np.ndarray, np.ndarray]:
    if wt.shape != variant.shape or wt.ndim != 2 or wt.shape[1] != 640:
        raise ValueError(f"incompatible hidden states: {wt.shape} vs {variant.shape}")
    if not 1 <= pos < wt.shape[0]:
        raise ValueError(f"position {pos} outside hidden state length {wt.shape[0] - 1}")
    delta = variant[1:] - wt[1:]
    norms = np.linalg.norm(delta, axis=1)
    idx = np.arange(1, len(norms) + 1)
    local = np.abs(idx - pos) <= RADIUS
    context = ~local
    distances = np.abs(idx - pos)
    return {
        "site_delta_norm": float(norms[pos - 1]),
        "local_radius_6_mean_delta_norm": float(norms[local].mean()),
        "non_site_context_mean_delta_norm": float(norms[context].mean()),
        "full_sequence_mean_delta_norm": float(norms.mean()),
        "sequence_length": int(len(norms)),
    }, distances, norms


def masked_hidden(sequence: str, tokenizer, model, device: str) -> np.ndarray:
    import torch
    encoded = tokenizer(sequence, return_tensors="pt", truncation=True, max_length=1022)
    if device == "cuda":
        encoded = {key: value.cuda() for key, value in encoded.items()}
    with torch.no_grad():
        return model(**encoded).last_hidden_state[0].cpu().float().numpy()


def load_pairs() -> tuple[dict, list[dict], dict[str, dict]]:
    data1 = json.loads((BRIDGE / "DATA1A.json").read_text(encoding="utf-8"))
    data2 = json.loads((BRIDGE / "DATA2X2.json").read_text(encoding="utf-8"))
    _, _, sequences = load_duongly()
    pair_table = json.loads((X0 / "X0_PAIR_TABLE.json").read_text(encoding="utf-8"))
    records = build_pair_records(pair_table, sequences)
    by_construct = {normalize_construct_name(record["construct"]): record for record in records}
    selected = []
    for pair_index in data2["covered_pair_indices"]:
        pair = data1["pairs"][pair_index]
        record = by_construct.get(normalize_construct_name(pair["var_label"]))
        if record is None:
            raise KeyError(f"missing verified sequence for {pair['var_label']}")
        if pair["pos"] != record["pos"]:
            raise ValueError(f"coordinate mismatch for {pair['var_label']}")
        selected.append({"pair_index": int(pair_index), "pair": pair, "record": record})
    return data2, selected, sequences


def run(output_dir: Path, device: str) -> dict:
    prereg = HERE / "PREREGISTRATION.md"
    data2, selected, _ = load_pairs()
    cache_path = X0C / "q1_esm_cache.npz"
    cache = np.load(cache_path, allow_pickle=True)
    aggregate = []
    max_distance = 0
    per_pair_curves = []
    for item in selected:
        pair, record = item["pair"], item["record"]
        wt_key = next(k for k in cache.files if k.startswith("wt:") and normalize_parent_name(k[3:]) == normalize_parent_name(pair["parent"]))
        mt_key = next(k for k in cache.files if k.startswith("mt:") and normalize_construct_name(k[3:]) == normalize_construct_name(pair["var_label"]))
        summary, distances, norms = summarize_delta(cache[wt_key], cache[mt_key], pair["pos"])
        max_distance = max(max_distance, int(distances.max()))
        aggregate.append({"pair": item["pair_index"], "parent": pair["parent"],
                          "mutation": pair["mutation"], **summary})
        per_pair_curves.append((item["pair_index"], distances, norms))

    curves = np.full((len(per_pair_curves), max_distance + 1), np.nan, dtype=np.float32)
    pair_indices = np.empty(len(per_pair_curves), dtype=np.int32)
    for row, (pair_index, distances, norms) in enumerate(per_pair_curves):
        pair_indices[row] = pair_index
        for distance in np.unique(distances):
            curves[row, distance] = norms[distances == distance].mean()
    np.savez_compressed(output_dir / "context_distance_curves.npz",
                        pair_indices=pair_indices, distances=np.arange(max_distance + 1),
                        mean_delta_norm_by_distance=curves)

    import torch
    from x0_i2 import load_esm
    tokenizer, model, actual_device = load_esm(device)
    torch.backends.cudnn.benchmark = False
    erasure = []
    for item in selected:
        pair, record = item["pair"], item["record"]
        wt_masked = erase_at(record["wt_seq"], pair["pos"])
        mt_masked = erase_at(record["mt_seq"], pair["pos"])
        if wt_masked != mt_masked:
            raise AssertionError(f"erasure inputs differ for {pair['var_label']}")
        hw = masked_hidden(wt_masked, tokenizer, model, actual_device)
        hm = masked_hidden(mt_masked, tokenizer, model, actual_device)
        max_abs = float(np.max(np.abs(hw - hm)))
        erasure.append({"pair": item["pair_index"], "parent": pair["parent"],
                        "mutation": pair["mutation"], "inputs_identical": True,
                        "max_absolute_embedding_delta": max_abs,
                        "within_tolerance": max_abs <= ERASURE_TOLERANCE})
    if not all(row["within_tolerance"] for row in erasure):
        raise AssertionError("mutation erasure did not eliminate ESM delta")

    curve_mean = np.nanmean(curves, axis=0)
    result = {
        "schema": "MetaSieve.CIIP.ContextPropagation.v1",
        "preregistration_sha256": sha256(prereg),
        "inputs": {"data2x2_sha256": sha256(BRIDGE / "DATA2X2.json"),
                   "esm_cache_sha256": sha256(cache_path), "model": "facebook/esm2_t30_150M_UR50D",
                   "device": actual_device},
        "scope": "49 ESM-covered verified single-mutation pairs; read-only representation audit; Duong-Ly functional percent inhibition is not modeled",
        "summary": {"n_pairs": len(aggregate), "radius": RADIUS,
                    "site_delta_norm_mean": float(np.mean([x["site_delta_norm"] for x in aggregate])),
                    "local_radius_6_mean_delta_norm_mean": float(np.mean([x["local_radius_6_mean_delta_norm"] for x in aggregate])),
                    "non_site_context_mean_delta_norm_mean": float(np.mean([x["non_site_context_mean_delta_norm"] for x in aggregate])),
                    "full_sequence_mean_delta_norm_mean": float(np.mean([x["full_sequence_mean_delta_norm"] for x in aggregate])),
                    "distance_curve_mean": curve_mean.tolist(),
                    "erasure_max_absolute_embedding_delta": max(x["max_absolute_embedding_delta"] for x in erasure),
                    "erasure_tolerance": ERASURE_TOLERANCE},
        "per_pair": aggregate,
        "mutation_erasure": erasure,
        "interpretation": {"random_window_can_carry_mutation_information_via_context": "measured by non-site delta; no predictive claim in this read-only stage",
                              "context_only_predictive_value": "NOT_EVALUATED: requires separately preregistered model training",
                              "site_specific_increment_beyond_context": "representation magnitude reported; predictive increment not estimated",
                              "successor_authorization": "NOT_AUTHORIZED"},
    }
    serialized = json.dumps(result, indent=2) + "\n"
    (output_dir / "CONTEXT_PROPAGATION_RESULT.json").write_text(serialized, encoding="utf-8")
    (output_dir / "RESULT.json").write_text(serialized, encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    print(json.dumps(run(HERE, args.device), indent=2))


if __name__ == "__main__":
    main()
