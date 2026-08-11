"""Run the frozen R0-B ceiling and power audit without fitting a model."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import torch

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM, MAX_ATOMS
from research.correspondence_router.r0_exact_distance import (
    additive_checkerboard_rps,
    component_macro,
    negative_log_likelihood,
    ranked_probability_score,
)
from research.e0_identifiability.run_tdir_pilot import _load_frozen_model
from scripts.cache_exact_structure_proteins import unpack_projected_row
from scripts.cache_r0b_exact_geometry import (
    FROZEN_LIGAND_BANK_SHA256,
    FROZEN_PANEL_SHA256,
    R0C_LIGAND_BANK_SHA256,
    R0C_PANEL_SHA256,
    unpack_exact_geometry_row,
    verify_exact_geometry_bank,
)
from scripts.data_contract import read_jsonl, write_jsonl
from scripts.structure_sources.rcsb import sha256_file


PREFIT_BOOTSTRAP_SEED = 1700
PREFIT_BOOTSTRAP_DRAWS = 10_000
FROZEN_CHECKPOINT_SHA256 = (
    "90b0010b81fa2758a2dbdd1a8dbe06adae2e05acbbc267ccb62ceee6ff6c4f37")


def bootstrap_mde80(component_scores: np.ndarray, *, seed: int,
                    draws: int = PREFIT_BOOTSTRAP_DRAWS) -> float:
    """Bootstrap MDE for a constant shift at one-sided alpha=.05, power=.80."""
    values = np.asarray(component_scores, dtype=np.float64)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("MDE80 needs at least two finite component scores")
    if draws < 1:
        raise ValueError("MDE80 bootstrap draws must be positive")
    centered = values - values.mean()
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(draws, len(values)))
    means = centered[indices].mean(axis=1)
    return float(
        np.quantile(means, 0.95, method="linear")
        - np.quantile(means, 0.20, method="linear"))


def _load_exact_proteins(root: Path, wanted: set[str]) -> dict[str, dict]:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    unresolved = set(wanted)
    result = {}
    for shard_info in manifest["shards"]:
        if not unresolved:
            break
        path = root / shard_info["path"]
        if sha256_file(path) != shard_info["sha256"]:
            raise ValueError(f"exact protein shard hash differs: {path}")
        with np.load(path, allow_pickle=False) as payload:
            keys = payload["keys"].tolist()
            for index, key in enumerate(keys):
                if key in unresolved:
                    result[key] = unpack_projected_row(payload, index)
                    unresolved.remove(key)
    if unresolved:
        raise KeyError(f"exact protein rows missing: {len(unresolved)}")
    return result


def _load_exact_geometry(root: Path, rows: list[dict]) -> dict[str, dict]:
    by_shard: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_shard[row["shard"]].append(row)
    result = {}
    for name, shard_rows in by_shard.items():
        with np.load(root / name, allow_pickle=False) as payload:
            for row in shard_rows:
                result[row["source_entry_id"]] = unpack_exact_geometry_row(
                    payload, int(row["row"]))
    return result


def _frozen_distance_prior(model, graph: dict, protein: dict,
                           device: str) -> np.ndarray:
    X = torch.zeros(1, MAX_ATOMS, ATOM_FEAT_DIM, device=device)
    A = torch.zeros(1, MAX_ATOMS, MAX_ATOMS, BOND_FEAT_DIM, device=device)
    mask = graph["mask"].float().unsqueeze(0).to(device)
    X[0] = graph["X"].float().to(device)
    edge = graph["edge_index"].long().to(device)
    A[0, edge[0], edge[1]] = graph["edge_attr"].float().to(device)
    residue = torch.from_numpy(protein["slot_projected"]).float().unsqueeze(0).to(device)
    residue_mask = torch.from_numpy(protein["slot_mask"]).float().unsqueeze(0).to(device)
    with torch.inference_mode():
        _, atom_state = model.ligand(X, A, mask)
        prediction = model.bridge(atom_state, mask, residue, residue_mask)
        probability = torch.softmax(prediction.distance_logits[0], dim=-1)
    return probability.cpu().numpy()


def run_prefit(panel_path: str | Path, geometry_dir: str | Path,
               protein_dir: str | Path, ligand_bank_path: str | Path,
               checkpoint_path: str | Path, output_path: str | Path, *,
               device: str = "cuda", contract: str = "r0b") -> dict:
    if not device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("R0-B frozen prior scoring is registered for CUDA")
    panel_path, geometry_dir = Path(panel_path), Path(geometry_dir)
    protein_dir, ligand_bank_path = Path(protein_dir), Path(ligand_bank_path)
    checkpoint_path, output_path = Path(checkpoint_path), Path(output_path)
    if contract == "r0b":
        expected_panel = FROZEN_PANEL_SHA256
        expected_ligand = FROZEN_LIGAND_BANK_SHA256
        heldout_split, expected_records, minimum_components = "heldout_a", 144, 30
    elif contract == "r0c":
        expected_panel = R0C_PANEL_SHA256
        expected_ligand = R0C_LIGAND_BANK_SHA256
        heldout_split, expected_records, minimum_components = "heldout_b", 219, 120
    else:
        raise ValueError(f"unknown prefit contract: {contract}")
    if output_path.exists():
        raise FileExistsError(f"{contract.upper()} prefit report already exists: {output_path}")
    if sha256_file(panel_path) != expected_panel:
        raise ValueError(f"prefit panel is not frozen {contract.upper()}")
    if sha256_file(ligand_bank_path) != expected_ligand:
        raise ValueError(f"prefit ligand bank is not frozen {contract.upper()} input")
    if sha256_file(checkpoint_path) != FROZEN_CHECKPOINT_SHA256:
        raise ValueError("prefit checkpoint is not frozen P1B")
    geometry_verification = verify_exact_geometry_bank(
        geometry_dir, contract=contract)

    protein_manifest = json.loads(
        (protein_dir / "manifest.json").read_text(encoding="utf-8"))
    if protein_manifest.get("records_sha256") != expected_panel or \
            protein_manifest.get("checkpoint_sha256") != FROZEN_CHECKPOINT_SHA256:
        raise ValueError("exact protein bank does not reference frozen R0-B inputs")

    panel = read_jsonl(panel_path)
    panel_by_entry = {row["source_entry_id"]: row for row in panel}
    heldout_index = [row for row in read_jsonl(geometry_dir / "index.jsonl")
                     if row["r0_split"] == heldout_split]
    if len(heldout_index) != expected_records:
        raise ValueError(f"{contract.upper()} heldout geometry count differs")
    geometry = _load_exact_geometry(geometry_dir, heldout_index)
    proteins = _load_exact_proteins(
        protein_dir, {row["sequence_sha256"] for row in heldout_index})
    ligands = torch.load(ligand_bank_path, map_location="cpu", weights_only=False)
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    model, _ = _load_frozen_model(
        checkpoint_path, int(checkpoint["protein_dim"]), device)

    prior_system, prior_nll_system, additive_system = {}, {}, {}
    component_of = {}
    bin_counts = np.zeros(5, dtype=np.int64)
    heldout_movable, heldout_residues = 0, 0
    system_rows = []
    for index_row in heldout_index:
        entry = index_row["source_entry_id"]
        record = panel_by_entry[entry]
        exact = geometry[entry]
        graph = ligands[record["ccd_sha256"]]
        protein = proteins[record["sequence_sha256"]]
        component = record.get(
            "homology_group_id", record.get("r0c_final_component_id"))
        if (len(protein["exact_projected"]) != int(index_row["residues"])
                or index_row["homology_group_id"] != component
                or index_row["ccd_sha256"] != record["ccd_sha256"]):
            raise ValueError(f"panel/index/protein metadata disagree: {entry}")
        slot_probability = _frozen_distance_prior(
            model, graph, protein, device)
        atoms = int(index_row["atoms"])
        slots = exact["slot_of_residue"].astype(np.int64)
        if not np.all(protein["slot_mask"][slots] > 0):
            raise ValueError(f"exact residues enter masked P1B slots: {entry}")
        probability = slot_probability[:atoms, slots]
        labels = exact["distance_bin"].astype(np.int64)
        if not np.isfinite(probability).all() or not np.allclose(
                probability.sum(-1), 1.0, atol=1e-6):
            raise ValueError(f"invalid lifted prior: {entry}")
        prior_system[entry] = float(ranked_probability_score(
            probability, labels).mean())
        prior_nll_system[entry] = float(negative_log_likelihood(
            probability, labels).mean())
        additive_system[entry] = additive_checkerboard_rps(labels, slots)
        component_of[entry] = component
        local_bins = np.bincount(labels.reshape(-1), minlength=5)[:5]
        bin_counts += local_bins
        multiplicity = np.bincount(slots, minlength=128)
        movable = int((multiplicity[slots] >= 2).sum())
        heldout_movable += movable
        heldout_residues += len(slots)
        system_rows.append({
            "source_entry_id": entry,
            "homology_group_id": component,
            "atoms": atoms,
            "residues": len(slots),
            "distance_bin_counts": local_bins.tolist(),
            "movable_residues": movable,
            "prior_rps": prior_system[entry],
            "prior_nll": prior_nll_system[entry],
            "additive_rps": additive_system[entry],
        })

    prior_component = component_macro(prior_system, component_of)
    nll_component = component_macro(prior_nll_system, component_of)
    additive_component = component_macro(additive_system, component_of)
    s_prior = float(np.mean(list(prior_component.values())))
    s_add_star = float(np.mean(list(additive_component.values())))
    delta_star = 0.05 * s_prior
    component_keys = sorted(prior_component)
    mde80 = bootstrap_mde80(
        np.asarray([prior_component[key] for key in component_keys]),
        seed=PREFIT_BOOTSTRAP_SEED)
    component_records = Counter(component_of.values())
    largest_share = max(component_records.values()) / len(component_of)
    geometry_manifest = json.loads(
        (geometry_dir / "manifest.json").read_text(encoding="utf-8"))
    heldout_movable_fraction = heldout_movable / heldout_residues
    gates = {
        "all_bins_supported": bool((bin_counts > 0).all()),
        "additive_pair_headroom": bool(s_add_star >= delta_star),
        "mde80_within_delta_star": bool(mde80 <= delta_star),
        "largest_component_below_0_20": bool(largest_share < 0.20),
        "at_least_required_components": len(prior_component) >= minimum_components,
        "movable_residue_fraction_at_least_0_50": (
            heldout_movable_fraction >= 0.50),
        "zero_geometry_exclusions": geometry_manifest["excluded_records"] == 0,
    }
    passed = all(gates.values())
    component_rows = [{
        "homology_group_id": key,
        "systems": int(component_records[key]),
        "prior_rps": prior_component[key],
        "prior_nll": nll_component[key],
        "additive_rps": additive_component[key],
    } for key in component_keys]
    systems_path = output_path.with_name(output_path.stem + ".systems.jsonl")
    components_path = output_path.with_name(output_path.stem + ".components.jsonl")
    write_jsonl(systems_path, system_rows)
    write_jsonl(components_path, component_rows)
    report = {
        "schema": ("MetaSieve.R0BPrefitAdmission.v1" if contract == "r0b"
                   else "MetaSieve.R0CPrefitAdmission.v1"),
        "contract": contract,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "device": device,
        "affinity_labels_used": False,
        "trainable_parameters": 0,
        "panel_sha256": sha256_file(panel_path),
        "geometry_manifest_sha256": sha256_file(geometry_dir / "manifest.json"),
        "protein_manifest_sha256": sha256_file(protein_dir / "manifest.json"),
        "ligand_bank_sha256": sha256_file(ligand_bank_path),
        "checkpoint_sha256": sha256_file(checkpoint_path),
        "geometry_verification": geometry_verification,
        "heldout_records": len(heldout_index),
        "heldout_components": len(prior_component),
        "largest_component_share": largest_share,
        "distance_bin_counts": bin_counts.tolist(),
        "heldout_movable_residue_fraction": heldout_movable_fraction,
        "s_prior": s_prior,
        "s_prior_nll": float(np.mean(list(nll_component.values()))),
        "s_exact_star": 0.0,
        "s_add_star": s_add_star,
        "delta_star": delta_star,
        "u_pair": s_add_star,
        "mde80": mde80,
        "bootstrap_seed": PREFIT_BOOTSTRAP_SEED,
        "bootstrap_draws": PREFIT_BOOTSTRAP_DRAWS,
        "systems_path": str(systems_path.resolve()),
        "systems_sha256": sha256_file(systems_path),
        "components_path": str(components_path.resolve()),
        "components_sha256": sha256_file(components_path),
        "runner_sha256": sha256_file(Path(__file__)),
        "gates": gates,
        "verdict": (
            ("R0B_PREFIT_ADMISSION_PASS" if passed else "R0B_NOT_RUN_FAIL_CLOSED")
            if contract == "r0b" else
            ("R0C_PREFIT_ADMISSION_PASS" if passed else "R0C_NOT_RUN_FAIL_CLOSED")
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("panel")
    parser.add_argument("geometry")
    parser.add_argument("proteins")
    parser.add_argument("ligands")
    parser.add_argument("checkpoint")
    parser.add_argument("output")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--contract", choices=("r0b", "r0c"), default="r0b")
    args = parser.parse_args()
    report = run_prefit(
        args.panel, args.geometry, args.proteins, args.ligands,
        args.checkpoint, args.output, device=args.device,
        contract=args.contract)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
