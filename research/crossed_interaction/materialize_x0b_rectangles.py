"""Materialize the original label-blind X0-B cell-disjoint rectangles."""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "research" / "crossed_interaction" / "recovered" / "eaff__x0_v1_cells.jsonl"
OUTPUT = ROOT / "research" / "crossed_interaction" / "artifacts" / "x1a_r_direct_dd"
CONTRACT = ROOT / "research" / "crossed_interaction" / "PREREG_X1A_R_DIRECT_DD_DEPENDENCE.md"
CELLS_SHA256 = "898df88235401a2be2341ae1ab222e6c5903202796c8312d8e9091cf76741562"
EXPECTED = {
    "Ki": {"rectangles": 11168, "clusters": 36, "cap": 32},
    "Kd": {"rectangles": 1041, "clusters": 12, "cap": 125},
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def pack_cell_disjoint(target_ligands: dict[str, set[str]]) -> list[tuple[str, str, str, str]]:
    available = set(target_ligands)
    packed = []
    while len(available) >= 2:
        targets = sorted(available)
        best = None
        for index, left in enumerate(targets):
            for right in targets[index + 1:]:
                shared = len(target_ligands[left] & target_ligands[right])
                if shared >= 2 and (best is None or shared > best[0]):
                    best = (shared, left, right)
        if best is None:
            break
        _, left, right = best
        available.remove(left)
        available.remove(right)
        common = sorted(target_ligands[left] & target_ligands[right])
        for offset in range(0, len(common) - 1, 2):
            packed.append((left, right, common[offset], common[offset + 1]))
    return packed


def rectangle_id(endpoint: str, panel: str, values: tuple[str, str, str, str]) -> str:
    payload = "|".join((endpoint, panel, *values)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def build_rows(cells: list[dict]) -> list[dict]:
    targets = defaultdict(lambda: defaultdict(set))
    cell_by_key = {}
    panel_meta = {}
    for cell in cells:
        panel = cell["panel_id"]
        target = cell["protein_sequence_sha256"]
        ligand = cell["ligand_connectivity_key"]
        targets[panel][target].add(ligand)
        cell_by_key[(panel, target, ligand)] = cell
        meta = (cell["endpoint_family"], cell["closure_component_id"])
        if panel in panel_meta and panel_meta[panel] != meta:
            raise RuntimeError(f"inconsistent panel metadata: {panel}")
        panel_meta[panel] = meta

    rows = []
    for panel in sorted(targets):
        endpoint, cluster = panel_meta[panel]
        for values in pack_cell_disjoint(targets[panel]):
            left, right, first, second = values
            keys = ((panel, left, first), (panel, left, second),
                    (panel, right, first), (panel, right, second))
            if len(set(keys)) != 4 or any(key not in cell_by_key for key in keys):
                raise RuntimeError("invalid packed rectangle")
            rows.append({
                "rectangle_id": rectangle_id(endpoint, panel, values),
                "endpoint": endpoint,
                "dependency_cluster": cluster,
                "panel_id": panel,
                "protein_a": left,
                "protein_b": right,
                "ligand_a": first,
                "ligand_b": second,
                "cells": [
                    {"sign": 1, **cell_by_key[keys[0]]},
                    {"sign": -1, **cell_by_key[keys[1]]},
                    {"sign": -1, **cell_by_key[keys[2]]},
                    {"sign": 1, **cell_by_key[keys[3]]},
                ],
            })

    by_cluster = defaultdict(list)
    for row in rows:
        by_cluster[(row["endpoint"], row["dependency_cluster"])].append(row)
    for (endpoint, _cluster), members in by_cluster.items():
        selected = {row["rectangle_id"] for row in
                    sorted(members, key=lambda row: row["rectangle_id"])
                    [:EXPECTED[endpoint]["cap"]]}
        for row in members:
            row["selected_at_frozen_cap"] = row["rectangle_id"] in selected
    return sorted(rows, key=lambda row: (row["endpoint"], row["rectangle_id"]))


def run(output: Path = OUTPUT) -> dict:
    if output.exists():
        raise FileExistsError(f"no-clobber output exists: {output}")
    if sha256_file(INPUT) != CELLS_SHA256:
        raise RuntimeError("X0 cells hash mismatch")
    rows = build_rows(list(read_jsonl(INPUT)))

    counts = {}
    for endpoint, expected in EXPECTED.items():
        endpoint_rows = [row for row in rows if row["endpoint"] == endpoint]
        clusters = {row["dependency_cluster"] for row in endpoint_rows}
        if len(endpoint_rows) != expected["rectangles"] or len(clusters) != expected["clusters"]:
            raise RuntimeError(f"X0-B reproduction failed for {endpoint}")
        counts[endpoint] = {
            "rectangles": len(endpoint_rows),
            "clusters": len(clusters),
            "cap": expected["cap"],
            "selected_rectangles": sum(row["selected_at_frozen_cap"] for row in endpoint_rows),
        }

    output.mkdir(parents=True)
    rectangles = output / "rectangles.jsonl"
    with rectangles.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    manifest = {
        "stage": "E-AFF-X1A-R_DIRECT_DD_DEPENDENCE",
        "label_blind": True,
        "input_cells_sha256": CELLS_SHA256,
        "contract_sha256": sha256_file(CONTRACT),
        "runner_sha256": sha256_file(Path(__file__)),
        "rectangles_sha256": sha256_file(rectangles),
        "counts": counts,
        "affinity_value_fields_selected": 0,
        "training_performed": False,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
