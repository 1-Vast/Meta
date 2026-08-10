"""Label-blind cycle-space census for the frozen ChEMBL crossed panels."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "research" / "crossed_interaction" / "recovered" / "eaff__x0_v1_cells.jsonl"
INPUT_SHA256 = "898df88235401a2be2341ae1ab222e6c5903202796c8312d8e9091cf76741562"
OUTPUT = ROOT / "report" / "crossed_interaction" / "cycle_quotient_feasibility"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def graph_stats(edges: set[tuple[str, str]]) -> dict:
    proteins = {protein for protein, _ in edges}
    ligands = {ligand for _, ligand in edges}
    adjacency = defaultdict(set)
    for protein, ligand in edges:
        left, right = ("P", protein), ("L", ligand)
        adjacency[left].add(right)
        adjacency[right].add(left)
    seen = set()
    components = 0
    for vertex in adjacency:
        if vertex in seen:
            continue
        components += 1
        stack = [vertex]
        seen.add(vertex)
        while stack:
            current = stack.pop()
            for neighbour in adjacency[current]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
    dimension = len(edges) - len(proteins) - len(ligands) + components
    if dimension < 0:
        raise RuntimeError("negative cycle-space dimension")
    return {
        "edges": len(edges),
        "proteins": len(proteins),
        "ligands": len(ligands),
        "graph_components": components,
        "cycle_dimension": dimension,
    }


def census(rows: list[dict]) -> dict:
    panels = defaultdict(set)
    panel_cluster = {}
    assays = defaultdict(set)
    for row in rows:
        endpoint = row["endpoint_family"]
        panel_key = (endpoint, row["panel_id"])
        edge = (row["protein_sequence_sha256"], row["ligand_connectivity_key"])
        panels[panel_key].add(edge)
        cluster = row["closure_component_id"]
        if panel_key in panel_cluster and panel_cluster[panel_key] != cluster:
            raise RuntimeError("panel spans multiple frozen dependency components")
        panel_cluster[panel_key] = cluster
        for assay_id in row["assay_activity_ids"]:
            assays[(endpoint, assay_id)].add(edge)

    result = {}
    for endpoint in ("Ki", "Kd"):
        panel_records = []
        cluster_dimensions = Counter()
        for (family, panel_id), edges in sorted(panels.items()):
            if family != endpoint:
                continue
            record = graph_stats(edges)
            record.update({"panel_id": panel_id, "dependency_component": panel_cluster[(family, panel_id)]})
            panel_records.append(record)
            cluster_dimensions[record["dependency_component"]] += record["cycle_dimension"]
        assay_records = [graph_stats(edges) for (family, _), edges in assays.items()
                         if family == endpoint]
        total_dimension = sum(record["cycle_dimension"] for record in panel_records)
        positive = [record for record in panel_records if record["cycle_dimension"] > 0]
        result[endpoint] = {
            "panel_graphs": len(panel_records),
            "cycle_positive_panels": len(positive),
            "raw_panel_cycle_dimension": total_dimension,
            "median_positive_panel_dimension": sorted(record["cycle_dimension"] for record in positive)[len(positive) // 2],
            "dependency_components": len(cluster_dimensions),
            "cycle_positive_dependency_components": sum(value > 0 for value in cluster_dimensions.values()),
            "largest_dependency_component_share": max(cluster_dimensions.values()) / total_dimension,
            "exact_assay_graphs": len(assay_records),
            "exact_assay_multi_target_graphs": sum(record["proteins"] > 1 for record in assay_records),
            "exact_assay_cycle_dimension": sum(record["cycle_dimension"] for record in assay_records),
            "top_dependency_component_dimensions": cluster_dimensions.most_common(5),
            "top_panel_dimensions": sorted(panel_records, key=lambda record: record["cycle_dimension"], reverse=True)[:5],
        }
    return result


def run(output: Path = OUTPUT) -> dict:
    if output.exists():
        raise FileExistsError(f"no-clobber output exists: {output}")
    if sha256_file(INPUT) != INPUT_SHA256:
        raise RuntimeError("frozen X0 cell hash mismatch")
    with INPUT.open(encoding="utf-8") as handle:
        rows = [json.loads(line) for line in handle if line.strip()]
    endpoints = census(rows)
    verdict = "CYCLE_QUOTIENT_ALGEBRAICALLY_AVAILABLE_BUT_DEPENDENCY_NOT_REPAIRED"
    output.mkdir(parents=True)
    report = {
        "stage": "E-AFF-CQ-R0_LABEL_BLIND_CYCLE_QUOTIENT_CENSUS",
        "status": "development_feasibility_audit",
        "terminal_verdict": verdict,
        "endpoints": endpoints,
        "input_sha256": INPUT_SHA256,
        "runner_sha256": sha256_file(Path(__file__)),
        "affinity_value_reads": 0,
        "training_performed": False,
        "gpu_used": False,
        "interpretation": (
            "Panel cycle spaces are large, but exact assays have zero crossed cycle dimension "
            "and panel coordinates remain concentrated in the same 36/12 dependency components."
        ),
        "training_authorized": False,
    }
    (output / "census.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                                        encoding="utf-8")
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    report = run(args.output)
    print(json.dumps({"terminal_verdict": report["terminal_verdict"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
