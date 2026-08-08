"""Check the L0 location-estimand identifiability contract.

Design metadata only: the audited value-free `x0_metadata.sql` projection plus
the label-blind governed row index. No affinity value is selected or read.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

from research.e0_identifiability.run_eaff_x0 import load_metadata
from research.e0_identifiability.x0_contract import stable_hash
from scripts.source_affinity.common import sha256_file


STAGE = "P1R2B-E-AFF-L0_LOCATION_ESTIMAND_IDENTIFIABILITY"
MIN_STRATA_WITH_TWO_PROTEINS = 30
MIN_PROTEINS_IN_TWO_STRATA = 30
MIN_LARGEST_COMPONENT_FRACTION = 0.50
# Target-dependent context keys that must never enter the stratum.
EXCLUDED_FROM_STRATUM = ("component_accession", "variant_id", "variant_mutation",
                         "document_chembl_id", "assay_chembl_id", "task_id")


def _read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def stratum_key(row: dict, parameters: list[dict]) -> str:
    payload = {
        "assay_organism": row.get("assay_organism"),
        "bao_format": row.get("bao_format"),
        "cell_id": row.get("cell_id"),
        "tissue_id": row.get("tissue_id"),
        "subcellular_fraction": row.get("subcellular_fraction"),
        "relationship_type": row.get("relationship_type"),
        "parameters": parameters,
    }
    if set(payload) & set(EXCLUDED_FROM_STRATUM):
        raise RuntimeError("stratum payload carries a target or identifier field")
    return stable_hash(payload)


def largest_component_fraction(edges: list[tuple[str, str]], proteins: set[str]) -> float:
    adjacency: dict[str, set[str]] = defaultdict(set)
    for protein, stratum in edges:
        adjacency["p:" + protein].add("s:" + stratum)
        adjacency["s:" + stratum].add("p:" + protein)
    seen: set[str] = set()
    best = 0
    for node in adjacency:
        if node in seen:
            continue
        stack, component = [node], set()
        seen.add(node)
        while stack:
            current = stack.pop()
            component.add(current)
            for neighbour in adjacency[current]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        best = max(best, sum(1 for value in component if value.startswith("p:")))
    return best / len(proteins) if proteins else 0.0


def run(args) -> dict:
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(f"output exists: {output}")
    input_root = Path(args.input)
    release_root = Path(args.release)

    rows = list(_read_jsonl(input_root / "rows.label_blind.jsonl"))
    by_activity = {int(row["activity_id"]): row for row in rows}

    release_manifest = json.loads(
        (release_root / "release_manifest.json").read_text(encoding="utf-8"))
    database = release_root / release_manifest["sqlite"]["path"]
    if not database.is_file() or database.stat().st_size != release_manifest["sqlite"]["bytes"]:
        raise RuntimeError("pinned SQLite disagrees with release manifest")
    metadata, parameters = load_metadata(database, Path(args.sql), set(by_activity))

    endpoints = {}
    strata_by_activity: dict[int, str] = {}
    for endpoint in ("Ki", "Kd"):
        edges: set[tuple[str, str]] = set()
        proteins_by_stratum: dict[str, set[str]] = defaultdict(set)
        strata_by_protein: dict[str, set[str]] = defaultdict(set)
        for activity, source in by_activity.items():
            if source["endpoint_family"] != endpoint:
                continue
            row = metadata[activity]
            key = stratum_key(row, parameters.get(int(row["assay_id"]), []))
            strata_by_activity[activity] = key
            protein = source["protein_sequence_sha256"]
            edges.add((protein, key))
            proteins_by_stratum[key].add(protein)
            strata_by_protein[protein].add(key)
        proteins = set(strata_by_protein)
        c1 = sum(1 for value in proteins_by_stratum.values() if len(value) >= 2)
        c2 = sum(1 for value in strata_by_protein.values() if len(value) >= 2)
        c3 = largest_component_fraction(sorted(edges), proteins)
        endpoints[endpoint] = {
            "proteins": len(proteins),
            "strata": len(proteins_by_stratum),
            "cells": len(edges),
            "C1_strata_with_two_or_more_proteins": c1,
            "C2_proteins_in_two_or_more_strata": c2,
            "C3_largest_component_protein_fraction": round(c3, 5),
            "C1_pass": c1 >= MIN_STRATA_WITH_TWO_PROTEINS,
            "C2_pass": c2 >= MIN_PROTEINS_IN_TWO_STRATA,
            "C3_pass": c3 >= MIN_LARGEST_COMPONENT_FRACTION,
        }
        endpoints[endpoint]["admitted"] = all(
            endpoints[endpoint][name] for name in ("C1_pass", "C2_pass", "C3_pass"))

    admitted = [name for name in ("Ki", "Kd") if endpoints[name]["admitted"]]
    verdict = ("L0_LOCATION_ESTIMAND_IDENTIFIED_" + "_".join(name.upper() for name in admitted)
               if admitted else "L0_NOT_RUN_LOCATION_ESTIMAND_NOT_IDENTIFIED")

    output.mkdir(parents=True)
    report = {
        "schema": "MetaSieve.EAffL0Estimand.v1",
        "stage": STAGE,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": verdict,
        "admitted_endpoints": admitted,
        "affinity_value_fields_selected": 0,
        "affinity_values_read": False,
        "davis_label_reads": 0,
        "recipient_label_reads": 0,
        "stratum_definition": (
            "target-independent context: assay_organism, bao_format, cell_id, "
            "tissue_id, subcellular_fraction, relationship_type, assay_parameters"),
        "excluded_from_stratum": list(EXCLUDED_FROM_STRATUM),
        "frozen_thresholds": {
            "C1_min_strata_with_two_or_more_proteins": MIN_STRATA_WITH_TWO_PROTEINS,
            "C2_min_proteins_in_two_or_more_strata": MIN_PROTEINS_IN_TWO_STRATA,
            "C3_min_largest_component_protein_fraction": MIN_LARGEST_COMPONENT_FRACTION,
        },
        "endpoints": endpoints,
        "input_rows": len(rows),
    }
    (output / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (output / "strata.jsonl").open("w", encoding="utf-8", newline="\n") as handle:
        for activity in sorted(strata_by_activity):
            handle.write(json.dumps(
                {"activity_id": activity, "stratum_id": strata_by_activity[activity]},
                sort_keys=True, separators=(",", ":")) + "\n")
    (output / "manifest.json").write_text(json.dumps({
        "stage": STAGE,
        "inputs": {
            "label_blind_rows": sha256_file(input_root / "rows.label_blind.jsonl"),
            "sql": sha256_file(Path(args.sql)),
            "data_contract": sha256_file(
                Path(__file__).with_name("EAFF_L0_DATA_CONTRACT.md")),
        },
        "pinned_sqlite": {"path": release_manifest["sqlite"]["path"],
                          "sha256": release_manifest["sqlite"]["sha256"]},
        "outputs": {name: sha256_file(output / name)
                    for name in ("report.json", "strata.jsonl")},
        "label_reads": {"affinity_values": 0, "davis": 0, "recipient": 0},
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="dataset/processed/source_affinity/e0_input_v1")
    parser.add_argument("--release", default="dataset/raw/source_affinity/chembl37_sqlite_v1")
    parser.add_argument("--sql", default="research/e0_identifiability/x0_metadata.sql")
    parser.add_argument("--output",
                        default="research/e0_identifiability/artifacts/eaff_l0_estimand_v1")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), indent=2, sort_keys=True))
