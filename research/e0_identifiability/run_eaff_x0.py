"""Run the label-blind E-AFF-X0 crossed target-ligand census."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import sqlite3

from research.e0_identifiability.x0_contract import (
    UnionFind,
    crossed_panel_payload,
    rectangle_count,
    stable_hash,
)
from scripts.source_affinity.common import sha256_file


STAGE = "P1R2B-E-AFF-X0_CROSSED_TARGET_LIGAND_CENSUS"
REQUIRED_COMPONENTS = 245
PROHIBITED_SQL = ("act.standard_value", "act.value", "pchembl_value")


def _read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_jsonl(path: Path, rows) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def load_metadata(database: Path, sql_path: Path,
                  wanted: set[int]) -> tuple[dict[int, dict], dict[int, list[dict]]]:
    sql = sql_path.read_text(encoding="utf-8")
    lowered = sql.lower()
    if any(value in lowered for value in PROHIBITED_SQL):
        raise RuntimeError("X0 SQL selects a prohibited affinity field")
    connection = sqlite3.connect(f"file:{database.as_posix()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        connection.execute("CREATE TEMP TABLE selected_activity(activity_id INTEGER PRIMARY KEY)")
        connection.executemany("INSERT INTO selected_activity VALUES (?)",
                               ((value,) for value in sorted(wanted)))
        rows = {int(row["activity_id"]): dict(row) for row in connection.execute(sql)}
        if set(rows) != wanted:
            raise RuntimeError(f"metadata mapping mismatch: {len(rows)} != {len(wanted)}")
        assays = sorted({int(row["assay_id"]) for row in rows.values()})
        connection.execute("CREATE TEMP TABLE selected_assay(assay_id INTEGER PRIMARY KEY)")
        connection.executemany("INSERT INTO selected_assay VALUES (?)", ((value,) for value in assays))
        parameters: dict[int, list[dict]] = defaultdict(list)
        query = """
            SELECT parameter.assay_id, parameter.standard_type,
                   parameter.standard_relation, parameter.standard_value,
                   parameter.standard_units, parameter.standard_text_value
            FROM assay_parameters AS parameter
            JOIN selected_assay AS selected ON selected.assay_id = parameter.assay_id
            ORDER BY parameter.assay_id, parameter.standard_type,
                     parameter.standard_value, parameter.standard_units,
                     parameter.standard_text_value
        """
        for row in connection.execute(query):
            parameters[int(row[0])].append({
                "type": row[1], "relation": row[2], "value": row[3],
                "units": row[4], "text_value": row[5],
            })
    finally:
        connection.close()
    return rows, parameters


def _dependency_components(panel_records: list[dict], endpoint: str,
                           replicate: bool) -> list[dict]:
    union = UnionFind()
    accepted = []
    rectangle_key = "replicate_rectangles" if replicate else "rectangles"
    closure_key = "replicate_closure_components" if replicate else "closure_components"
    for panel in panel_records:
        if panel["endpoint_family"] != endpoint or panel[rectangle_key] <= 0:
            continue
        node = "panel:" + panel["panel_id"]
        union.add(node)
        for closure in panel[closure_key]:
            union.union(node, "closure:" + closure)
        accepted.append(panel)
    grouped: dict[str, list[dict]] = defaultdict(list)
    for panel in accepted:
        grouped[union.find("panel:" + panel["panel_id"])].append(panel)
    result = []
    for panels in grouped.values():
        panel_ids = sorted(panel["panel_id"] for panel in panels)
        closures = sorted({value for panel in panels for value in panel[closure_key]})
        result.append({
            "endpoint_family": endpoint,
            "stratum": "replicate_supported" if replicate else "all_rectangles",
            "dependency_component_id": stable_hash([endpoint, replicate, panel_ids, closures]),
            "panels": panel_ids,
            "panel_count": len(panel_ids),
            "rectangle_count": int(sum(panel[rectangle_key] for panel in panels)),
            "closure_components": closures,
        })
    return sorted(result, key=lambda row: row["dependency_component_id"])


def run(args) -> dict:
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(f"output exists: {output}")
    output.mkdir(parents=True)
    input_root = Path(args.input)
    label_blind_rows = list(_read_jsonl(input_root / "rows.label_blind.jsonl"))
    by_activity = {int(row["activity_id"]): row for row in label_blind_rows}
    if len(by_activity) != len(label_blind_rows):
        raise RuntimeError("label-blind activity IDs are not unique")

    release_root = Path(args.release)
    release_manifest = json.loads((release_root / "release_manifest.json").read_text(encoding="utf-8"))
    database = release_root / release_manifest["sqlite"]["path"]
    if not database.is_file() or database.stat().st_size != release_manifest["sqlite"]["bytes"]:
        raise RuntimeError("pinned SQLite disagrees with release manifest")
    metadata, parameters = load_metadata(database, Path(args.sql), set(by_activity))

    panels: dict[str, dict] = {}
    excluded_variants = 0
    for activity, source in by_activity.items():
        row = metadata[activity]
        if row["variant_id"] is not None:
            excluded_variants += 1
            continue
        payload = crossed_panel_payload(row, parameters.get(int(row["assay_id"]), []))
        panel_id = stable_hash(payload)
        panel = panels.setdefault(panel_id, {
            "panel_id": panel_id,
            "endpoint_family": row["endpoint_family"],
            "payload": payload,
            "cells": defaultdict(lambda: {
                "activity_ids": [], "assay_ids": defaultdict(list),
                "closure_component_id": None,
            }),
            "assay_description_sha256": {},
        })
        cell_key = (source["protein_sequence_sha256"], source["ligand_connectivity_key"])
        cell = panel["cells"][cell_key]
        cell["activity_ids"].append(activity)
        cell["assay_ids"][str(row["assay_chembl_id"])].append(activity)
        cell["closure_component_id"] = source["closure_component_id"]
        panel["assay_description_sha256"][str(row["assay_chembl_id"])] = stable_hash(
            row.get("assay_description"))

    panel_records, cell_records = [], []
    for panel_id, panel in sorted(panels.items()):
        target_ligands: dict[str, set[str]] = defaultdict(set)
        replicate_target_ligands: dict[str, set[str]] = defaultdict(set)
        closures_by_target = {}
        for (protein, ligand), cell in panel["cells"].items():
            target_ligands[protein].add(ligand)
            closures_by_target[protein] = cell["closure_component_id"]
            replicate = max(map(len, cell["assay_ids"].values())) >= 2
            if replicate:
                replicate_target_ligands[protein].add(ligand)
            cell_records.append({
                "panel_id": panel_id,
                "endpoint_family": panel["endpoint_family"],
                "protein_sequence_sha256": protein,
                "closure_component_id": cell["closure_component_id"],
                "ligand_connectivity_key": ligand,
                "activity_ids": sorted(cell["activity_ids"]),
                "assay_activity_ids": {key: sorted(value)
                                       for key, value in sorted(cell["assay_ids"].items())},
                "exact_assay_replicate_supported": replicate,
            })
        target_pairs, rectangles = rectangle_count(target_ligands)
        replicate_pairs, replicate_rectangles = rectangle_count(replicate_target_ligands)
        if not rectangles:
            continue
        participating, replicate_participating = set(), set()
        targets = sorted(target_ligands)
        for index, left in enumerate(targets):
            for right in targets[index + 1:]:
                if len(target_ligands[left] & target_ligands[right]) >= 2:
                    participating.update((left, right))
                if len(replicate_target_ligands[left] & replicate_target_ligands[right]) >= 2:
                    replicate_participating.update((left, right))
        panel_records.append({
            "panel_id": panel_id,
            "endpoint_family": panel["endpoint_family"],
            "payload": panel["payload"],
            "target_count": len(target_ligands),
            "ligand_count": len({value for values in target_ligands.values() for value in values}),
            "cell_count": len(panel["cells"]),
            "target_pairs": target_pairs,
            "rectangles": rectangles,
            "replicate_target_pairs": replicate_pairs,
            "replicate_rectangles": replicate_rectangles,
            "closure_components": sorted({closures_by_target[value] for value in participating}),
            "replicate_closure_components": sorted(
                {closures_by_target[value] for value in replicate_participating}),
            "assay_description_sha256": panel["assay_description_sha256"],
        })

    dependencies = []
    endpoint_reports = {}
    for endpoint in ("Ki", "Kd"):
        total = _dependency_components(panel_records, endpoint, replicate=False)
        replicated = _dependency_components(panel_records, endpoint, replicate=True)
        dependencies.extend(total + replicated)
        endpoint_panels = [row for row in panel_records if row["endpoint_family"] == endpoint]
        endpoint_cells = [row for row in cell_records
                          if row["endpoint_family"] == endpoint
                          and row["panel_id"] in {panel["panel_id"] for panel in endpoint_panels}]
        total_rectangles = sum(row["rectangles"] for row in endpoint_panels)
        largest = max((row["rectangle_count"] for row in total), default=0)
        endpoint_reports[endpoint] = {
            "panels": len(endpoint_panels),
            "dependency_components": len(total),
            "replicate_supported_dependency_components": len(replicated),
            "targets": len({row["protein_sequence_sha256"] for row in endpoint_cells}),
            "ligands": len({row["ligand_connectivity_key"] for row in endpoint_cells}),
            "cells": len(endpoint_cells),
            "target_pairs": sum(row["target_pairs"] for row in endpoint_panels),
            "rectangles": total_rectangles,
            "replicate_supported_rectangles": sum(
                row["replicate_rectangles"] for row in endpoint_panels),
            "largest_dependency_rectangle_fraction": (
                largest / total_rectangles if total_rectangles else 0.0),
        }

    noise_supported = [name for name, values in endpoint_reports.items()
                       if values["replicate_supported_dependency_components"] >= REQUIRED_COMPONENTS]
    additive_supported = [name for name, values in endpoint_reports.items()
                          if values["dependency_components"] >= REQUIRED_COMPONENTS]
    if noise_supported:
        verdict = "X0_NOISE_AWARE_CROSSED_DATA_SUPPORTED_" + "_".join(name.upper() for name in noise_supported)
    elif additive_supported:
        verdict = "X0_ADDITIVE_NULL_ONLY_CROSSED_DATA_SUPPORTED_" + "_".join(
            name.upper() for name in additive_supported)
    else:
        verdict = "STOP_SOURCE_INTERACTION_UNDERDETERMINED"

    _write_jsonl(output / "panels.jsonl", panel_records)
    eligible_panel_ids = {row["panel_id"] for row in panel_records}
    _write_jsonl(output / "cells.jsonl",
                 (row for row in cell_records if row["panel_id"] in eligible_panel_ids))
    _write_jsonl(output / "dependency_components.jsonl", dependencies)
    report = {
        "schema": "MetaSieve.EAffX0.v1", "stage": STAGE,
        "created_utc": datetime.now(timezone.utc).isoformat(), "verdict": verdict,
        "label_blind": True, "affinity_value_fields_selected": 0,
        "affinity_values_read": False, "davis_label_reads": 0, "recipient_label_reads": 0,
        "input_rows": len(label_blind_rows), "variant_rows_excluded": excluded_variants,
        "panel_contract": "document_endpoint_target_independent_structured_context_no_variants",
        "required_effective_components": REQUIRED_COMPONENTS,
        "power_design": {"interaction_rms_over_assay_noise": 0.5,
                         "alpha_one_sided": 0.05, "power": 0.80,
                         "ideal_chi_square_required_n": REQUIRED_COMPONENTS},
        "endpoints": endpoint_reports,
        "interpretation_limits": [
            "same structured context does not prove identical wet-lab protocol",
            "rectangle counts are dependent and never used as IID sample size",
            "no-replicate panels cannot identify measurement-noise variance",
            "X0 does not access or test affinity interaction values",
        ],
    }
    _write_json(output / "report.json", report)
    manifest = {
        "stage": STAGE, "label_blind": True,
        "inputs": {
            "release_manifest": sha256_file(release_root / "release_manifest.json"),
            "label_blind_rows": sha256_file(input_root / "rows.label_blind.jsonl"),
            "input_manifest": sha256_file(input_root / "manifest.json"),
            "sql": sha256_file(Path(args.sql)),
            "preregistration": sha256_file(Path(__file__).with_name("EAFF_X0_PREREGISTRATION.md")),
        },
        "pinned_sqlite": {"path": release_manifest["sqlite"]["path"],
                          "sha256": release_manifest["sqlite"]["sha256"],
                          "bytes": release_manifest["sqlite"]["bytes"]},
        "outputs": {name: sha256_file(output / name) for name in (
            "panels.jsonl", "cells.jsonl", "dependency_components.jsonl", "report.json")},
        "label_reads": {"affinity_values": 0, "davis": 0, "recipient": 0},
    }
    _write_json(output / "manifest.json", manifest)
    return report


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="dataset/processed/source_affinity/e0_input_v1")
    parser.add_argument("--release", default="dataset/raw/source_affinity/chembl37_sqlite_v1")
    parser.add_argument("--sql", default="research/e0_identifiability/x0_metadata.sql")
    parser.add_argument("--output", default="research/e0_identifiability/artifacts/eaff_x0_v1")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), indent=2, sort_keys=True))
