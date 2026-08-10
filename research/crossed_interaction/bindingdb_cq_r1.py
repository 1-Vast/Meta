"""Build and audit exact BindingDB Ki/Kd quotient labels."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import re
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np

from research.crossed_interaction.bindingdb_cq_r0 import (
    _panel_keys,
    iter_projection,
    sha256_file,
)


VALUE_COLUMNS = {"Ki": "Ki (nM)", "Kd": "Kd (nM)"}
EXACT_NUMBER = re.compile(r"^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


def _zip_rows(path: Path):
    archive = zipfile.ZipFile(path)
    names = [name for name in archive.namelist() if not name.endswith("/")]
    if len(names) != 1:
        archive.close()
        raise ValueError(f"expected one file in {path}")
    binary = archive.open(names[0])
    text = io.TextIOWrapper(binary, encoding="utf-8-sig", errors="replace", newline="")
    return archive, text, csv.DictReader(text, delimiter="\t")


def parse_exact_nm(text: str) -> float | None:
    value = text.strip()
    if not EXACT_NUMBER.fullmatch(value):
        return None
    number = float(value)
    return number if math.isfinite(number) and number > 0 else None


def extract_labels(articles: Path, projection: Path, output: Path) -> dict:
    metadata = {row["source_row_id"]: row for row in iter_projection(projection)}
    archive, handle, reader = _zip_rows(articles)
    accepted = []
    excluded = defaultdict(int)
    try:
        for row in reader:
            source_row_id = row["BindingDB Reactant_set_id"].strip()
            record = metadata.get(source_row_id)
            if record is None:
                continue
            for endpoint, column in VALUE_COLUMNS.items():
                raw = row[column].strip()
                if not raw:
                    continue
                value_nm = parse_exact_nm(raw)
                if value_nm is None:
                    excluded[f"{endpoint}_non_exact"] += 1
                    continue
                panel_id = _panel_keys(record, endpoint)["document_protocol"][0]
                ligand = record.get("ligand_inchikey")
                target = record.get("target_sequence_sha256")
                if record.get("chain_count") != 1 or not ligand or not target:
                    excluded[f"{endpoint}_identity_incomplete"] += 1
                    continue
                accepted.append(
                    {
                        "source_row_id": source_row_id,
                        "document_id": record["document_id"],
                        "panel_id": panel_id,
                        "endpoint": endpoint,
                        "target_id": target,
                        "ligand_id": ligand,
                        "pK": 9.0 - math.log10(value_nm),
                    }
                )
    finally:
        handle.close()
        archive.close()
    accepted.sort(key=lambda item: (item["endpoint"], item["panel_id"], item["source_row_id"]))
    output.parent.mkdir(parents=True, exist_ok=True)
    raw_file = output.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw_file, mtime=0)
    writer = io.TextIOWrapper(compressed, encoding="utf-8", newline="\n")
    try:
        for row in accepted:
            writer.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    finally:
        writer.close()
        raw_file.close()
    return {
        "schema": "BindingDB.CQ-R1.ExactLabels.v1",
        "articles_sha256": sha256_file(articles),
        "projection_sha256": sha256_file(projection),
        "labels_sha256": sha256_file(output),
        "accepted_rows": len(accepted),
        "accepted_by_endpoint": {
            endpoint: sum(row["endpoint"] == endpoint for row in accepted)
            for endpoint in VALUE_COLUMNS
        },
        "excluded": dict(sorted(excluded.items())),
        "direction": "pK=9-log10(nM); higher_is_stronger",
    }


def iter_labels(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def panel_residual(rows: list[dict]) -> dict | None:
    cells: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        cells[(row["target_id"], row["ligand_id"])].append(float(row["pK"]))
    targets = sorted({key[0] for key in cells})
    ligands = sorted({key[1] for key in cells})
    edges = sorted(cells)
    if not edges:
        return None
    t_index = {value: index for index, value in enumerate(targets)}
    l_index = {value: index for index, value in enumerate(ligands)}
    design = np.zeros((len(edges), 1 + len(targets) + len(ligands)), dtype=np.float64)
    response = np.empty(len(edges), dtype=np.float64)
    replicate_mean_variances = []
    for index, (target, ligand) in enumerate(edges):
        values = np.asarray(cells[(target, ligand)], dtype=np.float64)
        response[index] = values.mean()
        design[index, 0] = 1.0
        design[index, 1 + t_index[target]] = 1.0
        design[index, 1 + len(targets) + l_index[ligand]] = 1.0
        if len(values) >= 2:
            replicate_mean_variances.append(float(values.var(ddof=1) / len(values)))
    fitted = design @ np.linalg.lstsq(design, response, rcond=None)[0]
    residual = response - fitted
    design_rank = int(np.linalg.matrix_rank(design))
    retained_rank = len(edges) - design_rank
    if retained_rank <= 0:
        return None
    return {
        "edges": len(edges),
        "targets": len(targets),
        "ligands": len(ligands),
        "retained_rank": retained_rank,
        "rank_normalized_mse": float(residual @ residual / retained_rank),
        "orthogonality_error": float(np.max(np.abs(design.T @ residual))),
        "replicate_cells": len(replicate_mean_variances),
        "replicate_mean_variance_sum": float(sum(replicate_mean_variances)),
    }


def bootstrap_rms(document_mse: np.ndarray, draws: int = 10000, seed: int = 20260810):
    if not len(document_mse):
        return [0.0, 0.0]
    rng = np.random.default_rng(seed)
    values = np.empty(draws, dtype=np.float64)
    for start in range(0, draws, 500):
        count = min(500, draws - start)
        indices = rng.integers(0, len(document_mse), size=(count, len(document_mse)))
        values[start : start + count] = np.sqrt(document_mse[indices].mean(axis=1))
    return [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))]


def audit(labels: Path) -> dict:
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in iter_labels(labels):
        grouped[(row["endpoint"], row["panel_id"])].append(row)
    endpoint_panels: dict[str, list[dict]] = defaultdict(list)
    for (endpoint, panel_id), rows in grouped.items():
        result = panel_residual(rows)
        if result is None:
            continue
        result["panel_id"] = panel_id
        result["document_id"] = rows[0]["document_id"]
        endpoint_panels[endpoint].append(result)

    summaries = {}
    for endpoint in VALUE_COLUMNS:
        panels = endpoint_panels.get(endpoint, [])
        by_document: dict[str, list[float]] = defaultdict(list)
        for panel in panels:
            by_document[panel["document_id"]].append(panel["rank_normalized_mse"])
        document_mse = np.asarray(
            [np.mean(values) for _, values in sorted(by_document.items())], dtype=np.float64
        )
        quotient_rms = float(np.sqrt(document_mse.mean())) if len(document_mse) else 0.0
        total_rank = sum(panel["retained_rank"] for panel in panels)
        max_orthogonality = max(
            (panel["orthogonality_error"] for panel in panels), default=0.0
        )
        conditions = {
            "cycle_positive_panels_ge_50": len(panels) >= 50,
            "total_retained_rank_ge_1000": total_rank >= 1000,
            "document_units_ge_30": len(by_document) >= 30,
            "quotient_rms_ge_0_25": quotient_rms >= 0.25,
            "bootstrap_lcb_positive": bootstrap_rms(document_mse)[0] > 0,
            "orthogonality_le_1e_7": max_orthogonality <= 1e-7,
        }
        summaries[endpoint] = {
            "cycle_positive_panels": len(panels),
            "document_units": len(by_document),
            "total_retained_rank": total_rank,
            "quotient_rms": quotient_rms,
            "quotient_rms_ci95": bootstrap_rms(document_mse),
            "max_orthogonality_error": max_orthogonality,
            "replicate_cells": sum(panel["replicate_cells"] for panel in panels),
            "replicate_mean_variance_sum": sum(
                panel["replicate_mean_variance_sum"] for panel in panels
            ),
            "conditions": conditions,
            "development_trainable": all(conditions.values()),
        }
    primary_pass = summaries["Ki"]["development_trainable"]
    verdict = (
        "CQ_R1_DEVELOPMENT_INTERACTION_OBSERVED"
        if primary_pass
        else "CQ_R1_DEVELOPMENT_INTERACTION_NOT_OBSERVED"
    )
    return {
        "schema": "BindingDB.CQ-R1.InteractionAudit.v1",
        "labels_sha256": sha256_file(labels),
        "endpoint_summaries": summaries,
        "primary_endpoint": "Ki",
        "verdict": verdict,
        "development_training_authorized": primary_pass,
        "biological_claim_authorized": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    extract = sub.add_parser("extract")
    extract.add_argument("--articles", type=Path, required=True)
    extract.add_argument("--projection", type=Path, required=True)
    extract.add_argument("--labels", type=Path, required=True)
    extract.add_argument("--manifest", type=Path, required=True)
    run = sub.add_parser("audit")
    run.add_argument("--labels", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "extract":
        result = extract_labels(args.articles, args.projection, args.labels)
        destination = args.manifest
    else:
        result = audit(args.labels)
        destination = args.output
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
