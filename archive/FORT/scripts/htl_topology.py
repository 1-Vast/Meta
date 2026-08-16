"""Read-only target-frequency topology audit for HTL-DTA."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


ROOT = Path("dataset/public/chembl_37/processed/dualcold")
REGISTRY = ROOT / "registry.parquet"
MANIFEST = ROOT / "manifest.json"
DEFAULT_OUTPUT = Path("dataset/processed/htl_target_topology.v1.json")
ENDPOINTS = ("pKi", "pKd")
THRESHOLD_GRID = (5, 10, 30, 100, 300)
K_GRID = (0, 1, 3, 5, 10, 20)
SOURCE_THRESHOLD = 100
RECIPIENT_THRESHOLD = 30
METADATA_COLUMNS = (
    "target",
    "conn",
    "endpoint",
    "n_records",
    "scaffold",
    "assays",
    "docs",
    "accession",
    "hcluster",
    "dual_cold_split",
)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tokens(value: object) -> tuple[str, ...]:
    """Return deterministic provenance tokens from the registry's pipe field."""

    values = {part.strip() for part in str(value).split("|") if part.strip()}
    if not values:
        raise ValueError("document and assay fields must contain at least one token")
    return tuple(sorted(values))


def token_key(value: object) -> str:
    return "|".join(tokens(value))


def quantiles(values: Iterable[int]) -> dict[str, float]:
    series = pd.Series(tuple(values), dtype="float64")
    if series.empty:
        return {name: 0.0 for name in ("min", "q25", "median", "q75", "q90", "max")}
    return {
        "min": float(series.min()),
        "q25": float(series.quantile(0.25)),
        "median": float(series.median()),
        "q75": float(series.quantile(0.75)),
        "q90": float(series.quantile(0.90)),
        "max": float(series.max()),
    }


def _target_record(group: pd.DataFrame) -> dict[str, Any]:
    target = str(group.target.iloc[0])
    endpoint = str(group.endpoint.iloc[0])
    if group.endpoint.nunique() != 1:
        raise ValueError(f"target {target} mixes endpoint strata")

    provenance_units = {
        (
            str(row.conn),
            token_key(row.docs),
            token_key(row.assays),
        )
        for row in group.itertuples(index=False)
    }
    documents = {document for value in group.docs for document in tokens(value)}
    assays = {assay for value in group.assays for assay in tokens(value)}
    components = sorted({str(value) for value in group.hcluster})
    if not components:
        raise ValueError(f"target {target} has no homology component")

    unique_scaffolds = int(group.scaffold.nunique())
    n_eff = len(provenance_units)
    return {
        "target": target,
        "endpoint": endpoint,
        "homology_components": components,
        "raw_registry_rows": int(len(group)),
        "replicate_records": int(group.n_records.sum()),
        "n_eff": int(n_eff),
        "unique_parent_connectivity": int(group.conn.nunique()),
        "unique_scaffolds": unique_scaffolds,
        "unique_documents": len(documents),
        "unique_assays": len(assays),
        "unique_accessions": int(group.accession.nunique()),
        "scaffold_closed_query_depth_upper_bound": {
            str(k): max(0, unique_scaffolds - k) for k in K_GRID
        },
        "provenance_closed_query_depth_upper_bound": {
            str(k): max(0, n_eff - k) for k in K_GRID
        },
    }


def _threshold_summary(records: list[dict[str, Any]], threshold: int) -> dict[str, Any]:
    tail = [item for item in records if item["n_eff"] < threshold]
    head = [item for item in records if item["n_eff"] >= threshold]
    analyzable = [item for item in tail if item["unique_scaffolds"] >= 6 and item["n_eff"] >= 6]
    return {
        "tail_rule": f"n_eff < {threshold}",
        "head_rule": f"n_eff >= {threshold}",
        "tail_targets": len(tail),
        "head_targets": len(head),
        "tail_homology_components": len({c for item in tail for c in item["homology_components"]}),
        "tail_k5_scaffold_and_provenance_analyzable": len(analyzable),
        "tail_k5_query_depth_upper_bound_median": float(
            pd.Series(
                [item["scaffold_closed_query_depth_upper_bound"]["5"] for item in analyzable],
                dtype="float64",
            ).median()
        )
        if analyzable
        else 0.0,
    }


def _transfer_summary(records: list[dict[str, Any]]) -> dict[str, Any]:
    """Describe candidate abundant-source/scarce-recipient topology only."""

    sources = [item for item in records if item["n_eff"] >= SOURCE_THRESHOLD]
    recipients = [item for item in records if item["n_eff"] < RECIPIENT_THRESHOLD]
    source_components = {
        component for item in sources for component in item["homology_components"]
    }
    same_component = [
        item
        for item in recipients
        if set(item["homology_components"]).intersection(source_components)
    ]
    analyzable = [
        item
        for item in recipients
        if item["unique_scaffolds"] >= 6 and item["n_eff"] >= 6
    ]
    return {
        "source_rule": f"n_eff >= {SOURCE_THRESHOLD}",
        "recipient_rule": f"n_eff < {RECIPIENT_THRESHOLD}",
        "source_targets": len(sources),
        "recipient_targets": len(recipients),
        "source_homology_components": len(source_components),
        "recipient_targets_with_same_homology_source": len(same_component),
        "recipient_targets_homology_cold_to_source": len(recipients) - len(same_component),
        "recipient_k5_scaffold_and_provenance_upper_bound_analyzable": len(analyzable),
    }


def build_audit(
    frame: pd.DataFrame,
    *,
    split: str,
    registry_sha256: str,
    manifest_sha256: str,
) -> dict[str, Any]:
    if split not in {"train", "development"}:
        raise ValueError("HTL topology audit only accepts train or development")
    missing = sorted(set(METADATA_COLUMNS).difference(frame.columns))
    if missing:
        raise ValueError(f"registry is missing metadata columns: {missing}")
    if frame.empty:
        raise ValueError(f"registry contains no rows for split={split}")
    if frame.dual_cold_split.astype(str).nunique() != 1:
        raise ValueError("filtered registry contains more than one split")
    if not bool(frame.endpoint.astype(str).isin(ENDPOINTS).all()):
        unexpected = sorted(set(frame.endpoint.astype(str)).difference(ENDPOINTS))
        raise ValueError(f"unexpected endpoint in exact-affinity audit: {unexpected}")

    target_records = [
        _target_record(group)
        for (_, _), group in frame.groupby(["endpoint", "target"], sort=True, dropna=False)
    ]
    by_endpoint: dict[str, Any] = {}
    for endpoint in ENDPOINTS:
        records = [item for item in target_records if item["endpoint"] == endpoint]
        by_endpoint[endpoint] = {
            "target_count": len(records),
            "registry_rows": int(sum(item["raw_registry_rows"] for item in records)),
            "n_eff_quantiles": quantiles(item["n_eff"] for item in records),
            "unique_scaffold_quantiles": quantiles(item["unique_scaffolds"] for item in records),
            "unique_document_quantiles": quantiles(item["unique_documents"] for item in records),
            "unique_assay_quantiles": quantiles(item["unique_assays"] for item in records),
            "homology_components": len({c for item in records for c in item["homology_components"]}),
            "candidate_thresholds": {
                str(threshold): _threshold_summary(records, threshold)
                for threshold in THRESHOLD_GRID
            },
            "candidate_transfer_topology": _transfer_summary(records),
            "targets": records,
        }

    document_counts = Counter(
        document
        for value in frame.docs
        for document in tokens(value)
    )
    return {
        "schema_version": "htl-target-topology-v1",
        "program": "HTL_DTA",
        "stage": "HTL_1_TARGET_FREQUENCY_TOPOLOGY",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "label_policy": "metadata-only; affinity, replicate_sd, and confirmation rows are not read",
        "input": {
            "registry": str(REGISTRY),
            "registry_sha256": registry_sha256,
            "manifest": str(MANIFEST),
            "manifest_sha256": manifest_sha256,
            "split": split,
            "rows": int(len(frame)),
            "columns_read": list(METADATA_COLUMNS),
        },
        "definitions": {
            "n_eff": "unique (target, endpoint, parent connectivity, canonical document-token set, canonical assay-token set) units",
            "scaffold_closed_query_depth_upper_bound": "max(0, unique scaffolds - support k); this is an upper bound before document/assay closure",
            "provenance_closed_query_depth_upper_bound": "max(0, n_eff - support k); this is an upper bound and not an episode admission",
            "candidate_threshold_grid": list(THRESHOLD_GRID),
            "candidate_transfer_source_threshold": SOURCE_THRESHOLD,
            "candidate_transfer_recipient_threshold": RECIPIENT_THRESHOLD,
            "thresholds_frozen": False,
            "natural_tail_status": "not admitted by this artifact",
            "source_lineage": "single ChEMBL-37 dual-cold registry; source-family identity is not inferred from document IDs",
        },
        "global": {
            "endpoints": list(ENDPOINTS),
            "by_endpoint": by_endpoint,
            "top_document_tokens": [
                {"document": document, "rows": count}
                for document, count in document_counts.most_common(10)
            ],
        },
    }


def run(*, registry: Path, manifest: Path, output: Path, split: str) -> dict[str, Any]:
    frame = pd.read_parquet(
        registry,
        columns=list(METADATA_COLUMNS),
        filters=[("dual_cold_split", "=", split)],
    )
    report = build_audit(
        frame,
        split=split,
        registry_sha256=file_sha256(registry),
        manifest_sha256=file_sha256(manifest),
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=REGISTRY)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--split", choices=("train", "development"), default="train")
    args = parser.parse_args(argv)
    report = run(
        registry=args.registry,
        manifest=args.manifest,
        output=args.output,
        split=args.split,
    )
    endpoint_summary = report["global"]["by_endpoint"]
    print(
        json.dumps(
            {
                "output": str(args.output),
                "split": args.split,
                "rows": report["input"]["rows"],
                "targets": {endpoint: endpoint_summary[endpoint]["target_count"] for endpoint in ENDPOINTS},
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
