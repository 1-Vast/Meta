"""BindingDB 202608 label-blind projection and cycle census.

The extractor is the only function allowed to open the monolithic BindingDB
article TSV. It records endpoint availability but never converts or emits an
affinity cell. The census consumes only the resulting metadata projection.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


AFFINITY_COLUMNS = ("Ki (nM)", "Kd (nM)")
METADATA_COLUMNS = (
    "BindingDB Reactant_set_id",
    "Ligand SMILES",
    "Ligand InChI Key",
    "Target Name",
    "Target Source Organism According to Curator or DataSource",
    "Curation/DataSource",
    "Article DOI",
    "PMID",
    "Date of publication",
    "Number of Protein Chains in Target (>1 implies a multichain complex)",
    "BindingDB Target Chain Sequence 1",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def normalize_text(value: str) -> str:
    value = re.sub(r"&(?:#\d+|[A-Za-z]+);", " ", value.lower())
    value = re.sub(r"\d+(?:\.\d+)?", "#", value)
    return " ".join(re.findall(r"[a-z]+|#", value))


def protocol_signature(name: str, description: str) -> str:
    return stable_hash(normalize_text(name) + "\n" + normalize_text(description))


def _zip_text(path: Path):
    archive = zipfile.ZipFile(path)
    names = [name for name in archive.namelist() if not name.endswith("/")]
    if len(names) != 1:
        archive.close()
        raise ValueError(f"expected one file in {path}, found {len(names)}")
    binary = archive.open(names[0])
    text = io.TextIOWrapper(binary, encoding="utf-8-sig", errors="replace", newline="")
    return archive, text


def read_assay_map(path: Path) -> dict[str, list[dict[str, str]]]:
    archive, handle = _zip_text(path)
    result: dict[str, list[dict[str, str]]] = defaultdict(list)
    try:
        for row in csv.DictReader(handle, delimiter="\t"):
            entry = row["ENTRYID"].strip()
            assay_id = row["ASSAYID"].strip()
            name = row["ASSAY_NAME"].strip()
            description = row["DESCRIPTION"].strip()
            result[entry].append(
                {
                    "assay_id": assay_id,
                    "assay_name_norm": normalize_text(name),
                    "protocol_sha256": protocol_signature(name, description),
                }
            )
    finally:
        handle.close()
        archive.close()
    return dict(result)


def _deterministic_gzip_text(path: Path):
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return raw, io.TextIOWrapper(compressed, encoding="utf-8", newline="\n")


def extract_projection(articles_zip: Path, assays_zip: Path, output: Path) -> dict:
    assay_map = read_assay_map(assays_zip)
    archive, handle = _zip_text(articles_zip)
    output.parent.mkdir(parents=True, exist_ok=True)
    raw_output, writer = _deterministic_gzip_text(output)
    rows = 0
    projected = 0
    missing_assay_map = 0
    endpoint_rows = defaultdict(int)
    try:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(METADATA_COLUMNS + AFFINITY_COLUMNS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"missing BindingDB columns: {sorted(missing)}")
        for row in reader:
            rows += 1
            endpoints = [name.split()[0] for name in AFFINITY_COLUMNS if row[name].strip()]
            if not endpoints:
                continue
            entry_id = row["BindingDB Reactant_set_id"].strip()
            assays = assay_map.get(entry_id, [])
            missing_assay_map += int(not assays)
            chain_count_text = row[
                "Number of Protein Chains in Target (>1 implies a multichain complex)"
            ].strip()
            chain_count = int(chain_count_text) if chain_count_text.isdigit() else None
            sequence = re.sub(r"\s+", "", row["BindingDB Target Chain Sequence 1"]).upper()
            doi = row["Article DOI"].strip().lower()
            pmid = row["PMID"].strip()
            document = "doi:" + doi if doi else "pmid:" + pmid if pmid else "entry:" + entry_id
            record = {
                "source_row_id": entry_id,
                "document_id": document,
                "publication_date": row["Date of publication"].strip(),
                "source": row["Curation/DataSource"].strip(),
                "endpoint_available": endpoints,
                "target_name": row["Target Name"].strip(),
                "organism": row[
                    "Target Source Organism According to Curator or DataSource"
                ].strip(),
                "chain_count": chain_count,
                "target_sequence": sequence,
                "target_sequence_sha256": stable_hash(sequence) if sequence else "",
                "ligand_inchikey": row["Ligand InChI Key"].strip(),
                "ligand_smiles": row["Ligand SMILES"].strip(),
                "assays": assays,
            }
            writer.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
            projected += 1
            for endpoint in endpoints:
                endpoint_rows[endpoint] += 1
    finally:
        writer.close()
        raw_output.close()
        handle.close()
        archive.close()
    return {
        "schema": "BindingDB.CQ-R0.MetadataProjection.v1",
        "article_rows_traversed": rows,
        "projected_rows": projected,
        "endpoint_rows": dict(sorted(endpoint_rows.items())),
        "missing_assay_mapping_rows": missing_assay_map,
        "affinity_bytes_traversed_by_trusted_extractor": True,
        "numeric_affinity_values_parsed": 0,
        "numeric_affinity_values_exposed": 0,
        "numeric_affinity_values_used": 0,
        "articles_sha256": sha256_file(articles_zip),
        "assays_sha256": sha256_file(assays_zip),
        "projection_sha256": sha256_file(output),
    }


@dataclass
class PanelGraph:
    targets: set[str]
    ligands: set[str]
    edges: set[tuple[str, str]]
    documents: set[str]

    @classmethod
    def empty(cls) -> "PanelGraph":
        return cls(set(), set(), set(), set())

    def add(self, target: str, ligand: str, document: str) -> None:
        self.targets.add(target)
        self.ligands.add(ligand)
        self.edges.add((target, ligand))
        self.documents.add(document)

    @property
    def cycle_rank(self) -> int:
        if not self.edges:
            return 0
        adjacency: dict[str, set[str]] = defaultdict(set)
        for target, ligand in self.edges:
            t_node = "t:" + target
            l_node = "l:" + ligand
            adjacency[t_node].add(l_node)
            adjacency[l_node].add(t_node)
        remaining = set(adjacency)
        components = 0
        while remaining:
            components += 1
            stack = [remaining.pop()]
            while stack:
                node = stack.pop()
                unseen = adjacency[node] & remaining
                remaining.difference_update(unseen)
                stack.extend(unseen)
        return max(0, len(self.edges) - len(self.targets) - len(self.ligands) + components)


def iter_projection(path: Path) -> Iterable[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def _panel_keys(record: dict, endpoint: str) -> dict[str, list[str]]:
    document = record["document_id"]
    assays = record.get("assays") or []
    protocol_hashes = sorted({item["protocol_sha256"] for item in assays})
    name_hashes = sorted({stable_hash(item["assay_name_norm"]) for item in assays})
    exact = [
        f"{record['source_row_id']}|{item['assay_id']}|{endpoint}"
        for item in assays
    ] or [f"{record['source_row_id']}|unmapped|{endpoint}"]
    return {
        "document_protocol": [
            f"{document}|{endpoint}|{'/'.join(protocol_hashes) or 'unmapped'}"
        ],
        "document_assay_name": [
            f"{document}|{endpoint}|{'/'.join(name_hashes) or 'unmapped'}"
        ],
        "document_endpoint": [f"{document}|{endpoint}"],
        "entry_assay_negative_control": exact,
    }


def _summarize(graphs: dict[str, PanelGraph]) -> dict:
    ranks = [graph.cycle_rank for graph in graphs.values()]
    positive = [(key, graph) for key, graph in graphs.items() if graph.cycle_rank > 0]
    positive_ranks = [graph.cycle_rank for _, graph in positive]
    return {
        "panels": len(graphs),
        "cycle_positive_panels": len(positive),
        "total_cycle_rank": sum(ranks),
        "max_cycle_rank": max(ranks, default=0),
        "median_positive_cycle_rank": (
            sorted(positive_ranks)[len(positive_ranks) // 2] if positive_ranks else 0
        ),
        "cycle_positive_edges": sum(len(graph.edges) for _, graph in positive),
        "cycle_positive_targets": len(
            {target for _, graph in positive for target in graph.targets}
        ),
        "cycle_positive_ligands": len(
            {ligand for _, graph in positive for ligand in graph.ligands}
        ),
    }


def run_census(projection: Path) -> dict:
    panel_graphs: dict[str, dict[str, PanelGraph]] = defaultdict(dict)
    total_rows = 0
    eligible_rows = 0
    mapped_construct_rows = 0
    missing_ligand_rows = 0
    for record in iter_projection(projection):
        total_rows += 1
        single_chain = record.get("chain_count") == 1
        sequence = record.get("target_sequence", "")
        ligand = record.get("ligand_inchikey") or stable_hash(record.get("ligand_smiles", ""))
        if not ligand:
            missing_ligand_rows += 1
            continue
        if single_chain:
            eligible_rows += 1
            mapped_construct_rows += int(bool(sequence))
        if not single_chain or not sequence:
            continue
        target = record["target_sequence_sha256"]
        for endpoint in record["endpoint_available"]:
            for definition, keys in _panel_keys(record, endpoint).items():
                for key in keys:
                    graph = panel_graphs[definition].setdefault(key, PanelGraph.empty())
                    graph.add(target, ligand, record["document_id"])
    summaries = {
        definition: _summarize(graphs)
        for definition, graphs in sorted(panel_graphs.items())
    }
    primary = summaries.get("document_protocol", {})
    traceable = eligible_rows > 0 and mapped_construct_rows / eligible_rows >= 0.95
    development_ready = bool(primary.get("total_cycle_rank", 0) > 0 and traceable)
    return {
        "schema": "BindingDB.CQ-R0.Census.v1",
        "projection_sha256": sha256_file(projection),
        "rows": total_rows,
        "single_chain_endpoint_rows": eligible_rows,
        "construct_mapping_coverage": (
            mapped_construct_rows / eligible_rows if eligible_rows else 0.0
        ),
        "missing_ligand_rows": missing_ligand_rows,
        "panel_definitions": summaries,
        "development_training_ready_preclosure": development_ready,
        "biological_claim_ready": False,
        "claim_blockers": [
            "homology/scaffold/document conflict partition not yet materialized",
            "component-level power simulation not yet executed",
            "numeric interaction-existence Gate not executed",
        ],
        "numeric_affinity_values_parsed": 0,
        "numeric_affinity_values_exposed": 0,
        "numeric_affinity_values_used": 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    extract = sub.add_parser("extract")
    extract.add_argument("--articles", type=Path, required=True)
    extract.add_argument("--assays", type=Path, required=True)
    extract.add_argument("--projection", type=Path, required=True)
    extract.add_argument("--manifest", type=Path, required=True)
    census = sub.add_parser("census")
    census.add_argument("--projection", type=Path, required=True)
    census.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "extract":
        result = extract_projection(args.articles, args.assays, args.projection)
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    else:
        result = run_census(args.projection)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
