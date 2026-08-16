"""Freeze strict, label-blind few-shot target-adaptation episodes."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from math import floor
from pathlib import Path
from typing import Any, Hashable, Iterable, Mapping
from uuid import uuid4

import numpy as np
import pandas as pd
import parasail
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import DataStructs
import torch

from .episode import selectsupport
from .preprocess import preparetable


ROOT = Path("dataset/public/chembl_37/processed/dualcold")
OUT = Path("dataset/processed")
FIELDS = [
    "target",
    "conn",
    "endpoint",
    "scaffold",
    "assays",
    "docs",
    "accession",
    "hcluster",
    "dual_cold_split",
]
ROSTERFIELDS = [
    "source_row",
    "episode",
    "role",
    "rank",
    "target",
    "endpoint",
    "homology_component",
    "chemical_component",
]
SOURCEFILES = (
    "registry.parquet",
    "target_sequences.json",
    "ligand_features.npz",
    "target_esm2.npz",
)
IDENTITY = 0.40
COVERAGE = 0.50
TANIMOTO = 0.90
MORGANBITS = 1024
MINQUERY = 5
GAPOPEN = 10
GAPEXTEND = 1
PROTEINMAP = "proteins.v1.json"


class UnionFind:
    """Small deterministic disjoint-set implementation used by both graphs."""

    def __init__(self, items: Iterable[Hashable]) -> None:
        self.parent = {item: item for item in items}
        self.size = {item: 1 for item in self.parent}

    def find(self, item: Hashable) -> Hashable:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            parent = self.parent[item]
            self.parent[item] = root
            item = parent
        return root

    def union(self, left: Hashable, right: Hashable) -> bool:
        leftroot = self.find(left)
        rightroot = self.find(right)
        if leftroot == rightroot:
            return False
        if self.size[leftroot] < self.size[rightroot]:
            leftroot, rightroot = rightroot, leftroot
        self.parent[rightroot] = leftroot
        self.size[leftroot] += self.size[rightroot]
        return True

    def groups(self) -> tuple[tuple[Hashable, ...], ...]:
        grouped: dict[Hashable, list[Hashable]] = {}
        for item in self.parent:
            grouped.setdefault(self.find(item), []).append(item)
        return tuple(tuple(items) for items in grouped.values())


@dataclass(frozen=True)
class AuditData:
    frame: pd.DataFrame
    feature: np.ndarray
    tokenkeys: frozenset[str]
    proteins: Mapping[str, str]
    proteinstats: Mapping[str, Any]
    sources: Mapping[str, Mapping[str, str]]


def canonical(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def contenthash(value: Any) -> str:
    return sha256(canonical(value).encode("utf-8")).hexdigest()


def filehash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def plain(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if pd.isna(value):
        return None
    return value


def framehash(frame: pd.DataFrame) -> str:
    columns = sorted(frame.columns)
    records: list[dict[str, Any]] = []
    for values in frame.loc[:, columns].itertuples(index=False, name=None):
        record = {column: plain(value) for column, value in zip(columns, values)}
        records.append(record)
    records.sort(key=canonical)
    return contenthash({"columns": columns, "records": records})


def tokens(value: object) -> frozenset[str]:
    return frozenset(part.strip() for part in str(value).split("|") if part.strip())


def cleansequence(value: object) -> str:
    alphabet = frozenset("ARNDCQEGHILKMFPSTWYV")
    sequence = "".join(str(value).split()).upper()
    result = "".join(letter if letter in alphabet else "X" for letter in sequence)
    if not result:
        raise ValueError("target sequences must be nonempty")
    return result


def tracecounts(left: str, right: str) -> tuple[int, int, int, int]:
    result = parasail.sw_trace_striped_32(
        left,
        right,
        GAPOPEN,
        GAPEXTEND,
        parasail.blosum62,
    )
    query = result.traceback.query
    reference = result.traceback.ref
    columns = len(query)
    if columns == 0:
        return 0, 0, 0, min(len(left), len(right))
    matches = sum(
        queryletter == referenceletter and queryletter != "-"
        for queryletter, referenceletter in zip(query, reference)
    )
    shorter = query if len(left) <= len(right) else reference
    covered = sum(letter != "-" for letter in shorter)
    return matches, columns, covered, min(len(left), len(right))


def alignmentcounts(left: str, right: str) -> tuple[int, int, int, int]:
    return tracecounts(cleansequence(left), cleansequence(right))


def alignmentmetrics(left: str, right: str) -> tuple[float, float]:
    matches, columns, covered, shorter = alignmentcounts(left, right)
    identity = matches / columns if columns else 0.0
    coverage = covered / shorter if shorter else 0.0
    return identity, coverage


def proteincomponents(
    sequences: Mapping[str, object],
) -> tuple[dict[str, str], dict[str, Any]]:
    keys = tuple(sorted(str(key) for key in sequences))
    normalized = {key: cleansequence(sequences[key]) for key in keys}
    graph = UnionFind(keys)
    edges = 0
    pairs = 0
    for position, left in enumerate(keys):
        for right in keys[position + 1 :]:
            pairs += 1
            matches, columns, covered, shorter = tracecounts(
                normalized[left], normalized[right]
            )
            identitypasses = columns > 0 and matches * 5 >= columns * 2
            coveragepasses = shorter > 0 and covered * 2 >= shorter
            if identitypasses and coveragepasses:
                graph.union(left, right)
                edges += 1
    labels: dict[str, str] = {}
    groups = graph.groups()
    for members in groups:
        label = f"P:{min(str(member) for member in members)}"
        labels.update({str(member): label for member in members})
    stats = {
        "targets": len(keys),
        "pairs": pairs,
        "edges": edges,
        "components": len(groups),
        "assignment_sha256": contenthash(sorted(labels.items())),
    }
    return labels, stats


def bitvector(packed: np.ndarray) -> Any:
    return DataStructs.CreateFromBinaryText(packed.tobytes())


def chemicalcomponents(
    frame: pd.DataFrame,
    feature: np.ndarray,
) -> tuple[dict[int, str], dict[str, Any]]:
    if frame.empty:
        edges = {
            name: {"qualifying_pairs": 0, "unions": 0}
            for name in ("parent", "scaffold", "identical_fingerprint", "tanimoto")
        }
        return {}, {
            "rows": 0,
            "fingerprints": 0,
            "pairs": 0,
            "edges": edges,
            "components": 0,
            "assignment_sha256": contenthash([]),
        }
    source = frame.source_row.to_numpy(dtype=np.int64)
    if feature.ndim != 2 or feature.shape[1] < MORGANBITS:
        raise ValueError("ligand features must contain Morgan-1024 bits")
    if source.min() < 0 or source.max() >= len(feature):
        raise ValueError("registry source rows do not index ligand features")
    rawbits = feature[source, :MORGANBITS]
    if not np.isfinite(rawbits).all() or not np.isin(rawbits, (0.0, 1.0)).all():
        raise ValueError("Morgan-1024 features must be finite binary values")
    bits = np.ascontiguousarray(rawbits, dtype=np.uint8)
    packed = np.packbits(bits, axis=1, bitorder="little")
    unique, first, inverse = np.unique(
        packed,
        axis=0,
        return_index=True,
        return_inverse=True,
    )
    graph = UnionFind(range(len(frame)))
    edgestats = {
        name: {"qualifying_pairs": 0, "unions": 0}
        for name in ("parent", "scaffold", "identical_fingerprint", "tanimoto")
    }
    for field in ("conn", "scaffold"):
        name = "parent" if field == "conn" else "scaffold"
        for indices in frame.groupby(field, sort=False).indices.values():
            indices = tuple(int(index) for index in indices)
            edgestats[name]["qualifying_pairs"] += len(indices) * (len(indices) - 1) // 2
            for index in indices[1:]:
                edgestats[name]["unions"] += int(graph.union(indices[0], index))
    fingerprintsizes = np.bincount(inverse, minlength=len(unique))
    edgestats["identical_fingerprint"]["qualifying_pairs"] = int(
        sum(size * (size - 1) // 2 for size in fingerprintsizes)
    )
    for index, fingerprint in enumerate(inverse):
        edgestats["identical_fingerprint"]["unions"] += int(
            graph.union(index, int(first[fingerprint]))
        )

    fingerprints = [bitvector(value) for value in unique]
    counts = np.fromiter(
        (fingerprint.GetNumOnBits() for fingerprint in fingerprints),
        dtype=np.int32,
        count=len(fingerprints),
    )
    order = np.argsort(counts, kind="stable")
    orderedcounts = counts[order]
    pairs = 0
    for position, uniqueindex in enumerate(order):
        count = int(counts[uniqueindex])
        if count == 0:
            continue
        maximum = floor(count / TANIMOTO + 1e-12)
        stop = int(np.searchsorted(orderedcounts, maximum, side="right"))
        candidates = order[position + 1 : stop]
        if not len(candidates):
            continue
        pairs += len(candidates)
        similarities = DataStructs.BulkTanimotoSimilarity(
            fingerprints[uniqueindex],
            [fingerprints[index] for index in candidates],
        )
        for candidate, similarity in zip(candidates, similarities):
            if similarity + 1e-12 >= TANIMOTO:
                edgestats["tanimoto"]["qualifying_pairs"] += 1
                edgestats["tanimoto"]["unions"] += int(
                    graph.union(int(first[uniqueindex]), int(first[candidate]))
                )

    labels: dict[int, str] = {}
    groups = graph.groups()
    for members in groups:
        rows = [int(source[int(member)]) for member in members]
        label = f"C:{min(rows)}"
        labels.update({row: label for row in rows})
    stats = {
        "rows": len(frame),
        "fingerprints": len(unique),
        "pairs": int(pairs),
        "edges": edgestats,
        "components": len(groups),
        "assignment_sha256": contenthash(sorted(labels.items())),
    }
    return labels, stats


def queryindices(
    rows: tuple,
    candidates: Iterable[int],
    support: tuple[int, ...],
    chemicals: Mapping[int, str],
    minimum: int = MINQUERY,
) -> tuple[int, ...]:
    supportset = set(support)
    supportchemicals = {chemicals[index] for index in support}
    supportscaffolds = {rows[index].scaffold_key for index in support}
    supportdocs = set().union(
        *(tokens(rows[index].document_or_provenance_key) for index in support)
    )
    supportassays = set().union(*(tokens(rows[index].assay_key) for index in support))
    query = tuple(
        index
        for index in sorted(candidates)
        if index not in supportset
        and chemicals[index] not in supportchemicals
        and rows[index].scaffold_key not in supportscaffolds
        and tokens(rows[index].document_or_provenance_key).isdisjoint(supportdocs)
        and tokens(rows[index].assay_key).isdisjoint(supportassays)
    )
    if len(query) < minimum:
        raise ValueError(f"strict closure leaves fewer than {minimum} query rows")
    return query


def protocolspec(support: int) -> dict[str, Any]:
    return {
        "schema": "strict-fewshot-roster.v1",
        "support_rows": support,
        "minimum_query_rows": MINQUERY,
        "protein_components": {
            "scope": "all local target sequences before split filtering",
            "algorithm": "Parasail Smith-Waterman traceback",
            "device": "cpu",
            "matrix": "BLOSUM62",
            "gap_open": GAPOPEN,
            "gap_extend": GAPEXTEND,
            "unknown_residue": "X",
            "identity": "exact nongap matches / alignment columns",
            "identity_threshold": IDENTITY,
            "shorter_coverage": "aligned residues from shorter input / shorter input length",
            "shorter_coverage_threshold": COVERAGE,
            "thresholds_inclusive": True,
            "closure": "global transitive union",
            "train_component_exclusion": True,
        },
        "chemical_components": {
            "scope": "protein-clean development rows within one endpoint",
            "edges": ["same parent", "same scaffold", "Morgan-1024 Tanimoto"],
            "morgan_radius": 2,
            "morgan_bits": MORGANBITS,
            "tanimoto_threshold": TANIMOTO,
            "threshold_inclusive": True,
            "closure": "global transitive union",
            "screening": "exact bit-count upper bound; no approximate neighbours",
            "device": "cpu",
        },
        "support_query_separation": ["chemical component", "document token", "assay token"],
        "support_selection": {
            "method": "label-blind greedy query-span",
            "design_columns": [0, 64],
            "chemical_component_cap": 1,
            "device": "cuda",
        },
        "rank_check_device": "cuda",
        "affinity_labels_read": False,
        "confirmation_or_sealed_affinity_labels_read": False,
    }


def proteinreport(data: AuditData) -> dict[str, Any]:
    """Bind every sequence target to the verified global homology graph."""

    assignments = {target: data.proteins[target] for target in sorted(data.proteins)}
    assignmenthash = contenthash(list(assignments.items()))
    stats = dict(data.proteinstats)
    registeredhash = stats.get("assignment_sha256")
    if registeredhash is not None and registeredhash != assignmenthash:
        raise ValueError("protein component statistics disagree with assignments")
    stats["assignment_sha256"] = assignmenthash
    protocol = protocolspec(5)["protein_components"]
    return {
        "schema": "protein-components.v1",
        "protocol": protocol,
        "protocol_sha256": contenthash(protocol),
        "sources": data.sources,
        "stats": stats,
        "assignments": assignments,
    }


def validateproteinreport(report: Mapping[str, Any]) -> None:
    """Fail closed when a serialized protein-component map is incomplete."""

    if report.get("schema") != "protein-components.v1":
        raise ValueError("protein component schema is invalid")
    protocol = report.get("protocol")
    if not isinstance(protocol, dict) or report.get("protocol_sha256") != contenthash(protocol):
        raise ValueError("protein component protocol hash is invalid")
    assignments = report.get("assignments")
    stats = report.get("stats")
    if not isinstance(assignments, dict) or not assignments:
        raise ValueError("protein component assignments are missing")
    if not isinstance(stats, dict):
        raise ValueError("protein component statistics are missing")
    normalized = {str(key): str(value) for key, value in assignments.items()}
    if assignments != normalized:
        raise ValueError("protein component assignments must be string mappings")
    assignmenthash = contenthash(list(sorted(normalized.items())))
    if stats.get("assignment_sha256") != assignmenthash:
        raise ValueError("protein component assignment hash is invalid")
    if int(stats.get("targets", -1)) != len(normalized):
        raise ValueError("protein component target count is invalid")
    if int(stats.get("components", -1)) != len(set(normalized.values())):
        raise ValueError("protein component count is invalid")


def writeproteinmap(path: Path, data: AuditData) -> dict[str, Any]:
    """Atomically persist and reread the complete component certificate."""

    report = proteinreport(data)
    validateproteinreport(report)
    atomicjson(path, report)
    reloaded = json.loads(path.read_text(encoding="utf-8"))
    validateproteinreport(reloaded)
    if reloaded != report:
        raise RuntimeError("protein component JSON round-trip changed content")
    return report


def sourcehashes(root: Path) -> dict[str, dict[str, str]]:
    return {
        name: {"path": name, "sha256": filehash(root / name)}
        for name in SOURCEFILES
    }


def connectionhash(frame: pd.DataFrame) -> str:
    values = pd.util.hash_pandas_object(frame["conn"], index=False).values
    return sha256(values.tobytes()).hexdigest()


def loadprotocol(root: Path) -> AuditData:
    frame = pd.read_parquet(root / "registry.parquet", columns=FIELDS).reset_index(
        names="source_row"
    )
    sequencepayload = json.loads((root / "target_sequences.json").read_text(encoding="utf-8"))
    sequences = sequencepayload.get("sequences")
    if not isinstance(sequences, dict):
        raise ValueError("target_sequences.json must contain a sequences mapping")
    sequencekeys = {str(key) for key in sequences}
    required = set(
        frame.loc[
            frame.dual_cold_split.isin(("train", "development")), "target"
        ].astype(str)
    )
    missing = sorted(required.difference(sequencekeys))
    if missing:
        raise ValueError(f"train/development targets lack sequences: {missing[:5]}")
    with np.load(root / "target_esm2.npz", allow_pickle=False) as archive:
        tokenkeys = frozenset(str(key) for key in archive["keys"])
    with np.load(root / "ligand_features.npz", allow_pickle=False) as archive:
        needed = {"feat", "keep", "conn_sha"}
        if not needed.issubset(archive.files):
            raise ValueError("ligand feature cache lacks feat/keep/conn_sha alignment fields")
        feature = archive["feat"]
        keep = archive["keep"]
        cachedhash = str(archive["conn_sha"].item())
    if len(feature) != len(frame):
        raise ValueError("ligand features and registry must have identical row counts")
    if keep.dtype != np.bool_ or keep.shape != (len(frame),) or not bool(keep.all()):
        raise ValueError("every registry row must have a valid cached ligand feature")
    if cachedhash != connectionhash(frame):
        raise ValueError("ligand feature cache is not aligned to registry connectivity rows")
    proteins, proteinstats = proteincomponents(sequences)
    return AuditData(
        frame=frame,
        feature=feature,
        tokenkeys=tokenkeys,
        proteins=proteins,
        proteinstats=proteinstats,
        sources=sourcehashes(root),
    )


def makerecords(frame: pd.DataFrame) -> tuple:
    return preparetable(frame[FIELDS].to_dict("records"), split="meta_test")


def episodepayload(
    frame: pd.DataFrame,
    rows: tuple,
    target: str,
    endpoint: str,
    homology: str,
    support: tuple[int, ...],
    query: tuple[int, ...],
    chemicals: Mapping[int, str],
    rank: int,
    protocolhash: str,
) -> dict[str, Any]:
    def source(indices: tuple[int, ...]) -> list[int]:
        return sorted(int(frame.source_row.iat[index]) for index in indices)

    def values(indices: tuple[int, ...], field: str) -> list[str]:
        return sorted({str(getattr(rows[index], field)) for index in indices})

    def metatokens(indices: tuple[int, ...], field: str) -> list[str]:
        return sorted(
            set().union(*(tokens(getattr(rows[index], field)) for index in indices))
        )

    payload = {
        "target": target,
        "endpoint": endpoint,
        "homology_component": homology,
        "support_source_rows": source(support),
        "query_source_rows": source(query),
        "support_chemical_components": sorted({chemicals[index] for index in support}),
        "query_chemical_components": sorted({chemicals[index] for index in query}),
        "support_ligand_parents": values(support, "ligand_parent_key"),
        "query_ligand_parents": values(query, "ligand_parent_key"),
        "support_scaffolds": values(support, "scaffold_key"),
        "query_scaffolds": values(query, "scaffold_key"),
        "support_documents": metatokens(support, "document_or_provenance_key"),
        "query_documents": metatokens(query, "document_or_provenance_key"),
        "support_assays": metatokens(support, "assay_key"),
        "query_assays": metatokens(query, "assay_key"),
        "support_count": len(support),
        "query_count": len(query),
        "rank": rank,
        "protocol_sha256": protocolhash,
    }
    return {"episode": contenthash(payload), **payload}


def makeroster(
    data: AuditData,
    endpoint: str,
    support: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if not torch.cuda.is_available():
        raise RuntimeError("the few-shot audit requires CUDA")
    if support <= 0:
        raise ValueError("support must be positive")
    frame = data.frame
    traincomponents = {
        data.proteins[target]
        for target in frame.loc[frame.dual_cold_split == "train", "target"].astype(str)
        if target in data.proteins
    }
    development = frame.loc[
        (frame.dual_cold_split == "development") & (frame.endpoint == endpoint)
    ].copy()
    skipped: dict[str, str] = {}
    for target in sorted(development.target.astype(str).unique()):
        if target not in data.proteins:
            skipped[target] = "missing_sequence"
        elif target not in data.tokenkeys:
            skipped[target] = "missing_embedding"
        elif data.proteins[target] in traincomponents:
            skipped[target] = "train_homology_component"
    development = development.loc[
        development.target.astype(str).map(
            lambda target: target in data.proteins
            and target in data.tokenkeys
            and data.proteins[target] not in traincomponents
        )
    ].reset_index(drop=True)
    chemicalbysource, chemicalstats = chemicalcomponents(development, data.feature)
    rows = makerecords(development)
    vectors = data.feature[
        development.source_row.to_numpy(dtype=np.int64), :64
    ].astype(np.float64, copy=False)
    chemicals = {
        index: chemicalbysource[int(development.source_row.iat[index])]
        for index in range(len(development))
    }
    protocol = protocolspec(support)
    protocolhash = contenthash(protocol)
    output: list[dict[str, Any]] = []
    episodes: list[dict[str, Any]] = []

    for target, group in development.groupby("target", sort=True):
        indices = tuple(int(index) for index in group.index)
        if len(indices) < support + MINQUERY:
            skipped[str(target)] = "row_depth"
            continue
        if group.scaffold.nunique() < support + 1:
            skipped[str(target)] = "scaffold_depth"
            continue
        design = {index: vectors[index] for index in indices}
        try:
            chosen = selectsupport(
                rows,
                indices,
                support,
                design,
                chemical_components=chemicals,
                chemical_component_cap=1,
            )
            query = queryindices(rows, indices, chosen, chemicals)
        except ValueError as error:
            skipped[str(target)] = str(error)
            continue
        rank = int(
            torch.linalg.matrix_rank(
                torch.as_tensor(vectors[list(chosen)], dtype=torch.float64, device="cuda")
            ).item()
        )
        homology = data.proteins[str(target)]
        payload = episodepayload(
            development,
            rows,
            str(target),
            endpoint,
            homology,
            chosen,
            query,
            chemicals,
            rank,
            protocolhash,
        )
        episodes.append(payload)
        for role, indiceset in (("support", chosen), ("query", query)):
            for index in indiceset:
                output.append(
                    {
                        "source_row": int(development.source_row.iat[index]),
                        "episode": payload["episode"],
                        "role": role,
                        "rank": rank,
                        "target": str(target),
                        "endpoint": endpoint,
                        "homology_component": homology,
                        "chemical_component": chemicals[index],
                    }
                )

    roster = pd.DataFrame(output, columns=ROSTERFIELDS)
    if not roster.empty:
        roleorder = roster.role.map({"support": 0, "query": 1})
        roster = (
            roster.assign(roleorder=roleorder)
            .sort_values(["episode", "roleorder", "source_row"], kind="stable")
            .drop(columns="roleorder")
            .reset_index(drop=True)
        )
    querycounts = [episode["query_count"] for episode in episodes]
    queryscaffolds = [len(episode["query_scaffolds"]) for episode in episodes]
    ranks = [episode["rank"] for episode in episodes]
    summary = {
        "endpoint": endpoint,
        "support": support,
        "candidate_targets": int(development.target.nunique()),
        "episodes": len(episodes),
        "independent_homology_components": len(
            {episode["homology_component"] for episode in episodes}
        ),
        "median_query_rows": float(np.median(querycounts)) if querycounts else 0.0,
        "median_query_scaffolds": float(np.median(queryscaffolds)) if queryscaffolds else 0.0,
        "minimum_support_rank": min(ranks) if ranks else 0,
        "skipped_targets": skipped,
        "cuda": torch.cuda.get_device_name(torch.cuda.current_device()),
    }
    report = {
        "schema": "strict-fewshot-roster.v1",
        "endpoint": endpoint,
        "protocol": protocol,
        "protocol_sha256": protocolhash,
        "labels_read": False,
        "sources": data.sources,
        "graphs": {
            "protein": dict(data.proteinstats),
            "chemical": chemicalstats,
        },
        "roster": {
            "columns": ROSTERFIELDS,
            "rows": len(roster),
            "content_sha256": framehash(roster),
        },
        "summary": summary,
        "episodes": episodes,
    }
    return roster, report


def buildroster(
    root: Path,
    endpoint: str,
    support: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    return makeroster(loadprotocol(root), endpoint, support)


def episodehash(episode: Mapping[str, Any]) -> str:
    return contenthash({key: value for key, value in episode.items() if key != "episode"})


def validatereport(roster: pd.DataFrame, report: Mapping[str, Any]) -> None:
    if list(roster.columns) != ROSTERFIELDS:
        raise ValueError("roster columns do not match the frozen schema")
    expected = framehash(roster)
    if report["roster"]["content_sha256"] != expected:
        raise ValueError("roster changed after its content hash was computed")
    if report["roster"]["columns"] != ROSTERFIELDS or report["roster"]["rows"] != len(roster):
        raise ValueError("roster schema or row count disagrees with sidecar metadata")
    protocolhash = contenthash(report["protocol"])
    if report["protocol_sha256"] != protocolhash:
        raise ValueError("protocol hash does not match protocol parameters")
    requiredseparation = {"chemical component", "document token", "assay token"}
    if set(report["protocol"]["support_query_separation"]) != requiredseparation:
        raise ValueError("protocol must separate chemical, document, and assay components")
    episodes = report["episodes"]
    identifiers = [episode["episode"] for episode in episodes]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("sidecar contains duplicate episode identifiers")
    if set(roster.episode) != set(identifiers):
        raise ValueError("roster and sidecar contain different episodes")
    minimum = int(report["protocol"]["minimum_query_rows"])
    supportsize = int(report["protocol"]["support_rows"])
    for episode in episodes:
        identifier = episode["episode"]
        if identifier != episodehash(episode):
            raise ValueError(f"episode hash is invalid: {identifier}")
        if episode["protocol_sha256"] != protocolhash:
            raise ValueError(f"episode protocol hash is invalid: {identifier}")
        group = roster.loc[roster.episode == identifier]
        if set(group.role) != {"support", "query"} or group.source_row.duplicated().any():
            raise ValueError(f"episode roles or source rows are invalid: {identifier}")
        supportrows = sorted(
            int(value) for value in group.loc[group.role == "support", "source_row"]
        )
        queryrows = sorted(
            int(value) for value in group.loc[group.role == "query", "source_row"]
        )
        if supportrows != episode["support_source_rows"]:
            raise ValueError(f"episode support rows disagree with roster: {identifier}")
        if queryrows != episode["query_source_rows"]:
            raise ValueError(f"episode query rows disagree with roster: {identifier}")
        if len(supportrows) != supportsize or len(queryrows) < minimum:
            raise ValueError(f"episode role counts violate protocol: {identifier}")
        if len(group) != len(supportrows) + len(queryrows):
            raise ValueError(f"episode contains unregistered role rows: {identifier}")
        if (
            episode["support_count"] != len(supportrows)
            or episode["query_count"] != len(queryrows)
        ):
            raise ValueError(f"episode sidecar counts are invalid: {identifier}")
        for column, key in (
            ("target", "target"),
            ("endpoint", "endpoint"),
            ("homology_component", "homology_component"),
            ("rank", "rank"),
        ):
            if set(group[column]) != {episode[key]}:
                raise ValueError(f"episode {column} disagrees with roster: {identifier}")
        for role, key in (
            ("support", "support_chemical_components"),
            ("query", "query_chemical_components"),
        ):
            values = sorted(set(group.loc[group.role == role, "chemical_component"]))
            if values != episode[key]:
                raise ValueError(f"episode chemical components disagree with roster: {identifier}")
        closures = (
            ("support_chemical_components", "query_chemical_components"),
            ("support_scaffolds", "query_scaffolds"),
            ("support_documents", "query_documents"),
            ("support_assays", "query_assays"),
        )
        for left, right in closures:
            if set(episode[left]).intersection(episode[right]):
                raise ValueError(f"episode support/query closure is invalid: {identifier}")


def atomicjson(path: Path, payload: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def writerooster(
    path: Path,
    roster: pd.DataFrame,
    report: dict[str, Any],
) -> tuple[dict[str, Any], Path]:
    validatereport(roster, report)
    expected = framehash(roster)
    metadata = {
        b"fort.protocol": canonical(report["protocol"]).encode("utf-8"),
        b"fort.sources": canonical(report["sources"]).encode("utf-8"),
        b"fort.roster_content_sha256": expected.encode("ascii"),
        b"fort.episodes": canonical(report["episodes"]).encode("utf-8"),
    }
    table = pa.Table.from_pandas(roster, preserve_index=False)
    table = table.replace_schema_metadata({**(table.schema.metadata or {}), **metadata})
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        pq.write_table(table, temporary)
        reloaded = pd.read_parquet(temporary)
        if framehash(reloaded) != expected:
            raise RuntimeError("Parquet round-trip changed the roster content")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    report["roster"].update({"path": path.name, "file_sha256": filehash(path)})
    sidecar = path.with_suffix(".json")
    atomicjson(sidecar, report)
    return report, sidecar


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--support", type=int, default=5)
    parser.add_argument("--endpoint", choices=("pKd", "pKi", "both"), default="both")
    args = parser.parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError("the few-shot audit requires CUDA")
    args.out.mkdir(parents=True, exist_ok=True)
    data = loadprotocol(args.root)
    proteinpath = args.out / PROTEINMAP
    proteinmap = writeproteinmap(proteinpath, data)
    results: list[dict[str, Any]] = []
    endpoints = ("pKd", "pKi") if args.endpoint == "both" else (args.endpoint,)
    for endpoint in endpoints:
        roster, report = makeroster(data, endpoint, args.support)
        report, sidecar = writerooster(
            args.out / f"episodes.{endpoint}.v1.parquet",
            roster,
            report,
        )
        results.append(
            {
                "endpoint": endpoint,
                "roster": report["roster"],
                "sidecar": {"path": sidecar.name, "sha256": filehash(sidecar)},
                "summary": report["summary"],
            }
        )
    result = {
        "schema": "strict-fewshot-roster-index.v1",
        "protocol": protocolspec(args.support),
        "sources": data.sources,
        "protein_components": {
            "path": proteinpath.name,
            "sha256": filehash(proteinpath),
            "assignment_sha256": proteinmap["stats"]["assignment_sha256"],
        },
        "results": results,
    }
    atomicjson(args.out / "episodes.v1.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
