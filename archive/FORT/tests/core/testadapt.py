"""GPU-focused contract tests for the active few-shot runtime."""

from __future__ import annotations

import ast
import inspect
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
import torch

from model.adapt import TargetAdapter
from model.gradadapt import (
    adapttarget,
    buildmatched,
    buildadapter,
    predictquery,
    countparams,
)
from model.ligandbase import LigandBaseline
from model.ligand import FingerprintEncoder
from model.likelihood import ObservationHeads
from model.interaction import InteractionEncoder
from model.posterior import JointPosterior
from model.protein import LandmarkAttention
from model.reorder import ReorderingModel, ReorderingPosterior, calibrationgeometry
from scripts.audit import (
    AuditData,
    FIELDS,
    ROOT,
    ROSTERFIELDS,
    UnionFind,
    canonical,
    chemicalcomponents,
    connectionhash,
    contenthash,
    filehash,
    framehash,
    loadprotocol,
    makeroster,
    proteincomponents,
    proteinreport,
    protocolspec,
    queryindices,
    validateproteinreport,
    writeproteinmap,
    writerooster,
)
from scripts.contract import AffinityRow, Episode
from scripts.episode import (
    buildregistry,
    buildwrong,
    selectsupport,
)
from scripts.guard import assertauthorized
from scripts.metric import evaluateprotocol, pairedcomponents
from scripts.preprocess import normalizeligands, preparerows, preparevectors
from scripts.train import (
    ROSTER,
    buildreordering,
    contrastloss,
    fitbase,
    loadroster,
    maketrainroster,
    verifyroster,
)


pytestmark = pytest.mark.skipif(not torch.cuda.is_available(), reason="core tests require CUDA")


def _rows() -> list[AffinityRow]:
    rows = []
    for target in ("T1", "T2"):
        for index in range(6):
            rows.append(
                AffinityRow(
                    target_key=target,
                    ligand_parent_key=f"L{index}",
                    scaffold_key=f"S{index}",
                    endpoint="pKd",
                    assay_key="A1",
                    document_or_provenance_key=f"P{target}",
                    affinity_value=float(index),
                    split_role="meta_test",
                )
            )
    return rows


def _registry() -> tuple[list[AffinityRow], tuple]:
    rows = _rows()
    vectors = {index: (float(index % 3), float((index + 1) % 3), 1.0) for index in range(len(rows))}
    chemical = {index: f"C{index % 6}" for index in range(len(rows))}
    episodes = buildregistry(
        rows,
        homology_components={"T1": "H1", "T2": "H2"},
        design_vectors=vectors,
        support_size=2,
        chemical_components=chemical,
    )
    return rows, episodes


def testgpuidentity() -> None:
    assert torch.cuda.is_available()
    rows, episodes = _registry()
    assert len(episodes) == 2
    for episode in episodes:
        episode.assertrows(rows)
        assert episode.episode_hash == episode.episode_hash
        assert not set(episode.support_scaffolds).intersection(episode.query_scaffolds)


def testsupportsafe() -> None:
    rows, _ = _registry()
    candidates = tuple(range(6))
    vectors = {index: (float(index), float(index % 2 + 1)) for index in candidates}
    first = selectsupport(rows, candidates, 2, vectors)
    relabeled = [
        AffinityRow(**{**row.__dict__, "affinity_value": 1000.0 - (row.affinity_value or 0.0)})
        for row in rows
    ]
    second = selectsupport(relabeled, candidates, 2, vectors)
    assert first == second


def testuniontransitivity() -> None:
    graph = UnionFind(("A", "B", "C", "D"))
    graph.union("A", "B")
    graph.union("B", "C")
    assert graph.find("A") == graph.find("C")
    assert graph.find("A") != graph.find("D")


def testproteintransitivity() -> None:
    components, stats = proteincomponents(
        {
            "A": "A" * 20,
            "B": "A" * 10 + "G" * 10,
            "C": "G" * 20,
            "D": "W" * 20,
        }
    )
    assert components["A"] == components["B"] == components["C"]
    assert components["A"] != components["D"]
    assert stats["targets"] == 4
    assert stats["pairs"] == 6
    assert stats["edges"] == 2
    assert stats["components"] == 2
    assert len(stats["assignment_sha256"]) == 64


def testproteinmapcertificate(tmp_path: Path) -> None:
    assignments, stats = proteincomponents(
        {"A": "A" * 20, "B": "A" * 10 + "G" * 10, "C": "G" * 20}
    )
    data = AuditData(
        frame=pd.DataFrame(),
        feature=np.empty((0, 1024), dtype=np.float32),
        tokenkeys=frozenset(),
        proteins=assignments,
        proteinstats=stats,
        sources={
            "target_sequences.json": {
                "path": "target_sequences.json",
                "sha256": "a" * 64,
            }
        },
    )
    path = tmp_path / "proteins.v1.json"

    report = writeproteinmap(path, data)

    assert report == json.loads(path.read_text(encoding="utf-8"))
    assert report["stats"]["components"] == 1
    assert len(report["assignments"]) == 3
    damaged = proteinreport(data)
    damaged["assignments"]["C"] = "P:C"
    with pytest.raises(ValueError, match="assignment hash"):
        validateproteinreport(damaged)


def testchemicaltransitivity() -> None:
    frame = pd.DataFrame(
        {
            "source_row": range(6),
            "conn": ("P0", "P1", "P2", "P2", "P4", "P5"),
            "scaffold": ("S0", "S1", "S2", "S3", "S3", "S5"),
        }
    )
    feature = np.zeros((6, 1024), dtype=np.float32)
    feature[0, :90] = 1
    feature[1, :100] = 1
    feature[2, 10:100] = 1
    feature[3, 200:220] = 1
    feature[4, 300:320] = 1
    feature[5, 500:520] = 1

    components, stats = chemicalcomponents(frame, feature)

    assert len({components[index] for index in range(5)}) == 1
    assert components[5] != components[0]
    assert stats["components"] == 2
    assert stats["edges"]["tanimoto"]["qualifying_pairs"] == 2
    assert stats["edges"]["parent"]["qualifying_pairs"] == 1
    assert stats["edges"]["scaffold"]["qualifying_pairs"] == 1
    assert len(stats["assignment_sha256"]) == 64


def teststrictqueryclosure() -> None:
    rows = []
    for index in range(10):
        document = "D0|Dshared" if index == 0 else f"D{index}"
        assay = "A0|Ashared" if index == 0 else f"A{index}"
        if index == 2:
            document = "Dshared|D2"
        if index == 3:
            assay = "Ashared|A3"
        rows.append(
            AffinityRow(
                target_key="T",
                ligand_parent_key=f"L{index}",
                scaffold_key="S0" if index in (0, 4) else f"S{index}",
                endpoint="pKi",
                assay_key=assay,
                document_or_provenance_key=document,
                affinity_value=None,
                split_role="meta_test",
            )
        )
    chemicals = {index: f"C{index}" for index in range(10)}
    chemicals[1] = chemicals[0]

    query = queryindices(tuple(rows), range(10), (0,), chemicals)

    assert query == (5, 6, 7, 8, 9)
    with pytest.raises(ValueError, match="fewer than 6"):
        queryindices(tuple(rows), range(10), (0,), chemicals, minimum=6)


def testcanonicalrosterhash() -> None:
    first = pd.DataFrame(
        [
            {"source_row": 2, "role": "query", "episode": "E"},
            {"source_row": 1, "role": "support", "episode": "E"},
        ]
    )
    reordered = first.iloc[::-1].reset_index(drop=True)
    changed = reordered.copy()
    changed.loc[0, "source_row"] = 3
    assert framehash(first) == framehash(reordered)
    assert framehash(first) != framehash(changed)
    assert framehash(pd.DataFrame(columns=("a",))) != framehash(
        pd.DataFrame(columns=("b",))
    )


def testsidecarmetadata(tmp_path: Path) -> None:
    protocol = protocolspec(1)
    protocolhash = contenthash(protocol)
    base = {
        "target": "T",
        "endpoint": "pKi",
        "homology_component": "P:T",
        "support_source_rows": [1],
        "query_source_rows": [2, 3, 4, 5, 6],
        "support_chemical_components": ["C:1"],
        "query_chemical_components": [f"C:{index}" for index in range(2, 7)],
        "support_ligand_parents": ["L1"],
        "query_ligand_parents": [f"L{index}" for index in range(2, 7)],
        "support_scaffolds": ["S1"],
        "query_scaffolds": [f"S{index}" for index in range(2, 7)],
        "support_documents": ["D1"],
        "query_documents": [f"D{index}" for index in range(2, 7)],
        "support_assays": ["A1"],
        "query_assays": [f"A{index}" for index in range(2, 7)],
        "support_count": 1,
        "query_count": 5,
        "rank": 1,
        "protocol_sha256": protocolhash,
    }
    episode = {"episode": contenthash(base), **base}
    rows = []
    for role, sources in (
        ("support", episode["support_source_rows"]),
        ("query", episode["query_source_rows"]),
    ):
        for source in sources:
            rows.append(
                {
                    "source_row": source,
                    "episode": episode["episode"],
                    "role": role,
                    "rank": 1,
                    "target": "T",
                    "endpoint": "pKi",
                    "homology_component": "P:T",
                    "chemical_component": f"C:{source}",
                }
            )
    roster = pd.DataFrame(rows, columns=ROSTERFIELDS)
    report = {
        "schema": "strict-fewshot-roster.v1",
        "endpoint": "pKi",
        "protocol": protocol,
        "protocol_sha256": protocolhash,
        "labels_read": False,
        "sources": {"registry.parquet": {"path": "registry.parquet", "sha256": "source"}},
        "graphs": {},
        "roster": {
            "columns": ROSTERFIELDS,
            "rows": len(roster),
            "content_sha256": framehash(roster),
        },
        "summary": {},
        "episodes": [episode],
    }
    path = tmp_path / "episodes.pKi.v1.parquet"

    written, sidecar = writerooster(path, roster, report)

    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    metadata = pq.read_schema(path).metadata
    assert payload["protocol"]["minimum_query_rows"] == 5
    assert payload["sources"]["registry.parquet"]["sha256"] == "source"
    assert payload["roster"]["content_sha256"] == framehash(roster)
    assert payload["episodes"][0]["support_source_rows"] == [1]
    assert written["roster"]["file_sha256"]
    assert metadata[b"fort.roster_content_sha256"].decode() == framehash(roster)
    assert json.loads(metadata[b"fort.protocol"])["chemical_components"]["morgan_bits"] == 1024

    report["protocol_sha256"] = "tampered"
    with pytest.raises(ValueError, match="protocol hash"):
        writerooster(path, roster, report)
    report["protocol"] = json.loads(json.dumps(protocol))
    report["protocol"]["support_query_separation"] = ["chemical component"]
    report["protocol_sha256"] = contenthash(report["protocol"])
    with pytest.raises(ValueError, match="must separate"):
        writerooster(path, roster, report)


def testauditdoesnotreadlabels(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    frame = pd.DataFrame(
        {
            **{field: ["T"] for field in FIELDS},
            "affinity": [9.0],
            "replicate_sd": [0.1],
        }
    )
    frame.loc[0, "endpoint"] = "pKi"
    frame.loc[0, "dual_cold_split"] = "development"
    frame.to_parquet(tmp_path / "registry.parquet", index=False)
    (tmp_path / "target_sequences.json").write_text(
        json.dumps({"sequences": {"T": "A" * 20}}), encoding="utf-8"
    )
    np.savez(tmp_path / "target_esm2.npz", keys=np.array(["T"]))
    np.savez(
        tmp_path / "ligand_features.npz",
        feat=np.zeros((1, 1024), dtype=np.float32),
        keep=np.ones(1, dtype=bool),
        conn_sha=np.array(connectionhash(frame)),
    )
    requested = []
    reader = pd.read_parquet

    def capture(*args, **kwargs):
        requested.append(tuple(kwargs.get("columns", ())))
        return reader(*args, **kwargs)

    monkeypatch.setattr(pd, "read_parquet", capture)
    data = loadprotocol(tmp_path)

    assert len(data.frame) == 1
    assert requested == [tuple(FIELDS)]
    assert {"affinity", "replicate_sd", "n_records"}.isdisjoint(requested[0])

    np.savez(
        tmp_path / "ligand_features.npz",
        feat=np.zeros((1, 1024), dtype=np.float32),
        keep=np.ones(1, dtype=bool),
        conn_sha=np.array("stale"),
    )
    with pytest.raises(ValueError, match="not aligned"):
        loadprotocol(tmp_path)
    np.savez(
        tmp_path / "ligand_features.npz",
        feat=np.zeros((1, 1024), dtype=np.float32),
        keep=np.zeros(1, dtype=bool),
        conn_sha=np.array(connectionhash(frame)),
    )
    with pytest.raises(ValueError, match="every registry row"):
        loadprotocol(tmp_path)
    assert all(columns == tuple(FIELDS) for columns in requested)


def testauditselectionandrankgpu() -> None:
    records = []
    for index in range(10):
        records.append(
            {
                "source_row": index,
                "target": "T",
                "conn": f"L{index}",
                "endpoint": "pKi",
                "scaffold": f"S{index}",
                "assays": f"A{index}",
                "docs": f"D{index}",
                "accession": "T",
                "hcluster": "old",
                "dual_cold_split": "development",
            }
        )
    feature = np.zeros((10, 1024), dtype=np.float32)
    for index in range(10):
        feature[index, index * 10 : index * 10 + 5] = 1
    data = AuditData(
        frame=pd.DataFrame(records),
        feature=feature,
        tokenkeys=frozenset({"T"}),
        proteins={"T": "P:T"},
        proteinstats={"targets": 1, "pairs": 0, "edges": 0, "components": 1},
        sources={},
    )

    roster, report = makeroster(data, "pKi", 5)

    assert report["summary"]["episodes"] == 1
    assert (roster.role == "support").sum() == 5
    assert (roster.role == "query").sum() == 5
    assert roster.loc[roster.role == "support", "rank"].min() > 0


def testtrainhomologybridgeexcludedgpu() -> None:
    proteins, stats = proteincomponents(
        {
            "TRAIN": "A" * 20,
            "BRIDGE": "A" * 10 + "G" * 10,
            "DEV": "G" * 20,
        }
    )
    records = [
        {
            "source_row": 0,
            "target": "TRAIN",
            "conn": "LT",
            "endpoint": "pKi",
            "scaffold": "ST",
            "assays": "AT",
            "docs": "DT",
            "accession": "TRAIN",
            "hcluster": "oldtrain",
            "dual_cold_split": "train",
        }
    ]
    for index in range(10):
        records.append(
            {
                "source_row": index + 1,
                "target": "DEV",
                "conn": f"L{index}",
                "endpoint": "pKi",
                "scaffold": f"S{index}",
                "assays": f"A{index}",
                "docs": f"D{index}",
                "accession": "DEV",
                "hcluster": "olddev",
                "dual_cold_split": "development",
            }
        )
    data = AuditData(
        frame=pd.DataFrame(records),
        feature=np.zeros((11, 1024), dtype=np.float32),
        tokenkeys=frozenset({"DEV"}),
        proteins=proteins,
        proteinstats=stats,
        sources={},
    )

    roster, report = makeroster(data, "pKi", 5)

    assert roster.empty
    assert list(roster.columns) == ROSTERFIELDS
    assert report["summary"]["candidate_targets"] == 0
    assert report["summary"]["skipped_targets"]["DEV"] == "train_homology_component"


def testtrainrosterlimitgpu() -> None:
    records = []
    source = 0
    for target, count in (("A", 2), ("B", 6), ("C", 6)):
        for index in range(count):
            records.append(
                {
                    "source_row": source,
                    "target": target,
                    "conn": f"L{target}{index}",
                    "endpoint": "pKi",
                    "scaffold": f"S{target}{index}",
                    "assays": "A1",
                    "docs": f"P{target}",
                    "accession": target,
                    "hcluster": target,
                    "dual_cold_split": "train",
                }
            )
            source += 1
    frame = pd.DataFrame(records)
    features = np.arange(source * 64, dtype=np.float32).reshape(source, 64)
    proteins = {
        target: torch.randn(3, 4, device="cuda") for target in ("A", "B", "C")
    }
    roster = maketrainroster(
        frame, proteins, features, targets=2, querycap=8, support=5
    )
    assert [item.target_key for item in roster] == ["B", "C"]


def testfrozenrosterloadgpu() -> None:
    roster, report, proof = verifyroster(ROSTER, ROOT, "pKi")
    frame = pd.read_parquet(ROOT / "registry.parquet", columns=FIELDS).reset_index(
        names="source_row"
    )
    frame = frame.loc[
        (frame.endpoint == "pKi") & (frame.dual_cold_split == "development")
    ].reset_index(drop=True)
    proteins = {str(item["target"]): torch.empty(1) for item in report["episodes"]}

    episodes, selected = loadroster(
        frame,
        proteins,
        ROSTER,
        targets=2,
        querycap=8,
        root=ROOT,
        endpoint="pKi",
    )
    frozen = {str(item["target"]): item["homology_component"] for item in report["episodes"]}

    assert len(roster) == proof["frozen_row_count"]
    assert len(report["episodes"]) == proof["frozen_episode_count"]
    assert set(proof["support_query_separation"]) == {
        "chemical component",
        "document token",
        "assay token",
    }
    assert all(item.homology_component == frozen[item.target_key] for item in episodes)
    assert all(
        item.homology_component
        != str(frame.loc[frame.target == item.target_key, "hcluster"].iat[0])
        for item in episodes
    )
    assert selected["target_cap"] == 2
    assert selected["query_cap"] == 8
    assert selected["selected_query_rows_after_cap"] == 16


def testfrozenrosterrejectsepisodeeditgpu(tmp_path: Path) -> None:
    path = tmp_path / ROSTER.name
    sidecar = path.with_suffix(".json")
    indexpath = path.with_name("episodes.v1.json")
    shutil.copy2(ROSTER, path)
    report = json.loads(ROSTER.with_suffix(".json").read_text(encoding="utf-8"))
    collection = json.loads(
        ROSTER.with_name("episodes.v1.json").read_text(encoding="utf-8")
    )
    report["episodes"][0]["target"] = "EDITED"

    table = pq.read_table(path)
    metadata = dict(table.schema.metadata or {})
    metadata[b"fort.episodes"] = canonical(report["episodes"]).encode("utf-8")
    pq.write_table(table.replace_schema_metadata(metadata), path)
    report["roster"]["file_sha256"] = filehash(path)
    sidecar.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    indexed = next(item for item in collection["results"] if item["endpoint"] == "pKi")
    indexed["roster"] = report["roster"]
    indexed["sidecar"] = {"path": sidecar.name, "sha256": filehash(sidecar)}
    indexpath.write_text(json.dumps(collection, indent=2, sort_keys=True), encoding="utf-8")

    with pytest.raises(ValueError, match="episode hash is invalid"):
        verifyroster(path, ROOT, "pKi")


def testclosedformbasegpu() -> None:
    feature = np.arange(36, dtype=np.float32).reshape(12, 3) / 10.0
    weight = np.array([0.5, -0.2, 0.3], dtype=np.float32)
    labels = 2.0 + feature @ weight
    train = pd.DataFrame({"source_row": np.arange(8), "affinity": labels[:8]})
    development = pd.DataFrame(
        {"source_row": np.arange(8, 12), "affinity": labels[8:]}
    )
    trainbase, developmentbase, variance = fitbase(
        train, development, feature, ridge=1e-4
    )
    assert np.sqrt(np.mean((trainbase - labels[:8]) ** 2)) < 1e-3
    assert np.sqrt(np.mean((developmentbase - labels[8:]) ** 2)) < 1e-3
    assert variance < 1e-6


@pytest.mark.parametrize("queries", (0, 1))
def testcontrastlossnullgpu(queries: int) -> None:
    prediction = torch.randn(queries, device="cuda", requires_grad=True)
    label = torch.randn(queries, device="cuda")

    loss = contrastloss(prediction, label)
    loss.backward()

    assert torch.equal(loss, torch.zeros((), device="cuda"))
    assert prediction.grad is not None
    assert torch.equal(prediction.grad, torch.zeros_like(prediction))


def testcontrastlosstranslationgpu() -> None:
    prediction = torch.tensor([-1.0, -0.25, 0.5, 1.25], device="cuda")
    label = torch.tensor([0.25, -0.5, 1.0, 1.75], device="cuda")

    original = contrastloss(prediction, label)
    translated = contrastloss(prediction + 8.0, label - 4.0)

    assert torch.equal(original, translated)


def testcontrastlossantisymmetrygpu() -> None:
    prediction = torch.tensor([-0.5, 0.25, 1.5], device="cuda")
    label = torch.tensor([1.0, -0.25, 0.5], device="cuda")
    left, right = torch.triu_indices(3, 3, offset=1, device="cuda")
    forward = prediction[left] - prediction[right]
    reverse = prediction[right] - prediction[left]
    expected = torch.nn.functional.huber_loss(
        forward,
        label[left] - label[right],
    )

    assert torch.equal(reverse, -forward)
    assert torch.equal(contrastloss(prediction, label), expected)
    assert torch.equal(
        contrastloss(prediction.flip(0), label.flip(0)),
        expected,
    )


def testcontrastlossgradientgpu() -> None:
    prediction = torch.tensor(
        [-0.4, 0.2, 0.9, 1.7], device="cuda", requires_grad=True
    )
    label = torch.tensor([1.1, -0.7, 0.3, 2.2], device="cuda")

    loss = contrastloss(prediction, label)
    loss.backward()

    assert torch.isfinite(loss)
    assert loss.item() > 0.0
    assert prediction.grad is not None
    assert torch.isfinite(prediction.grad).all()
    assert torch.count_nonzero(prediction.grad) > 0


def testbackboneinitializationgpu() -> None:
    hybrid = buildreordering(
        backbone="hybrid",
        gateconfig=(False, False),
        seed=17,
    )
    transformer = buildreordering(
        backbone="transformer",
        gateconfig=(False, False),
        seed=17,
    )
    hybridvalues = hybrid.state_dict()
    for name, value in transformer.state_dict().items():
        assert torch.equal(value, hybridvalues[name])
    assert sum(p.numel() for p in transformer.parameters()) == sum(
        p.numel() for p in hybrid.parameters()
    )
    assert hybrid.interaction is None
    assert hybrid.subspace is None
    assert hybrid.globalbasis is not None

    conditioned = buildreordering(
        backbone="hybrid",
        gateconfig=(False, False),
        seed=17,
        proteinconditioned=True,
    )
    assert conditioned.subspace is not None
    assert conditioned.globalbasis is None
    assert conditioned.interaction is not None


def testligandnormalizationgpu() -> None:
    fingerprints = np.arange(6 * 1024, dtype=np.float32).reshape(6, 1024) % 2
    descriptors = np.arange(60, dtype=np.float32).reshape(6, 10)
    features = np.concatenate((fingerprints, descriptors), axis=1)
    normalized, statistics = normalizeligands(features, trainrows=(0, 1, 2, 3))

    assert np.array_equal(normalized[:, :1024], fingerprints)
    assert np.allclose(normalized[:4, -10:].mean(axis=0), 0.0, atol=1e-6)
    assert np.allclose(normalized[:4, -10:].std(axis=0, ddof=1), 1.0, atol=1e-6)

    center = np.asarray(statistics["descriptor_center"], dtype=np.float32)
    scale = np.asarray(statistics["descriptor_scale"], dtype=np.float32)
    assert np.allclose(normalized[4:, -10:], (descriptors[4:] - center) / scale)


def testwrongtarget() -> None:
    rows, episodes = _registry()
    episode = next(item for item in episodes if item.target_key == "T1")
    wrong = buildwrong(episode, rows)
    assert len(wrong) == episode.support_size
    assert all(rows[index].target_key != "T1" for index in wrong)
    assert [rows[index].ligand_parent_key for index in wrong] == [
        rows[index].ligand_parent_key for index in episode.support_indices
    ]


def testzeronullgpu() -> None:
    torch.manual_seed(7)
    model = buildadapter(protein_dim=16, ligand_dim=12, backbone="hybrid").cuda()
    protein = torch.randn(18, 16, device="cuda")
    support_x = torch.randn(3, 12, device="cuda")
    support_y = torch.tensor([1.2, 0.8, 1.6], device="cuda")
    support_b0 = torch.tensor([1.0, 0.9, 1.1], device="cuda")
    query_x = torch.randn(4, 12, device="cuda")
    query_b0 = torch.tensor([0.7, 0.8, 1.0, 1.2], device="cuda")

    state = adapttarget(
        model,
        protein_tokens=protein,
        support_ligand=support_x,
        support_y=support_y,
        support_b0=support_b0,
    )
    output = predictquery(
        model, state, protein_tokens=protein, query_ligand=query_x, query_b0=query_b0
    )
    null_output = predictquery(
        model,
        state.zerocode(),
        protein_tokens=protein,
        query_ligand=query_x,
        query_b0=query_b0,
    )
    assert torch.allclose(null_output["residual"], torch.zeros_like(query_b0), atol=0.0, rtol=0.0)
    assert torch.allclose(null_output["prediction"], null_output["calibration"], atol=0.0, rtol=0.0)

    order = torch.tensor([2, 0, 1], device="cuda")
    permuted = adapttarget(
        model,
        protein_tokens=protein,
        support_ligand=support_x[order],
        support_y=support_y[order],
        support_b0=support_b0[order],
    )
    permuted_output = predictquery(
        model, permuted, protein_tokens=protein, query_ligand=query_x, query_b0=query_b0
    )
    assert torch.allclose(output["prediction"], permuted_output["prediction"], atol=1e-5, rtol=1e-5)


def testmatchedbudget() -> None:
    models = buildmatched(protein_dim=16, ligand_dim=12)
    counts = {name: countparams(model) for name, model in models.items()}
    assert set(counts) == {"transformer", "mamba", "hybrid"}
    assert all(count > 0 for count in counts.values())
    assert not any("padding" in name for model in models.values() for name, _ in model.named_parameters())
    with pytest.raises(PermissionError, match="Gate D0 is DATA_NOT_READY"):
        assertauthorized()


def testsupportencoderoutergradient() -> None:
    model = buildadapter(protein_dim=16, ligand_dim=12, d_model=64, task_dim=4, hybrid_stages=1).cuda()
    protein = torch.randn(10, 16, device="cuda")
    support_x = torch.randn(3, 12, device="cuda")
    support_y = torch.tensor([1.0, 1.5, 2.0], device="cuda")
    support_b0 = torch.tensor([0.8, 1.1, 1.3], device="cuda")
    query_x = torch.randn(4, 12, device="cuda")
    query_b0 = torch.tensor([0.7, 0.9, 1.2, 1.4], device="cuda")
    with torch.no_grad():
        model.readout.fill_(1e-3)
    state = adapttarget(
        model,
        protein_tokens=protein,
        support_ligand=support_x,
        support_y=support_y,
        support_b0=support_b0,
    )
    output = predictquery(model, state, protein_tokens=protein, query_ligand=query_x, query_b0=query_b0)
    output["prediction"].square().mean().backward()
    assert model.support.element[0].weight.grad is not None


def testreorderingfallbackgpu() -> None:
    posterior = ReorderingPosterior(ambientdim=6, rank=2).cuda().eval()
    query = torch.randn(4, 6, device="cuda")
    base = torch.randn(4, device="cuda")
    basis = torch.linalg.qr(torch.randn(6, 2, device="cuda"), mode="reduced").Q
    prior = torch.eye(2, device="cuda")
    empty = posterior(query, None, None, base, basis, prior)
    assert torch.equal(empty["pred"], base)
    assert torch.count_nonzero(empty["appliedranking"]) == 0
    assert torch.all(empty["totalvariance"] > 0)

    support = torch.randn(1, 6, device="cuda")
    residual = torch.tensor([0.5], device="cuda")
    one = posterior(query, support, residual, base, basis, prior)
    assert one["rankprobability"].item() == 0.0
    assert torch.count_nonzero(one["appliedranking"]) == 0


def testobservationvariancecountedoncegpu() -> None:
    posterior = ReorderingPosterior(ambientdim=6, rank=2).cuda().eval()
    query = torch.randn(4, 6, device="cuda")
    base = torch.randn(4, device="cuda")
    basis = torch.linalg.qr(torch.randn(6, 2, device="cuda"), mode="reduced").Q
    output = posterior(
        query,
        None,
        None,
        base,
        basis,
        torch.eye(2, device="cuda"),
        basevariance=0.25,
        observationvariance=0.75,
    )

    assert torch.allclose(output["latentvariance"], torch.full_like(base, 0.25))
    assert torch.allclose(output["observationvariance"], torch.full_like(base, 0.75))
    assert torch.allclose(output["totalvariance"], torch.ones_like(base))


def testgateandbmavariancegpu() -> None:
    posterior = ReorderingPosterior(ambientdim=6, rank=2).cuda().eval()
    probability, weight = posterior.gate(
        torch.zeros((), device="cuda"), 0.05, hardgate=True
    )
    assert 0.0 < probability.item() < 0.5
    assert weight.item() == 0.0

    soft = ReorderingPosterior(
        ambientdim=6,
        rank=2,
        rankinghardgate=False,
        calibrationhardgate=False,
    ).cuda().eval()
    softprobability, softweight = soft.gate(
        torch.zeros((), device="cuda"), 0.05, hardgate=False
    )
    assert torch.equal(softweight, softprobability)

    posterior.train()
    query = torch.randn(4, 6, device="cuda")
    support = torch.randn(5, 6, device="cuda")
    residual = torch.tensor([-1.0, -0.2, 0.1, 0.7, 1.4], device="cuda")
    base = torch.zeros(4, device="cuda")
    basis = torch.linalg.qr(torch.randn(6, 2, device="cuda"), mode="reduced").Q
    output = posterior(query, support, residual, base, basis, torch.eye(2, device="cuda"))
    expected = (
        output["rankweight"] * output["rawrankingvariance"]
        + output["rankweight"]
        * (1.0 - output["rankweight"])
        * output["rawrankingmean"].square()
    )
    assert torch.allclose(output["rankingvariance"], expected)
    assert torch.all(output["totalvariance"] >= 0)


def testhelmertinvariancegpu() -> None:
    posterior = ReorderingPosterior(ambientdim=6, rank=2).cuda().train()
    query = torch.randn(4, 6, device="cuda")
    support = torch.randn(5, 6, device="cuda")
    residual = torch.tensor([-1.0, -0.2, 0.1, 0.7, 1.4], device="cuda")
    base = torch.zeros(4, device="cuda")
    basis = torch.linalg.qr(torch.randn(6, 2, device="cuda"), mode="reduced").Q
    prior = torch.eye(2, device="cuda")
    original = posterior(query, support, residual, base, basis, prior)
    order = torch.tensor([3, 0, 4, 1, 2], device="cuda")
    permuted = posterior(query, support[order], residual[order], base, basis, prior)
    shifted = posterior(query, support, residual + 9.0, base, basis, prior)
    translation = torch.randn(1, 6, device="cuda")
    translated = posterior(
        query + translation,
        support + translation,
        residual,
        base,
        basis,
        prior,
    )
    anchored = posterior(
        support.mean(dim=0, keepdim=True),
        support,
        residual,
        torch.zeros(1, device="cuda"),
        basis,
        prior,
    )
    assert torch.allclose(original["rawrankingmean"], permuted["rawrankingmean"], atol=1e-5)
    assert torch.allclose(original["rawrankingmean"], shifted["rawrankingmean"], atol=1e-5)
    assert torch.allclose(original["rawrankingmean"], translated["rawrankingmean"], atol=1e-5)
    assert torch.allclose(
        original["rawrankingvariance"], translated["rawrankingvariance"], atol=1e-5
    )
    assert torch.allclose(
        anchored["rawrankingmean"], torch.zeros(1, device="cuda"), atol=1e-7
    )
    assert torch.allclose(
        anchored["rawrankingvariance"], torch.zeros(1, device="cuda"), atol=1e-7
    )
    assert not torch.allclose(original["rawcalibrationmean"], shifted["rawcalibrationmean"])


def testcalibrationnullspacegpu() -> None:
    supportbase = torch.tensor([-1.0, -0.2, 0.3, 0.9, 1.7], device="cuda")
    querybase = torch.tensor([-0.7, 0.1, 1.2], device="cuda")
    calibration, querycalibration, contrast = calibrationgeometry(
        supportbase, querybase
    )
    design = torch.stack((torch.ones_like(supportbase), supportbase), dim=1)

    assert calibration.shape == (5, 2)
    assert querycalibration.shape == (3, 2)
    assert contrast.shape == (3, 5)
    assert torch.allclose(
        contrast @ design, torch.zeros(3, 2, device="cuda"), atol=1e-6
    )

    posterior = ReorderingPosterior(ambientdim=6, rank=2).cuda().eval()
    support = torch.randn(5, 2, device="cuda")
    query = torch.randn(3, 2, device="cuda")
    affine = 0.4 - 1.3 * supportbase
    output = posterior.ranking(
        query,
        support,
        affine,
        torch.eye(2, device="cuda"),
        supportbase,
        querybase,
        mode="joint",
    )
    assert torch.allclose(output["mean"], torch.zeros(3, device="cuda"), atol=1e-5)


def testposteriorcovariancegpu() -> None:
    posterior = ReorderingPosterior(ambientdim=6, rank=2).cuda()
    query = torch.randn(4, 2, device="cuda")
    support = torch.randn(5, 2, device="cuda")
    residual = torch.randn(5, device="cuda")
    result = posterior.ranking(query, support, residual, torch.eye(2, device="cuda"))
    identity = torch.eye(2, device="cuda")
    covariance = torch.cholesky_solve(identity, result["precisionchol"])
    assert torch.linalg.eigvalsh(covariance).min() >= 0
    assert torch.all(result["variance"] >= 0)


def testbayesianlabelboundarygpu() -> None:
    model = ReorderingModel(
        proteindim=16,
        liganddim=12,
        dmodel=64,
        stages=1,
        backbone="transformer",
    ).cuda()
    protein = torch.randn(10, 16, device="cuda")
    ligand = torch.randn(5, 12, device="cuda")
    first = model.adapt(
        proteintokens=protein,
        supportligand=ligand,
    )
    second = model.adapt(
        proteintokens=protein,
        supportligand=ligand,
    )
    assert torch.equal(first.proteinfeature, second.proteinfeature)
    assert torch.equal(first.supportfeature, second.supportfeature)
    assert torch.equal(first.basis, second.basis)
    assert "supportlabel" not in inspect.signature(model.adapt).parameters
    assert "supportbase" not in inspect.signature(model.adapt).parameters

    other = model.adapt(
        proteintokens=torch.randn_like(protein),
        supportligand=ligand,
    )
    assert not torch.allclose(first.basis, other.basis)
    queryligand = torch.randn(4, 12, device="cuda")
    firstsupport, firstquery = model.rankfeatures(first, queryligand)
    othersupport, otherquery = model.rankfeatures(other, queryligand)
    assert not torch.allclose(firstsupport, othersupport)
    assert not torch.allclose(firstquery, otherquery)
    assert torch.allclose(firstsupport.mean(dim=0), torch.zeros(2, device="cuda"), atol=1e-6)
    assert model.interaction.pairfeature[0].in_features == 64
    assert model.interaction.pairfeature[0].bias is None


def testfeaturenormalizationgpu() -> None:
    encoder = InteractionEncoder(
        16, 12, dmodel=64, taskdim=8, stages=1, backbone="transformer", conditioned=False
    ).cuda()
    output = encoder(torch.randn(10, 16, device="cuda"), torch.randn(5, 12, device="cuda"))
    assert torch.allclose(output.mean(dim=1), torch.zeros(5, device="cuda"), atol=1e-5)
    variance = output.var(dim=1, unbiased=False)
    assert torch.all((variance > 0.99) & (variance <= 1.001))


def testlandmarkmaskedmeangpu() -> None:
    attention = LandmarkAttention(dmodel=16, landmarks=3).cuda()
    tokens = torch.ones(1, 5, 16, device="cuda")
    summary = attention.summarize(tokens)
    assert summary.shape == (1, 3, 16)
    assert torch.equal(summary, torch.ones_like(summary))


def testrankbound() -> None:
    with pytest.raises(ValueError):
        ReorderingModel(proteindim=16, liganddim=12, ambientdim=8, rank=5, primaryk=5)


def testlikelihoodgpu() -> None:
    heads = ObservationHeads(endpoints=2, sources=3).cuda()
    score = torch.randn(4, device="cuda", requires_grad=True)
    value = torch.randn(4, device="cuda")
    endpoint = torch.tensor([0, 1, 0, 1], device="cuda")
    source = torch.tensor([0, 1, 2, 0], device="cuda")
    loss = heads.exact_nll(score, value, endpoint, source)
    loss.backward()
    assert torch.isfinite(loss)
    assert score.grad is not None
    latent = heads.latent(score.detach())
    assert torch.equal(latent, score.detach())
    with torch.no_grad():
        heads.exact_source_bias.fill_(3.0)
    observed, _ = heads.exact_parameters(score.detach(), endpoint, source)
    assert not torch.equal(observed, latent)
    assert torch.equal(heads.latent(score.detach()), latent)


def testligandbasegpu() -> None:
    model = LigandBaseline(bits=12, rep=16, pdim=4).cuda()
    query = torch.randn(3, 12, device="cuda")
    output = model(query)
    assert torch.equal(output["pred"], output["base"])
    assert output["effective_k"] == 0
    default = FingerprintEncoder().cuda()
    with pytest.raises(ValueError):
        default(torch.randn(2, 1024, device="cuda"))


def _targetadapter(*, rankinghardgate: bool = False) -> TargetAdapter:
    ligand = LigandBaseline(bits=12, rep=16, pdim=4).cuda()
    increment = ReorderingModel(
        proteindim=4,
        liganddim=12,
        dmodel=32,
        ambientdim=4,
        rank=2,
        primaryk=5,
        stages=1,
        landmarks=2,
        backbone="transformer",
        rankinghardgate=rankinghardgate,
    ).cuda()
    return TargetAdapter(ligand, increment).cuda()


def _targettensors(k: int = 5) -> tuple[torch.Tensor, ...]:
    return (
        torch.randn(3, 4, device="cuda"),
        torch.randn(k, 12, device="cuda"),
        torch.linspace(0.5, 1.5, k, device="cuda"),
        torch.linspace(0.4, 1.2, k, device="cuda"),
        torch.randn(3, 12, device="cuda"),
        torch.tensor([0.7, 0.9, 1.1], device="cuda"),
    )


def testsharedligandfallbackgpu() -> None:
    torch.manual_seed(31)
    adapter = _targetadapter().eval()
    protein, support, label, supportb0, query, queryb0 = _targettensors()
    state = adapter.adapt(
        proteintokens=protein,
        supportligand=support,
        supportlabel=label,
        supportb0=supportb0,
    )
    expected = adapter.ligand(
        query,
        support,
        label,
        query_base=queryb0,
        support_base=supportb0,
    )
    output = adapter.predict(
        state, queryligand=query, queryb0=queryb0, useprotein=False
    )

    assert adapter.ligand is adapter._modules["ligand"]
    assert torch.equal(output["prediction"], expected["pred"])
    assert torch.equal(output["ligandvariance"], expected["epistemic"])
    assert torch.equal(output["observationvariance"], expected["aleatoric"])
    assert torch.equal(output["totalvariance"], expected["total_variance"])


def testjointposteriorgpu() -> None:
    torch.manual_seed(37)
    posterior = JointPosterior().cuda()
    support = torch.randn(5, 4, device="cuda")
    query = torch.randn(3, 4, device="cuda")
    residual = torch.randn(5, device="cuda")
    raw = torch.randn(4, 4, device="cuda")
    precision = raw @ raw.T + torch.eye(4, device="cuda")
    noise = torch.tensor(0.7, device="cuda")
    result = posterior.condition(query, support, residual, precision, noise)
    covariance = torch.linalg.inv(
        precision + support.T @ support / noise + 1e-6 * torch.eye(4, device="cuda")
    )
    weight = covariance @ support.T @ residual / noise
    assert torch.allclose(result["mean"], query @ weight, atol=2e-5)
    assert torch.allclose(
        result["variance"],
        torch.diagonal(query @ covariance @ query.T),
        atol=2e-5,
    )
    prior = torch.linalg.inv(precision)
    marginal = support @ prior @ support.T + noise * torch.eye(5, device="cuda")
    expected = torch.distributions.MultivariateNormal(
        torch.zeros(5, device="cuda"), covariance_matrix=marginal + 1e-6 * torch.eye(5, device="cuda")
    ).log_prob(residual)
    expected = expected + 0.5 * len(residual) * np.log(2.0 * np.pi)
    assert torch.allclose(
        posterior.logevidence(support, residual, prior, noise), expected, atol=2e-5
    )


def testfrozenligandgpu() -> None:
    torch.manual_seed(41)
    adapter = _targetadapter()
    before = {
        name: value.detach().clone() for name, value in adapter.ligand.state_dict().items()
    }
    adapter.train()
    assert not adapter.ligand.training
    assert all(not parameter.requires_grad for parameter in adapter.ligand.parameters())

    protein, support, label, supportb0, query, queryb0 = _targettensors()
    state = adapter.adapt(
        proteintokens=protein,
        supportligand=support,
        supportlabel=label,
        supportb0=supportb0,
    )
    output = adapter.predict(state, queryligand=query, queryb0=queryb0)
    optimizer = torch.optim.Adam(
        (parameter for parameter in adapter.parameters() if parameter.requires_grad),
        lr=1e-3,
    )
    optimizer.zero_grad(set_to_none=True)
    output["prediction"].square().mean().backward()
    optimizer.step()

    assert all(parameter.grad is None for parameter in adapter.ligand.parameters())
    for name, value in adapter.ligand.state_dict().items():
        assert torch.equal(value, before[name])


def testadaptervariancegpu() -> None:
    torch.manual_seed(43)
    adapter = _targetadapter().eval()
    protein, support, label, supportb0, query, queryb0 = _targettensors()
    state = adapter.adapt(
        proteintokens=protein,
        supportligand=support,
        supportlabel=label + 5.0,
        supportb0=supportb0,
    )
    base = adapter.ligand(
        query,
        support,
        label + 5.0,
        query_base=queryb0,
        support_base=supportb0,
    )
    joint = adapter.predict(state, queryligand=query, queryb0=queryb0, mode="joint")
    output = adapter.predict(state, queryligand=query, queryb0=queryb0, mode="soft")

    probability = output["proteinprobability"]
    expected = (
        output["calibrationvariance"]
        + probability * (joint["jointvariance"] - output["calibrationvariance"])
        + probability
        * (1.0 - probability)
        * (joint["jointprediction"] - output["calibrationprediction"]).square()
    )
    assert torch.equal(output["observationvariance"], base["aleatoric"])
    assert torch.allclose(output["latentvariance"], expected)
    assert torch.allclose(
        output["totalvariance"],
        output["latentvariance"] + base["aleatoric"],
    )


@pytest.mark.parametrize("k", (0, 1))
def testadapternullgpu(k: int) -> None:
    torch.manual_seed(47 + k)
    adapter = _targetadapter().eval()
    protein, support, label, supportb0, query, queryb0 = _targettensors(k)
    state = adapter.adapt(
        proteintokens=protein,
        supportligand=support,
        supportlabel=label,
        supportb0=supportb0,
    )
    expected = adapter.ligand(
        query,
        support,
        label,
        query_base=queryb0,
        support_base=supportb0,
    )
    output = adapter.predict(state, queryligand=query, queryb0=queryb0)

    assert torch.equal(output["prediction"], expected["pred"])
    assert torch.equal(output["totalvariance"], expected["total_variance"])
    assert torch.count_nonzero(output["appliedprotein"]) == 0
    assert torch.equal(output["jointvariance"], expected["epistemic"])


def testimportgraph() -> None:
    prohibited = {"ubse", "dcst", "pcic", "rdib", "pdmvr", "omut", "medip", "tau", "a1"}
    active_files = [
        *Path("model").glob("*.py"),
        *Path("scripts").glob("*.py"),
    ]
    for path in active_files:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name.lower() for alias in node.names)
            if isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module.lower())
        assert not any(part in name for name in imports for part in prohibited)


def testpreprocessgpu() -> None:
    records = [
        {
            "target_key": "T1",
            "ligand_parent_key": "L1",
            "scaffold_key": "S1",
            "endpoint": "pKd",
            "assay_key": "A1",
            "document_or_provenance_key": "P1",
            "split_role": "meta_test",
        }
    ]
    rows = preparerows(records)
    assert rows[0].affinity_value is None
    vectors = preparevectors([[1.0, 0.0], [0.0, 1.0]])
    assert vectors.is_cuda
    with pytest.raises(PermissionError):
        preparerows([{**records[0], "affinity_value": 1.0}])


def testmetricgpu() -> None:
    _, episodes = _registry()
    episode = episodes[0]
    indices = list(episode.query_indices)
    result = evaluateprotocol(
        predictions=[0.1 + index for index in range(len(indices))],
        labels=[0.2 + index for index in range(len(indices))],
        episodes=[episode],
        prediction_indices=indices,
        component_by_target={episode.target_key: episode.homology_component},
    )
    assert result["independent_components"] == 1.0
    tied = evaluateprotocol(
        predictions=[0.0 for _ in indices],
        labels=[float(index) for index in range(len(indices))],
        episodes=[episode],
        prediction_indices=indices,
        component_by_target={episode.target_key: episode.homology_component},
    )
    assert tied["pairwise_accuracy"] == 0.5
    assert tied["concordance_index"] == 0.5


def testcomponentbootstrapgpu() -> None:
    _, episodes = _registry()
    indices = [index for episode in episodes for index in episode.query_indices]
    labels = [float(position) for position in range(len(indices))]
    result = pairedcomponents(
        predictions={
            "bayesian": labels,
            "control": [value + 1.0 for value in labels],
            "constant": [0.0 for _ in labels],
        },
        labels=labels,
        episodes=episodes,
        prediction_indices=indices,
        component_by_target={item.target_key: item.homology_component for item in episodes},
        reference="bayesian",
        replicates=100,
    )
    assert result["control"]["rmse_gain"]["mean"] > 0
    assert result["control"]["mae_gain"]["probability_positive"] == 1.0
    assert result["constant"]["spearman_gain"]["components"] == 0.0
    assert np.isnan(result["constant"]["spearman_gain"]["mean"])
