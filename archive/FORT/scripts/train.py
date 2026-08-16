"""CUDA-only Bayesian kill test on frozen few-shot target roles."""

from __future__ import annotations

import argparse
import json
import subprocess
import threading
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import torch
from torch.nn import functional as F

from model.adapt import TargetAdapter
from model.gradadapt import adapttarget, buildadapter, predictquery
from model.ligandbase import LigandBaseline
from model.reorder import ReorderingModel

from .audit import (
    FIELDS,
    ROOT,
    ROSTERFIELDS,
    SOURCEFILES,
    canonical,
    contenthash,
    episodehash,
    filehash,
    framehash,
    tokens,
)
from .contract import Episode
from .episode import selectsupport
from .guard import assertauthorized
from .metric import evaluateprotocol, pairedcomponents
from .preprocess import normalizeligands, preparetable


REPORT = Path("reports/active")
ROSTER = Path("dataset/processed/strict/episodes.pKi.v1.parquet")


def telemetry(stop: threading.Event, samples: list[dict[str, float]]) -> None:
    """Sample CUDA utilization, board power, and allocated device memory."""

    while not stop.is_set():
        try:
            text = subprocess.check_output(
                [
                    "nvidia-smi",
                    "--query-gpu=utilization.gpu,power.draw,memory.used",
                    "--format=csv,noheader,nounits",
                ],
                text=True,
            ).strip().splitlines()[0]
            use, power, memory = (float(value.strip()) for value in text.split(","))
            samples.append({"utilization": use, "power": power, "memory": memory})
        except (OSError, subprocess.SubprocessError, ValueError):
            pass
        stop.wait(0.5)


def loadproteins(root: Path) -> dict[str, torch.Tensor]:
    data = np.load(root / "target_esm2.npz", allow_pickle=False)
    return {
        str(key): torch.as_tensor(value, device="cuda")
        for key, value in zip(data["keys"], data["segments"])
    }


def loadlabels(root: Path, endpoint: str) -> pd.DataFrame:
    """Load authorized train/development labels with split predicate pushdown."""

    path = Path(root) / "registry.parquet"
    metadata = pd.read_parquet(path, columns=FIELDS).reset_index(names="source_row")
    selected = metadata.loc[
        (metadata.endpoint == endpoint)
        & metadata.dual_cold_split.isin(["train", "development"])
    ].reset_index(drop=True)
    labeled = pd.read_parquet(
        path,
        columns=FIELDS + ["affinity"],
        filters=[
            ("endpoint", "==", endpoint),
            ("dual_cold_split", "in", ["train", "development"]),
        ],
    ).reset_index(drop=True)
    if not selected[FIELDS].equals(labeled[FIELDS]):
        raise ValueError("predicate-filtered labels do not align with registry source rows")
    selected["affinity"] = labeled.affinity.to_numpy()
    if not np.isfinite(selected.affinity.to_numpy(dtype=np.float64)).all():
        raise ValueError("authorized affinity labels must be finite")
    return selected


def maketrainroster(
    frame: pd.DataFrame,
    protein: dict[str, torch.Tensor],
    feature: np.ndarray,
    *,
    targets: int,
    querycap: int,
    support: int = 5,
) -> list[Episode]:
    """Apply the same label-free query-span selector used by the frozen roster."""

    rows = preparetable(frame[FIELDS].to_dict("records"), split="meta_train")
    available = [target for target in sorted(frame.target.unique()) if target in protein]
    active: list[Episode] = []
    for target in available:
        if targets and len(active) >= targets:
            break
        candidates = tuple(frame.index[frame.target == target].tolist())
        if len(candidates) < support + 1:
            continue
        design = {
            index: feature[int(frame.source_row.iat[index]), :64]
            for index in candidates
        }
        components = {index: str(frame.conn.iat[index]) for index in candidates}
        try:
            chosen = selectsupport(
                rows,
                candidates,
                support,
                design,
                chemical_components=components,
                chemical_component_cap=1,
            )
        except ValueError:
            continue
        supportscaffolds = {str(frame.scaffold.iat[index]) for index in chosen}
        supportcomponents = {components[index] for index in chosen}
        query = tuple(
            index
            for index in candidates
            if index not in chosen
            and str(frame.scaffold.iat[index]) not in supportscaffolds
            and components[index] not in supportcomponents
        )[:querycap]
        if not query:
            continue
        active.append(
            Episode.create(
                target_key=target,
                support_indices=chosen,
                query_indices=query,
                homology_component=str(frame.hcluster.iat[chosen[0]]),
                support_scaffolds=frame.scaffold.iloc[list(chosen)],
                query_scaffolds=frame.scaffold.iloc[list(query)],
                provenance_components=frame.docs.iloc[list(chosen + query)],
            )
        )
    if not active:
        raise RuntimeError("no strict k=5 train episodes are available")
    return active


def verifyroster(
    path: Path,
    root: Path,
    endpoint: str,
) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any]]:
    """Verify the complete frozen-roster certificate chain before label access."""

    path = Path(path)
    sidecar = path.with_suffix(".json")
    indexpath = path.with_name("episodes.v1.json")
    for required in (path, sidecar, indexpath):
        if not required.is_file():
            raise FileNotFoundError(f"required frozen-roster file is missing: {required}")

    try:
        report = json.loads(sidecar.read_text(encoding="utf-8"))
        collection = json.loads(indexpath.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError("frozen-roster JSON is unreadable") from error
    if not isinstance(report, dict) or not isinstance(collection, dict):
        raise ValueError("frozen-roster JSON must contain objects")
    if report.get("schema") != "strict-fewshot-roster.v1":
        raise ValueError("sidecar schema is not the frozen strict-roster schema")
    if collection.get("schema") != "strict-fewshot-roster-index.v1":
        raise ValueError("roster index schema is invalid")
    if report.get("endpoint") != endpoint:
        raise ValueError("roster endpoint does not match the requested endpoint")
    if report.get("labels_read") is not False:
        raise ValueError("roster preprocessing must be label blind")

    protocol = report.get("protocol")
    if not isinstance(protocol, dict):
        raise ValueError("sidecar protocol is missing")
    if protocol.get("schema") != "strict-fewshot-roster.v1":
        raise ValueError("protocol schema is invalid")
    if protocol.get("affinity_labels_read") is not False:
        raise ValueError("protocol opened development affinity labels")
    if protocol.get("confirmation_or_sealed_affinity_labels_read") is not False:
        raise ValueError("protocol opened forbidden confirmation labels")
    if int(protocol.get("support_rows", 0)) != 5:
        raise ValueError("active training requires a frozen k=5 roster")
    if int(protocol.get("minimum_query_rows", 0)) < 1:
        raise ValueError("protocol minimum query depth is invalid")
    separation = protocol.get("support_query_separation")
    allowedseparation = {"chemical component", "document token", "assay token"}
    if (
        not isinstance(separation, list)
        or len(separation) != len(set(separation))
        or "chemical component" not in separation
        or not set(separation).issubset(allowedseparation)
    ):
        raise ValueError("protocol support/query separation is invalid")
    protocolhash = contenthash(protocol)
    if report.get("protocol_sha256") != protocolhash:
        raise ValueError("protocol hash does not match its parameters")

    sources = report.get("sources")
    if not isinstance(sources, dict) or set(sources) != set(SOURCEFILES):
        raise ValueError("sidecar source manifest is incomplete")
    sourcehashes: dict[str, str] = {}
    for name in SOURCEFILES:
        source = sources.get(name)
        if not isinstance(source, dict) or source.get("path") != name:
            raise ValueError(f"source manifest entry is invalid: {name}")
        expectedhash = source.get("sha256")
        if not isinstance(expectedhash, str) or len(expectedhash) != 64:
            raise ValueError(f"source hash is invalid: {name}")
        sourcepath = Path(root) / name
        if not sourcepath.is_file():
            raise FileNotFoundError(f"registered source is missing: {sourcepath}")
        actualhash = filehash(sourcepath)
        if actualhash != expectedhash:
            raise ValueError(f"registered source hash changed: {name}")
        sourcehashes[name] = actualhash

    if collection.get("protocol") != protocol or collection.get("sources") != sources:
        raise ValueError("roster index disagrees with the sidecar protocol or sources")
    matches = [
        item
        for item in collection.get("results", [])
        if isinstance(item, dict) and item.get("endpoint") == endpoint
    ]
    if len(matches) != 1:
        raise ValueError("roster index must contain one requested endpoint")
    indexed = matches[0]
    if indexed.get("roster") != report.get("roster"):
        raise ValueError("roster index and sidecar metadata disagree")
    if indexed.get("summary") != report.get("summary"):
        raise ValueError("roster index and sidecar summary disagree")
    indexedsidecar = indexed.get("sidecar")
    if not isinstance(indexedsidecar, dict) or indexedsidecar.get("path") != sidecar.name:
        raise ValueError("roster index points to the wrong sidecar")
    sidecarhash = filehash(sidecar)
    if indexedsidecar.get("sha256") != sidecarhash:
        raise ValueError("sidecar content hash does not match the roster index")

    rosterinfo = report.get("roster")
    if not isinstance(rosterinfo, dict) or rosterinfo.get("path") != path.name:
        raise ValueError("sidecar points to the wrong roster")
    rosterhash = filehash(path)
    if rosterinfo.get("file_sha256") != rosterhash:
        raise ValueError("Parquet file hash does not match the sidecar")
    try:
        parquet = pq.ParquetFile(path)
        metadata = parquet.schema_arrow.metadata or {}
        roster = pd.read_parquet(path)
    except (OSError, ValueError) as error:
        raise ValueError("frozen roster Parquet is unreadable") from error
    expectedmetadata = {
        b"fort.protocol": canonical(protocol),
        b"fort.sources": canonical(sources),
        b"fort.roster_content_sha256": str(rosterinfo.get("content_sha256")),
        b"fort.episodes": canonical(report.get("episodes")),
    }
    for key, expected in expectedmetadata.items():
        value = metadata.get(key)
        try:
            decoded = value.decode("utf-8") if value is not None else None
        except UnicodeDecodeError as error:
            raise ValueError(f"Parquet certificate is not UTF-8: {key!r}") from error
        if decoded != expected:
            raise ValueError(f"Parquet certificate disagrees with sidecar: {key!r}")

    if list(roster.columns) != ROSTERFIELDS:
        raise ValueError("roster columns do not match the frozen schema")
    contentdigest = framehash(roster)
    if rosterinfo.get("content_sha256") != contentdigest:
        raise ValueError("roster content hash does not match the sidecar")
    if rosterinfo.get("columns") != ROSTERFIELDS or rosterinfo.get("rows") != len(roster):
        raise ValueError("roster shape disagrees with the sidecar")

    episodes = report.get("episodes")
    if not isinstance(episodes, list) or not episodes:
        raise ValueError("sidecar contains no episodes")
    identifiers = [item.get("episode") for item in episodes if isinstance(item, dict)]
    if len(identifiers) != len(episodes) or len(identifiers) != len(set(identifiers)):
        raise ValueError("sidecar episode identifiers are invalid")
    if set(roster.episode) != set(identifiers):
        raise ValueError("roster and sidecar contain different episodes")
    targets = [item.get("target") for item in episodes]
    if len(targets) != len(set(targets)):
        raise ValueError("frozen roster must contain one episode per target")

    supportsize = int(protocol["support_rows"])
    minimumquery = int(protocol["minimum_query_rows"])
    closurekeys = [("support_scaffolds", "query_scaffolds")]
    if "chemical component" in separation:
        closurekeys.append(
            ("support_chemical_components", "query_chemical_components")
        )
    if "document token" in separation:
        closurekeys.append(("support_documents", "query_documents"))
    if "assay token" in separation:
        closurekeys.append(("support_assays", "query_assays"))
    for item in episodes:
        identifier = item["episode"]
        if identifier != episodehash(item):
            raise ValueError(f"episode hash is invalid: {identifier}")
        if item.get("protocol_sha256") != protocolhash:
            raise ValueError(f"episode protocol hash is invalid: {identifier}")
        group = roster.loc[roster.episode == identifier]
        if set(group.role) != {"support", "query"} or group.source_row.duplicated().any():
            raise ValueError(f"episode roles or rows are invalid: {identifier}")
        supportrows = sorted(
            int(value) for value in group.loc[group.role == "support", "source_row"]
        )
        queryrows = sorted(
            int(value) for value in group.loc[group.role == "query", "source_row"]
        )
        if supportrows != item.get("support_source_rows"):
            raise ValueError(f"episode support rows disagree with roster: {identifier}")
        if queryrows != item.get("query_source_rows"):
            raise ValueError(f"episode query rows disagree with roster: {identifier}")
        if len(supportrows) != supportsize or len(queryrows) < minimumquery:
            raise ValueError(f"episode role counts violate protocol: {identifier}")
        if (
            item.get("support_count") != len(supportrows)
            or item.get("query_count") != len(queryrows)
        ):
            raise ValueError(f"episode sidecar counts are invalid: {identifier}")
        for column in ("target", "endpoint", "homology_component", "rank"):
            if set(group[column]) != {item.get(column)}:
                raise ValueError(f"episode {column} disagrees with roster: {identifier}")
        for role, key in (
            ("support", "support_chemical_components"),
            ("query", "query_chemical_components"),
        ):
            values = sorted(set(group.loc[group.role == role, "chemical_component"]))
            if values != item.get(key):
                raise ValueError(f"episode chemical components disagree: {identifier}")
        for left, right in closurekeys:
            if set(item.get(left, ())).intersection(item.get(right, ())):
                raise ValueError(f"episode support/query closure is invalid: {identifier}")

    summary = report.get("summary")
    if (
        not isinstance(summary, dict)
        or summary.get("endpoint") != endpoint
        or summary.get("support") != supportsize
        or summary.get("episodes") != len(episodes)
    ):
        raise ValueError("sidecar summary is inconsistent")
    proof = {
        "path": str(path.resolve()),
        "endpoint": endpoint,
        "sidecar_path": str(sidecar.resolve()),
        "index_path": str(indexpath.resolve()),
        "file_sha256": rosterhash,
        "content_sha256": contentdigest,
        "sidecar_sha256": sidecarhash,
        "protocol_sha256": protocolhash,
        "source_sha256": sourcehashes,
        "support_query_separation": list(separation),
        "frozen_episode_count": len(episodes),
        "frozen_row_count": len(roster),
        "support_rows": supportsize,
        "minimum_query_rows": minimumquery,
        "labels_read_during_preprocessing": False,
    }
    return roster, report, proof


def loadroster(
    frame: pd.DataFrame,
    protein: dict[str, torch.Tensor],
    path: Path,
    targets: int,
    querycap: int,
    *,
    root: Path,
    endpoint: str,
    verified: tuple[pd.DataFrame, dict[str, Any], dict[str, Any]] | None = None,
) -> tuple[list[Episode], dict[str, Any]]:
    """Load the frozen label-free development roster without reselecting rows."""

    if verified is None:
        roster, report, proof = verifyroster(path, root, endpoint)
    else:
        roster, report, proof = verified
        if proof.get("path") != str(Path(path).resolve()) or report.get("endpoint") != endpoint:
            raise ValueError("verified roster does not match the requested roster")
        proof = dict(proof)
    minimumquery = int(report["protocol"]["minimum_query_rows"])
    if targets < 0:
        raise ValueError("target cap must be non-negative")
    if querycap < minimumquery:
        raise ValueError(
            f"query cap must preserve the protocol minimum of {minimumquery} rows"
        )
    if frame.source_row.duplicated().any():
        raise ValueError("development frame contains duplicate source rows")
    position = {int(row): index for index, row in enumerate(frame.source_row)}
    payloads = sorted(
        report["episodes"], key=lambda item: (str(item["target"]), str(item["episode"]))
    )
    if targets:
        payloads = payloads[:targets]
    active: list[Episode] = []
    uncappedqueries = 0
    for payload in payloads:
        target = str(payload["target"])
        if target not in protein:
            raise ValueError(f"frozen roster target lacks a protein embedding: {target}")
        supportrows = tuple(int(row) for row in payload["support_source_rows"])
        queryrows = tuple(int(row) for row in payload["query_source_rows"])
        missing = [row for row in supportrows + queryrows if row not in position]
        if missing:
            raise ValueError(f"frozen roster rows are absent from development: {missing[:5]}")
        support = tuple(position[row] for row in supportrows)
        query = tuple(position[row] for row in queryrows[:querycap])
        allindices = support + tuple(position[row] for row in queryrows)
        if set(frame.target.iloc[list(allindices)].astype(str)) != {target}:
            raise ValueError(f"frozen roster target disagrees with registry rows: {target}")
        if set(frame.endpoint.iloc[list(allindices)].astype(str)) != {endpoint}:
            raise ValueError(f"frozen roster endpoint disagrees with registry rows: {target}")

        def rolevalues(indices: tuple[int, ...], column: str) -> list[str]:
            return sorted(set(frame[column].iloc[list(indices)].astype(str)))

        def rolemetatokens(indices: tuple[int, ...], column: str) -> list[str]:
            return sorted(
                set().union(*(tokens(value) for value in frame[column].iloc[list(indices)]))
            )

        fullquery = tuple(position[row] for row in queryrows)
        checks = (
            (rolevalues(support, "conn"), payload["support_ligand_parents"]),
            (rolevalues(fullquery, "conn"), payload["query_ligand_parents"]),
            (rolevalues(support, "scaffold"), payload["support_scaffolds"]),
            (rolevalues(fullquery, "scaffold"), payload["query_scaffolds"]),
            (rolemetatokens(support, "docs"), payload["support_documents"]),
            (rolemetatokens(fullquery, "docs"), payload["query_documents"]),
            (rolemetatokens(support, "assays"), payload["support_assays"]),
            (rolemetatokens(fullquery, "assays"), payload["query_assays"]),
        )
        if any(actual != expected for actual, expected in checks):
            raise ValueError(f"frozen episode metadata disagrees with registry: {target}")
        uncappedqueries += len(fullquery)
        active.append(
            Episode.create(
                target_key=target,
                support_indices=support,
                query_indices=query,
                homology_component=str(payload["homology_component"]),
                support_scaffolds=frame.scaffold.iloc[list(support)],
                query_scaffolds=frame.scaffold.iloc[list(query)],
                provenance_components=frame.docs.iloc[list(support + query)],
            )
        )
    if not active:
        raise RuntimeError("frozen roster contains no loadable development episodes")
    proof.update(
        {
            "target_cap": targets if targets else None,
            "query_cap": querycap,
            "selected_episode_count": len(active),
            "selected_query_rows_before_cap": uncappedqueries,
            "selected_query_rows_after_cap": sum(len(item.query_indices) for item in active),
        }
    )
    return active, proof


def episodetensors(
    frame: pd.DataFrame,
    feature: np.ndarray,
    item: Episode,
    base: np.ndarray,
) -> tuple[torch.Tensor, ...]:
    support = list(item.support_indices)
    query = list(item.query_indices)
    supportrows = frame.source_row.iloc[support].to_numpy()
    queryrows = frame.source_row.iloc[query].to_numpy()
    return (
        torch.as_tensor(feature[supportrows], device="cuda"),
        torch.as_tensor(frame.affinity.iloc[support].to_numpy(dtype=np.float32), device="cuda"),
        torch.as_tensor(feature[queryrows], device="cuda"),
        torch.as_tensor(frame.affinity.iloc[query].to_numpy(dtype=np.float32), device="cuda"),
        torch.as_tensor(base[support], device="cuda"),
        torch.as_tensor(base[query], device="cuda"),
    )


def fitbase(
    train: pd.DataFrame,
    development: pd.DataFrame,
    feature: np.ndarray,
    ridge: float = 10.0,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Fit a deterministic train-only ridge B0 and its residual variance."""

    if ridge <= 0:
        raise ValueError("base ridge must be positive")
    rows = train.source_row.to_numpy()
    values = torch.as_tensor(feature[rows], device="cuda")
    labels = torch.as_tensor(train.affinity.to_numpy(dtype=np.float32), device="cuda")
    center = values.mean(dim=0)
    scale = values.std(dim=0).clamp_min(1e-6)
    standardized = (values - center) / scale
    labelcenter = labels.mean()
    centeredlabels = labels - labelcenter
    gram = standardized.T @ standardized
    gram = gram + ridge * torch.eye(gram.shape[0], device="cuda", dtype=gram.dtype)
    weight = torch.linalg.solve(gram, standardized.T @ centeredlabels)
    batchsize = 8192
    with torch.no_grad():
        trainprediction = standardized @ weight + labelcenter
        trainbase = trainprediction.cpu().numpy().astype(np.float32)
        basevariance = float((trainprediction - labels).square().mean().item())
        developmentbase: list[np.ndarray] = []
        developmentrows = development.source_row.to_numpy()
        for start in range(0, len(developmentrows), batchsize):
            batch = torch.as_tensor(
                feature[developmentrows[start : start + batchsize]], device="cuda"
            )
            developmentbase.append(
                (((batch - center) / scale) @ weight + labelcenter).cpu().numpy()
            )
    return trainbase, np.concatenate(developmentbase).astype(np.float32), basevariance


def fitjoint(
    model: TargetAdapter,
    frame: pd.DataFrame,
    protein: dict[str, torch.Tensor],
    feature: np.ndarray,
    episodes: list[Episode],
    base: np.ndarray,
    epochs: int,
    nllweight: float,
    contrastweight: float,
) -> None:
    optimizer = torch.optim.AdamW(
        model.protein.parameters(), lr=2e-4, weight_decay=1e-4
    )
    for _ in range(epochs):
        model.train()
        for item in episodes:
            sx, sy, qx, qy, sb, qb = episodetensors(frame, feature, item, base)
            state = model.adapt(
                proteintokens=protein[item.target_key],
                supportligand=sx,
                supportlabel=sy,
                supportb0=sb,
            )
            output = model.predict(
                state,
                queryligand=qx,
                queryb0=qb,
                mode="joint",
            )
            variance = output["totalvariance"].clamp_min(1e-5)
            nll = 0.5 * (
                torch.log(variance) + (output["prediction"] - qy).square() / variance
            )
            loss = (
                F.huber_loss(output["prediction"], qy)
                + contrastweight * contrastloss(output["prediction"], qy)
                + nllweight * nll.mean()
            )
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.protein.parameters(), 1.0)
            optimizer.step()


def contrastloss(prediction: torch.Tensor, label: torch.Tensor) -> torch.Tensor:
    """Train one transitive score from antisymmetric query differences."""

    prediction = prediction.float().reshape(-1)
    label = label.float().reshape(-1)
    if prediction.numel() != label.numel():
        raise ValueError("contrast predictions and labels must align")
    if prediction.numel() <= 1:
        return prediction.sum() * 0.0
    left, right = torch.triu_indices(
        prediction.numel(), prediction.numel(), offset=1, device=prediction.device
    )
    return F.huber_loss(
        prediction[left] - prediction[right],
        label[left] - label[right],
    )


def fitligand(
    model: LigandBaseline,
    frame: pd.DataFrame,
    feature: np.ndarray,
    episodes: list[Episode],
    base: np.ndarray,
    epochs: int,
) -> None:
    optimizer = torch.optim.AdamW(model.parameters(), lr=2e-4, weight_decay=1e-4)
    for _ in range(epochs):
        model.train()
        for item in episodes:
            sx, sy, qx, qy, sb, qb = episodetensors(frame, feature, item, base)
            output = model(qx, sx, sy, query_base=qb, support_base=sb)
            variance = output["total_variance"].clamp_min(1e-5)
            nll = 0.5 * (torch.log(variance) + (output["pred"] - qy).square() / variance)
            loss = F.huber_loss(output["pred"], qy) + 0.05 * nll.mean()
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()


def fitgradient(
    model: torch.nn.Module,
    frame: pd.DataFrame,
    protein: dict[str, torch.Tensor],
    feature: np.ndarray,
    episodes: list[Episode],
    base: np.ndarray,
    epochs: int,
) -> None:
    optimizer = torch.optim.AdamW(model.parameters(), lr=2e-4, weight_decay=1e-4)
    for _ in range(epochs):
        model.train()
        for item in episodes:
            sx, sy, qx, qy, sb, qb = episodetensors(frame, feature, item, base)
            state = adapttarget(
                model,
                protein_tokens=protein[item.target_key],
                support_ligand=sx,
                support_y=sy,
                support_b0=sb,
            )
            output = predictquery(
                model,
                state,
                protein_tokens=protein[item.target_key],
                query_ligand=qx,
                query_b0=qb,
            )
            optimizer.zero_grad(set_to_none=True)
            F.huber_loss(output["prediction"], qy).backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()


def buildreordering(
    *,
    backbone: str,
    gateconfig: tuple[bool, bool],
    seed: int,
    proteinconditioned: bool = False,
) -> ReorderingModel:
    """Build backbone arms with identical initialization for shared parameters."""

    torch.manual_seed(seed)
    reference = ReorderingModel(
        proteindim=1280,
        liganddim=1034,
        backbone="hybrid",
        proteinconditioned=proteinconditioned,
        interactiononly=True,
        rankinghardgate=gateconfig[0],
        calibrationhardgate=gateconfig[1],
    )
    if backbone == "hybrid":
        return reference.cuda()
    candidate = ReorderingModel(
        proteindim=1280,
        liganddim=1034,
        backbone=backbone,
        proteinconditioned=proteinconditioned,
        interactiononly=True,
        rankinghardgate=gateconfig[0],
        calibrationhardgate=gateconfig[1],
    )
    candidatevalues = candidate.state_dict()
    shared = {
        name: value
        for name, value in reference.state_dict().items()
        if name in candidatevalues and value.shape == candidatevalues[name].shape
    }
    missing = set(candidatevalues) - set(shared)
    if missing:
        raise RuntimeError(f"unmatched {backbone} parameters: {sorted(missing)}")
    candidate.load_state_dict(shared)
    return candidate.cuda()


def evaluate(
    full: TargetAdapter,
    gradient: torch.nn.Module,
    frame: pd.DataFrame,
    protein: dict[str, torch.Tensor],
    feature: np.ndarray,
    episodes: list[Episode],
    base: np.ndarray,
    observationvariance: float,
) -> tuple[
    dict[str, dict[str, float]],
    dict[str, Any],
    dict[str, Any],
    list[dict[str, Any]],
]:
    for model in (full, gradient):
        model.eval()
    names = (
        "b0",
        "supportmean",
        "residualmean",
        "calibration",
        "ligand",
        "gradient",
        "bayesian",
        "joint",
        "exact",
        "soft",
        "wrongsupport",
        "wrongprotein",
        "permutedlabels",
        "proteinfree",
    )
    predictions = {name: [] for name in names}
    variances: dict[str, list[float]] = {
        name: []
        for name in (
            "b0",
            "calibration",
            "ligand",
            "bayesian",
            "joint",
            "exact",
            "soft",
            "wrongsupport",
            "wrongprotein",
            "permutedlabels",
            "proteinfree",
        )
    }
    labels: list[float] = []
    indices: list[int] = []
    ranklogbf: list[float] = []
    rankprobability: list[float] = []
    ledger: list[dict[str, Any]] = []

    for position, item in enumerate(episodes):
        sx, sy, qx, qy, sb, qb = episodetensors(frame, feature, item, base)
        wrongitem = episodes[(position + 1) % len(episodes)]
        for offset in range(1, len(episodes)):
            candidate = episodes[(position + offset) % len(episodes)]
            if candidate.homology_component != item.homology_component:
                wrongproteinitem = candidate
                break
        else:
            raise RuntimeError("wrong-protein control needs a distinct homology component")
        wsx, wsy, _, _, wsb, _ = episodetensors(frame, feature, wrongitem, base)
        permutation = torch.roll(torch.arange(len(sy), device="cuda"), shifts=1)

        fullstate = full.adapt(
            proteintokens=protein[item.target_key],
            supportligand=sx,
            supportlabel=sy,
            supportb0=sb,
        )
        gateoutputs = {
            name: full.predict(
                fullstate,
                queryligand=qx,
                queryb0=qb,
                mode=name,
            )
            for name in ("joint", "exact", "soft")
        }
        fulloutput = gateoutputs[full.mode]
        wrongproteinstate = full.adapt(
            proteintokens=protein[wrongproteinitem.target_key],
            supportligand=sx,
            supportlabel=sy,
            supportb0=sb,
        )
        wrongproteinoutput = full.predict(
            wrongproteinstate,
            queryligand=qx,
            queryb0=qb,
        )
        wrongstate = full.adapt(
            proteintokens=protein[item.target_key],
            supportligand=wsx,
            supportlabel=wsy,
            supportb0=wsb,
        )
        wrongoutput = full.predict(
            wrongstate,
            queryligand=qx,
            queryb0=qb,
        )
        permutedstate = full.adapt(
            proteintokens=protein[item.target_key],
            supportligand=sx,
            supportlabel=sy[permutation],
            supportb0=sb,
        )
        permutedoutput = full.predict(
            permutedstate,
            queryligand=qx,
            queryb0=qb,
        )
        freeoutput = full.predict(
            fullstate,
            queryligand=qx,
            queryb0=qb,
            useprotein=False,
        )
        ligandoutput = full.predict(
            fullstate, queryligand=qx, queryb0=qb, useprotein=False
        )
        calibration = full.protein.posterior.calibration(sy - sb, sb, qb)
        calibrationprediction = qb + calibration["weight"] * calibration["mean"]
        calibrationvariance = (
            calibration["weight"] * calibration["variance"]
            + calibration["weight"]
            * (1.0 - calibration["weight"])
            * calibration["mean"].square()
        )
        gradientstate = adapttarget(
            gradient,
            protein_tokens=protein[item.target_key],
            support_ligand=sx,
            support_y=sy,
            support_b0=sb,
        )
        gradientoutput = predictquery(
            gradient,
            gradientstate,
            protein_tokens=protein[item.target_key],
            query_ligand=qx,
            query_b0=qb,
        )

        armvalues = {
            "b0": qb,
            "supportmean": sy.mean().expand_as(qb),
            "residualmean": qb + (sy - sb).mean(),
            "calibration": calibrationprediction,
            "ligand": ligandoutput["prediction"],
            "gradient": gradientoutput["prediction"],
            "bayesian": fulloutput["prediction"],
            "joint": gateoutputs["joint"]["prediction"],
            "exact": gateoutputs["exact"]["prediction"],
            "soft": gateoutputs["soft"]["prediction"],
            "wrongsupport": wrongoutput["prediction"],
            "wrongprotein": wrongproteinoutput["prediction"],
            "permutedlabels": permutedoutput["prediction"],
            "proteinfree": freeoutput["prediction"],
        }
        armvariances = {
            "b0": torch.full_like(qb, observationvariance),
            "calibration": torch.full_like(qb, observationvariance)
            + calibrationvariance,
            "ligand": ligandoutput["totalvariance"],
            "bayesian": fulloutput["totalvariance"],
            "joint": gateoutputs["joint"]["totalvariance"],
            "exact": gateoutputs["exact"]["totalvariance"],
            "soft": gateoutputs["soft"]["totalvariance"],
            "wrongsupport": wrongoutput["totalvariance"],
            "wrongprotein": wrongproteinoutput["totalvariance"],
            "permutedlabels": permutedoutput["totalvariance"],
            "proteinfree": freeoutput["totalvariance"],
        }
        ledger.append(
            {
                "target": item.target_key,
                "component": item.homology_component,
                "support_rows": [
                    int(value)
                    for value in frame.source_row.iloc[list(item.support_indices)]
                ],
                "query_rows": [
                    int(value)
                    for value in frame.source_row.iloc[list(item.query_indices)]
                ],
                "query_depth": len(item.query_indices),
                "support_residual_rms": float(
                    torch.sqrt((sy - sb).square().mean())
                    .detach()
                    .cpu()
                ),
                "ranking_logbf": float(fulloutput["proteinlogbf"].detach().cpu()),
                "ranking_probability": float(
                    fulloutput["proteinprobability"].detach().cpu()
                ),
                "calibration_logbf": float(calibration["logbf"].detach().cpu()),
                "labels": qy.detach().cpu().tolist(),
                "predictions": {
                    name: values.detach().cpu().tolist()
                    for name, values in armvalues.items()
                },
                "variances": {
                    name: values.detach().cpu().tolist()
                    for name, values in armvariances.items()
                },
            }
        )
        for name, values in armvalues.items():
            predictions[name].extend(values.detach().cpu().tolist())
        for name, values in armvariances.items():
            variances[name].extend(values.detach().cpu().tolist())
        labels.extend(qy.detach().cpu().tolist())
        indices.extend(item.query_indices)
        ranklogbf.append(float(fulloutput["proteinlogbf"].detach().cpu()))
        rankprobability.append(float(fulloutput["proteinprobability"].detach().cpu()))

    components = {item.target_key: item.homology_component for item in episodes}
    metrics = {
        name: evaluateprotocol(
            predictions=values,
            labels=labels,
            episodes=episodes,
            prediction_indices=indices,
            component_by_target=components,
            variances=variances.get(name),
        )
        for name, values in predictions.items()
    }
    diagnostics = {
        "ranking_logbf": {
            "mean": float(np.mean(ranklogbf)),
            "median": float(np.median(ranklogbf)),
            "minimum": float(np.min(ranklogbf)),
            "maximum": float(np.max(ranklogbf)),
        },
        "ranking_probability": {
            "mean": float(np.mean(rankprobability)),
            "median": float(np.median(rankprobability)),
            "minimum": float(np.min(rankprobability)),
            "maximum": float(np.max(rankprobability)),
        },
    }
    attribution = pairedcomponents(
        predictions=predictions,
        labels=labels,
        episodes=episodes,
        prediction_indices=indices,
        component_by_target=components,
        reference="bayesian",
    )
    return metrics, diagnostics, attribution, ledger


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--endpoint", choices=("pKd", "pKi"), default="pKi")
    parser.add_argument("--epochs", type=int, default=1)
    parser.add_argument("--targets", type=int, default=64)
    parser.add_argument("--queries", type=int, default=64)
    parser.add_argument("--baseridge", type=float, default=10.0)
    parser.add_argument("--bayesnll", type=float, default=0.0)
    parser.add_argument("--contrast", type=float, default=0.0)
    parser.add_argument(
        "--gate",
        choices=("joint", "exact", "soft"),
        default="exact",
    )
    parser.add_argument(
        "--backbone", choices=("hybrid", "transformer"), default="hybrid"
    )
    parser.add_argument(
        "--protein-conditioned",
        action="store_true",
        help="use a protein-generated basis instead of the shared exact-null basis",
    )
    parser.add_argument("--roster", type=Path, default=ROSTER)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--out", type=Path, default=REPORT / "adaptjoint.v1.json")
    args = parser.parse_args()
    assertauthorized()
    if not torch.cuda.is_available():
        raise RuntimeError("training requires CUDA")
    if args.endpoint != "pKi":
        raise ValueError("the active adaptation experiment is registered for pKi only")
    if args.bayesnll < 0:
        raise ValueError("Bayesian NLL weight must be non-negative")
    if args.contrast < 0:
        raise ValueError("contrast weight must be non-negative")
    if args.targets < 0:
        raise ValueError("target cap must be non-negative")
    if args.queries <= 0:
        raise ValueError("query cap must be positive")

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    torch.cuda.reset_peak_memory_stats()
    verifiedroster = verifyroster(args.roster, args.root, args.endpoint)
    minimumquery = int(verifiedroster[1]["protocol"]["minimum_query_rows"])
    if args.queries < minimumquery:
        raise ValueError(
            f"query cap must preserve the protocol minimum of {minimumquery} rows"
        )
    raw = loadlabels(args.root, args.endpoint)
    train = raw.loc[raw.dual_cold_split == "train"].reset_index(drop=True)
    development = raw.loc[raw.dual_cold_split == "development"].reset_index(drop=True)
    rawfeature = np.load(args.root / "ligand_features.npz", allow_pickle=False)["feat"]
    if rawfeature.shape[1] != 1034:
        raise ValueError(f"expected 1034 ligand features, received {rawfeature.shape[1]}")
    feature, ligandnormalization = normalizeligands(
        rawfeature, train.source_row.to_numpy(), descriptors=10
    )
    protein = loadproteins(args.root)

    stop = threading.Event()
    samples: list[dict[str, float]] = []
    watcher = threading.Thread(target=telemetry, args=(stop, samples), daemon=True)
    watcher.start()
    started = time.perf_counter()
    try:
        trainepisodes = maketrainroster(
            train,
            protein,
            feature,
            targets=args.targets,
            querycap=args.queries,
        )
        developmentepisodes, rosterproof = loadroster(
            development,
            protein,
            args.roster,
            args.targets,
            args.queries,
            root=args.root,
            endpoint=args.endpoint,
            verified=verifiedroster,
        )
        trainbase, developmentbase, observationvariance = fitbase(
            train, development, feature, args.baseridge
        )
        gateconfig = {
            "joint": (False, False),
            "exact": (True, False),
            "soft": (False, False),
        }[args.gate]
        fullincrement = buildreordering(
            backbone=args.backbone,
            gateconfig=gateconfig,
            seed=args.seed,
            proteinconditioned=args.protein_conditioned,
        )
        torch.manual_seed(args.seed + 1)
        ligand = LigandBaseline(bits=1034, rep=128, pdim=8).cuda()
        fitligand(ligand, train, feature, trainepisodes, trainbase, args.epochs)
        full = TargetAdapter(ligand, fullincrement, mode=args.gate).cuda()
        torch.manual_seed(args.seed + 2)
        gradient = buildadapter(protein_dim=1280, ligand_dim=1034).cuda()

        fitjoint(
            full,
            train,
            protein,
            feature,
            trainepisodes,
            trainbase,
            args.epochs,
            args.bayesnll,
            args.contrast,
        )
        fitgradient(
            gradient,
            train,
            protein,
            feature,
            trainepisodes,
            trainbase,
            args.epochs,
        )
        metrics, diagnostics, attribution, ledger = evaluate(
            full,
            gradient,
            development,
            protein,
            feature,
            developmentepisodes,
            developmentbase,
            observationvariance,
        )
    finally:
        stop.set()
        watcher.join(timeout=2)

    elapsed = time.perf_counter() - started
    parametercounts = {
        "bayesian": sum(parameter.numel() for parameter in full.parameters()),
        "ligand": sum(parameter.numel() for parameter in ligand.parameters()),
        "gradient": sum(parameter.numel() for parameter in gradient.parameters()),
    }
    result = {
        "protocol": "LOO ligand residualized calibration-null contrast posterior",
        "seed": args.seed,
        "endpoint": args.endpoint,
        "support": 5,
        "ambient_dimension": 8,
        "posterior_rank": 2,
        "selector": "runtime label-free TRAIN selector; frozen certified development roster",
        "epochs": args.epochs,
        "base_ridge": args.baseridge,
        "bayesian_nll_weight": args.bayesnll,
        "contrast_weight": args.contrast,
        "gate": args.gate,
        "backbone": args.backbone,
        "protein_conditioned_basis": args.protein_conditioned,
        "ranking_inclusion_prior": full.inclusion,
        "ligand_base": "shared frozen TRAIN-only posterior",
        "support_likelihood": "leave-one-out ligand residuals projected off [1, ligand baseline]",
        "reordering_feature": (
            "calibration-null protein interaction rank feature"
            if args.protein_conditioned
            else "calibration-null frozen ligand rank feature"
        ),
        "train_episodes": len(trainepisodes),
        "development_episodes": len(developmentepisodes),
        "development_roster": rosterproof,
        "base_residual_variance": observationvariance,
        "ligand_normalization": ligandnormalization,
        "parameters": parametercounts,
        "metrics": metrics,
        "diagnostics": diagnostics,
        "paired_component_bootstrap": attribution,
        "episode_ledger": ledger,
        "family_dominance": "not evaluated: the frozen registry has no family field",
        "training_seconds": elapsed,
        "gpu": torch.cuda.get_device_name(0),
        "gpu_telemetry": {
            "samples": len(samples),
            "mean_utilization_percent": float(np.mean([item["utilization"] for item in samples])) if samples else None,
            "peak_utilization_percent": float(np.max([item["utilization"] for item in samples])) if samples else None,
            "mean_power_watts": float(np.mean([item["power"] for item in samples])) if samples else None,
            "peak_power_watts": float(np.max([item["power"] for item in samples])) if samples else None,
            "peak_nvidia_memory_mib": float(np.max([item["memory"] for item in samples])) if samples else None,
            "peak_torch_memory_mib": float(torch.cuda.max_memory_allocated() / 2**20),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    checkpoint = args.out.with_suffix(".pt")
    torch.save(
        {
            "state": full.state_dict(),
            "proteindim": 1280,
            "liganddim": 1034,
            "ambientdim": 8,
            "rank": 2,
            "backbone": args.backbone,
            "gate": args.gate,
            "proteininclusion": full.inclusion,
            "proteinconditioned": args.protein_conditioned,
            "architecture": "LOO ligand residualized calibration-null contrast posterior",
            "ligandnormalization": ligandnormalization,
            "developmentroster": rosterproof,
        },
        checkpoint,
    )
    result["checkpoint"] = str(checkpoint)
    args.out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
