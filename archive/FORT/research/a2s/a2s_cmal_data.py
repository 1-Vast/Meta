"""Build the label-blind, model-ready episode package for A2S-CMAL.

The formal ChEMBL corpus and the D0-R recipient roster are immutable inputs.
This builder adds only information that may be fixed before model fitting:

* homology-component source task splits;
* label-blind nested support/query episodes;
* frozen ESM-2 target embeddings; and
* target-mismatched counterfactual support mappings.

Affinity columns are deliberately never requested from parquet.  The emitted
package contains identifiers and covariates, not labels.  Training code joins
labels only after this package has been sealed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import torch


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset" / "formal_training" / "chembl37_pki_formal.v4"
ROSTER = ROOT / "dataset" / "formal_training" / "a2s_d0r_roster.v3"
ESM = ROOT / "dataset" / "public" / "chembl_37" / "processed" / "dualcold" / "target_esm2.npz"
DEFAULT_OUTPUT = ROOT / "dataset" / "formal_training" / "a2s_cmal_episodes.v3"

SEED = 1729
SUPPORT_K = (1, 3, 5)
SOURCE_DRAWS = 24
SOURCE_QUERY_MAX = 32
SOURCE_QUERY_MIN = 5
PROTOCOLS = ("ordered", "random")
SPLIT_NAMES = ("meta_train", "meta_validation", "meta_test")
FORBIDDEN_COLUMNS = {
    "pKi", "pKd", "pchembl_value", "standard_value", "label", "affinity",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def json_list(values: Iterable[str]) -> str:
    return canonical_json([str(value) for value in values])


def episode_id(
    *, target: str, protocol: str, draw: int, k: int, support: list[str], query: list[str],
    support_measurements: list[str] | None = None,
    query_measurements: list[str] | None = None,
) -> str:
    payload = canonical_json({
        "draw": draw,
        "k": k,
        "protocol": protocol,
        "query": query,
        "query_measurements": query_measurements,
        "support": support,
        "support_measurements": support_measurements,
        "target": target,
    })
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def assign_source_splits(sources: pd.DataFrame, seed: int = SEED) -> pd.DataFrame:
    """Allocate intact homology components in the published 8:1:1 ratio."""

    required = {"target_uid", "component_id"}
    if not required.issubset(sources.columns):
        raise ValueError(f"source table lacks {sorted(required - set(sources.columns))}")
    components = np.array(sorted(sources.component_id.unique()))
    if len(components) < 10:
        raise ValueError("8:1:1 source split needs at least ten homology components")
    rng = np.random.default_rng(seed)
    rng.shuffle(components)
    n_train = max(1, int(np.floor(0.8 * len(components))))
    n_validation = max(1, int(np.floor(0.1 * len(components))))
    if n_train + n_validation >= len(components):
        n_validation = 1
        n_train = len(components) - 2
    component_split = {
        int(component): (
            "meta_train" if position < n_train else
            "meta_validation" if position < n_train + n_validation else
            "meta_test"
        )
        for position, component in enumerate(components)
    }
    out = sources[["target_uid", "component_id"]].copy()
    out["role"] = "source"
    out["meta_split"] = out.component_id.map(component_split)
    if out.meta_split.isna().any():
        raise RuntimeError("a source homology component was not assigned")
    return out.sort_values("target_uid").reset_index(drop=True)


def load_target_features(
    targets: list[str], sequence_path: Path, esm_path: Path
) -> tuple[np.ndarray, dict[str, Any]]:
    """Return ESM-2 pooled embeddings in target UID order with full coverage."""

    sequence_meta = json.loads(sequence_path.read_text(encoding="utf-8"))
    with np.load(esm_path, allow_pickle=True) as archive:
        keys = [str(value) for value in archive["keys"]]
        pooled = archive["pooled"].astype(np.float32)
        plm = str(archive["plm"])
    by_chembl = {key: row for row, key in enumerate(keys)}
    rows: list[np.ndarray] = []
    missing: list[str] = []
    for target in targets:
        record = sequence_meta.get(target)
        chembl = None if record is None else str(record.get("target_chembl_id", ""))
        if not chembl or chembl not in by_chembl:
            missing.append(target)
            continue
        rows.append(pooled[by_chembl[chembl]])
    if missing:
        raise ValueError(f"ESM-2 coverage missing for {len(missing)} targets: {missing[:5]}")
    matrix = np.stack(rows).astype(np.float32)
    return matrix, {"model": plm, "dimension": int(matrix.shape[1]), "coverage": len(rows)}


def _first_observation(frame: pd.DataFrame) -> pd.DataFrame:
    tie_breakers = [
        column for column in ("assay_context_uid", "measurement_uid")
        if column in frame.columns
    ]
    return (
        frame.sort_values(
            ["document_year", "document_uid", *tie_breakers, "compound_parent_uid"],
            kind="stable",
        )
        .drop_duplicates(["target_uid", "compound_parent_uid"], keep="first")
        .reset_index(drop=True)
    )


def _sample_source_draw(
    frame: pd.DataFrame,
    protocol: str,
    rng: np.random.Generator,
    draw: int,
    k_max: int,
    query_min: int,
    query_max: int,
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    first = _first_observation(frame)
    if len(first) < k_max + query_min:
        return None

    if protocol == "ordered":
        valid: list[tuple[pd.DataFrame, pd.DataFrame]] = []
        for year in sorted(first.document_year.unique())[:-1]:
            support_pool = first[first.document_year <= year]
            support_docs = set(support_pool.document_uid)
            query_pool = first[
                (first.document_year > year)
                & ~first.document_uid.isin(support_docs)
            ]
            if len(support_pool) >= k_max and len(query_pool) >= query_min:
                valid.append((support_pool, query_pool))
        if not valid:
            return None
        support_pool, query_pool = valid[draw % len(valid)]
        support = support_pool.iloc[
            rng.choice(len(support_pool), size=k_max, replace=False)
        ]
    elif protocol == "random":
        support = first.iloc[rng.choice(len(first), size=k_max, replace=False)]
        support_parents = set(support.compound_parent_uid)
        support_docs = set(support.document_uid)
        query_pool = first[
            ~first.compound_parent_uid.isin(support_parents)
            & ~first.document_uid.isin(support_docs)
        ]
        if len(query_pool) < query_min:
            return None
    else:
        raise ValueError(f"unknown source protocol: {protocol}")

    query_size = min(query_max, len(query_pool))
    query = query_pool.iloc[rng.choice(len(query_pool), size=query_size, replace=False)]
    support = support.reset_index(drop=True)
    query = query.reset_index(drop=True)
    return support, query


def build_source_episodes(
    metadata: pd.DataFrame,
    target_splits: pd.DataFrame,
    *,
    seed: int = SEED,
    draws: int = SOURCE_DRAWS,
    protocols: tuple[str, ...] = PROTOCOLS,
    query_min: int = SOURCE_QUERY_MIN,
    query_max: int = SOURCE_QUERY_MAX,
) -> pd.DataFrame:
    """Create nested, label-blind source episodes from publication metadata."""

    if set(metadata.columns) & FORBIDDEN_COLUMNS:
        raise AssertionError("outcome firewall violated while building source episodes")
    rng = np.random.default_rng(seed)
    split_by_target = dict(zip(target_splits.target_uid, target_splits.meta_split))
    component_by_target = dict(zip(target_splits.target_uid, target_splits.component_id))
    output: list[dict[str, Any]] = []
    source_rows = metadata[metadata.target_uid.isin(split_by_target)]

    for protocol in protocols:
        for target, frame in source_rows.groupby("target_uid", sort=True):
            for draw in range(draws):
                sampled = _sample_source_draw(
                    frame, protocol, rng, draw, max(SUPPORT_K), query_min, query_max
                )
                if sampled is None:
                    continue
                support_frame, query_frame = sampled
                support_all = [str(value) for value in support_frame.compound_parent_uid]
                query = [str(value) for value in query_frame.compound_parent_uid]
                support_measurements_all = [
                    str(value) for value in support_frame.measurement_uid
                ]
                query_measurements = [str(value) for value in query_frame.measurement_uid]
                support_docs_all = [str(value) for value in support_frame.document_uid]
                query_docs = [str(value) for value in query_frame.document_uid]
                for k in SUPPORT_K:
                    support = support_all[:k]
                    support_measurements = support_measurements_all[:k]
                    support_docs = support_docs_all[:k]
                    output.append({
                        "episode_id": episode_id(
                            target=target, protocol=protocol, draw=draw, k=k,
                            support=support, query=query,
                            support_measurements=support_measurements,
                            query_measurements=query_measurements,
                        ),
                        "role": "source",
                        "meta_split": split_by_target[target],
                        "protocol": protocol,
                        "target_uid": target,
                        "component_id": int(component_by_target[target]),
                        "draw_id": draw,
                        "k": k,
                        "support_parent_uids": json_list(support),
                        "query_parent_uids": json_list(query),
                        "support_measurement_uids": json_list(support_measurements),
                        "query_measurement_uids": json_list(query_measurements),
                        "support_document_uids": json_list(support_docs),
                        "query_document_uids": json_list(query_docs),
                        "support_max_year": int(support_frame.document_year.iloc[:k].max()),
                        "query_min_year": int(query_frame.document_year.min()),
                    })
    if not output:
        raise RuntimeError("no source episodes were constructed")
    return pd.DataFrame(output)


def build_recipient_episodes(
    recipients: pd.DataFrame,
    draws: pd.DataFrame,
    query: pd.DataFrame,
) -> pd.DataFrame:
    """Project the already sealed D0-R support/query roster into episode rows."""

    if (set(draws.columns) | set(query.columns)) & FORBIDDEN_COLUMNS:
        raise AssertionError("outcome firewall violated while building recipient episodes")
    required = {"measurement_uid", "document_uid", "document_year"}
    if not required.issubset(draws.columns) or not required.issubset(query.columns):
        raise ValueError("sealed recipient roster lacks frozen measurement identity")
    query_by_target = {
        target: frame.sort_values(["document_year", "document_uid", "measurement_uid"])
        for target, frame in query.groupby("target_uid", sort=True)
    }
    recipient_by_target = recipients.set_index("target_uid")
    output: list[dict[str, Any]] = []
    for (target, draw_id, k), frame in draws.groupby(
        ["target_uid", "draw_id", "k"], sort=True
    ):
        support = [str(value) for value in frame.compound_parent_uid]
        support_measurements = [str(value) for value in frame.measurement_uid]
        query_frame = query_by_target.get(target)
        if query_frame is None:
            queries: list[str] = []
            query_measurements: list[str] = []
        else:
            queries = [str(value) for value in query_frame.compound_parent_uid]
            query_measurements = [str(value) for value in query_frame.measurement_uid]
        if len(support) != int(k) or not queries:
            raise ValueError(f"invalid sealed recipient episode {target}/{draw_id}/k={k}")
        support_docs = [str(value) for value in frame.document_uid]
        query_docs = [str(value) for value in query_frame.document_uid]
        if set(support) & set(queries):
            raise ValueError("recipient support/query parent overlap")
        if set(support_docs) & set(query_docs):
            raise ValueError("recipient support/query document overlap")
        output.append({
            "episode_id": episode_id(
                target=target, protocol="d0r", draw=int(draw_id), k=int(k),
                support=support, query=queries,
                support_measurements=support_measurements,
                query_measurements=query_measurements,
            ),
            "role": "recipient",
            "meta_split": "recipient_test",
            "protocol": "d0r",
            "target_uid": target,
            "component_id": int(recipient_by_target.loc[target, "component_id"]),
            "draw_id": int(draw_id),
            "k": int(k),
            "support_parent_uids": json_list(support),
            "query_parent_uids": json_list(queries),
            "support_measurement_uids": json_list(support_measurements),
            "query_measurement_uids": json_list(query_measurements),
            "support_document_uids": json_list(support_docs),
            "query_document_uids": json_list(query_docs),
            "support_max_year": int(frame.document_year.max()),
            "query_min_year": int(query_frame.document_year.min()),
        })
    return pd.DataFrame(output)


def _candidate_for_target(group: pd.DataFrame, target: str, draw_id: int) -> int:
    candidates = group.index[group.target_uid == target].to_numpy()
    if len(candidates) == 0:
        raise RuntimeError(f"counterfactual target {target} has no candidate episode")
    positions = group.loc[candidates, "draw_id"].to_numpy()
    distance = np.abs(positions - draw_id)
    return int(candidates[np.lexsort((candidates, distance))[0]])


def attach_counterfactuals(
    episodes: pd.DataFrame,
    *,
    parent_uids: list[str],
    ecfp4: np.ndarray,
    scaffolds: dict[str, str],
    target_uids: list[str],
    target_features: np.ndarray,
    seed: int = SEED,
    device: str | None = None,
) -> pd.DataFrame:
    """Attach random, protein-hard and support-chemistry-matched wrong targets.

    Chemical matching uses only the positive support compounds. It first
    maximizes support-set scaffold Jaccard and then breaks ties with ECFP4
    centroid cosine. It never reads a query compound or any affinity label.
    """

    output = episodes.copy().reset_index(drop=True)
    parent_index = {value: row for row, value in enumerate(parent_uids)}
    target_index = {value: row for row, value in enumerate(target_uids)}
    pnorm = target_features.astype(np.float32)
    pnorm /= np.linalg.norm(pnorm, axis=1, keepdims=True).clip(min=1e-12)
    rng = np.random.default_rng(seed)
    chosen_device = device or ("cuda" if torch.cuda.is_available() else "cpu")

    for column in (
        "random_negative_episode_id",
        "protein_hard_negative_episode_id",
        "chemical_match_negative_episode_id",
    ):
        output[column] = ""
    output["protein_hard_cosine"] = np.nan
    output["chemical_match_cosine"] = np.nan
    output["chemical_match_scaffold_jaccard"] = np.nan

    group_columns = ["role", "meta_split", "protocol", "k"]
    for _, group in output.groupby(group_columns, sort=True):
        group = group.sort_values(["target_uid", "draw_id", "episode_id"])
        indices = group.index.to_numpy()
        if group.target_uid.nunique() < 2:
            raise RuntimeError("counterfactual group contains fewer than two targets")

        # Random wrong-target support.  A shuffled cyclic search guarantees that
        # every selected candidate has a different target, unlike a rolled batch.
        order = indices.copy()
        rng.shuffle(order)
        order_targets = output.loc[order, "target_uid"].to_numpy()
        position = {int(index): row for row, index in enumerate(order)}
        random_choice: dict[int, int] = {}
        for index in indices:
            start = position[int(index)]
            own = output.at[index, "target_uid"]
            for offset in range(1, len(order) + 1):
                candidate_pos = (start + offset) % len(order)
                if order_targets[candidate_pos] != own:
                    random_choice[int(index)] = int(order[candidate_pos])
                    break

        # Protein-hard wrong target.  Protein similarity is used for mining only;
        # the adapter receives the recipient protein for both positive and wrong
        # support, so it cannot solve the task by observing a target-ID mismatch.
        present_targets = sorted(group.target_uid.unique())
        present_rows = np.array([target_index[target] for target in present_targets])
        similarity = pnorm[present_rows] @ pnorm[present_rows].T
        np.fill_diagonal(similarity, -np.inf)
        hard_target = {
            target: present_targets[int(np.argmax(similarity[row]))]
            for row, target in enumerate(present_targets)
        }
        hard_score = {
            target: float(np.max(similarity[row]))
            for row, target in enumerate(present_targets)
        }

        # Exact scaffold-distribution matching with ECFP4 tie-breaking. Both
        # representations are computed from support compounds only.
        centroids = []
        support_scaffold_sets: list[list[str]] = []
        for value in group.support_parent_uids:
            parents = json.loads(value)
            rows = [parent_index[parent] for parent in parents]
            centroids.append(ecfp4[rows].astype(np.float32).mean(axis=0))
            support_scaffold_sets.append(sorted({
                scaffolds.get(parent, "INVALID")
                for parent in parents
                if scaffolds.get(parent, "INVALID") not in {"", "INVALID"}
            }))
        centroid = np.stack(centroids)
        centroid /= np.linalg.norm(centroid, axis=1, keepdims=True).clip(min=1e-12)
        matrix = torch.tensor(centroid, dtype=torch.float32, device=chosen_device)
        scaffold_vocab = {
            scaffold: row
            for row, scaffold in enumerate(sorted({
                scaffold
                for values in support_scaffold_sets
                for scaffold in values
            }))
        }
        scaffold_width = max(1, max(map(len, support_scaffold_sets)))
        scaffold_ids = np.full((len(group), scaffold_width), -1, dtype=np.int64)
        scaffold_sizes = np.zeros(len(group), dtype=np.int64)
        for row, values in enumerate(support_scaffold_sets):
            encoded = [scaffold_vocab[value] for value in values]
            scaffold_ids[row, :len(encoded)] = encoded
            scaffold_sizes[row] = len(encoded)
        scaffold_tensor = torch.tensor(
            scaffold_ids, dtype=torch.long, device=chosen_device
        )
        scaffold_size_tensor = torch.tensor(
            scaffold_sizes, dtype=torch.float32, device=chosen_device
        )
        group_targets = group.target_uid.to_numpy()
        chemical_local = np.empty(len(group), dtype=np.int64)
        chemical_score = np.empty(len(group), dtype=np.float32)
        chemical_jaccard = np.empty(len(group), dtype=np.float32)
        for start in range(0, len(group), 512):
            stop = min(start + 512, len(group))
            score = matrix[start:stop] @ matrix.T
            same = torch.tensor(
                group_targets[start:stop, None] == group_targets[None, :],
                dtype=torch.bool,
                device=chosen_device,
            )
            left = scaffold_tensor[start:stop, :, None, None]
            right = scaffold_tensor[None, None, :, :]
            valid = (left >= 0) & (right >= 0)
            intersection = ((left == right) & valid).any(dim=3).sum(dim=1).float()
            union = (
                scaffold_size_tensor[start:stop, None]
                + scaffold_size_tensor[None, :]
                - intersection
            )
            jaccard = torch.where(
                union > 0, intersection / union, torch.zeros_like(union)
            ).masked_fill(same, -torch.inf)
            best_jaccard = jaccard.max(dim=1).values
            scaffold_tie = jaccard == best_jaccard[:, None]
            score = score.masked_fill(same | ~scaffold_tie, -torch.inf)
            values, selected = score.max(dim=1)
            chemical_local[start:stop] = selected.cpu().numpy()
            chemical_score[start:stop] = values.cpu().numpy()
            chemical_jaccard[start:stop] = jaccard.gather(
                1, selected[:, None]
            ).squeeze(1).cpu().numpy()
        del matrix, scaffold_tensor, scaffold_size_tensor

        for local, index in enumerate(indices):
            own = output.at[index, "target_uid"]
            random_index = random_choice[int(index)]
            protein_index = _candidate_for_target(group, hard_target[own], int(output.at[index, "draw_id"]))
            chemical_index = int(indices[chemical_local[local]])
            output.at[index, "random_negative_episode_id"] = output.at[random_index, "episode_id"]
            output.at[index, "protein_hard_negative_episode_id"] = output.at[protein_index, "episode_id"]
            output.at[index, "chemical_match_negative_episode_id"] = output.at[chemical_index, "episode_id"]
            output.at[index, "protein_hard_cosine"] = hard_score[own]
            output.at[index, "chemical_match_cosine"] = float(chemical_score[local])
            output.at[index, "chemical_match_scaffold_jaccard"] = float(
                chemical_jaccard[local]
            )

    id_target = dict(zip(output.episode_id, output.target_uid))
    for column in (
        "random_negative_episode_id",
        "protein_hard_negative_episode_id",
        "chemical_match_negative_episode_id",
    ):
        if output[column].eq("").any():
            raise RuntimeError(f"unassigned counterfactual in {column}")
        wrong = output[column].map(id_target)
        if (wrong == output.target_uid).any():
            raise RuntimeError(f"same-target counterfactual in {column}")
    return output


def audit_package(episodes: pd.DataFrame, target_splits: pd.DataFrame) -> dict[str, Any]:
    id_target = dict(zip(episodes.episode_id, episodes.target_uid))
    negative_columns = [
        "random_negative_episode_id",
        "protein_hard_negative_episode_id",
        "chemical_match_negative_episode_id",
    ]
    negative_violations = {
        column: int((episodes[column].map(id_target) == episodes.target_uid).sum())
        for column in negative_columns
    }
    parent_overlap = 0
    document_overlap = 0
    measurement_overlap = 0
    ordered_year_violations = 0
    for row in episodes.itertuples():
        support = set(json.loads(row.support_parent_uids))
        query = set(json.loads(row.query_parent_uids))
        support_measurements = json.loads(row.support_measurement_uids)
        query_measurements = json.loads(row.query_measurement_uids)
        support_docs = set(json.loads(row.support_document_uids))
        query_docs = set(json.loads(row.query_document_uids))
        parent_overlap += int(bool(support & query))
        document_overlap += int(bool(support_docs & query_docs))
        measurement_overlap += int(bool(set(support_measurements) & set(query_measurements)))
        if len(support_measurements) != row.k:
            raise ValueError(f"support measurement alignment failed for {row.episode_id}")
        if len(query_measurements) != len(query):
            raise ValueError(f"query measurement alignment failed for {row.episode_id}")
        if set(support_measurements) & set(query_measurements):
            raise ValueError(f"support/query measurement overlap for {row.episode_id}")
        if row.protocol in {"ordered", "d0r"}:
            ordered_year_violations += int(row.support_max_year >= row.query_min_year)

    nested_violations = 0
    source = episodes[episodes.role == "source"]
    for _, frame in source.groupby(["protocol", "target_uid", "draw_id"], sort=False):
        support = {
            int(row.k): set(json.loads(row.support_parent_uids))
            for row in frame.itertuples()
        }
        if set(SUPPORT_K).issubset(support):
            nested_violations += int(not (support[1] <= support[3] <= support[5]))

    component_overlap = {}
    source_splits = target_splits[target_splits.role == "source"]
    for left_pos, left in enumerate(SPLIT_NAMES):
        for right in SPLIT_NAMES[left_pos + 1:]:
            a = set(source_splits.loc[source_splits.meta_split == left, "component_id"])
            b = set(source_splits.loc[source_splits.meta_split == right, "component_id"])
            component_overlap[f"{left}__{right}"] = len(a & b)

    return {
        "label_columns_read": [],
        "label_blind": True,
        "frozen_label_key": "measurement_uid",
        "episode_counts": {
            "total": int(len(episodes)),
            "by_role": {str(key): int(value) for key, value in episodes.role.value_counts().items()},
            "by_split": {str(key): int(value) for key, value in episodes.meta_split.value_counts().items()},
            "by_protocol": {str(key): int(value) for key, value in episodes.protocol.value_counts().items()},
            "by_k": {str(key): int(value) for key, value in episodes.k.value_counts().sort_index().items()},
        },
        "target_counts": {
            "total": int(target_splits.target_uid.nunique()),
            "by_split": {
                str(key): int(value)
                for key, value in target_splits.groupby("meta_split").target_uid.nunique().items()
            },
        },
        "violations": {
            "support_query_parent": parent_overlap,
            "support_query_document": document_overlap,
            "support_query_measurement": measurement_overlap,
            "ordered_time": ordered_year_violations,
            "nested_support": nested_violations,
            "negative_same_target": negative_violations,
            "source_component_overlap": component_overlap,
        },
        "counterfactuals": {
            "random": "different target sampled within role/split/protocol/k",
            "protein_hard": "highest ESM-2 cosine different target; mining only",
            "chemical_match": "highest support-scaffold Jaccard, ECFP4 support-centroid cosine tie-break, different target; query-blind",
        },
    }


def build(output: Path = DEFAULT_OUTPUT, *, device: str | None = None) -> dict[str, Any]:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite model-ready package: {output}")

    source_manifest = json.loads((CORPUS / "manifest.json").read_text(encoding="utf-8"))
    roster_manifest = json.loads((ROSTER / "manifest.json").read_text(encoding="utf-8"))
    if source_manifest.get("status") != "FORMAL_PKI_CORPUS_READY_NATURAL_TAIL_BLOCKED":
        raise ValueError("unexpected formal corpus status")
    if roster_manifest.get("status") != "PASS":
        raise ValueError("D0-R roster has not passed")

    sources = pd.read_parquet(ROSTER / "sources.parquet")
    recipients = pd.read_parquet(ROSTER / "recipients.parquet")
    target_splits = assign_source_splits(sources)
    recipient_splits = recipients[["target_uid", "component_id", "homology_warm"]].copy()
    recipient_splits["role"] = "recipient"
    recipient_splits["meta_split"] = "recipient_test"
    target_splits["homology_warm"] = False
    target_splits = pd.concat([target_splits, recipient_splits], ignore_index=True)
    target_splits = target_splits[
        ["target_uid", "component_id", "role", "meta_split", "homology_warm"]
    ].sort_values(["role", "target_uid"]).reset_index(drop=True)

    closed = pd.read_parquet(
        ROSTER / "source_rows.parquet",
        columns=[
            "target_uid", "compound_parent_uid", "measurement_uid",
            "document_uid", "document_year", "assay_context_uid",
        ],
    ).drop_duplicates()
    if set(closed.columns) & FORBIDDEN_COLUMNS:
        raise AssertionError("outcome firewall violated")
    source_metadata = closed.copy()
    source_episodes = build_source_episodes(source_metadata, target_splits)
    recipient_episodes = build_recipient_episodes(
        recipients,
        pd.read_parquet(ROSTER / "support_draws.parquet"),
        pd.read_parquet(ROSTER / "query.parquet"),
    )
    episodes = pd.concat([source_episodes, recipient_episodes], ignore_index=True)
    episodes = episodes.sort_values(
        ["role", "meta_split", "protocol", "target_uid", "draw_id", "k"]
    ).reset_index(drop=True)

    targets = list(target_splits.target_uid)
    target_features, feature_info = load_target_features(
        targets, CORPUS / "features" / "target_sequences.json", ESM
    )
    with np.load(CORPUS / "features" / "ligand_features.npz", allow_pickle=True) as ligand:
        parent_uids = [str(value) for value in ligand["parent_uids"]]
        ecfp4 = ligand["ecfp4"]
    compounds = pd.read_parquet(
        CORPUS / "components" / "compounds.parquet",
        columns=["compound_parent_uid", "parent_canonical_smiles"],
    )
    from research.a2s.a2s_d0r import murcko_scaffolds
    scaffolds = murcko_scaffolds(dict(zip(
        compounds.compound_parent_uid, compounds.parent_canonical_smiles
    )))
    episodes = attach_counterfactuals(
        episodes,
        parent_uids=parent_uids,
        ecfp4=ecfp4,
        scaffolds=scaffolds,
        target_uids=targets,
        target_features=target_features,
        device=device,
    )
    audit = audit_package(episodes, target_splits)
    violations = audit["violations"]
    flat_violations = [
        violations["support_query_parent"],
        violations["support_query_document"],
        violations["support_query_measurement"],
        violations["ordered_time"],
        violations["nested_support"],
        *violations["negative_same_target"].values(),
        *violations["source_component_overlap"].values(),
    ]
    if any(flat_violations):
        raise RuntimeError(f"model-ready package audit failed: {violations}")

    output.mkdir(parents=True)
    episodes.to_parquet(output / "episodes.parquet", index=False)
    target_splits.to_parquet(output / "target_splits.parquet", index=False)
    np.savez_compressed(
        output / "target_features.npz",
        target_uids=np.asarray(targets),
        pooled=target_features,
        model=np.asarray(feature_info["model"]),
    )
    (output / "audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True), encoding="utf-8"
    )
    (output / "README.md").write_text(
        "# A2S-CMAL model-ready episodes v3\n\n"
        "This immutable package is label-blind. It contains source meta-task splits, "
        "nested support/query parent and measurement identifiers, frozen ESM-2 target "
        "features, and three "
        "wrong-target support mappings. Affinity labels remain in the formal ChEMBL "
        "corpus and are joined only by the external trainer. `recipient_test` is the "
        "sealed D0-R target-disjoint evaluation role; local runs are smoke tests only.\n",
        encoding="utf-8",
    )

    inputs = [
        CORPUS / "manifest.json",
        CORPUS / "canonical" / "pki_measurements_context_main.parquet",
        CORPUS / "features" / "ligand_features.npz",
        CORPUS / "features" / "target_sequences.json",
        ROSTER / "manifest.json",
        ROSTER / "sources.parquet",
        ROSTER / "recipients.parquet",
        ROSTER / "source_rows.parquet",
        ROSTER / "support_draws.parquet",
        ROSTER / "query.parquet",
        ESM,
    ]
    output_files = [
        output / "episodes.parquet",
        output / "target_splits.parquet",
        output / "target_features.npz",
        output / "audit.json",
        output / "README.md",
    ]
    manifest: dict[str, Any] = {
        "schema": "a2s-cmal-model-ready-v3",
        "status": "READY_FOR_EXTERNAL_FORMAL_TRAINING",
        "seed": SEED,
        "source_package": str(CORPUS.relative_to(ROOT)).replace("\\", "/"),
        "recipient_roster": str(ROSTER.relative_to(ROOT)).replace("\\", "/"),
        "endpoint": "pKi",
        "support_budgets": list(SUPPORT_K),
        "source_protocols": list(PROTOCOLS),
        "source_draws_requested_per_target": SOURCE_DRAWS,
        "source_query_bounds": [SOURCE_QUERY_MIN, SOURCE_QUERY_MAX],
        "task_split": {
            "method": "intact 40%-identity homology components",
            "ratio": [8, 1, 1],
            "paper_protocol_source": "AdaMBind, Nature Communications 2026",
        },
        "protein_features": feature_info,
        "label_firewall": {
            "labels_in_package": False,
            "label_columns_read_during_build": [],
            "counterfactual_mining_query_blind": True,
            "frozen_external_label_key": "measurement_uid",
            "label_table": "canonical/pki_measurements_context_main.parquet",
        },
        "audit": audit,
        "inputs": {
            str(path.relative_to(ROOT)).replace("\\", "/"): {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in inputs
        },
        "files": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in output_files
        },
    }
    manifest["content_sha256"] = hashlib.sha256(
        canonical_json(manifest).encode("utf-8")
    ).hexdigest()
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="prepare label-blind A2S-CMAL episodes")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--device", choices=("cpu", "cuda"), default=None)
    args = parser.parse_args()
    result = build(args.output, device=args.device)
    print(json.dumps({
        "status": result["status"],
        "output": str(args.output),
        "episodes": result["audit"]["episode_counts"],
        "targets": result["audit"]["target_counts"],
        "content_sha256": result["content_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()
