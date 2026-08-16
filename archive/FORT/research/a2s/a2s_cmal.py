"""Protein-conditioned Counterfactual Meta-Adaptation Learning for A2S-DTA.

The final predictor has two deliberately separated parts:

1. a support-free DTA base ``f_theta(protein, ligand)``; and
2. a learned, non-closed-form operator ``A_phi(protein, support, query)``.

The base is trained first and frozen.  The adapter then learns from abundant
source-target episodes.  For every positive episode it sees three target-
mismatched support sets (random target, protein-hard target and support-
chemistry-matched target), while the recipient protein and query remain fixed.
The contrastive score is post-adaptation query ranking loss, never embedding
similarity.  Thus the objective directly rewards a support set only when it
causes the right query-specific ranking change for that protein.

Formal multi-seed training is intentionally guarded for an external CUDA host.
Local execution is limited to ``--smoke`` and never reads recipient labels.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

from research.a2s.a2s_bir import (
    METRIC_NAMES,
    PooledPrior,
    component_bootstrap,
    episode_metrics,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset" / "formal_training" / "chembl37_pki_formal.v4"
ROSTER = ROOT / "dataset" / "formal_training" / "a2s_d0r_roster.v3"
EPISODES = ROOT / "dataset" / "formal_training" / "a2s_cmal_episodes.v3"
DEFAULT_OUT = ROOT / "reports" / "active"

SUPPORT_K = (1, 3, 5)
NEGATIVE_COLUMNS = (
    "random_negative_episode_id",
    "protein_hard_negative_episode_id",
    "chemical_match_negative_episode_id",
)
EPISODE_NEGATIVE_NAMES = ("random", "protein_hard", "chemical_match")
NEGATIVE_NAMES = (*EPISODE_NEGATIVE_NAMES, "label_swap")
RANK_WEIGHT = 0.3
COUNTERFACTUAL_WEIGHT = 1.0
COUNTERFACTUAL_TEMPERATURE = 0.2
N_BOOTSTRAP = 5000
FORMAL_SEEDS = (1729, 1731, 1733, 1741, 1753)
SOURCE_PARTITION_SEED = 1729


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_model_ready_package(path: Path) -> dict[str, Any]:
    manifest_path = path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "READY_FOR_EXTERNAL_FORMAL_TRAINING":
        raise ValueError("CMAL episode package is not training-ready")
    if manifest.get("label_firewall", {}).get("labels_in_package") is not False:
        raise ValueError("CMAL package label firewall is not closed")
    for name, record in manifest["files"].items():
        file = path / name
        if not file.is_file() or file.stat().st_size != record["bytes"]:
            raise ValueError(f"CMAL package file size mismatch: {file}")
        if sha256_file(file) != record["sha256"]:
            raise ValueError(f"CMAL package file hash mismatch: {file}")
    violations = manifest["audit"]["violations"]
    scalar = [
        violations["support_query_parent"],
        violations["support_query_document"],
        violations["support_query_measurement"],
        violations["ordered_time"],
        violations["nested_support"],
        *violations["negative_same_target"].values(),
        *violations["source_component_overlap"].values(),
    ]
    if any(scalar):
        raise ValueError("CMAL episode package contains a leakage violation")
    return manifest


@dataclass
class CmalData:
    ligand_features: np.ndarray
    parent_uids: list[str]
    parent_index: dict[str, int]
    target_features: np.ndarray
    target_uids: list[str]
    target_index: dict[str, int]
    target_splits: pd.DataFrame
    episodes: pd.DataFrame
    labels: dict[str, float]
    components: dict[str, int]
    manifest: dict[str, Any]


def verify_episode_measurement_identity(
    episodes: pd.DataFrame, observations: pd.DataFrame
) -> dict[str, int]:
    required = {"target_uid", "compound_parent_uid", "measurement_uid"}
    if not required.issubset(observations.columns):
        raise ValueError("measurement identity table is incomplete")
    if observations.measurement_uid.duplicated().any():
        raise ValueError("measurement_uid is not unique in the admitted label table")
    identity = {
        str(row.measurement_uid): (str(row.target_uid), str(row.compound_parent_uid))
        for row in observations.itertuples()
    }
    used: set[str] = set()
    occurrences = 0
    for row in episodes.itertuples():
        for parent_column, measurement_column in (
            ("support_parent_uids", "support_measurement_uids"),
            ("query_parent_uids", "query_measurement_uids"),
        ):
            parents = json.loads(getattr(row, parent_column))
            measurements = json.loads(getattr(row, measurement_column))
            if len(parents) != len(measurements):
                raise ValueError(f"measurement alignment failed for {row.episode_id}")
            for parent, measurement in zip(parents, measurements):
                if identity.get(str(measurement)) != (str(row.target_uid), str(parent)):
                    raise ValueError(
                        f"measurement identity mismatch in episode {row.episode_id}"
                    )
                used.add(str(measurement))
                occurrences += 1
    return {"unique_measurements": len(used), "episode_occurrences": occurrences}


def _read_measurement_labels(targets: list[str]) -> pd.DataFrame:
    """Read only requested target labels from the admitted context table."""

    columns = ["target_uid", "measurement_uid", "pKi"]
    return pd.read_parquet(
        CORPUS / "canonical" / "pki_measurements_context_main.parquet",
        columns=columns,
        filters=[("target_uid", "in", targets)],
    )


def load_data(
    episode_root: Path = EPISODES,
    *,
    include_recipient_labels: bool = False,
) -> CmalData:
    manifest = verify_model_ready_package(episode_root)
    episodes = pd.read_parquet(episode_root / "episodes.parquet")
    target_splits = pd.read_parquet(episode_root / "target_splits.parquet")
    observation_identity = pd.read_parquet(
        CORPUS / "canonical" / "pki_measurements_context_main.parquet",
        columns=["target_uid", "compound_parent_uid", "measurement_uid"],
        filters=[("target_uid", "in", list(target_splits.target_uid))],
    )
    verify_episode_measurement_identity(episodes, observation_identity)
    with np.load(episode_root / "target_features.npz", allow_pickle=False) as target_archive:
        target_uids = [str(value) for value in target_archive["target_uids"]]
        target_features = target_archive["pooled"].astype(np.float32)
    with np.load(CORPUS / "features" / "ligand_features.npz", allow_pickle=True) as ligand:
        parent_uids = [str(value) for value in ligand["parent_uids"]]
        ecfp = ligand["ecfp4"].astype(np.float32)
        descriptors = np.nan_to_num(
            ligand["descriptors"].astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0
        )
    ligand_features = np.hstack([ecfp, descriptors]).astype(np.float32)

    source_targets = list(
        target_splits.loc[target_splits.role == "source", "target_uid"]
    )
    source_values = _read_measurement_labels(source_targets)
    closed_measurements = set(pd.read_parquet(
        ROSTER / "source_rows.parquet",
        columns=["measurement_uid"],
    ).measurement_uid.astype(str))
    source_values = source_values[
        source_values.measurement_uid.astype(str).isin(closed_measurements)
    ]
    if source_values.measurement_uid.duplicated().any():
        raise ValueError("source measurement_uid is not unique in context-main labels")
    labels = {
        str(row.measurement_uid): float(row.pKi)
        for row in source_values.itertuples()
    }

    target_index = {target: row for row, target in enumerate(target_uids)}
    parent_index = {parent: row for row, parent in enumerate(parent_uids)}
    if set(target_splits.target_uid) != set(target_uids):
        raise ValueError("target split and target feature identities differ")
    components = dict(zip(target_splits.target_uid, target_splits.component_id))
    data = CmalData(
        ligand_features=ligand_features,
        parent_uids=parent_uids,
        parent_index=parent_index,
        target_features=target_features,
        target_uids=target_uids,
        target_index=target_index,
        target_splits=target_splits,
        episodes=episodes,
        labels=labels,
        components=components,
        manifest=manifest,
    )
    if include_recipient_labels:
        add_recipient_labels(data)
    return data


def add_recipient_labels(data: CmalData) -> int:
    """Open the recipient-label firewall only after the source mechanism gate."""

    recipient_targets = list(
        data.target_splits.loc[data.target_splits.role == "recipient", "target_uid"]
    )
    recipient_values = _read_measurement_labels(recipient_targets)
    recipient_measurements: set[str] = set()
    recipient_episodes = data.episodes[data.episodes.role == "recipient"]
    for row in recipient_episodes.itertuples():
        recipient_measurements.update(json.loads(row.support_measurement_uids))
        recipient_measurements.update(json.loads(row.query_measurement_uids))
    recipient_labels = recipient_values[
        recipient_values.measurement_uid.astype(str).isin(recipient_measurements)
    ]
    if recipient_labels.measurement_uid.duplicated().any():
        raise ValueError("recipient measurement_uid is not unique in context-main labels")
    found = set(recipient_labels.measurement_uid.astype(str))
    missing = recipient_measurements - found
    if missing:
        raise KeyError(f"recipient labels missing for {len(missing)} measurements")
    data.labels.update({
        str(row.measurement_uid): float(row.pKi)
        for row in recipient_labels.itertuples()
    })
    return len(recipient_labels)


def _episode_observations(frame: pd.DataFrame) -> list[tuple[str, str]]:
    observations: set[tuple[str, str]] = set()
    for row in frame.itertuples():
        observations.update(zip(
            json.loads(row.support_parent_uids),
            json.loads(row.support_measurement_uids),
        ))
        observations.update(zip(
            json.loads(row.query_parent_uids),
            json.loads(row.query_measurement_uids),
        ))
    return sorted(observations)


def fit_source_prior(
    data: CmalData, targets: set[str] | None = None
) -> PooledPrior:
    train = data.episodes[
        (data.episodes.meta_split == "meta_train")
        & (data.episodes.protocol == "ordered")
        & (data.episodes.k == max(SUPPORT_K))
    ]
    if targets is not None:
        train = train[train.target_uid.isin(targets)]
    if train.empty:
        raise ValueError("source prior has no eligible base-pretraining episodes")
    observations = _episode_observations(train)
    missing = [measurement for _, measurement in observations if measurement not in data.labels]
    if missing:
        raise KeyError(f"closed source labels missing for {len(missing)} measurements")
    parent_rows = np.array([data.parent_index[parent] for parent, _ in observations])
    measurement_to_target: dict[str, str] = {}
    for row in train.itertuples():
        for measurement in json.loads(row.support_measurement_uids):
            measurement_to_target[measurement] = row.target_uid
        for measurement in json.loads(row.query_measurement_uids):
            measurement_to_target[measurement] = row.target_uid
    target_ids = np.array([measurement_to_target[measurement] for _, measurement in observations])
    y = np.array([data.labels[measurement] for _, measurement in observations], dtype=np.float32)
    return PooledPrior().fit(data.ligand_features[parent_rows], y, target_ids)


def component_disjoint_source_partition(
    target_splits: pd.DataFrame,
    seed: int = SOURCE_PARTITION_SEED,
) -> tuple[set[str], set[str], dict[str, Any]]:
    """Freeze disjoint base-pretraining and meta-adapter source components."""

    train = target_splits[
        (target_splits.role == "source")
        & (target_splits.meta_split == "meta_train")
    ]
    components = np.asarray(sorted(train.component_id.unique()), dtype=np.int64)
    if len(components) < 2:
        raise ValueError("component-disjoint source training needs two components")
    rng = np.random.default_rng(seed)
    rng.shuffle(components)
    split = (len(components) + 1) // 2
    base_components = set(int(value) for value in components[:split])
    adapter_components = set(int(value) for value in components[split:])
    base_targets = set(
        str(value)
        for value in train.loc[
            train.component_id.isin(base_components), "target_uid"
        ]
    )
    adapter_targets = set(
        str(value)
        for value in train.loc[
            train.component_id.isin(adapter_components), "target_uid"
        ]
    )
    if not base_targets or not adapter_targets or base_targets & adapter_targets:
        raise RuntimeError("source training partition is empty or target-overlapping")
    return base_targets, adapter_targets, {
        "seed": seed,
        "method": "fixed 1:1 split of intact meta-train homology components",
        "base_pretrain_components": len(base_components),
        "meta_adapter_components": len(adapter_components),
        "base_pretrain_targets": len(base_targets),
        "meta_adapter_targets": len(adapter_targets),
        "target_overlap": len(base_targets & adapter_targets),
        "component_overlap": len(base_components & adapter_components),
    }


@dataclass
class EpisodeTensorGroup:
    target: torch.Tensor
    support_parent: torch.Tensor
    support_y: torch.Tensor
    query_parent: torch.Tensor
    query_y: torch.Tensor
    query_mask: torch.Tensor
    negative_support_parent: torch.Tensor
    negative_support_y: torch.Tensor
    target_uids: list[str]
    draw_ids: np.ndarray
    episode_ids: list[str]
    sampling_weight: torch.Tensor

    @property
    def size(self) -> int:
        return int(self.target.shape[0])


class GpuEpisodeStore:
    """All model inputs and episode tensors resident on one device."""

    def __init__(self, data: CmalData, prior: PooledPrior, device: torch.device) -> None:
        self.data = data
        self.device = device
        mean = torch.tensor(prior.mean, dtype=torch.float32, device=device)
        scale = torch.tensor(prior.scale, dtype=torch.float32, device=device)
        ligand = torch.tensor(data.ligand_features, dtype=torch.float32, device=device)
        self.ligand = (ligand - mean) / scale
        beta = prior.beta.to(device=device, dtype=torch.float32)
        self.f0 = self.ligand @ beta + float(prior.intercept)
        self.target = torch.tensor(data.target_features, dtype=torch.float32, device=device)
        self.ligand_dim = int(self.ligand.shape[1])
        self.protein_dim = int(self.target.shape[1])
        self._episode_by_id = data.episodes.set_index("episode_id", drop=False)

    def materialize(
        self,
        *,
        meta_split: str,
        protocol: str,
        k: int,
        targets: set[str] | None = None,
    ) -> EpisodeTensorGroup:
        frame = self.data.episodes[
            (self.data.episodes.meta_split == meta_split)
            & (self.data.episodes.protocol == protocol)
            & (self.data.episodes.k == k)
        ].sort_values(["target_uid", "draw_id", "episode_id"])
        if targets is not None:
            frame = frame[frame.target_uid.isin(targets)]
        if frame.empty:
            raise ValueError(f"no episodes for {meta_split}/{protocol}/k={k}")
        max_query = max(len(json.loads(value)) for value in frame.query_parent_uids)

        target_rows = []
        support_rows = []
        support_labels = []
        query_rows = []
        query_labels = []
        query_masks = []
        negative_rows = []
        negative_labels = []
        targets: list[str] = []
        draws: list[int] = []
        ids: list[str] = []

        for row in frame.itertuples():
            support = json.loads(row.support_parent_uids)
            support_measurements = json.loads(row.support_measurement_uids)
            query = json.loads(row.query_parent_uids)
            query_measurements = json.loads(row.query_measurement_uids)
            if len(support) != k:
                raise ValueError(f"episode {row.episode_id} support size mismatch")
            if len(support_measurements) != len(support) or len(query_measurements) != len(query):
                raise ValueError(f"episode {row.episode_id} measurement alignment mismatch")
            missing = [
                measurement for measurement in support_measurements + query_measurements
                if measurement not in self.data.labels
            ]
            if missing:
                raise KeyError(f"labels unavailable for {missing[:3]}")
            q_parent = [self.data.parent_index[parent] for parent in query]
            q_y = [self.data.labels[measurement] for measurement in query_measurements]
            pad_parent = q_parent + [q_parent[0]] * (max_query - len(q_parent))
            pad_y = q_y + [q_y[0]] * (max_query - len(q_y))
            mask = [1.0] * len(q_parent) + [0.0] * (max_query - len(q_parent))

            row_negative_parent = []
            row_negative_y = []
            for column in NEGATIVE_COLUMNS:
                negative = self._episode_by_id.loc[getattr(row, column)]
                if negative.target_uid == row.target_uid:
                    raise ValueError("same-target counterfactual reached the trainer")
                parents = json.loads(negative.support_parent_uids)
                measurements = json.loads(negative.support_measurement_uids)
                if len(parents) != k:
                    raise ValueError("counterfactual support budget mismatch")
                if len(measurements) != len(parents):
                    raise ValueError("counterfactual measurement alignment mismatch")
                if any(measurement not in self.data.labels for measurement in measurements):
                    raise KeyError(f"counterfactual labels unavailable for {negative.episode_id}")
                row_negative_parent.append([self.data.parent_index[parent] for parent in parents])
                row_negative_y.append([self.data.labels[measurement] for measurement in measurements])

            target_rows.append(self.data.target_index[row.target_uid])
            support_rows.append([self.data.parent_index[parent] for parent in support])
            support_labels.append([self.data.labels[measurement] for measurement in support_measurements])
            query_rows.append(pad_parent)
            query_labels.append(pad_y)
            query_masks.append(mask)
            negative_rows.append(row_negative_parent)
            negative_labels.append(row_negative_y)
            targets.append(row.target_uid)
            draws.append(int(row.draw_id))
            ids.append(row.episode_id)

        tensor = lambda value, dtype: torch.tensor(value, dtype=dtype, device=self.device)
        target_array = np.asarray(target_rows, dtype=np.int64)
        _, inverse, counts = np.unique(
            target_array, return_inverse=True, return_counts=True
        )
        sampling_weight = 1.0 / counts[inverse].astype(np.float32)
        return EpisodeTensorGroup(
            target=tensor(target_rows, torch.long),
            support_parent=tensor(support_rows, torch.long),
            support_y=tensor(support_labels, torch.float32),
            query_parent=tensor(query_rows, torch.long),
            query_y=tensor(query_labels, torch.float32),
            query_mask=tensor(query_masks, torch.float32),
            negative_support_parent=tensor(negative_rows, torch.long),
            negative_support_y=tensor(negative_labels, torch.float32),
            target_uids=targets,
            draw_ids=np.asarray(draws, dtype=np.int64),
            episode_ids=ids,
            sampling_weight=tensor(sampling_weight, torch.float32),
        )


@dataclass
class AdaptationState:
    protein: torch.Tensor
    support_tokens: torch.Tensor
    support_residual: torch.Tensor


class ProteinConditionedMetaAdapter(nn.Module):
    """Learned support-to-state-to-query operator; contains no analytic solve."""

    def __init__(
        self,
        ligand_dim: int,
        protein_dim: int,
        hidden: int = 128,
        heads: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if hidden % heads:
            raise ValueError("hidden width must be divisible by attention heads")
        self.ligand_encoder = nn.Sequential(
            nn.Linear(ligand_dim, 2 * hidden),
            nn.GELU(),
            nn.LayerNorm(2 * hidden),
            nn.Linear(2 * hidden, hidden),
            nn.GELU(),
            nn.LayerNorm(hidden),
        )
        self.protein_encoder = nn.Sequential(
            nn.LayerNorm(protein_dim),
            nn.Linear(protein_dim, 2 * hidden),
            nn.GELU(),
            nn.Linear(2 * hidden, hidden),
            nn.LayerNorm(hidden),
        )
        self.pair_encoder = nn.Sequential(
            nn.Linear(4 * hidden, 2 * hidden),
            nn.GELU(),
            nn.LayerNorm(2 * hidden),
            nn.Linear(2 * hidden, hidden),
            nn.GELU(),
            nn.LayerNorm(hidden),
        )
        self.base_head = nn.Sequential(
            nn.Linear(hidden, hidden), nn.GELU(), nn.Linear(hidden, 1)
        )
        self.residual_encoder = nn.Sequential(
            nn.Linear(1, hidden), nn.Tanh(), nn.Linear(hidden, hidden)
        )
        self.support_attention = nn.MultiheadAttention(
            hidden, heads, dropout=dropout, batch_first=True
        )
        self.support_norm = nn.LayerNorm(hidden)
        self.query_attention = nn.MultiheadAttention(
            hidden, heads, dropout=dropout, batch_first=True
        )
        self.query_norm = nn.LayerNorm(hidden)
        self.delta_head = nn.Sequential(
            nn.Linear(4 * hidden, 2 * hidden),
            nn.GELU(),
            nn.LayerNorm(2 * hidden),
            nn.Linear(2 * hidden, hidden),
            nn.GELU(),
            nn.Linear(hidden, 1),
        )
        nn.init.zeros_(self.base_head[-1].weight)
        nn.init.zeros_(self.base_head[-1].bias)
        nn.init.zeros_(self.delta_head[-1].weight)
        nn.init.zeros_(self.delta_head[-1].bias)

    def base_modules(self) -> tuple[nn.Module, ...]:
        return (
            self.ligand_encoder,
            self.protein_encoder,
            self.pair_encoder,
            self.base_head,
        )

    def adapter_modules(self) -> tuple[nn.Module, ...]:
        return (
            self.residual_encoder,
            self.support_attention,
            self.support_norm,
            self.query_attention,
            self.query_norm,
            self.delta_head,
        )

    @staticmethod
    def _collect_parameters(modules: Iterable[nn.Module]) -> list[nn.Parameter]:
        return [parameter for module in modules for parameter in module.parameters()]

    def base_parameters(self) -> list[nn.Parameter]:
        return self._collect_parameters(self.base_modules())

    def adapter_parameters(self) -> list[nn.Parameter]:
        return self._collect_parameters(self.adapter_modules())

    def freeze_base(self) -> None:
        for module in self.base_modules():
            module.requires_grad_(False)
            module.eval()

    def train(self, mode: bool = True) -> "ProteinConditionedMetaAdapter":
        super().train(mode)
        if not any(parameter.requires_grad for parameter in self.base_parameters()):
            for module in self.base_modules():
                module.eval()
        return self

    def encode_pair(
        self,
        protein: torch.Tensor,
        ligand: torch.Tensor,
        *,
        use_protein: bool = True,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        if protein.ndim != 2 or ligand.ndim != 3:
            raise ValueError("protein [B,P] and ligand [B,N,D] tensors are required")
        b, n, _ = ligand.shape
        ligand_hidden = self.ligand_encoder(ligand.reshape(b * n, -1)).reshape(b, n, -1)
        protein_hidden = self.protein_encoder(protein)
        if not use_protein:
            protein_hidden = torch.zeros_like(protein_hidden)
        expanded = protein_hidden.unsqueeze(1).expand(-1, n, -1)
        pair = torch.cat(
            [ligand_hidden, expanded, ligand_hidden * expanded, (ligand_hidden - expanded).abs()],
            dim=-1,
        )
        return self.pair_encoder(pair), protein_hidden

    def base_predict(
        self,
        protein: torch.Tensor,
        ligand: torch.Tensor,
        f0: torch.Tensor,
        *,
        use_protein: bool = True,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        pair, protein_hidden = self.encode_pair(protein, ligand, use_protein=use_protein)
        return f0 + self.base_head(pair).squeeze(-1), pair, protein_hidden

    def adapt(
        self,
        protein: torch.Tensor,
        support_ligand: torch.Tensor,
        support_y: torch.Tensor,
        support_f0: torch.Tensor,
        *,
        use_protein: bool = True,
    ) -> AdaptationState:
        base, pair, protein_hidden = self.base_predict(
            protein, support_ligand, support_f0, use_protein=use_protein
        )
        residual = support_y - base
        tokens = pair + self.residual_encoder(residual.unsqueeze(-1))
        attended, _ = self.support_attention(tokens, tokens, tokens, need_weights=False)
        tokens = self.support_norm(tokens + attended)
        return AdaptationState(
            protein=protein_hidden,
            support_tokens=tokens,
            support_residual=residual,
        )

    def predict(
        self,
        state: AdaptationState,
        protein: torch.Tensor,
        query_ligand: torch.Tensor,
        query_f0: torch.Tensor,
        *,
        use_protein: bool = True,
        use_support: bool = True,
    ) -> dict[str, torch.Tensor]:
        base, query_pair, _ = self.base_predict(
            protein, query_ligand, query_f0, use_protein=use_protein
        )
        if not use_support:
            return {"prediction": base, "base_prediction": base, "delta": torch.zeros_like(base)}
        context, support_weights = self.query_attention(
            query_pair,
            state.support_tokens,
            state.support_tokens,
            need_weights=True,
            average_attn_weights=True,
        )
        context = self.query_norm(query_pair + context)
        expanded_protein = state.protein.unsqueeze(1).expand_as(query_pair)
        delta_input = torch.cat(
            [query_pair, context, query_pair * context, expanded_protein], dim=-1
        )
        delta_scale = self.delta_head(delta_input).squeeze(-1)
        measured_residual = torch.einsum(
            "bqk,bk->bq", support_weights, state.support_residual
        )
        delta = delta_scale * measured_residual
        return {"prediction": base + delta, "base_prediction": base, "delta": delta}

    def forward(
        self,
        protein: torch.Tensor,
        support_ligand: torch.Tensor,
        support_y: torch.Tensor,
        support_f0: torch.Tensor,
        query_ligand: torch.Tensor,
        query_f0: torch.Tensor,
        *,
        use_protein: bool = True,
        use_support: bool = True,
    ) -> dict[str, torch.Tensor]:
        state = self.adapt(
            protein,
            support_ligand,
            support_y,
            support_f0,
            use_protein=use_protein,
        )
        return self.predict(
            state,
            protein,
            query_ligand,
            query_f0,
            use_protein=use_protein,
            use_support=use_support,
        )


def masked_mse_per_episode(
    prediction: torch.Tensor, target: torch.Tensor, mask: torch.Tensor
) -> torch.Tensor:
    return ((prediction - target).square() * mask).sum(1) / mask.sum(1).clamp_min(1.0)


def ranking_loss_per_episode(
    prediction: torch.Tensor, target: torch.Tensor, mask: torch.Tensor
) -> torch.Tensor:
    dp = prediction.unsqueeze(2) - prediction.unsqueeze(1)
    dy = target.unsqueeze(2) - target.unsqueeze(1)
    pair_mask = mask.unsqueeze(2) * mask.unsqueeze(1) * (dy.abs() > 1e-6)
    loss = F.softplus(-torch.sign(dy) * dp) * pair_mask
    return loss.sum((1, 2)) / pair_mask.sum((1, 2)).clamp_min(1.0)


def counterfactual_meta_loss(
    positive_prediction: torch.Tensor,
    negative_predictions: torch.Tensor,
    base_prediction: torch.Tensor,
    target: torch.Tensor,
    mask: torch.Tensor,
    temperature: float = COUNTERFACTUAL_TEMPERATURE,
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    """Base-anchored InfoNCE over post-adaptation query ranking gains."""

    if temperature <= 0:
        raise ValueError("counterfactual temperature must be positive")
    if negative_predictions.ndim != 3:
        raise ValueError("negative predictions must have shape [B,N,Q]")
    base = ranking_loss_per_episode(base_prediction.detach(), target, mask)
    positive = ranking_loss_per_episode(positive_prediction, target, mask)
    negative = torch.stack(
        [
            ranking_loss_per_episode(negative_predictions[:, row], target, mask)
            for row in range(negative_predictions.shape[1])
        ],
        dim=1,
    )
    positive_gain = base - positive
    # A wrong support is a competitor only while it improves over the frozen
    # base. Once it is no better than base, stop rewarding further damage.
    negative_gain = F.relu(base.unsqueeze(1) - negative)
    logits = torch.cat([positive_gain.unsqueeze(1), negative_gain], dim=1) / temperature
    labels = torch.zeros(logits.shape[0], dtype=torch.long, device=logits.device)
    loss = F.cross_entropy(logits, labels)
    return loss, {
        "base_ranking_loss": base,
        "positive_ranking_loss": positive,
        "negative_ranking_loss": negative,
        "positive_ranking_gain": positive_gain,
        "negative_ranking_gain": negative_gain,
        "wrong_support_harm": F.relu(negative - base.unsqueeze(1)),
        "ranking_gap": negative.mean(1) - positive,
    }


class GpuTelemetry:
    """Low-frequency NVIDIA telemetry outside the training synchronization path."""

    def __init__(self, enabled: bool, interval: float = 0.5) -> None:
        self.enabled = enabled
        self.interval = interval
        self.samples: list[dict[str, float]] = []
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def _sample(self) -> None:
        flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        while not self._stop.wait(self.interval):
            try:
                result = subprocess.run(
                    [
                        "nvidia-smi",
                        "--query-gpu=utilization.gpu,power.draw,memory.used",
                        "--format=csv,noheader,nounits",
                    ],
                    capture_output=True,
                    text=True,
                    timeout=3,
                    check=True,
                    creationflags=flags,
                )
                values = [float(value.strip()) for value in result.stdout.splitlines()[0].split(",")]
                self.samples.append({
                    "utilization_percent": values[0],
                    "power_watts": values[1],
                    "memory_mib": values[2],
                })
            except (OSError, subprocess.SubprocessError, ValueError, IndexError):
                return

    def __enter__(self) -> "GpuTelemetry":
        if self.enabled:
            self._thread = threading.Thread(target=self._sample, daemon=True)
            self._thread.start()
        return self

    def __exit__(self, *_: object) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2)

    def summary(self) -> dict[str, Any]:
        if not self.samples:
            return {"samples": 0}
        util = np.array([row["utilization_percent"] for row in self.samples])
        power = np.array([row["power_watts"] for row in self.samples])
        memory = np.array([row["memory_mib"] for row in self.samples])
        return {
            "samples": len(self.samples),
            "utilization_mean_percent": float(util.mean()),
            "utilization_p10_percent": float(np.quantile(util, 0.10)),
            "utilization_p90_percent": float(np.quantile(util, 0.90)),
            "busy_fraction_ge_40_percent": float((util >= 40).mean()),
            "power_mean_watts": float(power.mean()),
            "power_peak_watts": float(power.max()),
            "memory_peak_mib": float(memory.max()),
        }


def _phase_timer(device: torch.device) -> float:
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    return time.perf_counter()


def _batch(group: EpisodeTensorGroup, selected: torch.Tensor, store: GpuEpisodeStore) -> dict[str, torch.Tensor]:
    support_parent = group.support_parent[selected]
    query_parent = group.query_parent[selected]
    negative_parent = group.negative_support_parent[selected]
    return {
        "protein": store.target[group.target[selected]],
        "support_ligand": store.ligand[support_parent],
        "support_y": group.support_y[selected],
        "support_f0": store.f0[support_parent],
        "query_ligand": store.ligand[query_parent],
        "query_y": group.query_y[selected],
        "query_f0": store.f0[query_parent],
        "query_mask": group.query_mask[selected],
        "negative_support_ligand": store.ligand[negative_parent],
        "negative_support_y": group.negative_support_y[selected],
        "negative_support_f0": store.f0[negative_parent],
    }


def sample_target_balanced(
    group: EpisodeTensorGroup, batch_size: int
) -> torch.Tensor:
    return torch.multinomial(group.sampling_weight, batch_size, replacement=True)


@torch.no_grad()
def label_swap_support_y(
    model: ProteinConditionedMetaAdapter, batch: dict[str, torch.Tensor]
) -> torch.Tensor:
    """Move a wrong-target residual onto the correct support compounds."""

    correct_base, _, _ = model.base_predict(
        batch["protein"], batch["support_ligand"], batch["support_f0"]
    )
    donor_base, _, _ = model.base_predict(
        batch["protein"],
        batch["negative_support_ligand"][:, 0],
        batch["negative_support_f0"][:, 0],
    )
    donor_residual = batch["negative_support_y"][:, 0] - donor_base
    return correct_base + donor_residual


def parameter_inventory(model: ProteinConditionedMetaAdapter) -> list[dict[str, Any]]:
    base = {id(parameter) for parameter in model.base_parameters()}
    adapter = {id(parameter) for parameter in model.adapter_parameters()}
    inventory = []
    for name, parameter in model.named_parameters():
        parameter_id = id(parameter)
        if parameter_id in base:
            optimizer_group = "base_phase"
        elif parameter_id in adapter:
            optimizer_group = "adapter_phase"
        else:
            optimizer_group = "none"
        inventory.append({
            "module": name.split(".", 1)[0],
            "parameter_name": name,
            "shape": list(parameter.shape),
            "numel": parameter.numel(),
            "requires_grad_at_initialization": bool(parameter.requires_grad),
            "optimizer_group": optimizer_group,
        })
    return inventory


def _module_snapshot(
    model: ProteinConditionedMetaAdapter, modules: tuple[nn.Module, ...]
) -> dict[str, torch.Tensor]:
    module_ids = {id(parameter) for module in modules for parameter in module.parameters()}
    return {
        name: parameter.detach().clone()
        for name, parameter in model.named_parameters()
        if id(parameter) in module_ids
    }


def _module_gradient_norms(
    model: ProteinConditionedMetaAdapter, modules: tuple[nn.Module, ...]
) -> dict[str, float]:
    allowed = {id(parameter) for module in modules for parameter in module.parameters()}
    squared: dict[str, torch.Tensor] = {}
    for name, parameter in model.named_parameters():
        if id(parameter) not in allowed or parameter.grad is None:
            continue
        module = name.split(".", 1)[0]
        value = parameter.grad.detach().float().square().sum()
        squared[module] = value if module not in squared else squared[module] + value
    return {module: float(value.sqrt()) for module, value in squared.items()}


def _module_relative_delta(
    model: ProteinConditionedMetaAdapter,
    before: dict[str, torch.Tensor],
) -> dict[str, float]:
    numerator: dict[str, torch.Tensor] = {}
    denominator: dict[str, torch.Tensor] = {}
    for name, parameter in model.named_parameters():
        if name not in before:
            continue
        module = name.split(".", 1)[0]
        delta = (parameter.detach() - before[name]).float().square().sum()
        base = before[name].float().square().sum()
        numerator[module] = delta if module not in numerator else numerator[module] + delta
        denominator[module] = base if module not in denominator else denominator[module] + base
    return {
        module: float(numerator[module].sqrt() / (denominator[module].sqrt() + 1e-12))
        for module in numerator
    }


def _audit_step_indices(steps: int) -> set[int]:
    if steps <= 0:
        return set()
    fixed = {0, 99, 299, 999, 2999, steps - 1}
    fixed.update(round(position * (steps - 1) / 4) for position in range(5))
    return {step for step in fixed if 0 <= step < steps}


def train_base(
    model: ProteinConditionedMetaAdapter,
    store: GpuEpisodeStore,
    groups: dict[int, EpisodeTensorGroup],
    *,
    steps: int,
    batch_size: int,
    lr: float,
    log: list[dict[str, Any]],
) -> dict[str, Any]:
    model.train()
    optimizer = torch.optim.AdamW(model.base_parameters(), lr=lr, weight_decay=1e-4)
    scaler = torch.amp.GradScaler("cuda", enabled=store.device.type == "cuda")
    before = _module_snapshot(model, model.base_modules())
    audit_steps = _audit_step_indices(steps)
    started = _phase_timer(store.device)
    final = {"loss": float("nan"), "mse": float("nan"), "rank": float("nan")}
    group = groups[max(SUPPORT_K)]
    for step in range(steps):
        selected = sample_target_balanced(group, batch_size)
        batch = _batch(group, selected, store)
        with torch.amp.autocast("cuda", enabled=store.device.type == "cuda"):
            prediction, _, _ = model.base_predict(
                batch["protein"], batch["query_ligand"], batch["query_f0"]
            )
            mse = masked_mse_per_episode(
                prediction, batch["query_y"], batch["query_mask"]
            ).mean()
            rank = ranking_loss_per_episode(
                prediction, batch["query_y"], batch["query_mask"]
            ).mean()
            loss = mse + RANK_WEIGHT * rank
        optimizer.zero_grad(set_to_none=True)
        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        torch.nn.utils.clip_grad_norm_(model.base_parameters(), 5.0)
        scaler.step(optimizer)
        scaler.update()
        if step in audit_steps:
            final = {"loss": float(loss), "mse": float(mse), "rank": float(rank)}
            log.append({
                "phase": "base",
                "completed_steps": step + 1,
                **final,
                "gradient_l2_by_module": _module_gradient_norms(
                    model, model.base_modules()
                ),
            })
    elapsed = _phase_timer(store.device) - started
    return {
        **final,
        "seconds": elapsed,
        "steps_per_second": steps / max(elapsed, 1e-9),
        "episodes_per_second": steps * batch_size / max(elapsed, 1e-9),
        "relative_parameter_delta_by_module": _module_relative_delta(model, before),
    }


def train_adapter(
    model: ProteinConditionedMetaAdapter,
    store: GpuEpisodeStore,
    groups: dict[int, EpisodeTensorGroup],
    *,
    steps: int,
    batch_size: int,
    lr: float,
    counterfactual_weight: float,
    temperature: float,
    log: list[dict[str, Any]],
) -> dict[str, Any]:
    model.freeze_base()
    model.train()
    optimizer = torch.optim.AdamW(model.adapter_parameters(), lr=lr, weight_decay=1e-4)
    scaler = torch.amp.GradScaler("cuda", enabled=store.device.type == "cuda")
    before = _module_snapshot(model, model.adapter_modules())
    audit_steps = _audit_step_indices(steps)
    started = _phase_timer(store.device)
    final = {
        "loss": float("nan"), "mse": float("nan"), "rank": float("nan"),
        "counterfactual": float("nan"), "positive_ranking_gain": float("nan"),
        "wrong_support_harm": float("nan"), "ranking_gap": float("nan"),
    }
    n_arms = 1 + len(NEGATIVE_NAMES)
    for step in range(steps):
        k = SUPPORT_K[step % len(SUPPORT_K)]
        group = groups[k]
        selected = sample_target_balanced(group, batch_size)
        batch = _batch(group, selected, store)
        swapped_y = label_swap_support_y(model, batch)

        # One fused call for the positive, three wrong-target supports, and a
        # same-compound label-swap arm that cannot be identified by chemistry.
        support_ligand = torch.cat([
            batch["support_ligand"],
            *[
                batch["negative_support_ligand"][:, row]
                for row in range(len(EPISODE_NEGATIVE_NAMES))
            ],
            batch["support_ligand"],
        ])
        support_y = torch.cat([
            batch["support_y"],
            *[
                batch["negative_support_y"][:, row]
                for row in range(len(EPISODE_NEGATIVE_NAMES))
            ],
            swapped_y,
        ])
        support_f0 = torch.cat([
            batch["support_f0"],
            *[
                batch["negative_support_f0"][:, row]
                for row in range(len(EPISODE_NEGATIVE_NAMES))
            ],
            batch["support_f0"],
        ])
        protein = batch["protein"].repeat(n_arms, 1)
        query_ligand = batch["query_ligand"].repeat(n_arms, 1, 1)
        query_f0 = batch["query_f0"].repeat(n_arms, 1)

        with torch.amp.autocast("cuda", enabled=store.device.type == "cuda"):
            model_output = model(
                protein, support_ligand, support_y, support_f0, query_ligand, query_f0
            )
            prediction = model_output["prediction"].reshape(n_arms, batch_size, -1)
            base_prediction = model_output["base_prediction"].reshape(
                n_arms, batch_size, -1
            )[0]
            positive = prediction[0]
            negative = prediction[1:].transpose(0, 1)
            mse = masked_mse_per_episode(
                positive, batch["query_y"], batch["query_mask"]
            ).mean()
            rank = ranking_loss_per_episode(
                positive, batch["query_y"], batch["query_mask"]
            ).mean()
            contrast, detail = counterfactual_meta_loss(
                positive,
                negative,
                base_prediction,
                batch["query_y"],
                batch["query_mask"],
                temperature,
            )
            loss = mse + RANK_WEIGHT * rank + counterfactual_weight * contrast
        optimizer.zero_grad(set_to_none=True)
        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        torch.nn.utils.clip_grad_norm_(model.adapter_parameters(), 5.0)
        scaler.step(optimizer)
        scaler.update()
        if step in audit_steps:
            final = {
                "loss": float(loss),
                "mse": float(mse),
                "rank": float(rank),
                "counterfactual": float(contrast),
                "positive_ranking_gain": float(
                    detail["positive_ranking_gain"].mean()
                ),
                "wrong_support_harm": float(detail["wrong_support_harm"].mean()),
                "ranking_gap": float(detail["ranking_gap"].mean()),
            }
            log.append({
                "phase": "adapter",
                "completed_steps": step + 1,
                **final,
                "gradient_l2_by_module": _module_gradient_norms(
                    model, model.adapter_modules()
                ),
            })
    elapsed = _phase_timer(store.device) - started
    return {
        **final,
        "seconds": elapsed,
        "steps_per_second": steps / max(elapsed, 1e-9),
        "episodes_per_second": steps * batch_size / max(elapsed, 1e-9),
        "relative_parameter_delta_by_module": _module_relative_delta(model, before),
    }


@torch.no_grad()
def predict_group(
    model: ProteinConditionedMetaAdapter,
    store: GpuEpisodeStore,
    group: EpisodeTensorGroup,
    *,
    batch_size: int,
) -> dict[str, np.ndarray]:
    model.eval()
    arms: dict[str, list[np.ndarray]] = {
        "base": [],
        "adapted": [],
        **{f"adapted_{name}": [] for name in NEGATIVE_NAMES},
        "adapted_permuted_labels": [],
    }
    labels: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    n_arms = 1 + len(NEGATIVE_NAMES)
    for start in range(0, group.size, batch_size):
        selected = torch.arange(
            start, min(start + batch_size, group.size), device=store.device
        )
        batch = _batch(group, selected, store)
        swapped_y = label_swap_support_y(model, batch)
        b = len(selected)
        support_ligand = torch.cat([
            batch["support_ligand"],
            *[
                batch["negative_support_ligand"][:, row]
                for row in range(len(EPISODE_NEGATIVE_NAMES))
            ],
            batch["support_ligand"],
        ])
        support_y = torch.cat([
            batch["support_y"],
            *[
                batch["negative_support_y"][:, row]
                for row in range(len(EPISODE_NEGATIVE_NAMES))
            ],
            swapped_y,
        ])
        support_f0 = torch.cat([
            batch["support_f0"],
            *[
                batch["negative_support_f0"][:, row]
                for row in range(len(EPISODE_NEGATIVE_NAMES))
            ],
            batch["support_f0"],
        ])
        prediction = model(
            batch["protein"].repeat(n_arms, 1),
            support_ligand,
            support_y,
            support_f0,
            batch["query_ligand"].repeat(n_arms, 1, 1),
            batch["query_f0"].repeat(n_arms, 1),
        )
        block = prediction["prediction"].reshape(n_arms, b, -1)
        base = prediction["base_prediction"].reshape(n_arms, b, -1)[0]
        arms["base"].append(base.float().cpu().numpy())
        arms["adapted"].append(block[0].float().cpu().numpy())
        for row, name in enumerate(NEGATIVE_NAMES, start=1):
            arms[f"adapted_{name}"].append(block[row].float().cpu().numpy())

        if group.support_y.shape[1] > 1:
            permuted = batch["support_y"].roll(1, dims=1)
            permuted_prediction = model(
                batch["protein"],
                batch["support_ligand"],
                permuted,
                batch["support_f0"],
                batch["query_ligand"],
                batch["query_f0"],
            )["prediction"]
        else:
            permuted_prediction = torch.full_like(block[0], torch.nan)
        arms["adapted_permuted_labels"].append(
            permuted_prediction.float().cpu().numpy()
        )
        labels.append(batch["query_y"].cpu().numpy())
        masks.append(batch["query_mask"].cpu().numpy())
    return {
        **{name: np.concatenate(values) for name, values in arms.items()},
        "labels": np.concatenate(labels),
        "mask": np.concatenate(masks).astype(bool),
    }


def summarize_group(
    prediction: dict[str, np.ndarray],
    group: EpisodeTensorGroup,
    components: dict[str, int],
    *,
    seed: int,
    bootstrap: int,
) -> dict[str, Any]:
    def finite_mean(values: Iterable[float]) -> float:
        finite = np.asarray(list(values), dtype=np.float64)
        finite = finite[np.isfinite(finite)]
        return float(finite.mean()) if finite.size else float("nan")

    arm_names = [name for name in prediction if name not in {"labels", "mask"}]
    per_episode: dict[str, list[dict[str, float]]] = {name: [] for name in arm_names}
    for row in range(group.size):
        mask = prediction["mask"][row]
        truth = prediction["labels"][row][mask]
        for name in arm_names:
            values = prediction[name][row][mask]
            if np.isnan(values).all():
                per_episode[name].append({metric: float("nan") for metric in METRIC_NAMES})
            else:
                per_episode[name].append(episode_metrics(truth, values))

    per_target: dict[str, dict[str, dict[str, float]]] = {}
    for target in sorted(set(group.target_uids)):
        rows = [row for row, value in enumerate(group.target_uids) if value == target]
        per_target[target] = {
            arm: {
                metric: finite_mean(per_episode[arm][row][metric] for row in rows)
                for metric in METRIC_NAMES
            }
            for arm in arm_names
        }

    summary = {
        arm: {
            metric: finite_mean(value[arm][metric] for value in per_target.values())
            for metric in METRIC_NAMES
        }
        for arm in arm_names
    }
    rng = np.random.default_rng(seed)
    specificity: dict[str, Any] = {}
    for wrong in (f"adapted_{name}" for name in NEGATIVE_NAMES):
        entry = {}
        for metric in ("ci", "spearman", "ndcg10", "rmse"):
            sign = -1.0 if metric == "rmse" else 1.0
            delta = {
                target: sign * (values["adapted"][metric] - values[wrong][metric])
                for target, values in per_target.items()
            }
            entry[metric] = component_bootstrap(delta, components, rng, bootstrap)
        specificity[wrong] = entry
    if not np.isnan(summary["adapted_permuted_labels"]["ci"]):
        specificity["adapted_permuted_labels"] = {
            metric: component_bootstrap(
                {
                    target: values["adapted"][metric]
                    - values["adapted_permuted_labels"][metric]
                    for target, values in per_target.items()
                },
                components,
                rng,
                bootstrap,
            )
            for metric in ("ci", "spearman", "ndcg10")
        }
    return {
        "n_episodes": group.size,
        "n_targets": len(per_target),
        "summary": summary,
        "target_specific_ranking": specificity,
    }


def mechanism_snapshot(
    model: ProteinConditionedMetaAdapter,
    store: GpuEpisodeStore,
    group: EpisodeTensorGroup,
    *,
    batch_size: int,
    seed: int,
    split: str = "meta_validation",
) -> dict[str, Any]:
    predicted = predict_group(model, store, group, batch_size=batch_size)
    result = summarize_group(
        predicted,
        group,
        store.data.components,
        seed=seed,
        bootstrap=200,
    )
    metrics = ("rmse", "ci", "spearman", "ndcg10")
    return {
        "split": split,
        "k": int(group.support_y.shape[1]),
        "n_targets": result["n_targets"],
        "absolute": {
            arm: {metric: result["summary"][arm][metric] for metric in metrics}
            for arm in ("base", "adapted")
        },
        "correct_support_advantage": {
            arm: {
                metric: result["target_specific_ranking"][arm][metric]
                for metric in metrics
            }
            for arm in (f"adapted_{name}" for name in NEGATIVE_NAMES)
        },
    }


def mechanism_change(
    before: dict[str, Any], after: dict[str, Any]
) -> dict[str, Any]:
    return {
        "adapted_absolute_after_minus_before": {
            metric: after["absolute"]["adapted"][metric]
            - before["absolute"]["adapted"][metric]
            for metric in after["absolute"]["adapted"]
        },
        "specificity_mean_after_minus_before": {
            arm: {
                metric: after["correct_support_advantage"][arm][metric]["mean"]
                - before["correct_support_advantage"][arm][metric]["mean"]
                for metric in after["correct_support_advantage"][arm]
            }
            for arm in after["correct_support_advantage"]
        },
    }


def source_mechanism_gate(
    before: dict[str, Any], after: dict[str, Any]
) -> dict[str, Any]:
    primary = ("ci", "spearman", "ndcg10")
    absolute_gain = {
        metric: after["absolute"]["adapted"][metric]
        - before["absolute"]["adapted"][metric]
        for metric in primary
    }
    expected_arms = tuple(f"adapted_{name}" for name in NEGATIVE_NAMES)
    missing = set(expected_arms) - set(after["correct_support_advantage"])
    if missing:
        raise ValueError(f"source mechanism gate lacks counterfactual arms: {sorted(missing)}")
    specificity = {
        arm: {
            metric: after["correct_support_advantage"][arm][metric]["mean"]
            for metric in primary
        }
        for arm in expected_arms
    }
    values = [*absolute_gain.values()]
    values.extend(
        value for arm in specificity.values() for value in arm.values()
    )
    passed = all(np.isfinite(value) and value > 0.0 for value in values)
    return {
        "status": "PASS" if passed else "FAIL",
        "rule": (
            "after-adapter source meta-validation CI, Spearman and NDCG@10 must "
            "all improve over the frozen-base/no-adapter state, and correct support "
            "must beat each random, protein-hard, scaffold+ECFP chemical-matched, "
            "and same-compound label-swap support on all three metrics"
        ),
        "uses_recipient_labels": False,
        "absolute_ranking_gain": absolute_gain,
        "correct_support_advantage": specificity,
    }


def run(
    *,
    seed: int,
    formal: bool,
    protocol: str,
    base_steps: int,
    meta_steps: int,
    batch_size: int,
    lr: float,
    hidden: int,
    episode_root: Path,
    output: Path,
) -> dict[str, Any]:
    if formal and os.environ.get("A2S_FORMAL_EXTERNAL") != "1":
        raise RuntimeError(
            "formal training is external-only; set A2S_FORMAL_EXTERNAL=1 on the designated host"
        )
    if not torch.cuda.is_available():
        raise RuntimeError("A2S-CMAL training requires CUDA")
    if seed not in FORMAL_SEEDS:
        raise ValueError(f"seed must be one of the preregistered values {FORMAL_SEEDS}")

    device = torch.device("cuda")
    torch.manual_seed(seed)
    np.random.seed(seed)
    torch.cuda.reset_peak_memory_stats(device)
    # Recipient labels stay sealed until the source-only mechanism gate passes.
    data = load_data(episode_root, include_recipient_labels=False)
    base_targets, adapter_targets, source_partition = (
        component_disjoint_source_partition(data.target_splits)
    )
    prior = fit_source_prior(data, base_targets)
    store = GpuEpisodeStore(data, prior, device)
    base_groups = {
        max(SUPPORT_K): store.materialize(
            meta_split="meta_train",
            protocol=protocol,
            k=max(SUPPORT_K),
            targets=base_targets,
        )
    }
    adapter_groups = {
        k: store.materialize(
            meta_split="meta_train",
            protocol=protocol,
            k=k,
            targets=adapter_targets,
        )
        for k in SUPPORT_K
    }
    model = ProteinConditionedMetaAdapter(
        store.ligand_dim, store.protein_dim, hidden=hidden
    ).to(device)
    inventory = parameter_inventory(model)
    audit_group = store.materialize(
        meta_split="meta_validation", protocol=protocol, k=max(SUPPORT_K)
    )
    source_holdout_group = store.materialize(
        meta_split="meta_test", protocol=protocol, k=max(SUPPORT_K)
    )
    initial_snapshot = mechanism_snapshot(
        model, store, audit_group, batch_size=batch_size, seed=seed + 100
    )
    log: list[dict[str, Any]] = []
    started = time.perf_counter()
    with GpuTelemetry(enabled=True) as base_telemetry:
        base_result = train_base(
            model,
            store,
            base_groups,
            steps=base_steps,
            batch_size=batch_size,
            lr=lr,
            log=log,
        )
    pre_adapter_snapshot = mechanism_snapshot(
        model, store, audit_group, batch_size=batch_size, seed=seed + 101
    )
    pre_adapter_holdout_snapshot = mechanism_snapshot(
        model,
        store,
        source_holdout_group,
        batch_size=batch_size,
        seed=seed + 201,
        split="meta_test",
    )
    with GpuTelemetry(enabled=True) as adapter_telemetry:
        adapter_result = train_adapter(
            model,
            store,
            adapter_groups,
            steps=meta_steps,
            batch_size=batch_size,
            lr=lr,
            counterfactual_weight=COUNTERFACTUAL_WEIGHT,
            temperature=COUNTERFACTUAL_TEMPERATURE,
            log=log,
        )
    post_adapter_snapshot = mechanism_snapshot(
        model, store, audit_group, batch_size=batch_size, seed=seed + 102
    )
    post_adapter_holdout_snapshot = mechanism_snapshot(
        model,
        store,
        source_holdout_group,
        batch_size=batch_size,
        seed=seed + 202,
        split="meta_test",
    )
    mechanism_gate = source_mechanism_gate(
        pre_adapter_snapshot, post_adapter_snapshot
    )
    holdout_gate = source_mechanism_gate(
        pre_adapter_holdout_snapshot, post_adapter_holdout_snapshot
    )
    elapsed = time.perf_counter() - started
    combined_telemetry = GpuTelemetry(enabled=False)
    combined_telemetry.samples = base_telemetry.samples + adapter_telemetry.samples

    recipient_labels_read = False
    recipient_label_count = 0
    gate_stopped_formal = formal and (
        mechanism_gate["status"] != "PASS" or holdout_gate["status"] != "PASS"
    )
    if formal and not gate_stopped_formal:
        recipient_label_count = add_recipient_labels(data)
        recipient_labels_read = True

    evaluation_split = (
        None if gate_stopped_formal else
        "recipient_test" if formal else "meta_validation"
    )
    evaluation_protocol = "d0r" if formal else protocol
    budgets = () if gate_stopped_formal else SUPPORT_K if formal else (5,)
    evaluation = {}
    for k in budgets:
        group = store.materialize(
            meta_split=evaluation_split, protocol=evaluation_protocol, k=k
        )
        predicted = predict_group(model, store, group, batch_size=batch_size)
        evaluation[str(k)] = summarize_group(
            predicted,
            group,
            data.components,
            seed=seed + k,
            bootstrap=N_BOOTSTRAP if formal else 200,
        )

    report: dict[str, Any] = {
        "schema": "a2s-cmal-run-v2",
        "status": (
            "SOURCE_VALIDATION_MECHANISM_GATE_FAIL"
            if gate_stopped_formal else
            "FORMAL_RESULT" if formal else
            "MECHANISM_SMOKE_ONLY"
        ),
        "model": "A2S-CMAL protein-conditioned counterfactual meta-adaptation operator",
        "seed": seed,
        "formal": formal,
        "device": str(device),
        "gpu": torch.cuda.get_device_name(device),
        "data_package": {
            "path": str(episode_root),
            "content_sha256": data.manifest["content_sha256"],
            "recipient_labels_read": recipient_labels_read,
            "recipient_label_count": recipient_label_count,
        },
        "architecture": {
            "support_free_base": True,
            "protein_conditioned": True,
            "query_dependent": True,
            "closed_form_inner_update": False,
            "budget_gate": False,
            "deep_kernel": False,
            "parameters": sum(parameter.numel() for parameter in model.parameters()),
            "hidden": hidden,
        },
        "objective": {
            "base": "masked query MSE + pairwise ranking loss",
            "adaptation": "positive query fit + InfoNCE over post-adaptation query ranking loss",
            "counterfactuals": list(NEGATIVE_NAMES),
            "counterfactual_weight": COUNTERFACTUAL_WEIGHT,
            "temperature": COUNTERFACTUAL_TEMPERATURE,
        },
        "training": {
            "protocol": protocol,
            "source_partition": source_partition,
            "base_steps": base_steps,
            "meta_steps": meta_steps,
            "batch_episodes": batch_size,
            "base_phase": base_result,
            "adapter_phase": adapter_result,
            "wall_time_seconds": elapsed,
            "log": log,
        },
        "gradient_parameter_audit": {
            "inventory": inventory,
            "total_parameters": sum(row["numel"] for row in inventory),
            "base_optimizer_parameters": sum(
                row["numel"] for row in inventory
                if row["optimizer_group"] == "base_phase"
            ),
            "adapter_optimizer_parameters": sum(
                row["numel"] for row in inventory
                if row["optimizer_group"] == "adapter_phase"
            ),
            "base_frozen_during_adapter_phase": True,
            "snapshots": {
                "random_initialization": initial_snapshot,
                "after_base_before_meta_adapter": pre_adapter_snapshot,
                "after_meta_adapter": post_adapter_snapshot,
                "source_holdout_after_base_before_meta_adapter": (
                    pre_adapter_holdout_snapshot
                ),
                "source_holdout_after_meta_adapter": post_adapter_holdout_snapshot,
            },
            "trained_vs_untrained": mechanism_change(
                pre_adapter_snapshot, post_adapter_snapshot
            ),
            "source_mechanism_gate": mechanism_gate,
            "source_holdout_gate": holdout_gate,
        },
        "evaluation": {
            "split": evaluation_split,
            "statistical_unit": "homology component" if formal else "source validation target",
            "results": evaluation,
        },
        "gpu_pipeline": {
            "all_features_labels_episode_indices_resident": True,
            "pandas_in_training_loop": False,
            "per_step_host_transfer": False,
            "counterfactual_arms_fused": True,
            "telemetry": combined_telemetry.summary(),
            "telemetry_by_phase": {
                "base": base_telemetry.summary(),
                "adapter": adapter_telemetry.summary(),
            },
            "peak_torch_memory_mib": torch.cuda.max_memory_allocated(device) / 2**20,
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    if formal and not gate_stopped_formal:
        torch.save(
            {
                "state_dict": model.state_dict(),
                "seed": seed,
                "data_content_sha256": data.manifest["content_sha256"],
                "architecture": report["architecture"],
            },
            output.with_suffix(".pt"),
        )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="A2S-CMAL external trainer / local smoke")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--smoke", action="store_true")
    mode.add_argument("--formal", action="store_true")
    parser.add_argument("--seed", type=int, default=FORMAL_SEEDS[0])
    parser.add_argument("--protocol", choices=("ordered", "random"), default="ordered")
    parser.add_argument("--base-steps", type=int, default=None)
    parser.add_argument("--meta-steps", type=int, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--lr", type=float, default=3e-4)
    parser.add_argument("--hidden", type=int, default=128)
    parser.add_argument("--episodes", type=Path, default=EPISODES)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    formal = bool(args.formal)
    base_steps = args.base_steps if args.base_steps is not None else (3000 if formal else 20)
    meta_steps = args.meta_steps if args.meta_steps is not None else (3000 if formal else 20)
    batch_size = args.batch_size if args.batch_size is not None else (64 if formal else 16)
    default_name = (
        f"a2s_cmal_{args.protocol}_seed{args.seed}.json"
        if formal else f"a2s_cmal_smoke_{args.protocol}_seed{args.seed}.json"
    )
    output = args.out or (DEFAULT_OUT / default_name)
    report = run(
        seed=args.seed,
        formal=formal,
        protocol=args.protocol,
        base_steps=base_steps,
        meta_steps=meta_steps,
        batch_size=batch_size,
        lr=args.lr,
        hidden=args.hidden,
        episode_root=args.episodes,
        output=output,
    )
    print(json.dumps({
        "status": report["status"],
        "output": str(output),
        "base_phase": report["training"]["base_phase"],
        "adapter_phase": report["training"]["adapter_phase"],
        "gpu_pipeline": report["gpu_pipeline"],
    }, indent=2))


if __name__ == "__main__":
    main()
