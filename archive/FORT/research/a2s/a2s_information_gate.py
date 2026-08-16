"""Run the source-only support-information admission gate.

The script consumes the metadata-only lock emitted by ``a2s_source_lock.py``.
Only ``fit`` and ``probe`` source targets are opened.  The ``locked`` role is
never requested from the labeled Parquet read.  A bounded linear probe is used
as an information diagnostic, not as a proposed final adaptation mechanism.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "dataset" / "public" / "chembl_37" / "processed" / "dualcold"
REGISTRY = DATA / "registry.parquet"
FEATURES = DATA / "ligand_features.npz"
PROTEINS = DATA / "target_esm2.npz"
DEFAULT_LOCK = ROOT / "reports" / "active" / "a2s_source_information_gate_lock_2026-08-01.json"
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_source_information_gate_2026-08-01.json"
DEFAULT_OOF = ROOT / "reports" / "active" / "a2s_source_information_gate_oof_2026-08-01.npz"

SEED = 20260801
K_VALUES = (1, 3, 5)
MIN_QUERY = 8
MAX_QUERY = 32
N_FOLDS = 5
RIDGE = 10.0
PROBE_RIDGE = 10.0
SYNTHETIC_RIDGE = 0.25
PROTEIN_PROJECTION_DIM = 16
MORGAN_BITS = 1024
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
KEY_COLUMNS = [
    "target",
    "conn",
    "endpoint",
    "scaffold",
    "assays",
    "docs",
    "hcluster",
    "dual_cold_split",
]
LABEL_COLUMNS = [*KEY_COLUMNS, "affinity"]


@dataclass(frozen=True)
class Episode:
    episode_id: str
    target: str
    component: str
    role: str
    k: int
    support_rows: tuple[int, ...]
    query_rows: tuple[int, ...]


@dataclass
class RidgeFit:
    mean: np.ndarray
    scale: np.ndarray
    weight: np.ndarray
    label_mean: float

    def predict(self, values: np.ndarray) -> np.ndarray:
        return ((values - self.mean) / self.scale) @ self.weight + self.label_mean


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_lock(path: Path) -> dict[str, object]:
    lock = json.loads(path.read_text(encoding="utf-8"))
    expected = str(lock.get("content_sha256", ""))
    payload = dict(lock)
    payload.pop("content_sha256", None)
    actual = sha256(canonical(payload).encode("utf-8")).hexdigest()
    if not expected or expected != actual:
        raise ValueError("source lock content hash is invalid")
    if lock.get("status") != "METADATA_PREFLIGHT_ONLY":
        raise ValueError("unexpected source lock status")
    if lock.get("model_family") not in {
        "a2s_information_gate_v1",
        "a2s_information_gate_v2",
    }:
        raise ValueError("unexpected model family in source lock")
    if lock.get("source_policy", {}).get("labels_read") is not False:
        raise ValueError("source lock is not label-blind")
    return lock


def tokens(value: object) -> set[str]:
    return {part.strip() for part in str(value).split("|") if part.strip()}


def load_metadata(lock: dict[str, object]) -> pd.DataFrame:
    frame = pd.read_parquet(REGISTRY, columns=KEY_COLUMNS).reset_index(names="source_row")
    frame = frame[
        (frame.endpoint.astype(str) == "pKi")
        & (frame.dual_cold_split.astype(str) == "train")
    ].copy()
    retained_rows = lock.get("row_policy", {}).get("retained_source_rows")
    if retained_rows is not None:
        expected_hash = lock["row_policy"].get("retained_source_rows_sha256")
        actual_hash = sha256(canonical(retained_rows).encode("utf-8")).hexdigest()
        if not expected_hash or expected_hash != actual_hash:
            raise ValueError("retained source-row hash is invalid")
        frame = frame.loc[frame.source_row.isin(set(retained_rows))].copy()
    assignment = lock["components"]["assignment"]
    role_by_component = lock["components"]["role_by_component"]
    frame["component"] = frame.target.astype(str).map(assignment)
    frame["role"] = frame.component.map(role_by_component)
    if frame.component.isna().any() or frame.role.isna().any():
        raise RuntimeError("metadata rows are missing locked component roles")
    if set(frame.role) != {"fit", "probe", "locked"}:
        raise RuntimeError("source roles are incomplete")
    oof_assignment = lock.get("oof", {}).get("role_by_component", {})
    frame["oof_fold"] = frame.component.map(oof_assignment)
    if oof_assignment and frame.loc[frame.role == "fit", "oof_fold"].isna().any():
        raise RuntimeError("retained fit rows are missing locked OOF folds")
    return frame


def load_labeled_fit_probe(metadata: pd.DataFrame) -> pd.DataFrame:
    allowed = sorted(metadata.loc[metadata.role.isin(("fit", "probe")), "target"].unique())
    if not allowed:
        raise RuntimeError("fit/probe roles are empty")
    # The affinity read is predicate-filtered by target; locked targets are not
    # requested from Parquet.  The global source_row is recovered by joining to
    # the label-free full-registry identity table.
    labels = pd.read_parquet(
        REGISTRY,
        filters=[
            ("endpoint", "==", "pKi"),
            ("dual_cold_split", "==", "train"),
            ("target", "in", allowed),
        ],
        columns=LABEL_COLUMNS,
    )
    if set(labels.target.astype(str)) - set(allowed):
        raise RuntimeError("filtered label read returned a locked target")
    identity = metadata.loc[metadata.role.isin(("fit", "probe")), ["source_row", *KEY_COLUMNS]]
    merged = identity.merge(labels, on=KEY_COLUMNS, how="left", validate="one_to_one")
    if merged.affinity.isna().any():
        raise RuntimeError("retained source rows could not be aligned to labels")
    merged["source_row"] = merged.source_row.astype(np.int64)
    merged["target"] = merged.target.astype(str)
    merged["component"] = merged.target.map(
        metadata.set_index("target").component.to_dict()
    )
    merged["role"] = merged.target.map(metadata.set_index("target").role.to_dict())
    merged["oof_fold"] = merged.target.map(
        metadata.set_index("target").oof_fold.to_dict()
    )
    if set(merged.role) != {"fit", "probe"}:
        raise RuntimeError("label read contains a role outside fit/probe")
    if len(merged) != len(metadata.loc[metadata.role.isin(("fit", "probe"))]):
        raise RuntimeError("label-filtered rows do not cover fit/probe metadata")
    if not np.isfinite(merged.affinity.to_numpy(dtype=np.float64)).all():
        raise ValueError("source affinity values must be finite")
    return merged.sort_values("source_row", kind="stable").reset_index(drop=True)


def load_design(metadata: pd.DataFrame, labeled: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    archive = np.load(FEATURES, allow_pickle=False)
    feature = np.asarray(archive["feat"], dtype=np.float32)
    if len(feature) != len(pd.read_parquet(REGISTRY, columns=["target"])):
        raise RuntimeError("feature cache is not aligned to the registry row count")
    if feature.ndim != 2 or feature.shape[1] < MORGAN_BITS + 10:
        raise ValueError("ligand feature cache lacks the expected Morgan/descriptor layout")
    protein_archive = np.load(PROTEINS, allow_pickle=False)
    keys = [str(value) for value in protein_archive["keys"]]
    pooled = np.asarray(protein_archive["pooled"], dtype=np.float32)
    target_to_protein = {key: row for row, key in enumerate(keys)}
    missing = sorted(set(labeled.target.astype(str)) - set(target_to_protein))
    if missing:
        raise RuntimeError(f"protein cache lacks source targets: {missing[:5]}")

    # A fixed, label-free projection keeps the diagnostic base bounded while
    # retaining target conditioning.  The supervised model still receives the
    # complete ligand design vector and is cross-fitted by component.
    rng = np.random.default_rng(SEED)
    projection = rng.normal(0.0, 1.0 / np.sqrt(pooled.shape[1]), size=(pooled.shape[1], PROTEIN_PROJECTION_DIM)).astype(np.float32)
    protein = pooled @ projection
    target_rows = np.asarray([target_to_protein[str(target)] for target in labeled.target], dtype=np.int64)
    ligand = feature[labeled.source_row.to_numpy(dtype=np.int64)].copy()
    fit_rows = labeled.loc[labeled.role == "fit", "source_row"].to_numpy(dtype=np.int64)
    if len(fit_rows) == 0:
        raise RuntimeError("fit role has no rows")
    desc_center = feature[fit_rows, -10:].mean(axis=0)
    desc_scale = feature[fit_rows, -10:].std(axis=0).clip(min=1e-6)
    ligand[:, -10:] = (ligand[:, -10:] - desc_center) / desc_scale
    return np.hstack([ligand, protein[target_rows]]).astype(np.float32), feature, target_to_protein


def target_fold(component: str) -> int:
    digest = sha256(f"{SEED}:oof:{component}".encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % N_FOLDS


def oof_folds(labeled: pd.DataFrame) -> np.ndarray:
    fit_mask = labeled.role.to_numpy() == "fit"
    if "oof_fold" in labeled and labeled.loc[fit_mask, "oof_fold"].notna().all():
        return labeled.oof_fold.to_numpy(dtype=float)
    return labeled.component.astype(str).map(target_fold).to_numpy(dtype=float)


def fit_ridge(values: np.ndarray, labels: np.ndarray, targets: Iterable[str], ridge: float) -> RidgeFit:
    if not torch.cuda.is_available():
        raise RuntimeError("the source admission gate requires the drug CUDA environment")
    values = np.asarray(values, dtype=np.float32)
    labels = np.asarray(labels, dtype=np.float32)
    target_array = np.asarray(list(targets), dtype=str)
    unique, inverse, counts = np.unique(target_array, return_inverse=True, return_counts=True)
    weights = (1.0 / counts[inverse]).astype(np.float32)
    weights /= weights.mean()
    mean = np.average(values, axis=0, weights=weights).astype(np.float32)
    variance = np.average((values - mean) ** 2, axis=0, weights=weights)
    scale = np.sqrt(np.maximum(variance, 1e-8)).astype(np.float32)
    x = torch.as_tensor((values - mean) / scale, device=DEVICE, dtype=torch.float32)
    y = torch.as_tensor(labels, device=DEVICE, dtype=torch.float32)
    w = torch.as_tensor(weights, device=DEVICE, dtype=torch.float32)
    sqrt_w = torch.sqrt(w).unsqueeze(1)
    design = x * sqrt_w
    gram = design.T @ design
    center = (w * y).sum() / w.sum()
    rhs = design.T @ (sqrt_w[:, 0] * (y - center))
    eye = torch.eye(gram.shape[0], device=DEVICE, dtype=torch.float32)
    beta = torch.linalg.solve(gram + ridge * eye, rhs).cpu().numpy()
    return RidgeFit(mean=mean, scale=scale, weight=beta, label_mean=float(center.cpu()))


def component_oof(metadata: pd.DataFrame, labeled: pd.DataFrame, values: np.ndarray) -> tuple[np.ndarray, dict[str, object]]:
    base = np.full(len(labeled), np.nan, dtype=np.float32)
    fit_mask = labeled.role.to_numpy() == "fit"
    probe_mask = labeled.role.to_numpy() == "probe"
    components = labeled.component.astype(str)
    fold_assignment = oof_folds(labeled)
    fold_stats: list[dict[str, object]] = []
    for fold in range(N_FOLDS):
        hold = fit_mask & (fold_assignment == fold)
        train = fit_mask & ~hold
        if not hold.any() or not train.any():
            raise RuntimeError(f"OOF fold {fold} has an empty train/holdout role")
        model = fit_ridge(
            values[train],
            labeled.loc[train, "affinity"].to_numpy(np.float32),
            labeled.loc[train, "target"].astype(str),
            RIDGE,
        )
        base[hold] = model.predict(values[hold])
        fold_stats.append({"fold": fold, "train_rows": int(train.sum()), "holdout_rows": int(hold.sum()), "holdout_components": sorted(set(components[hold]))})
    fit_model = fit_ridge(
        values[fit_mask],
        labeled.loc[fit_mask, "affinity"].to_numpy(np.float32),
        labeled.loc[fit_mask, "target"].astype(str),
        RIDGE,
    )
    base[probe_mask] = fit_model.predict(values[probe_mask])
    if not np.isfinite(base).all():
        raise RuntimeError("component OOF base left missing predictions")
    return base, {"folds": fold_stats, "fit_rows": int(fit_mask.sum()), "probe_rows": int(probe_mask.sum()), "model": "target-balanced ligand+fixed-protein-projection ridge"}


def select_support_query(group: pd.DataFrame, k: int) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    ordered = group.sort_values(["scaffold", "conn", "docs", "assays", "source_row"], kind="stable")
    support: list[int] = []
    used = {field: set() for field in ("conn", "scaffold", "docs", "assays")}
    for row in ordered.itertuples(index=False):
        if any(str(getattr(row, field)) in used[field] for field in used):
            continue
        support.append(int(row.source_row))
        for field in used:
            used[field].add(str(getattr(row, field)))
        if len(support) == k:
            break
    if len(support) < k:
        return None
    support_set = set(support)
    query: list[int] = []
    for row in ordered.itertuples(index=False):
        source_row = int(row.source_row)
        if source_row in support_set:
            continue
        if any(str(getattr(row, field)) in used[field] for field in used):
            continue
        query.append(source_row)
        if len(query) == MAX_QUERY:
            break
    if len(query) < MIN_QUERY:
        return None
    return tuple(support), tuple(query)


def make_episodes(labeled: pd.DataFrame) -> list[Episode]:
    episodes: list[Episode] = []
    for (role, target), group in labeled.groupby(["role", "target"], sort=True):
        component = str(group.component.iloc[0])
        for k in K_VALUES:
            selected = select_support_query(group, k)
            if selected is None:
                continue
            support, query = selected
            episodes.append(Episode(
                episode_id=sha256(canonical([role, target, k, support, query]).encode("utf-8")).hexdigest()[:16],
                target=str(target), component=component, role=str(role), k=k,
                support_rows=support, query_rows=query,
            ))
    return episodes


def tanimoto(query_bits: np.ndarray, support_bits: np.ndarray) -> np.ndarray:
    intersection = support_bits @ query_bits.astype(np.float32)
    support_count = support_bits.sum(axis=1)
    query_count = float(query_bits.sum())
    union = support_count + query_count - intersection
    return intersection / np.maximum(union, 1.0)


def row_index(labeled: pd.DataFrame) -> dict[int, int]:
    return {int(value): position for position, value in enumerate(labeled.source_row)}


def episode_features(
    episode: Episode,
    labeled: pd.DataFrame,
    values: np.ndarray,
    raw_features: np.ndarray,
    base: np.ndarray,
    *,
    position: dict[int, int] | None = None,
    residuals: dict[int, float] | None = None,
    label_mask: bool = False,
    derange: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if position is None:
        position = row_index(labeled)
    support = list(episode.support_rows)
    query = list(episode.query_rows)
    support_indices = [position[row] for row in support]
    query_indices = [position[row] for row in query]
    support_bits = raw_features[np.asarray(support, dtype=np.int64), :MORGAN_BITS]
    support_residual = np.asarray(
        [residuals[row] if residuals is not None else float(labeled.iloc[position[row]].affinity - base[position[row]]) for row in support],
        dtype=np.float32,
    )
    if derange and len(support_residual) >= 3:
        support_residual = np.roll(support_residual, 1)
    if label_mask:
        support_residual = np.zeros_like(support_residual)
    rows: list[np.ndarray] = []
    labels: list[float] = []
    base_values: list[float] = []
    for query_row, query_index in zip(query, query_indices):
        similarity = tanimoto(raw_features[query_row, :MORGAN_BITS], support_bits)
        weight = similarity / max(float(np.abs(similarity).sum()), 1e-6)
        label_features = np.asarray([
            float(weight @ support_residual),
            float(np.abs(weight) @ np.abs(support_residual)),
            float(support_residual.mean()),
            float(support_residual.std()),
            float(np.clip(values[query_index, MORGAN_BITS], -3.0, 3.0) * support_residual.mean()),
        ], dtype=np.float32)
        if label_mask:
            label_features.fill(0.0)
        chemistry_features = np.asarray([
            float(similarity.mean()),
            float(similarity.max(initial=0.0)),
            float(similarity.std()),
            float(np.mean([base[position[row]] for row in support])),
            float(np.std([base[position[row]] for row in support])),
            float(episode.k) / 5.0,
        ], dtype=np.float32)
        rows.append(np.concatenate((values[query_index], chemistry_features, label_features)))
        labels.append(float(labeled.iloc[query_index].affinity))
        base_values.append(float(base[query_index]))
    return np.asarray(rows, dtype=np.float32), np.asarray(labels, dtype=np.float32), np.asarray(base_values, dtype=np.float32), np.asarray([episode.component] * len(rows))


def build_real_matrix(
    episodes: list[Episode],
    labeled: pd.DataFrame,
    values: np.ndarray,
    raw_features: np.ndarray,
    base: np.ndarray,
    k: int,
    *,
    role: str,
    mask: bool = False,
    derange: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], list[str], np.ndarray]:
    selected = [episode for episode in episodes if episode.k == k and episode.role == role]
    if not selected:
        raise RuntimeError(f"no {role} episodes for k={k}")
    position = row_index(labeled)
    matrices: list[np.ndarray] = []
    labels: list[np.ndarray] = []
    bases: list[np.ndarray] = []
    components: list[str] = []
    targets: list[str] = []
    episode_sizes: list[int] = []
    for episode in selected:
        x, y, b, component = episode_features(
            episode, labeled, values, raw_features, base,
            position=position, label_mask=mask, derange=derange,
        )
        matrices.append(x)
        labels.append(y)
        bases.append(b)
        components.extend(component.tolist())
        targets.extend([episode.target] * len(y))
        episode_sizes.append(len(y))
    return np.vstack(matrices), np.concatenate(labels), np.concatenate(bases), components, targets, np.asarray(episode_sizes, dtype=np.int64)


def standardize_train_eval(train: np.ndarray, eval_values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mean = train.mean(axis=0).astype(np.float32)
    scale = train.std(axis=0).clip(min=1e-6).astype(np.float32)
    return (train - mean) / scale, (eval_values - mean) / scale, mean, scale


def ranks(values: np.ndarray) -> np.ndarray:
    return pd.Series(values).rank(method="average").to_numpy(dtype=np.float64)


def metric_loss(label: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    if len(label) < 2:
        return {"rmse": float("nan"), "ci": float("nan"), "spearman": float("nan"), "ndcg10": float("nan"), "rank_loss": float("nan")}
    error = label - prediction
    left, right = np.triu_indices(len(label), k=1)
    truth = np.sign(label[left] - label[right])
    active = truth != 0
    pred = np.sign(prediction[left] - prediction[right])
    ci = float((pred[active] == truth[active]).mean() + 0.5 * (pred[active] == 0).mean()) if active.any() else float("nan")
    label_rank = ranks(label)
    prediction_rank = ranks(prediction)
    spearman = float(np.corrcoef(label_rank, prediction_rank)[0, 1]) if np.std(label_rank) and np.std(prediction_rank) else float("nan")
    cutoff = min(10, len(label))
    order = np.argsort(-prediction, kind="stable")[:cutoff]
    ideal = np.argsort(-label, kind="stable")[:cutoff]
    gain = label - label.min() + 1e-6
    discounts = 1.0 / np.log2(np.arange(2, cutoff + 2))
    dcg = float(np.sum(gain[order] * discounts))
    ideal_dcg = float(np.sum(gain[ideal] * discounts))
    ndcg = dcg / ideal_dcg if ideal_dcg > 0 else float("nan")
    return {"rmse": float(np.sqrt(np.mean(error**2))), "ci": ci, "spearman": spearman, "ndcg10": ndcg, "rank_loss": 1.0 - ndcg}


def component_bootstrap(values: pd.DataFrame, value_column: str, seed: int = SEED, draws: int = 2000) -> dict[str, float]:
    grouped = values.groupby("component", sort=True)[value_column].mean().to_numpy(dtype=np.float64)
    grouped = grouped[np.isfinite(grouped)]
    if len(grouped) == 0:
        return {"components": 0, "mean": float("nan"), "lower95": float("nan"), "upper95": float("nan")}
    rng = np.random.default_rng(seed)
    sample = rng.integers(0, len(grouped), size=(draws, len(grouped)))
    bootstrap = grouped[sample].mean(axis=1)
    return {"components": int(len(grouped)), "mean": float(grouped.mean()), "lower95": float(np.quantile(bootstrap, 0.025)), "upper95": float(np.quantile(bootstrap, 0.975))}


def high_data_oracle(
    episodes: list[Episode], labeled: pd.DataFrame, values: np.ndarray, base: np.ndarray, role: str, k: int
) -> dict[str, np.ndarray]:
    """Fit a target-specific ridge on all non-query source rows for headroom only."""
    position = row_index(labeled)
    outputs: dict[str, np.ndarray] = {}
    for episode in episodes:
        if episode.role != role or episode.k != k or episode.target in outputs:
            continue
        target_mask = labeled.target.astype(str).to_numpy() == episode.target
        query_set = set(episode.query_rows)
        train_mask = target_mask & ~labeled.source_row.isin(query_set).to_numpy()
        if train_mask.sum() < 8:
            continue
        model = fit_ridge(values[train_mask], labeled.loc[train_mask, "affinity"].to_numpy(np.float32), labeled.loc[train_mask, "target"].astype(str), RIDGE)
        query_indices = [position[row] for row in episode.query_rows]
        outputs[episode.target] = model.predict(values[query_indices])
    return outputs


def synthetic_matrix(
    episodes: list[Episode], labeled: pd.DataFrame, values: np.ndarray, raw_features: np.ndarray, base: np.ndarray, k: int, role: str, *, mask: bool = False, derange: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], list[str]]:
    position = row_index(labeled)
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    bs: list[np.ndarray] = []
    components: list[str] = []
    targets: list[str] = []
    for episode in episodes:
        if episode.role != role or episode.k != k:
            continue
        rng = np.random.default_rng(int(sha256(f"{SEED}:synthetic:{episode.episode_id}".encode()).hexdigest()[:8], 16))
        support_residual = rng.normal(0.0, 1.0, size=episode.k).astype(np.float32)
        if derange and episode.k >= 3:
            support_residual = np.roll(support_residual, 1)
        if mask:
            support_residual = np.zeros_like(support_residual)
        support_bits = raw_features[np.asarray(episode.support_rows), :MORGAN_BITS]
        residual_map = {
            row: float(value)
            for row, value in zip(episode.support_rows, support_residual)
        }
        episode_x_all, _, _, _ = episode_features(
            episode,
            labeled,
            values,
            raw_features,
            base,
            position=position,
            residuals=residual_map,
            label_mask=mask,
            derange=derange,
        )
        episode_x: list[np.ndarray] = []
        episode_y: list[float] = []
        episode_b: list[float] = []
        for query_row in episode.query_rows:
            query_index = position[query_row]
            similarity = tanimoto(raw_features[query_row, :MORGAN_BITS], support_bits)
            weights = similarity / max(float(np.abs(similarity).sum()), 1e-6)
            query_modulation = float(np.clip(values[query_index, MORGAN_BITS], -3.0, 3.0))
            latent_delta = float(2.0 * query_modulation * support_residual.mean())
            offset = list(episode.query_rows).index(query_row)
            episode_x.append(episode_x_all[offset])
            episode_y.append(float(base[query_index] + latent_delta))
            episode_b.append(float(base[query_index]))
        xs.append(np.asarray(episode_x, dtype=np.float32))
        ys.append(np.asarray(episode_y, dtype=np.float32))
        bs.append(np.asarray(episode_b, dtype=np.float32))
        components.extend([episode.component] * len(episode_y))
        targets.extend([episode.target] * len(episode_y))
    if not xs:
        raise RuntimeError(f"no synthetic {role} episodes for k={k}")
    return np.vstack(xs), np.concatenate(ys), np.concatenate(bs), components, targets


def evaluate_probe(
    train_x: np.ndarray, train_y: np.ndarray, train_targets: list[str], eval_x: np.ndarray, eval_y: np.ndarray, eval_base: np.ndarray, eval_components: list[str], eval_targets: list[str],
) -> tuple[dict[str, object], RidgeFit]:
    train_xs, eval_xs, mean, scale = standardize_train_eval(train_x, eval_x)
    model = fit_ridge(train_xs, train_y, train_targets, PROBE_RIDGE)
    prediction = model.predict(eval_xs)
    records: list[dict[str, object]] = []
    for component, target in sorted(set(zip(eval_components, eval_targets))):
        active = np.asarray([(c == component and t == target) for c, t in zip(eval_components, eval_targets)], dtype=bool)
        metrics = metric_loss(eval_y[active], prediction[active])
        base_metrics = metric_loss(eval_y[active], eval_base[active])
        records.append({"component": component, "target": target, "n": int(active.sum()), "metrics": metrics, "base_metrics": base_metrics})
    table = pd.DataFrame([
        {
            "component": row["component"],
            "target": row["target"],
            "n": row["n"],
            **{f"adapt_{key}": value for key, value in row["metrics"].items()},
            **{f"base_{key}": value for key, value in row["base_metrics"].items()},
        }
        for row in records
    ])
    summary = {
        "rows": int(len(eval_y)),
        "targets": int(table.target.nunique()),
        "components": int(table.component.nunique()),
        "adaptation": {key: float(np.nanmean(table[f"adapt_{key}"].to_numpy())) for key in ("rmse", "ci", "spearman", "ndcg10", "rank_loss")},
        "base": {key: float(np.nanmean(table[f"base_{key}"].to_numpy())) for key in ("rmse", "ci", "spearman", "ndcg10", "rank_loss")},
        "component_bootstrap_adapt_minus_base": {
            key: component_bootstrap(pd.DataFrame({"component": table.component, "value": table[f"adapt_{key}"] - table[f"base_{key}"]}), "value")
            for key in ("rmse", "ci", "spearman", "ndcg10", "rank_loss")
        },
    }
    return {"summary": summary, "records": records, "prediction": prediction, "mean": mean.tolist(), "scale": scale.tolist()}, model


def run(lock_path: Path, output: Path, oof_output: Path) -> dict[str, object]:
    if DEVICE.type != "cuda":
        raise RuntimeError("run this gate with D:\\anaconda\\envs\\drug\\python.exe")
    lock = verify_lock(lock_path)
    metadata = load_metadata(lock)
    labeled = load_labeled_fit_probe(metadata)
    values, raw_features, _ = load_design(metadata, labeled)
    try:
        cached = np.load(oof_output, allow_pickle=False) if oof_output.exists() else None
        if cached is None:
            raise ValueError("no cached OOF")
        if not np.array_equal(cached["source_row"], labeled.source_row.to_numpy(np.int64)):
            raise RuntimeError("cached OOF source rows disagree with the current label lock")
        base = np.asarray(cached["base"], dtype=np.float32)
        if not np.isfinite(base).all():
            raise RuntimeError("cached OOF base contains non-finite predictions")
        fit_mask = labeled.role.to_numpy() == "fit"
        components = labeled.component.astype(str)
        fold_assignment = oof_folds(labeled)
        fold_stats = []
        for fold in range(N_FOLDS):
            hold = fit_mask & (fold_assignment == fold)
            train = fit_mask & ~hold
            fold_stats.append({
                "fold": fold,
                "train_rows": int(train.sum()),
                "holdout_rows": int(hold.sum()),
                "holdout_components": sorted(set(components[hold])),
            })
        oof_stats = {
            "reused": True,
            "fit_rows": int(fit_mask.sum()),
            "probe_rows": int((labeled.role == "probe").sum()),
            "folds": fold_stats,
            "model": "target-balanced ligand+fixed-protein-projection ridge",
        }
    except (OSError, KeyError, ValueError):
        base, oof_stats = component_oof(metadata, labeled, values)
    residual = labeled.affinity.to_numpy(np.float32) - base
    np.savez_compressed(
        oof_output,
        source_row=labeled.source_row.to_numpy(np.int64),
        base=base,
        residual=residual,
    )
    episodes = make_episodes(labeled)
    if not episodes:
        raise RuntimeError("no source episodes could be constructed")
    result: dict[str, object] = {
        "schema": (
            "a2s-source-information-gate-v2"
            if lock["model_family"] == "a2s_information_gate_v2"
            else "a2s-source-information-gate-v1"
        ),
        "status": "SOURCE_ONLY_DIAGNOSTIC",
        "decision": "NO_GO_INFORMATION_NOT_ADMITTED",
        "lock": {"path": str(lock_path), "content_sha256": lock["content_sha256"]},
        "labels": {"opened_roles": ["fit", "probe"], "locked_labels_requested": False, "recipient_labels_requested": False},
        "device": torch.cuda.get_device_name(0),
        "oof": oof_stats,
        "episodes": {
            "total": len(episodes),
            "by_role": {role: {str(k): sum(e.role == role and e.k == k for e in episodes) for k in K_VALUES} for role in ("fit", "probe")},
        },
        "gates": {},
        "synthetic_positive_control": {},
    }
    fold_rows = [int(item["holdout_rows"]) for item in oof_stats.get("folds", [])]
    max_holdout_share = (
        float(max(fold_rows) / sum(fold_rows))
        if fold_rows and sum(fold_rows)
        else None
    )
    oof_balance_pass = max_holdout_share is not None and max_holdout_share <= 0.25
    result["quality"] = {
        # Probe labels are development-only. Admission requires a separately
        # frozen protocol and a one-time locked-role evaluation.
        "admission_grade": False,
        "oof_fold_holdout_rows": fold_rows,
        "oof_max_holdout_share": max_holdout_share,
        "oof_balance_pass": oof_balance_pass,
        "oof_warning": (
            "OOF row balance failed; interpret real gates as non-confirmatory"
            if not oof_balance_pass
            else "OOF row balance passed; the result remains development-only until the protocol is frozen and the locked role is opened once"
        ),
    }
    for k in K_VALUES:
        train_x, train_y, _, _, train_targets, _ = build_real_matrix(episodes, labeled, values, raw_features, base, k, role="fit", mask=False)
        masked_train_x, _, _, _, masked_train_targets, _ = build_real_matrix(episodes, labeled, values, raw_features, base, k, role="fit", mask=True)
        eval_x, eval_y, eval_base, eval_components, eval_targets, _ = build_real_matrix(episodes, labeled, values, raw_features, base, k, role="probe", mask=False)
        masked_x, _, _, _, _, _ = build_real_matrix(episodes, labeled, values, raw_features, base, k, role="probe", mask=True)
        deranged_x, _, _, _, _, _ = build_real_matrix(episodes, labeled, values, raw_features, base, k, role="probe", mask=False, derange=True)
        train_xs, eval_xs, mean, scale = standardize_train_eval(train_x, eval_x)
        _, deranged_xs, _, _ = standardize_train_eval(train_x, deranged_x)
        masked_train_xs, masked_xs, _, _ = standardize_train_eval(masked_train_x, masked_x)
        # G0 and G1 have the same design width, ridge penalty and training
        # budget. Their only difference is whether the label channel is masked
        # in both training and evaluation.
        g0_model = fit_ridge(masked_train_xs, train_y, masked_train_targets, PROBE_RIDGE)
        g1_model = fit_ridge(train_xs, train_y, train_targets, PROBE_RIDGE)
        g0_pred = g0_model.predict(masked_xs)
        g1_pred = g1_model.predict(eval_xs)
        deranged_pred = g1_model.predict(deranged_xs)
        records: list[dict[str, object]] = []
        for component, target in sorted(set(zip(eval_components, eval_targets))):
            active = np.asarray([(c == component and t == target) for c, t in zip(eval_components, eval_targets)], dtype=bool)
            g0 = metric_loss(eval_y[active], g0_pred[active])
            g1 = metric_loss(eval_y[active], g1_pred[active])
            deranged = metric_loss(eval_y[active], deranged_pred[active])
            base_metrics = metric_loss(eval_y[active], eval_base[active])
            records.append({"component": component, "target": target, "n": int(active.sum()), "g0": g0, "g1": g1, "deranged": deranged, "base": base_metrics})
        gate_table = pd.DataFrame([
            {
                "component": row["component"],
                "target": row["target"],
                "delta_label_rank_loss": row["g0"]["rank_loss"] - row["g1"]["rank_loss"],
                "delta_assign_rank_loss": row["deranged"]["rank_loss"] - row["g1"]["rank_loss"],
                "delta_label_rmse": row["g0"]["rmse"] - row["g1"]["rmse"],
                "delta_assign_rmse": row["deranged"]["rmse"] - row["g1"]["rmse"],
                "base_rank_loss": row["base"]["rank_loss"],
                "g1_rank_loss": row["g1"]["rank_loss"],
            }
            for row in records
        ])
        gate_result: dict[str, object] = {
            "rows": int(len(eval_y)),
            "targets": int(len(records)),
            "components": int(gate_table.component.nunique()),
            "mean": {column: float(np.nanmean(gate_table[column].to_numpy())) for column in gate_table.columns if column not in ("component", "target")},
            "component_bootstrap": {
                column: component_bootstrap(gate_table[["component", column]].rename(columns={column: "value"}), "value")
                for column in gate_table.columns if column.startswith("delta_")
            },
            "record_count": len(records),
            "k1_assignment_control": "not applicable for k>=3",
        }
        if k == 1:
            # A true within-support assignment permutation does not exist at k=1.
            # The result retains this arm as a sensitivity placeholder only.
            gate_result["k1_assignment_control"] = "undefined; no within-support permutation was run"
        oracle = high_data_oracle(episodes, labeled, values, base, "probe", k)
        oracle_rows: list[dict[str, object]] = []
        position = row_index(labeled)
        for episode in episodes:
            if episode.role != "probe" or episode.k != k or episode.target not in oracle:
                continue
            query_indices = [position[row] for row in episode.query_rows]
            oracle_metrics = metric_loss(labeled.loc[query_indices, "affinity"].to_numpy(np.float32), oracle[episode.target])
            base_metrics = metric_loss(labeled.loc[query_indices, "affinity"].to_numpy(np.float32), base[query_indices])
            oracle_rows.append({
                "component": episode.component,
                "target": episode.target,
                "oracle_rank_loss": oracle_metrics["rank_loss"],
                "base_rank_loss": base_metrics["rank_loss"],
                "oracle_rmse": oracle_metrics["rmse"],
                "base_rmse": base_metrics["rmse"],
            })
        oracle_table = pd.DataFrame(oracle_rows)
        gate_result["headroom"] = {
            "oracle_targets": int(len(oracle_table)),
            "delta_headroom_rank_loss": component_bootstrap(pd.DataFrame({"component": oracle_table.component, "value": oracle_table.base_rank_loss - oracle_table.oracle_rank_loss}), "value") if len(oracle_table) else None,
            "delta_headroom_rmse": component_bootstrap(pd.DataFrame({"component": oracle_table.component, "value": oracle_table.base_rmse - oracle_table.oracle_rmse}), "value") if len(oracle_table) else None,
        }
        result["gates"][f"k{k}"] = gate_result

        syn_train_x, syn_train_y, _, _, syn_train_targets = synthetic_matrix(episodes, labeled, values, raw_features, base, k, "fit")
        syn_g0_train_x, _, _, _, syn_g0_train_targets = synthetic_matrix(episodes, labeled, values, raw_features, base, k, "fit", mask=True)
        syn_eval_x, syn_eval_y, syn_base, syn_components, syn_targets = synthetic_matrix(episodes, labeled, values, raw_features, base, k, "probe")
        syn_g0_x, _, _, _, _ = synthetic_matrix(episodes, labeled, values, raw_features, base, k, "probe", mask=True)
        # The positive control isolates the five label-channel coordinates.
        # This prevents the deliberately random synthetic support code from
        # being drowned by the 1,050-dimensional support-free design vector.
        syn_train_x = syn_train_x[:, -5:]
        syn_g0_train_x = syn_g0_train_x[:, -5:]
        syn_eval_x = syn_eval_x[:, -5:]
        syn_g0_x = syn_g0_x[:, -5:]
        syn_train_xs, syn_eval_xs, _, _ = standardize_train_eval(syn_train_x, syn_eval_x)
        syn_g0_train_xs, syn_g0_xs, _, _ = standardize_train_eval(syn_g0_train_x, syn_g0_x)
        syn_g1_model = fit_ridge(syn_train_xs, syn_train_y, syn_train_targets, SYNTHETIC_RIDGE)
        syn_g0_model = fit_ridge(syn_g0_train_xs, syn_train_y, syn_g0_train_targets, SYNTHETIC_RIDGE)
        syn_g1 = syn_g1_model.predict(syn_eval_xs)
        syn_g0 = syn_g0_model.predict(syn_g0_xs)
        result["synthetic_positive_control"][f"k{k}"] = {
            "g0_rank_loss": float(metric_loss(syn_eval_y, syn_g0)["rank_loss"]),
            "g1_rank_loss": float(metric_loss(syn_eval_y, syn_g1)["rank_loss"]),
            "delta_label_rank_loss": float(metric_loss(syn_eval_y, syn_g0)["rank_loss"] - metric_loss(syn_eval_y, syn_g1)["rank_loss"]),
            "rows": int(len(syn_eval_y)),
            "components": int(len(set(syn_components))),
        }
    result["interpretation"] = {
        "fact": "This artifact uses only TRAIN source labels from fit/probe roles; no locked or recipient labels were requested.",
        "inference": "A positive G0/G1 or assignment gap is evidence of exploitable information within this bounded probe class, not proof of a transferable biological mechanism.",
        "hypothesis": "If real gates are positive on a new locked role and survive the registered controls, a discrete learned adaptation policy may be worth implementing.",
    }
    result["content_sha256"] = sha256(canonical(result).encode("utf-8")).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True, default=float) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--oof-out", type=Path, default=DEFAULT_OOF)
    args = parser.parse_args()
    result = run(args.lock.resolve(), args.out.resolve(), args.oof_out.resolve())
    print(json.dumps({
        "status": result["status"],
        "output": str(args.out.resolve()),
        "oof_output": str(args.oof_out.resolve()),
        "episodes": result["episodes"],
        "gates": {key: {"mean": value["mean"], "component_bootstrap": value["component_bootstrap"]} for key, value in result["gates"].items()},
        "synthetic_positive_control": result["synthetic_positive_control"],
    }, indent=2, sort_keys=True, default=float))


if __name__ == "__main__":
    main()
