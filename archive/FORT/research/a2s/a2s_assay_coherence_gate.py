"""Source-only assay-coherence information gate for passive A2S-DTA.

The gate asks whether k labels measured in one exact ChEMBL assay carry
target-specific ranking information for chemically distant compounds measured
in that same assay.  It does not train a meta-adaptation model.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from research.a2s.a2s_hotspot_falsification import standardise
from research.a2s.a2s_mode_gates import build_basis, source_heads
from research.a2s.a2s_mode_generalization import SourceSubspace, build_subspace, empirical_bayes_head
from research.a2s.a2s_trace import (
    DEFAULT_LOCK,
    DEFAULT_OOF,
    DEVICE,
    Substrate,
    analytic_delta,
    load_substrate,
    tanimoto,
)
from research.a2s.a2s_trace_stratum import metric_loss, paired_bootstrap


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_assay_coherence_gate_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "a2s_assay_coherence_gate_records_2026-08-02.parquet"

SEED = 20260802
K_VALUES = (1, 3, 5)
DRAWS = 8
MAX_QUERY = 32
MIN_QUERY = 8
KRR_RIDGE = 0.03
KRR_SCALE = 1.5
BOOTSTRAP_DRAWS = 2000
MDE = 0.005
SIMILARITY_EDGES = (0.0, 0.20, 0.35, 0.55, 1.0001)
SIMILARITY_NAMES = ("t00_20", "t20_35", "t35_55", "t55_100")
METHODS = ("base", "level", "krr", "scaled_krr", "eb_desc", "eb_original")
ARMS = ("correct", "permuted", "wrong_target")


@dataclass(frozen=True)
class Episode:
    episode_id: str
    target: str
    component: str
    assay: str
    role: str
    draw: int
    support5: tuple[int, ...]
    query: tuple[int, ...]


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def stable_seed(*values: object) -> int:
    token = ":".join(str(value) for value in values)
    return int(sha256(token.encode()).hexdigest()[:8], 16)


def similarity_bin(values: np.ndarray) -> np.ndarray:
    positions = np.digitize(values, np.asarray(SIMILARITY_EDGES[1:-1]), right=False)
    return np.asarray([SIMILARITY_NAMES[int(position)] for position in positions], dtype=object)


def build_episodes(substrate: Substrate, role: str) -> list[Episode]:
    frame = substrate.labeled.loc[substrate.labeled.role == role]
    episodes: list[Episode] = []
    for (target, assay), group in frame.groupby(["target", "assays"], sort=True):
        if len(group) < max(K_VALUES) + MIN_QUERY:
            continue
        rows = group.sort_values(["scaffold", "conn", "source_row"]).index.to_numpy()
        for draw in range(DRAWS):
            rng = np.random.default_rng(stable_seed(SEED, role, target, assay, draw))
            ordered = rows[rng.permutation(len(rows))]
            support = tuple(int(value) for value in ordered[: max(K_VALUES)])
            query = tuple(int(value) for value in ordered[max(K_VALUES) : max(K_VALUES) + MAX_QUERY])
            if len(query) < MIN_QUERY:
                continue
            episodes.append(
                Episode(
                    episode_id=sha256(
                        canonical([role, str(target), str(assay), draw, support, query]).encode()
                    ).hexdigest()[:20],
                    target=str(target),
                    component=str(group.component.iloc[0]),
                    assay=str(assay),
                    role=role,
                    draw=draw,
                    support5=support,
                    query=query,
                )
            )
    return episodes


def donor_indices(episodes: list[Episode]) -> list[int]:
    donors: list[int] = []
    for index, episode in enumerate(episodes):
        donor = None
        for offset in range(1, len(episodes)):
            candidate = (index + offset) % len(episodes)
            if episodes[candidate].target != episode.target:
                donor = candidate
                break
        if donor is None:
            raise RuntimeError("no wrong-target donor episode exists")
        donors.append(donor)
    return donors


def arm_rows(
    episode: Episode,
    donor: Episode,
    k: int,
    arm: str,
    residual: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    source = donor if arm == "wrong_target" else episode
    rows = np.asarray(source.support5[:k], dtype=np.int64)
    values = residual[rows].copy()
    if arm == "permuted":
        if k == 1:
            values = -values
        else:
            rng = np.random.default_rng(stable_seed(SEED, episode.episode_id, k, arm))
            permutation = rng.permutation(k)
            if np.array_equal(permutation, np.arange(k)):
                permutation = np.roll(permutation, 1)
            values = values[permutation]
    return rows, values


def nearest_tanimoto(substrate: Substrate, support: np.ndarray, query: np.ndarray) -> np.ndarray:
    support_bits = substrate.bits[torch.as_tensor(support, device=DEVICE)].unsqueeze(0)
    query_bits = substrate.bits[torch.as_tensor(query, device=DEVICE)].unsqueeze(0)
    return tanimoto(query_bits, support_bits).amax(-1).squeeze(0).cpu().numpy()


def eb_delta(
    basis: np.ndarray,
    subspace: SourceSubspace,
    support: np.ndarray,
    support_residual: np.ndarray,
    query: np.ndarray,
) -> np.ndarray:
    head = empirical_bayes_head(
        subspace,
        basis[support],
        support_residual,
        rank=basis.shape[1],
    )
    # The source mean is a support-free prior and was already shown not to
    # improve the base.  Subtracting it gives an exact no-support no-op.
    return basis[query] @ (head - subspace.mean_head)


def predict(
    substrate: Substrate,
    bases: dict[str, np.ndarray],
    subspaces: dict[str, SourceSubspace],
    support: np.ndarray,
    support_residual: np.ndarray,
    query: np.ndarray,
) -> dict[str, np.ndarray]:
    query_tensor = torch.as_tensor(query, device=DEVICE)
    support_tensor = torch.as_tensor(support, device=DEVICE)
    query_bits = substrate.bits[query_tensor].unsqueeze(0)
    support_bits = substrate.bits[support_tensor].unsqueeze(0)
    residual_tensor = torch.as_tensor(support_residual, device=DEVICE, dtype=torch.float32).unsqueeze(0)
    krr = analytic_delta(
        query_bits,
        support_bits,
        residual_tensor,
        estimator="krr",
        ridge=KRR_RIDGE,
    ).squeeze(0).cpu().numpy()
    return {
        "level": np.full(len(query), float(np.mean(support_residual))),
        "krr": krr,
        "scaled_krr": KRR_SCALE * krr,
        "eb_desc": eb_delta(bases["desc"], subspaces["desc"], support, support_residual, query),
        "eb_original": eb_delta(
            bases["original"], subspaces["original"], support, support_residual, query
        ),
    }


def evaluate(
    substrate: Substrate,
    bases: dict[str, np.ndarray],
    subspaces: dict[str, SourceSubspace],
    episodes: list[Episode],
) -> pd.DataFrame:
    residual = (substrate.affinity - substrate.base).cpu().numpy().astype(np.float64)
    affinity = substrate.affinity.cpu().numpy().astype(np.float64)
    base_all = substrate.base.cpu().numpy().astype(np.float64)
    donors = donor_indices(episodes)
    records: list[dict[str, object]] = []
    for index, episode in enumerate(episodes):
        donor = episodes[donors[index]]
        query = np.asarray(episode.query, dtype=np.int64)
        labels = affinity[query]
        base = base_all[query]
        if float(np.std(labels)) < 1e-9:
            continue
        for k in K_VALUES:
            correct_support = np.asarray(episode.support5[:k], dtype=np.int64)
            similarity = nearest_tanimoto(substrate, correct_support, query)
            bins = similarity_bin(similarity)
            predictions: dict[str, np.ndarray] = {"base__correct": base}
            for arm in ARMS:
                support, support_residual = arm_rows(episode, donor, k, arm, residual)
                deltas = predict(
                    substrate, bases, subspaces, support, support_residual, query
                )
                for method, delta in deltas.items():
                    predictions[f"{method}__{arm}"] = base + delta
            for stratum in (*SIMILARITY_NAMES, "all"):
                active = np.ones(len(query), dtype=bool) if stratum == "all" else bins == stratum
                if int(active.sum()) < MIN_QUERY:
                    continue
                truth = labels[active]
                if float(np.std(truth)) < 1e-9:
                    continue
                entry: dict[str, object] = {
                    "episode_id": episode.episode_id,
                    "target": episode.target,
                    "component": episode.component,
                    "assay": episode.assay,
                    "role": episode.role,
                    "draw": episode.draw,
                    "k": k,
                    "stratum": stratum,
                    "n_query": int(active.sum()),
                    "nearest_tanimoto_mean": float(similarity[active].mean()),
                }
                for name, values in predictions.items():
                    for metric, value in metric_loss(truth, values[active]).items():
                        entry[f"{name}__{metric}"] = float(value)
                records.append(entry)
    return pd.DataFrame.from_records(records)


def contrast(frame: pd.DataFrame, left: str, right: str, metric: str) -> dict[str, float]:
    left_column = f"{left}__{metric}"
    right_column = f"{right}__{metric}"
    if left_column not in frame or right_column not in frame:
        return {}
    sign = -1.0 if metric == "rmse" else 1.0
    working = frame[["component", "target", left_column, right_column]].copy()
    working["value"] = sign * (working[left_column] - working[right_column])
    working = working.groupby(["component", "target"], as_index=False).value.mean()
    return paired_bootstrap(working, "value", draws=BOOTSTRAP_DRAWS)


def summarise(records: pd.DataFrame) -> dict[str, object]:
    contrasts = {
        "eb_desc_minus_base": ("eb_desc__correct", "base__correct"),
        "eb_original_minus_base": ("eb_original__correct", "base__correct"),
        "scaled_krr_minus_base": ("scaled_krr__correct", "base__correct"),
        "eb_desc_correct_minus_permuted": ("eb_desc__correct", "eb_desc__permuted"),
        "eb_desc_correct_minus_wrong": ("eb_desc__correct", "eb_desc__wrong_target"),
        "eb_original_correct_minus_permuted": ("eb_original__correct", "eb_original__permuted"),
        "eb_original_correct_minus_wrong": ("eb_original__correct", "eb_original__wrong_target"),
    }
    summary: dict[str, object] = {}
    for role in sorted(records.role.unique()):
        for k in K_VALUES:
            for stratum in (*SIMILARITY_NAMES, "all"):
                frame = records.loc[
                    (records.role == role) & (records.k == k) & (records.stratum == stratum)
                ]
                if len(frame) < 5:
                    continue
                cell: dict[str, object] = {
                    "targets": int(frame.target.nunique()),
                    "components": int(frame.component.nunique()),
                    "assays": int(frame.assay.nunique()),
                    "episodes": int(frame.episode_id.nunique()),
                    "nearest_tanimoto_mean": float(frame.nearest_tanimoto_mean.mean()),
                    "contrasts": {},
                }
                for name, (left, right) in contrasts.items():
                    cell["contrasts"][name] = {
                        metric: contrast(frame, left, right, metric)
                        for metric in ("ci", "ndcg10", "rmse")
                    }
                summary.setdefault(role, {}).setdefault(f"k{k}", {})[stratum] = cell
    return summary


def noise_audit(substrate: Substrate) -> dict[str, float]:
    residual = (substrate.affinity - substrate.base).cpu().numpy().astype(np.float64)
    frame = substrate.labeled.assign(residual=residual)
    probe = frame.loc[frame.role == "probe"]
    assay_sd = probe.groupby(["target", "assays"]).residual.std().dropna()
    target_sd = probe.groupby("target").residual.std().dropna()
    return {
        "median_within_assay_residual_sd": float(assay_sd.median()),
        "median_within_target_residual_sd": float(target_sd.median()),
        "ratio": float(assay_sd.median() / target_sd.median()),
        "assay_groups": int(len(assay_sd)),
        "targets": int(probe.target.nunique()),
    }


def synthetic_control(
    substrate: Substrate, basis: np.ndarray, episodes: list[Episode]
) -> dict[str, object]:
    design = basis[:, :2]
    fit_targets = sorted(substrate.labeled.loc[substrate.labeled.role == "fit", "target"].unique())
    fit_heads = []
    for target in fit_targets:
        rng = np.random.default_rng(stable_seed("synthetic-head", target))
        fit_heads.append(rng.choice([-1.0, 1.0], size=2) * rng.uniform(0.8, 1.2, size=2))
    subspace = build_subspace(np.asarray(fit_heads), sigma=0.1)
    rows: list[dict[str, object]] = []
    donors = donor_indices(episodes)
    for index, episode in enumerate(episodes):
        if episode.draw > 1:
            continue
        rng = np.random.default_rng(stable_seed("synthetic-head", episode.target))
        truth_head = rng.choice([-1.0, 1.0], size=2) * rng.uniform(0.8, 1.2, size=2)
        donor = episodes[donors[index]]
        donor_rng = np.random.default_rng(stable_seed("synthetic-head", donor.target))
        donor_head = donor_rng.choice([-1.0, 1.0], size=2) * donor_rng.uniform(0.8, 1.2, size=2)
        query = np.asarray(episode.query, dtype=np.int64)
        labels = design[query] @ truth_head
        for k in (3, 5):
            support = np.asarray(episode.support5[:k], dtype=np.int64)
            donor_support = np.asarray(donor.support5[:k], dtype=np.int64)
            correct_residual = design[support] @ truth_head
            wrong_residual = design[donor_support] @ donor_head
            correct = design[query] @ empirical_bayes_head(
                subspace, design[support], correct_residual, rank=2
            )
            wrong = design[query] @ empirical_bayes_head(
                subspace, design[donor_support], wrong_residual, rank=2
            )
            rows.append(
                {
                    "component": episode.component,
                    "target": episode.target,
                    "k": k,
                    "correct_ci": metric_loss(labels, correct)["ci"],
                    "wrong_ci": metric_loss(labels, wrong)["ci"],
                }
            )
    frame = pd.DataFrame(rows)
    output: dict[str, object] = {}
    for k in (3, 5):
        cell = frame.loc[frame.k == k].copy()
        cell["value"] = cell.correct_ci - cell.wrong_ci
        cell = cell.groupby(["component", "target"], as_index=False).value.mean()
        output[f"k{k}"] = paired_bootstrap(cell, "value", draws=BOOTSTRAP_DRAWS)
    output["pass"] = bool(all(value["lower95"] > MDE for value in output.values()))
    return output


def decide(summary: dict[str, object], synthetic: dict[str, object]) -> dict[str, object]:
    gates: list[dict[str, object]] = []
    for k in (3, 5):
        for stratum in ("t00_20", "t20_35"):
            cell = summary.get("probe", {}).get(f"k{k}", {}).get(stratum)
            if not cell:
                continue
            for method in ("eb_desc", "eb_original"):
                gain = cell["contrasts"][f"{method}_minus_base"]["ci"]
                assignment = cell["contrasts"][f"{method}_correct_minus_permuted"]["ci"]
                wrong = cell["contrasts"][f"{method}_correct_minus_wrong"]["ci"]
                gates.append(
                    {
                        "k": k,
                        "stratum": stratum,
                        "method": method,
                        "components": cell["components"],
                        "gain_lower95": gain.get("lower95"),
                        "assignment_lower95": assignment.get("lower95"),
                        "wrong_target_lower95": wrong.get("lower95"),
                        "pass": bool(
                            gain.get("lower95", -1.0) > MDE
                            and assignment.get("lower95", -1.0) > 0.0
                            and wrong.get("lower95", -1.0) > 0.0
                        ),
                    }
                )
    k3 = any(row["pass"] for row in gates if row["k"] == 3)
    k5 = any(row["pass"] for row in gates if row["k"] == 5)
    admitted = bool(synthetic.get("pass") and k3 and k5)
    return {
        "synthetic_pass": bool(synthetic.get("pass")),
        "k3_low_similarity_pass": k3,
        "k5_low_similarity_pass": k5,
        "cells": gates,
        "verdict": (
            "ASSAY_COHERENT_GLOBAL_ADAPTATION_INFORMATION_ADMITTED"
            if admitted
            else "ASSAY_COHERENT_GLOBAL_ADAPTATION_INFORMATION_NOT_ADMITTED"
        ),
    }


def run(lock_path: Path, output: Path, records_path: Path, oof_cache: Path) -> dict[str, object]:
    substrate, context = load_substrate(lock_path, oof_cache)
    if set(substrate.labeled.role.unique()) - {"fit", "probe"}:
        raise AssertionError("an unauthorized role entered the assay-coherence gate")
    original, original_stats = build_basis(substrate)
    fit_mask = (substrate.labeled.role == "fit").to_numpy()
    desc = standardise(substrate.desc.cpu().numpy().astype(np.float64), fit_mask)
    bases = {"desc": desc, "original": original}
    subspaces = {}
    for name, basis in bases.items():
        heads, _, sigma, _ = source_heads(substrate, basis)
        subspaces[name] = build_subspace(heads, sigma)
    fit_episodes = build_episodes(substrate, "fit")
    probe_episodes = build_episodes(substrate, "probe")
    records = pd.concat(
        [
            evaluate(substrate, bases, subspaces, fit_episodes),
            evaluate(substrate, bases, subspaces, probe_episodes),
        ],
        ignore_index=True,
    )
    records_path.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    summary = summarise(records)
    synthetic = synthetic_control(substrate, desc, probe_episodes)
    result = {
        "schema": "a2s-assay-coherence-gate-v1",
        "status": "SOURCE_ONLY",
        "protocol": {
            "context": "exact ChEMBL assay token string",
            "support_policy": "nested random-within-exact-assay",
            "k_values": list(K_VALUES),
            "draws": DRAWS,
            "krr_ridge": KRR_RIDGE,
            "scaled_krr_scale": KRR_SCALE,
            "similarity_edges": list(SIMILARITY_EDGES),
            "basis": original_stats,
            "aggregation": "draw/assay mean within target, then component bootstrap",
        },
        "data": {
            "fit_episodes": len(fit_episodes),
            "probe_episodes": len(probe_episodes),
            "fit_targets": len({episode.target for episode in fit_episodes}),
            "probe_targets": len({episode.target for episode in probe_episodes}),
            "probe_components": len({episode.component for episode in probe_episodes}),
            "roles_opened": ["fit", "probe"],
            "locked_labels_requested": False,
            "recipient_labels_requested": False,
            "source_context": context,
        },
        "noise_audit": noise_audit(substrate),
        "synthetic_positive_control": synthetic,
        "summary": summary,
        "decision": decide(summary, synthetic),
        "records": str(records_path.relative_to(ROOT)).replace("\\", "/"),
    }
    payload = canonical(result)
    result["content_sha256"] = sha256(payload.encode()).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--oof-cache", type=Path, default=DEFAULT_OOF)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    args = parser.parse_args()
    result = run(args.lock, args.output, args.records, args.oof_cache)
    print(json.dumps(result["decision"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
