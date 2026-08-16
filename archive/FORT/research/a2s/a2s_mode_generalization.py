"""A2S-MODE Gates G1-G4: is the adaptation object generalizable, and is it worth building?

Gates A0/A1 established that a per-target head on a compact label-free ligand
basis carries ranking headroom in every relation stratum, including the ones
where every support-transport operator measures zero.  Two things were left
unexamined, and both decide whether A2S-MODE is worth implementing.

**The confound.**  A0's head is fitted on ~130 of the target's own rows and
evaluated on queries drawn from the same target.  Those queries were stratified
by their distance to the *support set*, never by their distance to the head's own
training rows.  If the headroom disappears once the evaluation compounds are
scaffold-cold with respect to the head's training rows, then A0 is the same
locality phenomenon measured at k ~ 130 instead of k ~ 5, and A2S-MODE inherits
the limit it was designed to escape.

**The value question.**  Even a genuinely global object is worthless if it cannot
be estimated from few labels.  The decisive number is not "does a discrete mode
help" but "how many labels does this object actually take", measured directly.

Four gates:

``G1`` locality audit  - within-target **scaffold-disjoint** head evaluation, and
                         evaluation rows stratified by their similarity to the
                         head's training rows.
``G2`` intrinsic dimension - spectrum of the source-target head matrix, and the
                         retained gain when a target's own head is projected onto
                         the top-``r`` **source** subspace.
``G3`` zero-shot prior - predict the head from the protein embedding alone, with
                         no support labels at all.  If this already pays, the
                         few-shot framing is not where the value is.
``G4`` label learning curve - empirical-Bayes head estimated from ``k`` support
                         rows *inside the rank-``r`` source subspace*, swept over
                         ``k`` and ``r``.  This is simultaneously the value
                         assessment for A2S-MODE and a direct re-test of the
                         A2S-IDA rank-``m`` code question on a better-conditioned
                         basis.

Only ``fit`` and ``probe`` roles are opened.  ``locked`` and the recipient roster
are never requested.
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

from research.a2s.a2s_information_gate import FEATURES, PROTEINS, REGISTRY, canonical, sha256_file
from research.a2s.a2s_trace import DEVICE, Substrate, load_substrate, tanimoto
from research.a2s.a2s_trace_stratum import DEFAULT_LOCK, DEFAULT_OOF, metric_loss, paired_bootstrap
from research.a2s.a2s_mode_gates import build_basis, fit_head, source_heads


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_mode_generalization_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "a2s_mode_generalization_records_2026-08-02.parquet"

SEED = 20260802
MIN_TARGET_ROWS = 40
MIN_HEAD_TRAIN = 20
MIN_EVAL_ROWS = 8
HEAD_RIDGE = 1.0
K_SWEEP = (1, 3, 5, 10, 20, 40)
RANK_SWEEP = (1, 2, 3, 5, 26)
DRAWS = 8
BOOTSTRAP_DRAWS = 2000
MDE = 0.005
SIMILARITY_EDGES = (0.0, 0.30, 0.45, 0.60, 1.0001)
SIMILARITY_NAMES = ("s00_30", "s30_45", "s45_60", "s60_100")


@dataclass
class TargetSplit:
    target: str
    component: str
    train_rows: np.ndarray
    eval_rows: np.ndarray
    split: str                 # "scaffold_disjoint" | "random"


# --------------------------------------------------------------------------- #
# Within-target splits
# --------------------------------------------------------------------------- #


def target_splits(substrate: Substrate, role: str = "probe") -> list[TargetSplit]:
    """Two within-target splits per target: scaffold-disjoint and random.

    The scaffold-disjoint split is the honest hard case: the head never sees a
    scaffold it is later scored on, so a gain there cannot be within-target
    memorisation of a chemical series.
    """

    frame = substrate.labeled
    splits: list[TargetSplit] = []
    for target, group in frame.loc[frame.role == role].groupby("target", sort=True):
        if len(group) < MIN_TARGET_ROWS:
            continue
        component = str(group.component.iloc[0])
        rows = group.index.to_numpy()
        scaffolds = group.scaffold.astype(str).to_numpy()
        digest = int(sha256(f"{SEED}:{target}".encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(digest)

        unique = sorted(set(scaffolds))
        order = rng.permutation(len(unique))
        held: set[str] = set()
        for index in order:
            held.add(unique[index])
            held_mask = np.isin(scaffolds, list(held))
            if held_mask.sum() >= max(MIN_EVAL_ROWS, int(0.3 * len(rows))):
                break
        held_mask = np.isin(scaffolds, list(held))
        if held_mask.sum() >= MIN_EVAL_ROWS and (~held_mask).sum() >= MIN_HEAD_TRAIN:
            splits.append(
                TargetSplit(str(target), component, rows[~held_mask], rows[held_mask], "scaffold_disjoint")
            )

        permutation = rng.permutation(len(rows))
        cut = max(MIN_HEAD_TRAIN, int(0.7 * len(rows)))
        if len(rows) - cut >= MIN_EVAL_ROWS:
            splits.append(
                TargetSplit(
                    str(target), component, rows[permutation[:cut]], rows[permutation[cut:]], "random"
                )
            )
    return splits


def similarity_to_train(
    substrate: Substrate, train_rows: np.ndarray, eval_rows: np.ndarray
) -> np.ndarray:
    """Max Tanimoto from each evaluation compound to the head's training rows."""

    train_bits = substrate.bits[torch.as_tensor(train_rows, device=DEVICE)].unsqueeze(0)
    eval_bits = substrate.bits[torch.as_tensor(eval_rows, device=DEVICE)].unsqueeze(0)
    return tanimoto(eval_bits, train_bits).amax(-1).squeeze(0).cpu().numpy()


def similarity_bin(values: np.ndarray) -> np.ndarray:
    index = np.digitize(values, np.asarray(SIMILARITY_EDGES[1:-1]), right=False)
    return np.asarray([SIMILARITY_NAMES[int(position)] for position in index], dtype=object)


# --------------------------------------------------------------------------- #
# Source subspace and empirical-Bayes head estimation
# --------------------------------------------------------------------------- #


@dataclass
class SourceSubspace:
    mean_head: np.ndarray        # (d,)
    directions: np.ndarray       # (d, d) columns are principal head directions
    variances: np.ndarray        # (d,) source variance along each direction
    spectrum: np.ndarray         # (d,) singular values
    sigma: float
    n_targets: int


def build_subspace(heads: np.ndarray, sigma: float) -> SourceSubspace:
    mean_head = heads.mean(axis=0)
    centred = heads - mean_head
    _, singular, right = np.linalg.svd(centred, full_matrices=False)
    variances = (singular**2) / max(len(heads) - 1, 1)
    return SourceSubspace(
        mean_head=mean_head,
        directions=right.T,
        variances=variances,
        spectrum=singular,
        sigma=sigma,
        n_targets=int(len(heads)),
    )


def empirical_bayes_head(
    subspace: SourceSubspace,
    design: np.ndarray,
    residual: np.ndarray,
    rank: int,
) -> np.ndarray:
    """Head estimated from `k` labels inside the top-`rank` source subspace.

    The prior is the source distribution of heads along each direction, so the
    ridge is not tuned: `lambda_j = sigma^2 / tau_j^2` with `tau_j` measured on
    source targets.  The level is profiled out (it is rank-null).
    """

    directions = subspace.directions[:, :rank]
    prior_prediction = design @ subspace.mean_head
    target = residual - prior_prediction
    projected = design @ directions
    centre = projected.mean(axis=0)
    projected = projected - centre
    target = target - target.mean()
    penalty = np.diag(subspace.sigma**2 / np.maximum(subspace.variances[:rank], 1e-8))
    gram = projected.T @ projected + penalty
    coefficients = np.linalg.solve(gram, projected.T @ target)
    return subspace.mean_head + directions @ coefficients


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #


def evaluate(
    substrate: Substrate,
    basis: np.ndarray,
    subspace: SourceSubspace,
    protein_head: np.ndarray | None,
    splits: list[TargetSplit],
) -> pd.DataFrame:
    residual_all = (substrate.affinity - substrate.base).cpu().numpy().astype(np.float64)
    affinity_all = substrate.affinity.cpu().numpy().astype(np.float64)
    base_all = substrate.base.cpu().numpy().astype(np.float64)
    target_rows = substrate.target_row.cpu().numpy()

    records: list[dict[str, object]] = []
    for split in splits:
        train_design = basis[split.train_rows]
        train_residual = residual_all[split.train_rows]
        eval_design = basis[split.eval_rows]
        labels = affinity_all[split.eval_rows]
        base = base_all[split.eval_rows]
        if float(np.std(labels)) < 1e-9:
            continue
        similarity = similarity_to_train(substrate, split.train_rows, split.eval_rows)
        bins = similarity_bin(similarity)

        predictions: dict[str, np.ndarray] = {"base": base}
        # G1: the full per-target oracle head.
        oracle_head, oracle_level = fit_head(train_design, train_residual, HEAD_RIDGE)
        predictions["target_head"] = base + eval_design @ oracle_head + oracle_level
        # The source mean head: a target-independent control.
        predictions["source_mean_head"] = base + eval_design @ subspace.mean_head
        # G3: zero-shot protein-predicted head, no support labels at all.
        if protein_head is not None:
            row = int(target_rows[split.train_rows[0]])
            predictions["protein_zero_shot"] = base + eval_design @ protein_head[row]
        # G2: the target's own head projected onto the top-r source subspace.
        centred = oracle_head - subspace.mean_head
        for rank in RANK_SWEEP:
            directions = subspace.directions[:, :rank]
            projected = subspace.mean_head + directions @ (directions.T @ centred)
            predictions[f"oracle_rank{rank}"] = base + eval_design @ projected + oracle_level
        # G4: empirical-Bayes head from k labels inside the rank-r source subspace.
        for k in K_SWEEP:
            if k > len(split.train_rows):
                continue
            for rank in RANK_SWEEP:
                estimates = []
                for draw in range(DRAWS):
                    digest = int(
                        sha256(f"{SEED}:{split.target}:{split.split}:{k}:{rank}:{draw}".encode()).hexdigest()[:8],
                        16,
                    )
                    rng = np.random.default_rng(digest)
                    chosen = rng.choice(len(split.train_rows), size=k, replace=False)
                    head = empirical_bayes_head(
                        subspace, train_design[chosen], train_residual[chosen], rank
                    )
                    estimates.append(eval_design @ head)
                predictions[f"eb_k{k}_rank{rank}"] = base + np.mean(estimates, axis=0)

        for stratum in (*SIMILARITY_NAMES, "all"):
            active = np.ones(len(labels), dtype=bool) if stratum == "all" else (bins == stratum)
            if int(active.sum()) < MIN_EVAL_ROWS:
                continue
            truth = labels[active]
            if float(np.std(truth)) < 1e-9:
                continue
            entry: dict[str, object] = {
                "target": split.target,
                "component": split.component,
                "split": split.split,
                "stratum": stratum,
                "n_eval": int(active.sum()),
                "n_train": int(len(split.train_rows)),
                "similarity_to_train_mean": float(similarity[active].mean()),
            }
            for name, prediction in predictions.items():
                for metric, value in metric_loss(truth, prediction[active]).items():
                    entry[f"{name}__{metric}"] = float(value)
            records.append(entry)
    return pd.DataFrame.from_records(records)


def protein_predicted_heads(
    substrate: Substrate, basis: np.ndarray, heads: np.ndarray, names: list[str], ridge: float = 100.0
) -> np.ndarray:
    """G3: ridge from the pooled protein embedding to head coefficients.

    Fitted on `fit` targets only; applied to every protein row so probe targets
    receive a head without any label of theirs being read.
    """

    frame = substrate.labeled
    target_rows = substrate.target_row.cpu().numpy()
    protein = substrate.protein.cpu().numpy().astype(np.float64)
    lookup = {
        str(target): int(target_rows[group.index.to_numpy()[0]])
        for target, group in frame.groupby("target", sort=True)
    }
    rows = np.asarray([lookup[name] for name in names], dtype=np.int64)
    design = protein[rows]
    centre = design.mean(axis=0)
    scale = design.std(axis=0).clip(min=1e-6)
    design = (design - centre) / scale
    gram = design.T @ design + ridge * np.eye(design.shape[1])
    weight = np.linalg.solve(gram, design.T @ (heads - heads.mean(axis=0)))
    all_design = (protein - centre) / scale
    return heads.mean(axis=0) + all_design @ weight


# --------------------------------------------------------------------------- #
# Summaries
# --------------------------------------------------------------------------- #


def contrast(frame: pd.DataFrame, left: str, right: str, metric: str) -> dict[str, float]:
    left_column, right_column = f"{left}__{metric}", f"{right}__{metric}"
    if left_column not in frame or right_column not in frame:
        return {}
    sign = -1.0 if metric == "rmse" else 1.0
    working = frame[["component", "target", left_column, right_column]].copy()
    working["value"] = sign * (working[left_column] - working[right_column])
    return paired_bootstrap(working, "value", draws=BOOTSTRAP_DRAWS)


def summarise(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for split in sorted(records.split.unique()):
        for stratum in sorted(records.stratum.unique()):
            frame = records.loc[(records.split == split) & (records.stratum == stratum)]
            if len(frame) < 5:
                continue
            cell: dict[str, object] = {
                "targets": int(frame.target.nunique()),
                "components": int(frame.component.nunique()),
                "similarity_to_train_mean": float(frame.similarity_to_train_mean.mean()),
                "absolute_ci": {},
                "vs_base": {},
            }
            names = sorted({name.rsplit("__", 1)[0] for name in frame.columns if name.endswith("__ci")})
            for name in names:
                cell["absolute_ci"][name] = float(np.nanmean(frame[f"{name}__ci"].to_numpy()))
                if name != "base":
                    cell["vs_base"][name] = contrast(frame, name, "base", "ci")
            summary.setdefault(split, {})[stratum] = cell
    return summary


def decide(summary: dict[str, object], subspace: SourceSubspace) -> dict[str, object]:
    def get(split: str, stratum: str, name: str) -> dict[str, float] | None:
        return summary.get(split, {}).get(stratum, {}).get("vs_base", {}).get(name)

    gates: dict[str, object] = {}

    # G1: does the per-target head survive a within-target scaffold-cold split?
    g1 = {
        stratum: get("scaffold_disjoint", stratum, "target_head")
        for stratum in (*SIMILARITY_NAMES, "all")
    }
    g1 = {key: value for key, value in g1.items() if value}
    gates["G1"] = {
        "records": g1,
        "pass": bool(g1.get("all", {}).get("lower95", -1.0) > MDE),
        "low_similarity_pass": bool(g1.get("s00_30", {}).get("lower95", -1.0) > MDE),
        "criterion": "per-target head beats base on scaffold-disjoint within-target evaluation",
    }

    # G2: how many source directions are needed?
    variance = subspace.spectrum**2
    fraction = float(variance[:3].sum() / variance.sum())
    ranks = {
        f"rank{rank}": get("scaffold_disjoint", "all", f"oracle_rank{rank}") for rank in RANK_SWEEP
    }
    ranks = {key: value for key, value in ranks.items() if value}
    full = ranks.get("rank26", {}).get("mean", float("nan"))
    gates["G2"] = {
        "head_variance_in_top3": fraction,
        "ranks": ranks,
        "rank2_retained_fraction": (
            float(ranks.get("rank2", {}).get("mean", float("nan")) / full) if full and full == full else None
        ),
        "criterion": "fraction of the full-head gain retained by a low-rank source projection",
    }

    # G3: is the object predictable from protein alone, with no labels?
    zero_shot = get("scaffold_disjoint", "all", "protein_zero_shot")
    gates["G3"] = {
        "records": zero_shot,
        "pass": bool(zero_shot and zero_shot.get("lower95", -1.0) > MDE),
        "criterion": "protein-predicted head beats base with no support labels",
    }

    # G4: the label learning curve.
    curve: dict[str, object] = {}
    for k in K_SWEEP:
        for rank in RANK_SWEEP:
            value = get("scaffold_disjoint", "all", f"eb_k{k}_rank{rank}")
            if value:
                curve[f"k{k}_rank{rank}"] = value
    best_small = max(
        (value["lower95"] for key, value in curve.items() if key.startswith(("k1_", "k3_", "k5_"))),
        default=-1.0,
    )
    gates["G4"] = {
        "curve": curve,
        "best_lower95_at_k_le_5": best_small,
        "pass": bool(best_small > MDE),
        "criterion": "an empirical-Bayes head from k<=5 labels beats the frozen base",
    }

    if gates["G1"]["pass"] and gates["G4"]["pass"]:
        verdict = "GENERALIZABLE_AND_FEW_SHOT_REACHABLE"
    elif gates["G1"]["pass"]:
        verdict = "GENERALIZABLE_BUT_NOT_YET_FEW_SHOT_REACHABLE"
    else:
        verdict = "NOT_GENERALIZABLE_OBJECT_IS_LOCAL"
    return {"gates": gates, "verdict": verdict}


# --------------------------------------------------------------------------- #
# Runner
# --------------------------------------------------------------------------- #


def run(lock_path: Path, output: Path, records_path: Path, oof_cache: Path) -> dict[str, object]:
    if DEVICE.type != "cuda":
        raise RuntimeError("run these gates with D:\\anaconda\\envs\\drug\\python.exe")
    substrate, context = load_substrate(lock_path, oof_cache)
    basis, basis_stats = build_basis(substrate)
    heads, names, sigma, level_sd = source_heads(substrate, basis)
    subspace = build_subspace(heads, sigma)
    protein_head = protein_predicted_heads(substrate, basis, heads, names)
    splits = target_splits(substrate, "probe")
    records = evaluate(substrate, basis, subspace, protein_head, splits)
    records_path.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    summary = summarise(records)
    verdict = decide(summary, subspace)

    result: dict[str, object] = {
        "schema": "a2s-mode-generalization-v1",
        "status": "SOURCE_ONLY_DIAGNOSTIC",
        "question": "G1-G4: is the target adaptation object global, low-dimensional, protein-predictable, and few-shot reachable?",
        "labels": {
            "opened_roles": ["fit", "probe"],
            "locked_labels_requested": False,
            "recipient_labels_requested": False,
        },
        "device": torch.cuda.get_device_name(0),
        "lock": {"path": str(lock_path), "content_sha256": context["lock"]["content_sha256"]},
        "inputs": {
            "registry_sha256": sha256_file(REGISTRY),
            "features_sha256": sha256_file(FEATURES),
            "proteins_sha256": sha256_file(PROTEINS),
        },
        "protocol": {
            "basis": basis_stats,
            "source_targets": subspace.n_targets,
            "sigma": subspace.sigma,
            "level_sd": level_sd,
            "head_spectrum": subspace.spectrum.tolist(),
            "head_variance_fraction": (subspace.spectrum**2 / (subspace.spectrum**2).sum()).tolist(),
            "k_sweep": list(K_SWEEP),
            "rank_sweep": list(RANK_SWEEP),
            "draws": DRAWS,
            "similarity_edges": list(SIMILARITY_EDGES),
            "mde": MDE,
            "prior": "empirical Bayes; lambda_j = sigma^2 / tau_j^2 measured on source targets, never tuned on probe",
        },
        "splits": {
            "targets": int(len({split.target for split in splits})),
            "scaffold_disjoint": int(sum(1 for split in splits if split.split == "scaffold_disjoint")),
            "random": int(sum(1 for split in splits if split.split == "random")),
        },
        "summary": summary,
        "decision": verdict,
    }
    result["content_sha256"] = sha256(canonical(result).encode("utf-8")).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True, default=float) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--oof-cache", type=Path, default=DEFAULT_OOF)
    args = parser.parse_args()
    result = run(
        args.lock.resolve(), args.out.resolve(), args.records.resolve(), args.oof_cache.resolve()
    )
    print(
        json.dumps(
            {
                "verdict": result["decision"]["verdict"],
                "gates": {
                    name: gate.get("pass") for name, gate in result["decision"]["gates"].items()
                },
                "splits": result["splits"],
            },
            indent=2,
            sort_keys=True,
            default=float,
        )
    )


if __name__ == "__main__":
    main()
