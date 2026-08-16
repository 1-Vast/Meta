"""A2S-RIP Gate R0: can a certified subset of ranking interventions beat the wholesale head?

Gate G4 measured that at ``k <= 5`` the empirical-Bayes target head is estimable
but too noisy to apply wholesale (+0.011 CI, interval crossing zero; the knee of
the learning curve sits at ``k ~ 10``).  A2S-RIP's premise is that the noisy head
should be applied only where it is certifiably right.

This gate measures the ceiling of that idea before any policy is trained.  It
adds one object to the G4 machinery: the **posterior covariance** of the head, so
each compound carries an honest per-compound uncertainty and therefore an
evidence margin ``z_q``.

* ``R0a`` selection ceiling - risk-coverage curve under a greedy *oracle* choice
  of which compounds to intervene on.  If its peak sits at or below the wholesale
  head, selection buys nothing and A2S-RIP dies here.
* ``R0b`` is the margin informative at all - AUC of ``|z_q|`` for predicting
  whether the proposed edit points the right way, on unseen components.  Chance
  is 0.5.  This is the existence test for the transferable object and it needs no
  training.  Compared against ``|delta_q|`` alone, which tests whether the
  posterior covariance is load-bearing.
* ``R0c`` does the threshold transfer - a harm-rate threshold fitted on ``fit``
  targets, applied unchanged to ``probe`` targets.  This pre-tests the
  meta-conformal premise without a model.
* ``R0d`` the magnitude confound - every selective curve is compared against a
  wholesale head rescaled to the *same mean edit magnitude*.  Correction
  magnitude alone was measured to be worth +0.009 CI, so without this control a
  selection gain is a shrinkage artefact.

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
from research.a2s.a2s_trace import DEVICE, Substrate, load_substrate
from research.a2s.a2s_trace_stratum import DEFAULT_LOCK, DEFAULT_OOF, metric_loss, paired_bootstrap
from research.a2s.a2s_mode_gates import build_basis, source_heads
from research.a2s.a2s_mode_generalization import (
    SourceSubspace,
    TargetSplit,
    build_subspace,
    target_splits,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_rip_gate_r0_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "a2s_rip_gate_r0_records_2026-08-02.parquet"

SEED = 20260802
K_SWEEP = (3, 5, 10)
DRAWS = 8
COVERAGE_GRID = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
RULES = ("oracle", "margin", "abs_delta", "random")
HARM_ALPHAS = (0.20, 0.30, 0.40)
MIN_EVAL_ROWS = 10
BOOTSTRAP_DRAWS = 2000
MDE = 0.005


@dataclass
class Posterior:
    mean: np.ndarray            # (d,) posterior mean head
    covariance: np.ndarray      # (d, d) posterior covariance of the head


def head_posterior(
    subspace: SourceSubspace, design: np.ndarray, residual: np.ndarray
) -> Posterior:
    """Closed-form empirical-Bayes posterior over the target head.

    Full rank by construction: Gate G2 measured that projecting onto the dominant
    source directions destroys the signal, so nothing is compressed here.  The
    prior is the measured source dispersion along each direction, so there is no
    free parameter to tune on the recipient.
    """

    directions = subspace.directions
    variances = np.maximum(subspace.variances, 1e-8)
    sigma2 = subspace.sigma**2
    centre = design.mean(axis=0)
    centred = design - centre
    projected = centred @ directions
    target = residual - residual.mean() - centred @ subspace.mean_head
    precision = projected.T @ projected / sigma2 + np.diag(1.0 / variances)
    covariance_c = np.linalg.inv(precision)
    coefficients = covariance_c @ projected.T @ target / sigma2
    return Posterior(
        mean=subspace.mean_head + directions @ coefficients,
        covariance=directions @ covariance_c @ directions.T,
    )


def marginal_concordance_gain(
    base: np.ndarray, labels: np.ndarray, delta: np.ndarray
) -> np.ndarray:
    """Change in this compound's own pairwise concordance if it alone is edited.

    Used only to define the hindsight oracle selection order.  The reported CI is
    always recomputed after applying the whole selected set, so the number is
    real; only the *choice* is oracular.
    """

    gaps_before = base[:, None] - base[None, :]
    gaps_after = (base + delta)[:, None] - base[None, :]
    truth = np.sign(labels[:, None] - labels[None, :])
    active = truth != 0
    np.fill_diagonal(active, False)

    def score(gaps: np.ndarray) -> np.ndarray:
        agree = (np.sign(gaps) == truth).astype(np.float64)
        tied = (np.sign(gaps) == 0).astype(np.float64)
        return ((agree + 0.5 * tied) * active).sum(axis=1)

    return score(gaps_after) - score(gaps_before)


def select(rule: str, coverage: float, statistics: dict[str, np.ndarray], rng: np.random.Generator) -> np.ndarray:
    count = int(round(coverage * len(statistics["margin"])))
    chosen = np.zeros(len(statistics["margin"]), dtype=bool)
    if count <= 0:
        return chosen
    if rule == "random":
        chosen[rng.choice(len(chosen), size=count, replace=False)] = True
        return chosen
    order = np.argsort(-statistics[rule], kind="stable")[:count]
    chosen[order] = True
    return chosen


def evaluate_split(
    split: TargetSplit,
    role: str,
    basis: np.ndarray,
    subspace: SourceSubspace,
    residual_all: np.ndarray,
    affinity_all: np.ndarray,
    base_all: np.ndarray,
) -> list[dict[str, object]]:
    train_design = basis[split.train_rows]
    train_residual = residual_all[split.train_rows]
    eval_design = basis[split.eval_rows]
    labels = affinity_all[split.eval_rows]
    base = base_all[split.eval_rows]
    if len(labels) < MIN_EVAL_ROWS or float(np.std(labels)) < 1e-9:
        return []
    # Rank-null reference: the centred true residual is what an edit should match.
    truth_direction = np.sign(residual_all[split.eval_rows] - residual_all[split.eval_rows].mean())

    records: list[dict[str, object]] = []
    for k in K_SWEEP:
        if k > len(split.train_rows):
            continue
        for draw in range(DRAWS):
            digest = int(
                sha256(f"{SEED}:{split.target}:{split.split}:{k}:{draw}".encode()).hexdigest()[:8], 16
            )
            rng = np.random.default_rng(digest)
            chosen_rows = rng.choice(len(split.train_rows), size=k, replace=False)
            posterior = head_posterior(
                subspace, train_design[chosen_rows], train_residual[chosen_rows]
            )
            raw = eval_design @ posterior.mean
            # Support-centred, so the rule stays a per-compound decision: no
            # candidate-set statistic may enter a deployable score.
            support_centre = float((train_design[chosen_rows] @ posterior.mean).mean())
            delta = raw - support_centre
            variance = np.einsum("ij,jk,ik->i", eval_design, posterior.covariance, eval_design)
            margin = np.abs(delta) / np.sqrt(np.maximum(variance, 0.0) + subspace.sigma**2)
            correct = (np.sign(delta) == truth_direction) & (truth_direction != 0)
            defined = truth_direction != 0

            statistics = {
                "margin": margin,
                "abs_delta": np.abs(delta),
                "oracle": marginal_concordance_gain(base, labels, delta),
                "random": np.zeros_like(margin),
            }
            wholesale = metric_loss(labels, base + delta)
            base_metrics = metric_loss(labels, base)
            mean_abs_full = float(np.mean(np.abs(delta)))

            for rule in RULES:
                for coverage in COVERAGE_GRID:
                    picked = select(rule, coverage, statistics, rng)
                    applied = np.where(picked, delta, 0.0)
                    metrics = metric_loss(labels, base + applied)
                    mean_abs = float(np.mean(np.abs(applied)))
                    # R0d: a wholesale head rescaled to the same mean edit
                    # magnitude.  If this matches the selective curve, the effect
                    # is the global magnitude scalar, not selection.
                    scale = mean_abs / mean_abs_full if mean_abs_full > 1e-12 else 0.0
                    matched = metric_loss(labels, base + scale * delta)
                    harm = (
                        float((~correct[picked & defined]).mean())
                        if int((picked & defined).sum()) > 0
                        else float("nan")
                    )
                    records.append(
                        {
                            "target": split.target,
                            "component": split.component,
                            "split": split.split,
                            "role": role,
                            "k": k,
                            "draw": draw,
                            "rule": rule,
                            "coverage": coverage,
                            "n_eval": int(len(labels)),
                            "mean_abs_delta": mean_abs,
                            "harm_rate": harm,
                            "base_ci": base_metrics["ci"],
                            "wholesale_ci": wholesale["ci"],
                            "selective_ci": metrics["ci"],
                            "matched_ci": matched["ci"],
                            "selective_ndcg10": metrics["ndcg10"],
                            "wholesale_ndcg10": wholesale["ndcg10"],
                            "base_ndcg10": base_metrics["ndcg10"],
                            "margin_auc": float("nan"),
                        }
                    )
            # R0b: is the margin informative at all?  One row per (k, draw).
            for name, statistic in (("margin", margin), ("abs_delta", np.abs(delta))):
                auc = rank_auc(statistic[defined], correct[defined])
                records.append(
                    {
                        "target": split.target,
                        "component": split.component,
                        "split": split.split,
                        "role": role,
                        "k": k,
                        "draw": draw,
                        "rule": f"auc_{name}",
                        "coverage": float("nan"),
                        "n_eval": int(defined.sum()),
                        "mean_abs_delta": mean_abs_full,
                        "harm_rate": float((~correct[defined]).mean()) if int(defined.sum()) else float("nan"),
                        "base_ci": base_metrics["ci"],
                        "wholesale_ci": wholesale["ci"],
                        "selective_ci": float("nan"),
                        "matched_ci": float("nan"),
                        "selective_ndcg10": float("nan"),
                        "wholesale_ndcg10": wholesale["ndcg10"],
                        "base_ndcg10": base_metrics["ndcg10"],
                        "margin_auc": auc,
                    }
                )
    return records


def rank_auc(scores: np.ndarray, positive: np.ndarray) -> float:
    """AUC of `scores` for separating correct edits from incorrect ones."""

    positive = positive.astype(bool)
    if positive.all() or (~positive).all() or len(scores) < 2:
        return float("nan")
    order = np.argsort(scores, kind="stable")
    ranks = np.empty(len(scores), dtype=np.float64)
    ranks[order] = np.arange(1, len(scores) + 1)
    n_pos = int(positive.sum())
    n_neg = int((~positive).sum())
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def collect(substrate: Substrate, basis: np.ndarray, subspace: SourceSubspace) -> pd.DataFrame:
    residual_all = (substrate.affinity - substrate.base).cpu().numpy().astype(np.float64)
    affinity_all = substrate.affinity.cpu().numpy().astype(np.float64)
    base_all = substrate.base.cpu().numpy().astype(np.float64)
    records: list[dict[str, object]] = []
    for role in ("fit", "probe"):
        for split in target_splits(substrate, role):
            if split.split != "scaffold_disjoint":
                continue
            records.extend(
                evaluate_split(split, role, basis, subspace, residual_all, affinity_all, base_all)
            )
    return pd.DataFrame.from_records(records)


def summarise(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    curves = records.loc[records.rule.isin(RULES)]
    for role in sorted(curves.role.unique()):
        for k in sorted(curves.k.unique()):
            for rule in RULES:
                frame = curves.loc[
                    (curves.role == role) & (curves.k == k) & (curves.rule == rule)
                ]
                if frame.empty:
                    continue
                points: dict[str, object] = {}
                for coverage in COVERAGE_GRID:
                    cell = frame.loc[frame.coverage == coverage].copy()
                    if cell.empty:
                        continue
                    cell["value"] = cell.selective_ci - cell.base_ci
                    versus_base = paired_bootstrap(cell, "value", draws=BOOTSTRAP_DRAWS)
                    cell["value"] = cell.selective_ci - cell.wholesale_ci
                    versus_wholesale = paired_bootstrap(cell, "value", draws=BOOTSTRAP_DRAWS)
                    cell["value"] = cell.selective_ci - cell.matched_ci
                    versus_matched = paired_bootstrap(cell, "value", draws=BOOTSTRAP_DRAWS)
                    points[f"c{coverage:.1f}"] = {
                        "vs_base": versus_base,
                        "vs_wholesale": versus_wholesale,
                        "vs_magnitude_matched": versus_matched,
                        "mean_harm_rate": float(np.nanmean(cell.harm_rate.to_numpy())),
                        "mean_abs_delta": float(cell.mean_abs_delta.mean()),
                    }
                summary.setdefault(role, {}).setdefault(f"k{int(k)}", {})[rule] = points

    auc_rows = records.loc[records.rule.str.startswith("auc_")]
    for role in sorted(auc_rows.role.unique()):
        for k in sorted(auc_rows.k.unique()):
            entry: dict[str, object] = {}
            for name in ("auc_margin", "auc_abs_delta"):
                cell = auc_rows.loc[
                    (auc_rows.role == role) & (auc_rows.k == k) & (auc_rows.rule == name)
                ].copy()
                if cell.empty:
                    continue
                cell["value"] = cell.margin_auc
                entry[name] = paired_bootstrap(cell.dropna(subset=["value"]), "value", draws=BOOTSTRAP_DRAWS)
                entry[f"{name}_harm_rate"] = float(np.nanmean(cell.harm_rate.to_numpy()))
            summary.setdefault(role, {}).setdefault(f"k{int(k)}", {})["auc"] = entry
    return summary


def threshold_transfer(records: pd.DataFrame) -> dict[str, object]:
    """R0c: fit a margin threshold on `fit` targets, apply it unchanged to `probe`."""

    curves = records.loc[records.rule == "margin"]
    output: dict[str, object] = {}
    for k in sorted(curves.k.unique()):
        fit_frame = curves.loc[(curves.role == "fit") & (curves.k == k)]
        probe_frame = curves.loc[(curves.role == "probe") & (curves.k == k)]
        if fit_frame.empty or probe_frame.empty:
            continue
        entry: dict[str, object] = {}
        for alpha in HARM_ALPHAS:
            eligible = [
                coverage
                for coverage in COVERAGE_GRID
                if coverage > 0.0
                and float(np.nanmean(fit_frame.loc[fit_frame.coverage == coverage].harm_rate.to_numpy())) <= alpha
            ]
            if not eligible:
                entry[f"alpha{alpha:.2f}"] = {"selected_coverage": None, "note": "no coverage met the harm budget on fit"}
                continue
            coverage = max(eligible)
            probe_cell = probe_frame.loc[probe_frame.coverage == coverage].copy()
            probe_cell["value"] = probe_cell.selective_ci - probe_cell.base_ci
            entry[f"alpha{alpha:.2f}"] = {
                "selected_coverage": coverage,
                "fit_harm_rate": float(np.nanmean(fit_frame.loc[fit_frame.coverage == coverage].harm_rate.to_numpy())),
                "probe_harm_rate": float(np.nanmean(probe_cell.harm_rate.to_numpy())),
                "probe_gain_vs_base": paired_bootstrap(probe_cell, "value", draws=BOOTSTRAP_DRAWS),
                "valid": bool(float(np.nanmean(probe_cell.harm_rate.to_numpy())) <= alpha + 0.02),
            }
        output[f"k{int(k)}"] = entry
    return output


def decide(summary: dict[str, object], transfer: dict[str, object]) -> dict[str, object]:
    gates: dict[str, object] = {}

    # R0a: does an oracle-selected subset beat the wholesale head?
    peaks: list[dict[str, object]] = []
    for k in ("k3", "k5"):
        points = summary.get("probe", {}).get(k, {}).get("oracle", {})
        best = None
        for label, value in points.items():
            candidate = value["vs_wholesale"]
            if best is None or candidate["mean"] > best["mean"]:
                best = {"coverage": label, **candidate}
        if best:
            peaks.append({"k": k, **best})
    gates["R0a"] = {
        "oracle_peak_vs_wholesale": peaks,
        "pass": any(item["lower95"] > MDE for item in peaks),
        "criterion": "a hindsight-selected subset beats applying the whole head",
    }

    # R0b: is the margin informative on unseen components?
    aucs: list[dict[str, object]] = []
    for k in ("k3", "k5"):
        entry = summary.get("probe", {}).get(k, {}).get("auc", {})
        if "auc_margin" in entry:
            aucs.append({"k": k, "statistic": "margin", **entry["auc_margin"]})
        if "auc_abs_delta" in entry:
            aucs.append({"k": k, "statistic": "abs_delta", **entry["auc_abs_delta"]})
    gates["R0b"] = {
        "records": aucs,
        "pass": any(
            item["statistic"] == "margin" and item["lower95"] > 0.5 for item in aucs
        ),
        "posterior_load_bearing": any(
            item["statistic"] == "margin" and item["mean"] > 0.0 for item in aucs
        ),
        "criterion": "AUC of the evidence margin for edit correctness exceeds chance on unseen components",
    }

    # R0c: does a fit-fitted threshold stay valid on probe?
    gates["R0c"] = {
        "records": transfer,
        "pass": any(
            isinstance(cell, dict) and cell.get("valid")
            for by_k in transfer.values()
            for cell in by_k.values()
        ),
        "criterion": "a harm-rate threshold fitted on fit targets remains valid on probe",
    }

    # R0d: is any selective gain distinguishable from a magnitude rescale?
    matched: list[dict[str, object]] = []
    for k in ("k3", "k5"):
        for rule in ("oracle", "margin"):
            points = summary.get("probe", {}).get(k, {}).get(rule, {})
            best = None
            for label, value in points.items():
                candidate = value["vs_magnitude_matched"]
                if best is None or candidate["mean"] > best["mean"]:
                    best = {"coverage": label, **candidate}
            if best:
                matched.append({"k": k, "rule": rule, **best})
    gates["R0d"] = {
        "records": matched,
        "pass": any(
            item["rule"] == "margin" and item["lower95"] > MDE for item in matched
        ),
        "oracle_pass": any(
            item["rule"] == "oracle" and item["lower95"] > MDE for item in matched
        ),
        "criterion": "selection beats a wholesale head rescaled to the same mean edit magnitude",
    }

    if gates["R0a"]["pass"] and gates["R0b"]["pass"] and gates["R0d"]["oracle_pass"]:
        verdict = "RIP_CEILING_ADMITTED"
    elif gates["R0a"]["pass"]:
        verdict = "CEILING_EXISTS_BUT_NOT_REACHABLE_BY_AN_OBSERVABLE_MARGIN"
    else:
        verdict = "NO_SELECTION_CEILING"
    return {"gates": gates, "verdict": verdict}


def run(lock_path: Path, output: Path, records_path: Path, oof_cache: Path) -> dict[str, object]:
    if DEVICE.type != "cuda":
        raise RuntimeError("run this gate with D:\\anaconda\\envs\\drug\\python.exe")
    substrate, context = load_substrate(lock_path, oof_cache)
    basis, basis_stats = build_basis(substrate)
    heads, names, sigma, level_sd = source_heads(substrate, basis)
    subspace = build_subspace(heads, sigma)
    records = collect(substrate, basis, subspace)
    records_path.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    summary = summarise(records)
    transfer = threshold_transfer(records)
    verdict = decide(summary, transfer)

    result: dict[str, object] = {
        "schema": "a2s-rip-gate-r0-v1",
        "status": "SOURCE_ONLY_DIAGNOSTIC",
        "question": "R0: is there a selection ceiling above the wholesale head, and is it reachable from an observable margin?",
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
            "k_sweep": list(K_SWEEP),
            "draws": DRAWS,
            "coverage_grid": list(COVERAGE_GRID),
            "rules": list(RULES),
            "harm_alphas": list(HARM_ALPHAS),
            "split": "within-target Murcko-scaffold-disjoint",
            "head": "full-rank empirical-Bayes posterior; prior from source head dispersion, never tuned on probe",
            "score_centring": "support-set mean, so the rule stays a per-compound decision",
            "mde": MDE,
        },
        "summary": summary,
        "threshold_transfer": transfer,
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
                "gates": {name: gate.get("pass") for name, gate in result["decision"]["gates"].items()},
            },
            indent=2,
            sort_keys=True,
            default=float,
        )
    )


if __name__ == "__main__":
    main()
