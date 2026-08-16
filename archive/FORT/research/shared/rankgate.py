"""TRAIN-only gate for support-anchored protein-conditioned reordering."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import threading
import time

import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F

from research.shared.pairprior import PairPrior
from research.shared.priorgate import (
    LIGANDS,
    PROTEINS,
    REGISTRY,
    descriptorstats,
    filesha,
    limitrows,
    loadarrays,
    loadframe,
    ridgebase,
    splitcomponents,
    telemetry,
    telemetrysummary,
    wrongtargets,
)
from scripts.metric import evaluateprotocol, pairedcomponents
from scripts.train import contrastloss, maketrainroster


REPORT = Path("reports/active/rankgate.v1.json")


def jsonsafe(value):
    if isinstance(value, dict):
        return {key: jsonsafe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [jsonsafe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def proteins(
    frame: pd.DataFrame,
    targetindex: dict[str, int],
    segments: torch.Tensor,
) -> dict[str, torch.Tensor]:
    return {
        target: segments[targetindex[target]]
        for target in sorted(frame.target.unique())
        if target in targetindex
    }


def tensors(
    frame: pd.DataFrame,
    feature: np.ndarray,
    base: np.ndarray,
    support: tuple[int, ...],
    query: tuple[int, ...],
    center: torch.Tensor,
    scale: torch.Tensor,
) -> tuple[torch.Tensor, ...]:
    def ligand(indices: tuple[int, ...]) -> torch.Tensor:
        source = frame.source_row.iloc[list(indices)].to_numpy(dtype=np.int64)
        value = torch.as_tensor(feature[source], device="cuda", dtype=torch.float32)
        value[:, -10:] = (value[:, -10:] - center) / scale
        return value

    def values(indices: tuple[int, ...], column: str) -> torch.Tensor:
        if column == "base":
            source = frame.source_row.iloc[list(indices)].to_numpy(dtype=np.int64)
            value = base[source]
        else:
            value = frame[column].iloc[list(indices)].to_numpy(dtype=np.float32)
        return torch.as_tensor(value, device="cuda", dtype=torch.float32)

    return (
        ligand(support),
        values(support, "affinity"),
        values(support, "base"),
        ligand(query),
        values(query, "affinity"),
        values(query, "base"),
    )


def anchored(
    model: PairPrior,
    protein: torch.Tensor,
    support: torch.Tensor,
    query: torch.Tensor,
    supportlabel: torch.Tensor,
    supportbase: torch.Tensor,
    querybase: torch.Tensor,
) -> torch.Tensor:
    ligand = torch.cat((support, query))
    repeated = protein.unsqueeze(0).expand(len(ligand), -1, -1)
    score = model.predict(repeated, ligand)["prediction"]
    supportscore = score[: len(support)]
    queryscore = score[len(support) :]
    calibration = (supportlabel - supportbase).mean()
    return querybase + calibration + (queryscore - supportscore.mean())


def episodes(
    frame: pd.DataFrame,
    protein: dict[str, torch.Tensor],
    feature: np.ndarray,
    targets: int,
    queries: int,
) -> tuple[pd.DataFrame, list]:
    local = frame.reset_index(drop=True).copy()
    roster = maketrainroster(
        local,
        protein,
        feature,
        targets=targets,
        querycap=queries,
        support=5,
    )
    return local, roster


def fit(
    model: PairPrior,
    frame: pd.DataFrame,
    roster: list,
    protein: dict[str, torch.Tensor],
    feature: np.ndarray,
    base: np.ndarray,
    center: torch.Tensor,
    scale: torch.Tensor,
    epochs: int,
    rate: float,
    rankweight: float,
    specificityweight: float,
    seed: int,
) -> list[dict[str, float]]:
    optimizer = torch.optim.AdamW(model.parameters(), lr=rate, weight_decay=1e-4)
    wrong = wrongtargets(frame, seed + 11)
    history = []
    generator = np.random.default_rng(seed)
    for epoch in range(epochs):
        model.train()
        totals = torch.zeros(4, device="cuda")
        order = generator.permutation(len(roster))
        for position in order:
            item = roster[int(position)]
            sx, sy, sb, qx, qy, qb = tensors(
                frame,
                feature,
                base,
                item.support_indices,
                item.query_indices,
                center,
                scale,
            )
            correct = anchored(
                model, protein[item.target_key], sx, qx, sy, sb, qb
            )
            wrongprediction = anchored(
                model, protein[wrong[item.target_key]], sx, qx, sy, sb, qb
            )
            absolute = F.huber_loss(correct, qy)
            ranking = contrastloss(correct, qy)
            wrongloss = F.huber_loss(wrongprediction, qy)
            specificity = F.softplus(0.05 + absolute - wrongloss)
            loss = absolute + rankweight * ranking + specificityweight * specificity
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            totals += torch.stack(
                (loss.detach(), absolute.detach(), ranking.detach(), specificity.detach())
            )
        means = totals / len(roster)
        record = {
            "epoch": epoch + 1,
            "loss": float(means[0].cpu()),
            "absolute": float(means[1].cpu()),
            "ranking": float(means[2].cpu()),
            "specificity": float(means[3].cpu()),
        }
        history.append(record)
        print(
            f"epoch {epoch + 1}/{epochs} loss={record['loss']:.4f} "
            f"absolute={record['absolute']:.4f} ranking={record['ranking']:.4f}",
            flush=True,
        )
    return history


@torch.no_grad()
def evaluate(
    model: PairPrior,
    frame: pd.DataFrame,
    roster: list,
    protein: dict[str, torch.Tensor],
    feature: np.ndarray,
    base: np.ndarray,
    center: torch.Tensor,
    scale: torch.Tensor,
    seed: int,
) -> dict[str, object]:
    model.eval()
    wrong = wrongtargets(frame, seed + 19)
    prediction = {name: [] for name in ("correct", "calibration", "wrongprotein")}
    labels: list[float] = []
    indices: list[int] = []
    for item in roster:
        sx, sy, sb, qx, qy, qb = tensors(
            frame,
            feature,
            base,
            item.support_indices,
            item.query_indices,
            center,
            scale,
        )
        correct = anchored(model, protein[item.target_key], sx, qx, sy, sb, qb)
        wrongvalue = anchored(
            model, protein[wrong[item.target_key]], sx, qx, sy, sb, qb
        )
        calibration = qb + (sy - sb).mean()
        prediction["correct"].extend(correct.cpu().tolist())
        prediction["calibration"].extend(calibration.cpu().tolist())
        prediction["wrongprotein"].extend(wrongvalue.cpu().tolist())
        labels.extend(qy.cpu().tolist())
        indices.extend(item.query_indices)
    components = {item.target_key: item.homology_component for item in roster}
    metrics = {
        name: evaluateprotocol(
            predictions=value,
            labels=labels,
            episodes=roster,
            prediction_indices=indices,
            component_by_target=components,
        )
        for name, value in prediction.items()
    }
    paired = pairedcomponents(
        predictions=prediction,
        labels=labels,
        episodes=roster,
        prediction_indices=indices,
        component_by_target=components,
        reference="correct",
        seed=seed,
    )
    required = []
    for control in ("calibration", "wrongprotein"):
        for metric in ("rmse_gain", "spearman_gain", "pairwise_gain"):
            required.append(paired[control][metric]["ci95"][0] > 0.0)
    return {
        "gate": "PASS" if all(required) else "STOP",
        "metrics": metrics,
        "paired_component_bootstrap": paired,
        "episodes": len(roster),
        "components": len(set(components.values())),
        "queries": len(labels),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=8)
    parser.add_argument("--width", type=int, default=192)
    parser.add_argument("--fit-targets", type=int, default=0)
    parser.add_argument("--gate-targets", type=int, default=0)
    parser.add_argument("--queries", type=int, default=64)
    parser.add_argument("--rows-per-target", type=int, default=0)
    parser.add_argument("--learning-rate", type=float, default=2e-4)
    parser.add_argument("--rank-weight", type=float, default=0.5)
    parser.add_argument("--specificity-weight", type=float, default=0.15)
    parser.add_argument("--ridge", type=float, default=10.0)
    parser.add_argument("--holdout", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--out", type=Path, default=REPORT)
    args = parser.parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError("the support-anchored rank gate requires CUDA")
    if args.epochs < 1 or args.queries < 1:
        raise ValueError("epochs and queries must be positive")
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    torch.set_float32_matmul_precision("high")
    torch.cuda.reset_peak_memory_stats()
    started = time.perf_counter()
    full = loadframe("pKi")
    feature, targetindex, segments = loadarrays(full)
    fitframe, gateframe = splitcomponents(full, args.holdout, args.seed)
    fitframe = limitrows(
        fitframe, args.fit_targets, args.rows_per_target, args.seed + 1
    )
    gateframe = limitrows(
        gateframe, args.gate_targets, args.rows_per_target, args.seed + 2
    )
    base, baseinfo = ridgebase(
        fitframe, gateframe, feature, args.ridge, batchsize=512
    )
    center, scale = descriptorstats(fitframe, feature, batchsize=512)
    targetprotein = proteins(
        pd.concat((fitframe, gateframe)), targetindex, segments
    )
    fitframe, fitroster = episodes(
        fitframe, targetprotein, feature, args.fit_targets, args.queries
    )
    gateframe, gateroster = episodes(
        gateframe, targetprotein, feature, args.gate_targets, args.queries
    )
    model = PairPrior(width=args.width).cuda()
    with torch.no_grad():
        model.head[-1].weight[0].mul_(1e-3)
        model.head[-1].bias[0].zero_()
    samples: list[dict[str, float]] = []
    stop = threading.Event()
    monitor = threading.Thread(target=telemetry, args=(stop, samples), daemon=True)
    monitor.start()
    try:
        history = fit(
            model,
            fitframe,
            fitroster,
            targetprotein,
            feature,
            base,
            center,
            scale,
            args.epochs,
            args.learning_rate,
            args.rank_weight,
            args.specificity_weight,
            args.seed,
        )
        result = evaluate(
            model,
            gateframe,
            gateroster,
            targetprotein,
            feature,
            base,
            center,
            scale,
            args.seed,
        )
    finally:
        stop.set()
        monitor.join(timeout=2.0)
    result |= {
        "protocol": "TRAIN-only support-anchored protein reordering gate",
        "roles_read": ["train"],
        "support": 5,
        "seed": args.seed,
        "fit_episodes": len(fitroster),
        "gate_episodes": len(gateroster),
        "epochs": args.epochs,
        "queries_per_episode": args.queries,
        "parameters": sum(parameter.numel() for parameter in model.parameters()),
        "base": baseinfo,
        "history": history,
        "training_seconds": time.perf_counter() - started,
        "gpu": torch.cuda.get_device_name(),
        "gpu_telemetry": telemetrysummary(
            samples, torch.cuda.max_memory_allocated() / (1024 ** 2)
        ),
        "sources": {
            "registry_sha256": filesha(REGISTRY),
            "ligand_sha256": filesha(LIGANDS),
            "protein_sha256": filesha(PROTEINS),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    serialized = jsonsafe(result)
    args.out.write_text(
        json.dumps(serialized, indent=2, allow_nan=False), encoding="utf-8"
    )
    torch.save(model.state_dict(), args.out.with_suffix(".pt"))
    print(
        json.dumps(
            {
                "gate": serialized["gate"],
                "metrics": serialized["metrics"],
                "paired": serialized["paired_component_bootstrap"],
                "out": str(args.out),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
