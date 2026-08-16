"""TRAIN-only gate for a sequence-conditioned global pair prior."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import threading
import time
from typing import Iterator

import numpy as np
import pandas as pd
import torch

from .pairprior import PairPrior, evidencecontrast, gaussiannll
from scripts.audit import FIELDS


ROOT = Path("dataset/public/chembl_37/processed/dualcold")
REGISTRY = ROOT / "registry.parquet"
LIGANDS = ROOT / "ligand_features.npz"
PROTEINS = ROOT / "target_esm2.npz"
REPORT = Path("reports/active/priorgate.v1.json")


def stablekey(value: object, seed: int) -> str:
    return hashlib.sha256(f"{seed}:{value}".encode("utf-8")).hexdigest()


def splitcomponents(
    frame: pd.DataFrame, holdout: float, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Freeze a homology-component split using metadata only."""

    if not 0.0 < holdout < 0.5:
        raise ValueError("holdout fraction must lie in (0, 0.5)")
    components = sorted(
        frame.hcluster.unique(), key=lambda value: stablekey(value, seed)
    )
    count = max(1, int(math.ceil(len(components) * holdout)))
    gatecomponents = set(components[:count])
    gate = frame[frame.hcluster.isin(gatecomponents)].copy()
    fit = frame[~frame.hcluster.isin(gatecomponents)].copy()
    if fit.empty or gate.empty:
        raise RuntimeError("component split produced an empty role")
    if set(fit.hcluster).intersection(gate.hcluster):
        raise RuntimeError("homology components crossed the TRAIN gate")
    return fit, gate


def limitrows(
    frame: pd.DataFrame,
    targets: int,
    rowspertarget: int,
    seed: int,
) -> pd.DataFrame:
    """Deterministic diagnostic limiter; zero preserves the complete role."""

    selected = frame
    if targets:
        names = sorted(frame.target.unique(), key=lambda value: stablekey(value, seed))
        selected = frame[frame.target.isin(names[:targets])]
    if rowspertarget:
        parts = []
        for target, group in selected.groupby("target", sort=True):
            order = sorted(
                group.index,
                key=lambda value: stablekey(f"{target}:{int(value)}", seed),
            )
            parts.append(group.loc[order[:rowspertarget]])
        selected = pd.concat(parts, ignore_index=False)
    return selected.sort_index().copy()


def wrongtargets(
    frame: pd.DataFrame, seed: int
) -> dict[str, str]:
    """Map every target to a deterministic target in another homology component."""

    component = frame.groupby("target", sort=True).hcluster.first().to_dict()
    targets = sorted(component, key=lambda value: stablekey(value, seed))
    if len({component[target] for target in targets}) < 2:
        raise ValueError("wrong-protein control requires at least two components")
    mapping: dict[str, str] = {}
    for position, target in enumerate(targets):
        for offset in range(1, len(targets) + 1):
            candidate = targets[(position + offset) % len(targets)]
            if component[candidate] != component[target]:
                mapping[target] = candidate
                break
    return mapping


def balancedorder(
    frame: pd.DataFrame, seed: int
) -> np.ndarray:
    """Visit every row once in a shuffled target-round-robin order."""

    generator = np.random.default_rng(seed)
    groups: list[np.ndarray] = []
    for _, group in frame.groupby("target", sort=True):
        values = group.index.to_numpy(dtype=np.int64, copy=True)
        generator.shuffle(values)
        groups.append(values)
    rounds: list[np.ndarray] = []
    maximum = max(map(len, groups))
    for position in range(maximum):
        active = np.asarray(
            [values[position] for values in groups if position < len(values)],
            dtype=np.int64,
        )
        generator.shuffle(active)
        rounds.append(active)
    return np.concatenate(rounds)


def batches(order: np.ndarray, batchsize: int) -> Iterator[np.ndarray]:
    if batchsize < 2:
        raise ValueError("batch size must be at least two")
    for start in range(0, len(order), batchsize):
        yield order[start : start + batchsize]


def telemetry(stop: threading.Event, samples: list[dict[str, float]]) -> None:
    while not stop.is_set():
        try:
            line = subprocess.check_output(
                [
                    "nvidia-smi",
                    "--query-gpu=utilization.gpu,power.draw,memory.used",
                    "--format=csv,noheader,nounits",
                ],
                text=True,
            ).strip().splitlines()[0]
            use, power, memory = (float(value.strip()) for value in line.split(","))
            samples.append({"utilization": use, "power": power, "memory": memory})
        except (OSError, subprocess.SubprocessError, ValueError):
            pass
        stop.wait(0.5)


def loadframe(endpoint: str) -> pd.DataFrame:
    """Read affinity values only after the Parquet TRAIN predicate is applied."""

    role = pd.read_parquet(REGISTRY, columns=["dual_cold_split"])
    sourcerows = np.flatnonzero(role.dual_cold_split.to_numpy() == "train")
    frame = pd.read_parquet(
        REGISTRY,
        filters=[("dual_cold_split", "==", "train")],
        columns=[
            *FIELDS,
            "affinity",
        ],
    )
    if len(frame) != len(sourcerows):
        raise RuntimeError("TRAIN row alignment is inconsistent")
    frame["source"] = sourcerows
    frame["source_row"] = sourcerows
    frame = frame[frame.endpoint == endpoint].copy()
    if frame.empty or set(frame.dual_cold_split) != {"train"}:
        raise RuntimeError("the research gate may read TRAIN affinity rows only")
    if not np.isfinite(frame.affinity.to_numpy()).all():
        raise ValueError("TRAIN affinity values must be finite")
    return frame


def loadarrays(
    frame: pd.DataFrame,
) -> tuple[np.ndarray, dict[str, int], torch.Tensor]:
    ligand = np.load(LIGANDS, allow_pickle=False)
    feature = ligand["feat"]
    keep = ligand["keep"]
    frame.drop(frame.index[~keep[frame.source.to_numpy()]], inplace=True)
    protein = np.load(PROTEINS, allow_pickle=False)
    targetindex = {str(key): position for position, key in enumerate(protein["keys"])}
    frame.drop(frame.index[~frame.target.isin(targetindex)], inplace=True)
    segments = torch.as_tensor(protein["segments"], device="cuda")
    return feature, targetindex, segments


def descriptorstats(
    frame: pd.DataFrame, feature: np.ndarray, batchsize: int
) -> tuple[torch.Tensor, torch.Tensor]:
    total = torch.zeros(10, device="cuda", dtype=torch.float64)
    square = torch.zeros_like(total)
    count = 0
    rows = frame.source.to_numpy(dtype=np.int64)
    for start in range(0, len(rows), batchsize * 8):
        values = torch.as_tensor(
            feature[rows[start : start + batchsize * 8], -10:],
            device="cuda",
            dtype=torch.float64,
        )
        total += values.sum(dim=0)
        square += values.square().sum(dim=0)
        count += values.shape[0]
    center = total / count
    variance = (square - count * center.square()) / max(count - 1, 1)
    return center.float(), variance.clamp_min(1e-12).sqrt().float()


def ridgebase(
    fit: pd.DataFrame,
    gate: pd.DataFrame,
    feature: np.ndarray,
    ridge: float,
    batchsize: int,
) -> tuple[np.ndarray, dict[str, object]]:
    """Fit a frozen ligand-only ridge prior on the fitting role using CUDA."""

    if ridge <= 0:
        raise ValueError("ridge must be positive")
    rows = fit.source.to_numpy(dtype=np.int64)
    values = torch.as_tensor(feature[rows], device="cuda", dtype=torch.float32)
    labels = torch.as_tensor(
        fit.affinity.to_numpy(dtype=np.float32), device="cuda"
    )
    center = values.mean(dim=0)
    scale = values.std(dim=0).clamp_min(1e-6)
    standardized = (values - center) / scale
    labelcenter = labels.mean()
    gram = standardized.T @ standardized
    weight = torch.linalg.solve(
        gram + ridge * torch.eye(gram.shape[0], device="cuda"),
        standardized.T @ (labels - labelcenter),
    )
    base = np.full(len(feature), np.nan, dtype=np.float32)
    for frame in (fit, gate):
        source = frame.source.to_numpy(dtype=np.int64)
        output = []
        for start in range(0, len(source), batchsize * 8):
            batch = torch.as_tensor(
                feature[source[start : start + batchsize * 8]],
                device="cuda",
                dtype=torch.float32,
            )
            output.append(((batch - center) / scale @ weight + labelcenter).cpu())
        base[source] = torch.cat(output).numpy()
    fitbase = torch.as_tensor(base[rows], device="cuda")
    variance = (fitbase - labels).square().mean().clamp_min(1e-6)
    return base, {
        "kind": "ligand ridge",
        "ridge": ridge,
        "fit_residual_variance": float(variance.cpu()),
        "label_center": float(labelcenter.cpu()),
        "feature_center": center.cpu().tolist(),
        "feature_scale": scale.cpu().tolist(),
        "weight": weight.cpu().tolist(),
    }


def meanbase(
    fit: pd.DataFrame,
    gate: pd.DataFrame,
    rows: int,
) -> tuple[np.ndarray, dict[str, object]]:
    label = torch.as_tensor(
        fit.affinity.to_numpy(dtype=np.float32), device="cuda"
    )
    center = label.mean()
    variance = label.var().clamp_min(1e-6)
    base = np.full(rows, np.nan, dtype=np.float32)
    for frame in (fit, gate):
        base[frame.source.to_numpy(dtype=np.int64)] = float(center.cpu())
    return base, {
        "kind": "global mean",
        "fit_residual_variance": float(variance.cpu()),
        "label_center": float(center.cpu()),
    }


def makebatch(
    frame: pd.DataFrame,
    rowindex: np.ndarray,
    feature: np.ndarray,
    targetindex: dict[str, int],
    segments: torch.Tensor,
    wrong: dict[str, str],
    center: torch.Tensor,
    scale: torch.Tensor,
    counts: dict[str, int],
    base: np.ndarray,
) -> tuple[torch.Tensor, ...]:
    selected = frame.loc[rowindex]
    source = selected.source.to_numpy(dtype=np.int64)
    ligand = torch.as_tensor(feature[source], device="cuda").float()
    ligand[:, -10:] = (ligand[:, -10:] - center) / scale
    names = selected.target.astype(str).tolist()
    correctindex = torch.tensor(
        [targetindex[name] for name in names], device="cuda", dtype=torch.long
    )
    wrongindex = torch.tensor(
        [targetindex[wrong[name]] for name in names], device="cuda", dtype=torch.long
    )
    label = torch.as_tensor(
        selected.affinity.to_numpy(dtype=np.float32), device="cuda"
    )
    weight = torch.tensor(
        [1.0 / counts[name] for name in names], device="cuda", dtype=torch.float32
    )
    baseline = torch.as_tensor(base[source], device="cuda", dtype=torch.float32)
    return (
        segments[correctindex],
        segments[wrongindex],
        ligand,
        label,
        weight,
        baseline,
    )


def applybase(
    output: dict[str, torch.Tensor], baseline: torch.Tensor
) -> dict[str, torch.Tensor]:
    return output | {"prediction": baseline + output["prediction"]}


def fitprior(
    model: PairPrior,
    frame: pd.DataFrame,
    feature: np.ndarray,
    targetindex: dict[str, int],
    segments: torch.Tensor,
    center: torch.Tensor,
    scale: torch.Tensor,
    epochs: int,
    batchsize: int,
    learningrate: float,
    contrastweight: float,
    seed: int,
    useamp: bool,
    base: np.ndarray,
) -> list[dict[str, float]]:
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=learningrate, weight_decay=1e-4
    )
    scaler = torch.amp.GradScaler("cuda", enabled=useamp)
    wrong = wrongtargets(frame, seed + 11)
    counts = frame.target.value_counts().astype(int).to_dict()
    history: list[dict[str, float]] = []
    for epoch in range(epochs):
        model.train()
        totals = torch.zeros(4, device="cuda")
        seen = 0
        order = balancedorder(frame, seed + epoch)
        for rowindex in batches(order, batchsize):
            correctprotein, wrongprotein, ligand, label, weight, baseline = makebatch(
                frame,
                rowindex,
                feature,
                targetindex,
                segments,
                wrong,
                center,
                scale,
                counts,
                base,
            )
            with torch.autocast("cuda", dtype=torch.float16, enabled=useamp):
                correct = applybase(
                    model.predict(correctprotein, ligand), baseline
                )
                wrongoutput = applybase(
                    model.predict(wrongprotein, ligand), baseline
                )
                rownll = gaussiannll(correct, label)
                meanerror = (correct["prediction"] - label).square()
                nll = (rownll * weight).sum() / weight.sum()
                meanloss = (meanerror * weight).sum() / weight.sum()
                contrast = evidencecontrast(
                    correct, wrongoutput, label, weight=weight
                )
                loss = nll + 0.25 * meanloss + contrastweight * contrast
            optimizer.zero_grad(set_to_none=True)
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            scaler.step(optimizer)
            scaler.update()
            totals += torch.stack(
                (loss.detach(), nll.detach(), meanloss.detach(), contrast.detach())
            )
            seen += 1
        history.append(
            {
                "epoch": epoch + 1,
                "loss": float((totals[0] / seen).cpu()),
                "nll": float((totals[1] / seen).cpu()),
                "mean_square": float((totals[2] / seen).cpu()),
                "contrast": float((totals[3] / seen).cpu()),
            }
        )
        print(
            f"epoch {epoch + 1}/{epochs} "
            f"loss={history[-1]['loss']:.4f} nll={history[-1]['nll']:.4f} "
            f"mse={history[-1]['mean_square']:.4f} "
            f"contrast={history[-1]['contrast']:.4f}",
            flush=True,
        )
    return history


def initializehead(
    model: PairPrior,
    frame: pd.DataFrame,
    base: np.ndarray,
) -> dict[str, float]:
    """Initialize the correction and variance from fitting residuals on CUDA."""

    source = frame.source.to_numpy(dtype=np.int64)
    label = torch.as_tensor(frame.affinity.to_numpy(dtype=np.float32), device="cuda")
    baseline = torch.as_tensor(base[source], device="cuda")
    residual = label - baseline
    center = residual.mean()
    variance = residual.var().clamp_min(1e-4)
    with torch.no_grad():
        final = model.head[-1]
        final.weight.mul_(1e-3)
        final.bias[0] = center
        final.bias[1] = variance.log()
    return {
        "residual_center": float(center.cpu()),
        "residual_variance": float(variance.cpu()),
    }


@torch.no_grad()
def evaluate(
    model: PairPrior,
    frame: pd.DataFrame,
    feature: np.ndarray,
    targetindex: dict[str, int],
    segments: torch.Tensor,
    center: torch.Tensor,
    scale: torch.Tensor,
    batchsize: int,
    seed: int,
    base: np.ndarray,
    baselinevariance: float,
) -> dict[str, object]:
    model.eval()
    wrong = wrongtargets(frame, seed + 19)
    counts = frame.target.value_counts().astype(int).to_dict()
    correctprediction: list[torch.Tensor] = []
    wrongprediction: list[torch.Tensor] = []
    correctnll: list[torch.Tensor] = []
    wrongnll: list[torch.Tensor] = []
    labels: list[torch.Tensor] = []
    baselines: list[torch.Tensor] = []
    targetids: list[torch.Tensor] = []
    targetmap = {name: position for position, name in enumerate(sorted(counts))}
    order = frame.index.to_numpy(dtype=np.int64)
    for rowindex in batches(order, batchsize):
        correctprotein, wrongprotein, ligand, label, _, baseline = makebatch(
            frame,
            rowindex,
            feature,
            targetindex,
            segments,
            wrong,
            center,
            scale,
            counts,
            base,
        )
        correct = applybase(model.predict(correctprotein, ligand), baseline)
        wrongoutput = applybase(model.predict(wrongprotein, ligand), baseline)
        names = frame.loc[rowindex].target.astype(str).tolist()
        correctprediction.append(correct["prediction"])
        wrongprediction.append(wrongoutput["prediction"])
        correctnll.append(gaussiannll(correct, label))
        wrongnll.append(gaussiannll(wrongoutput, label))
        labels.append(label)
        baselines.append(baseline)
        targetids.append(
            torch.tensor(
                [targetmap[name] for name in names], device="cuda", dtype=torch.long
            )
        )
    correctvalue = torch.cat(correctprediction)
    wrongvalue = torch.cat(wrongprediction)
    labelvalue = torch.cat(labels)
    baselinevalue = torch.cat(baselines)
    correctnllvalue = torch.cat(correctnll)
    wrongnllvalue = torch.cat(wrongnll)
    targetvalue = torch.cat(targetids)
    correctrmse = []
    wrongrmse = []
    correctmae = []
    wrongmae = []
    correcttargetnll = []
    wrongtargetnll = []
    baselinermse = []
    baselinemae = []
    baselinetargetnll = []
    correctspearman = []
    wrongspearman = []
    baselinespearman = []
    variance = torch.tensor(baselinevariance, device="cuda").clamp_min(1e-6)

    def spearman(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
        if len(left) < 2 or left.unique().numel() < 2 or right.unique().numel() < 2:
            return left.new_tensor(float("nan"))
        leftrank = torch.argsort(torch.argsort(left, stable=True), stable=True).float()
        rightrank = torch.argsort(torch.argsort(right, stable=True), stable=True).float()
        return torch.corrcoef(torch.stack((leftrank, rightrank)))[0, 1]

    for target in range(len(targetmap)):
        active = targetvalue == target
        correcterror = correctvalue[active] - labelvalue[active]
        wrongerror = wrongvalue[active] - labelvalue[active]
        correctrmse.append(correcterror.square().mean().sqrt())
        wrongrmse.append(wrongerror.square().mean().sqrt())
        correctmae.append(correcterror.abs().mean())
        wrongmae.append(wrongerror.abs().mean())
        correcttargetnll.append(correctnllvalue[active].mean())
        wrongtargetnll.append(wrongnllvalue[active].mean())
        baselineerror = baselinevalue[active] - labelvalue[active]
        baselinermse.append(baselineerror.square().mean().sqrt())
        baselinemae.append(baselineerror.abs().mean())
        baselinetargetnll.append(
            0.5
            * (
                torch.log(2.0 * torch.pi * variance)
                + baselineerror.square() / variance
            ).mean()
        )
        correctspearman.append(spearman(labelvalue[active], correctvalue[active]))
        wrongspearman.append(spearman(labelvalue[active], wrongvalue[active]))
        baselinespearman.append(spearman(labelvalue[active], baselinevalue[active]))
    metric = {
        "correct": {
            "target_macro_rmse": float(torch.stack(correctrmse).mean().cpu()),
            "target_macro_mae": float(torch.stack(correctmae).mean().cpu()),
            "target_macro_nll": float(torch.stack(correcttargetnll).mean().cpu()),
            "within_target_spearman": float(torch.stack(correctspearman).nanmean().cpu()),
        },
        "wrongprotein": {
            "target_macro_rmse": float(torch.stack(wrongrmse).mean().cpu()),
            "target_macro_mae": float(torch.stack(wrongmae).mean().cpu()),
            "target_macro_nll": float(torch.stack(wrongtargetnll).mean().cpu()),
            "within_target_spearman": float(torch.stack(wrongspearman).nanmean().cpu()),
        },
        "baseline": {
            "target_macro_rmse": float(torch.stack(baselinermse).mean().cpu()),
            "target_macro_mae": float(torch.stack(baselinemae).mean().cpu()),
            "target_macro_nll": float(torch.stack(baselinetargetnll).mean().cpu()),
            "within_target_spearman": float(torch.stack(baselinespearman).nanmean().cpu()),
        },
    }
    rmsegain = torch.stack(wrongrmse) - torch.stack(correctrmse)
    nllgain = torch.stack(wrongtargetnll) - torch.stack(correcttargetnll)
    basermsegain = torch.stack(baselinermse) - torch.stack(correctrmse)
    basenllgain = torch.stack(baselinetargetnll) - torch.stack(correcttargetnll)
    wrongspearmangain = torch.stack(correctspearman) - torch.stack(wrongspearman)
    basespearmangain = torch.stack(correctspearman) - torch.stack(baselinespearman)
    generator = torch.Generator(device="cuda").manual_seed(seed + 29)
    samples = torch.randint(
        len(rmsegain),
        (2000, len(rmsegain)),
        device="cuda",
        generator=generator,
    )
    rmsebootstrap = rmsegain[samples].mean(dim=1)
    nllbootstrap = nllgain[samples].mean(dim=1)
    basermsebootstrap = basermsegain[samples].mean(dim=1)
    basenllbootstrap = basenllgain[samples].mean(dim=1)
    wrongspearmanbootstrap = wrongspearmangain[samples].nanmean(dim=1)
    basespearmanbootstrap = basespearmangain[samples].nanmean(dim=1)
    change = (correctvalue - wrongvalue).abs()
    comparison = {
        "rmse_gain": float(rmsegain.mean().cpu()),
        "rmse_ci95": [
            float(torch.quantile(rmsebootstrap, 0.025).cpu()),
            float(torch.quantile(rmsebootstrap, 0.975).cpu()),
        ],
        "nll_gain": float(nllgain.mean().cpu()),
        "nll_ci95": [
            float(torch.quantile(nllbootstrap, 0.025).cpu()),
            float(torch.quantile(nllbootstrap, 0.975).cpu()),
        ],
        "mean_absolute_protein_effect": float(change.mean().cpu()),
        "changed_fraction": float((change > 1e-6).float().mean().cpu()),
        "baseline_rmse_gain": float(basermsegain.mean().cpu()),
        "baseline_rmse_ci95": [
            float(torch.quantile(basermsebootstrap, 0.025).cpu()),
            float(torch.quantile(basermsebootstrap, 0.975).cpu()),
        ],
        "baseline_nll_gain": float(basenllgain.mean().cpu()),
        "baseline_nll_ci95": [
            float(torch.quantile(basenllbootstrap, 0.025).cpu()),
            float(torch.quantile(basenllbootstrap, 0.975).cpu()),
        ],
        "spearman_gain": float(wrongspearmangain.nanmean().cpu()),
        "spearman_ci95": [
            float(torch.quantile(wrongspearmanbootstrap, 0.025).cpu()),
            float(torch.quantile(wrongspearmanbootstrap, 0.975).cpu()),
        ],
        "baseline_spearman_gain": float(basespearmangain.nanmean().cpu()),
        "baseline_spearman_ci95": [
            float(torch.quantile(basespearmanbootstrap, 0.025).cpu()),
            float(torch.quantile(basespearmanbootstrap, 0.975).cpu()),
        ],
    }
    attributionpassed = (
        comparison["rmse_ci95"][0] > 0.0
        and comparison["nll_ci95"][0] > 0.0
        and comparison["spearman_ci95"][0] > 0.0
        and comparison["mean_absolute_protein_effect"] > 1e-6
    )
    performancepassed = (
        comparison["baseline_rmse_ci95"][0] > 0.0
        and comparison["baseline_nll_ci95"][0] > 0.0
        and comparison["baseline_spearman_ci95"][0] > 0.0
    )
    return {
        "metrics": metric,
        "comparison": comparison,
        "attribution_gate": "PASS" if attributionpassed else "STOP",
        "performance_gate": "PASS" if performancepassed else "STOP",
        "gate": "PASS" if attributionpassed and performancepassed else "STOP",
        "targets": len(targetmap),
        "queries": len(frame),
        "components": int(frame.hcluster.nunique()),
        "wrong_target_hash": hashlib.sha256(
            json.dumps(wrong, sort_keys=True).encode("utf-8")
        ).hexdigest(),
    }


def filesha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while block := source.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def telemetrysummary(
    samples: list[dict[str, float]], peakmemory: float
) -> dict[str, float | int | None]:
    if not samples:
        return {
            "samples": 0,
            "mean_utilization_percent": None,
            "peak_utilization_percent": None,
            "mean_power_watts": None,
            "peak_power_watts": None,
            "peak_nvidia_memory_mib": None,
            "peak_torch_memory_mib": peakmemory,
        }
    return {
        "samples": len(samples),
        "mean_utilization_percent": float(
            np.mean([sample["utilization"] for sample in samples])
        ),
        "peak_utilization_percent": max(sample["utilization"] for sample in samples),
        "mean_power_watts": float(
            np.mean([sample["power"] for sample in samples])
        ),
        "peak_power_watts": max(sample["power"] for sample in samples),
        "peak_nvidia_memory_mib": max(sample["memory"] for sample in samples),
        "peak_torch_memory_mib": peakmemory,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoint", choices=("pKi", "pKd"), default="pKi")
    parser.add_argument("--epochs", type=int, default=8)
    parser.add_argument("--batch", type=int, default=512)
    parser.add_argument("--width", type=int, default=192)
    parser.add_argument("--holdout", type=float, default=0.2)
    parser.add_argument("--learning-rate", type=float, default=2e-4)
    parser.add_argument("--contrast-weight", type=float, default=0.15)
    parser.add_argument("--base", choices=("mean", "ridge"), default="ridge")
    parser.add_argument("--ridge", type=float, default=10.0)
    parser.add_argument("--fit-targets", type=int, default=0)
    parser.add_argument("--gate-targets", type=int, default=0)
    parser.add_argument("--rows-per-target", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--no-amp", action="store_true")
    parser.add_argument("--out", type=Path, default=REPORT)
    args = parser.parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError("the TRAIN-only prior gate requires CUDA")
    if args.epochs < 1:
        raise ValueError("epochs must be positive")
    if args.batch < 2:
        raise ValueError("batch size must be at least two")
    if args.ridge <= 0:
        raise ValueError("ridge must be positive")
    if min(args.fit_targets, args.gate_targets, args.rows_per_target) < 0:
        raise ValueError("diagnostic limits must be non-negative")
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    torch.set_float32_matmul_precision("high")
    torch.backends.cudnn.benchmark = True
    started = time.perf_counter()
    frame = loadframe(args.endpoint)
    feature, targetindex, proteinsegments = loadarrays(frame)
    fit, gate = splitcomponents(frame, args.holdout, args.seed)
    fit = limitrows(
        fit, args.fit_targets, args.rows_per_target, args.seed + 1
    )
    gate = limitrows(
        gate, args.gate_targets, args.rows_per_target, args.seed + 2
    )
    if set(fit.hcluster).intersection(gate.hcluster):
        raise RuntimeError("diagnostic limiting broke homology closure")
    if args.base == "ridge":
        base, baseinfo = ridgebase(fit, gate, feature, args.ridge, args.batch)
    else:
        base, baseinfo = meanbase(fit, gate, len(feature))
    center, scale = descriptorstats(fit, feature, args.batch)
    model = PairPrior(width=args.width).cuda()
    initialization = initializehead(model, fit, base)
    torch.cuda.reset_peak_memory_stats()
    samples: list[dict[str, float]] = []
    stop = threading.Event()
    monitor = threading.Thread(target=telemetry, args=(stop, samples), daemon=True)
    monitor.start()
    trainingstarted = time.perf_counter()
    try:
        history = fitprior(
            model,
            fit,
            feature,
            targetindex,
            proteinsegments,
            center,
            scale,
            args.epochs,
            args.batch,
            args.learning_rate,
            args.contrast_weight,
            args.seed,
            not args.no_amp,
            base,
        )
        trainseconds = time.perf_counter() - trainingstarted
        evaluationstarted = time.perf_counter()
        evaluation = evaluate(
            model,
            gate,
            feature,
            targetindex,
            proteinsegments,
            center,
            scale,
            args.batch,
            args.seed,
            base,
            float(baseinfo["fit_residual_variance"]),
        )
        evaluationseconds = time.perf_counter() - evaluationstarted
    finally:
        stop.set()
        monitor.join(timeout=2.0)
    peakmemory = torch.cuda.max_memory_allocated() / (1024 ** 2)
    checkpoint = args.out.with_suffix(".pt")
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "state": model.state_dict(),
            "descriptor_center": center.cpu(),
            "descriptor_scale": scale.cpu(),
            "seed": args.seed,
            "endpoint": args.endpoint,
        },
        checkpoint,
    )
    result = {
        "protocol": "TRAIN-only homology-held-out sequence-conditioned k=0 pair prior",
        "roles_read": ["train"],
        "endpoint": args.endpoint,
        "support": 0,
        "seed": args.seed,
        "fit_rows": len(fit),
        "fit_targets": int(fit.target.nunique()),
        "fit_components": int(fit.hcluster.nunique()),
        "gate_rows": len(gate),
        "gate_targets": int(gate.target.nunique()),
        "gate_components": int(gate.hcluster.nunique()),
        "target_balancing": (
            "all fit rows once per epoch in target round-robin order "
            "with inverse-frequency loss"
        ),
        "wrong_protein": "deterministic cross-homology target rotation",
        "base": baseinfo,
        "parameters": sum(parameter.numel() for parameter in model.parameters()),
        "epochs": args.epochs,
        "batch": args.batch,
        "width": args.width,
        "learning_rate": args.learning_rate,
        "contrast_weight": args.contrast_weight,
        "amp": not args.no_amp,
        "initialization": initialization,
        "mean_square_weight": 0.25,
        "history": history,
        **evaluation,
        "training_seconds": trainseconds,
        "evaluation_seconds": evaluationseconds,
        "total_seconds": time.perf_counter() - started,
        "gpu": torch.cuda.get_device_name(),
        "gpu_telemetry": telemetrysummary(samples, peakmemory),
        "sources": {
            "registry_sha256": filesha(REGISTRY),
            "ligand_sha256": filesha(LIGANDS),
            "protein_sha256": filesha(PROTEINS),
        },
        "checkpoint": str(checkpoint),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(result, indent=2, allow_nan=False), encoding="utf-8"
    )
    print(
        json.dumps({
            "gate": result["gate"],
            "attribution_gate": result["attribution_gate"],
            "performance_gate": result["performance_gate"],
            "correct": result["metrics"]["correct"],
            "wrongprotein": result["metrics"]["wrongprotein"],
            "baseline": result["metrics"]["baseline"],
            "comparison": result["comparison"],
            "out": str(args.out),
        }, indent=2),
        flush=True,
    )


if __name__ == "__main__":
    main()
