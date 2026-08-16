"""TRAIN-only target-function memory gate for five-shot adaptation."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
import hashlib
import json
import math
from pathlib import Path
import threading
import time
from typing import Iterable

import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F

from research.shared.priorgate import (
    LIGANDS,
    PROTEINS,
    REGISTRY,
    filesha,
    loadarrays,
    loadframe,
    ridgebase,
    telemetry,
    telemetrysummary,
)
from scripts.audit import ROOT, proteincomponents
from scripts.metric import evaluateprotocol, pairedcomponents
from scripts.train import maketrainroster


REPORT = Path("reports/active/memorygate.v1.json")
SEQUENCES = ROOT / "target_sequences.json"
BITS = 1024
DESCRIPTORS = 10


def stablekey(value: object, seed: int) -> str:
    return hashlib.sha256(f"{seed}:{value}".encode("utf-8")).hexdigest()


def strictsplit(
    frame: pd.DataFrame, holdout: float, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split complete strict homology components without reading labels."""

    if not 0.0 < holdout < 0.5:
        raise ValueError("holdout fraction must lie in (0, 0.5)")
    if "component" not in frame or frame.component.isna().any():
        raise ValueError("every target must have a strict homology component")
    components = sorted(
        frame.component.unique(), key=lambda value: stablekey(value, seed)
    )
    count = max(1, math.ceil(len(components) * holdout))
    held = set(components[:count])
    gate = frame[frame.component.isin(held)].copy()
    fit = frame[~frame.component.isin(held)].copy()
    if fit.empty or gate.empty:
        raise RuntimeError("strict component split produced an empty role")
    if set(fit.component).intersection(gate.component):
        raise RuntimeError("strict homology components crossed the gate")
    return fit, gate


def targetlimit(frame: pd.DataFrame, limit: int, seed: int) -> pd.DataFrame:
    """Deterministic diagnostic target cap; zero preserves the complete role."""

    if limit < 0:
        raise ValueError("target limit cannot be negative")
    if not limit:
        return frame.copy()
    targets = sorted(frame.target.unique(), key=lambda value: stablekey(value, seed))
    return frame[frame.target.isin(targets[:limit])].copy()


def wrongtargets(frame: pd.DataFrame, seed: int) -> dict[str, str]:
    """Choose a deterministic protein from another strict component."""

    component = frame.groupby("target", sort=True).component.first().to_dict()
    targets = sorted(component, key=lambda value: stablekey(value, seed))
    if len(set(component.values())) < 2:
        raise ValueError("wrong-protein control needs two strict components")
    output: dict[str, str] = {}
    for position, target in enumerate(targets):
        for offset in range(1, len(targets) + 1):
            candidate = targets[(position + offset) % len(targets)]
            if component[candidate] != component[target]:
                output[target] = candidate
                break
    return output


@dataclass(frozen=True)
class FeatureMap:
    """Frozen low-dimensional ligand coordinates selected on fitting labels."""

    columns: torch.Tensor
    center: torch.Tensor
    scale: torch.Tensor

    def transform(self, feature: np.ndarray, rows: Iterable[int]) -> torch.Tensor:
        source = np.asarray(tuple(rows), dtype=np.int64)
        values = torch.as_tensor(
            feature[source][:, self.columns.cpu().numpy()],
            device="cuda",
            dtype=torch.float32,
        )
        return (values - self.center) / self.scale


@torch.no_grad()
def fitfeatures(
    frame: pd.DataFrame,
    feature: np.ndarray,
    base: np.ndarray,
    latent: int,
) -> FeatureMap:
    """Select fragments with reproducible within-target affinity covariance."""

    if not 1 <= latent <= BITS:
        raise ValueError("latent ligand dimension is outside Morgan-1024")
    score = torch.zeros(BITS, device="cuda", dtype=torch.float64)
    for _, group in frame.groupby("target", sort=True):
        rows = group.source.to_numpy(dtype=np.int64)
        ligand = torch.as_tensor(
            feature[rows, :BITS], device="cuda", dtype=torch.float64
        )
        residual = torch.as_tensor(
            group.affinity.to_numpy(dtype=np.float64) - base[rows], device="cuda"
        )
        residual = residual - residual.mean()
        ligand = ligand - ligand.mean(dim=0)
        covariance = ligand.T @ residual / max(len(rows) - 1, 1)
        score += covariance.square() * math.sqrt(len(rows))
    fragments = torch.topk(score, latent, sorted=True).indices
    descriptors = torch.arange(
        BITS, BITS + DESCRIPTORS, device="cuda", dtype=torch.long
    )
    columns = torch.cat((fragments, descriptors))
    rows = frame.source.to_numpy(dtype=np.int64)
    values = torch.as_tensor(
        feature[rows][:, columns.cpu().numpy()], device="cuda", dtype=torch.float64
    )
    center = values.mean(dim=0)
    scale = values.std(dim=0).clamp_min(1e-4)
    return FeatureMap(columns, center.float(), scale.float())


@dataclass(frozen=True)
class ExpertBank:
    """Target-specific ligand functions and their label-free protein keys."""

    targets: tuple[str, ...]
    coefficient: torch.Tensor
    protein: torch.Tensor
    proteincenter: torch.Tensor
    mapalpha: torch.Tensor
    fitnoise: float

    def proteinvalue(self, pooled: torch.Tensor) -> torch.Tensor:
        return F.normalize(pooled.float() - self.proteincenter, dim=-1, eps=1e-8)

    def mappedfunction(self, pooled: torch.Tensor) -> torch.Tensor:
        value = self.proteinvalue(pooled)
        return (value @ self.protein.T) @ self.mapalpha


@torch.no_grad()
def fitexperts(
    frame: pd.DataFrame,
    feature: np.ndarray,
    base: np.ndarray,
    mapping: FeatureMap,
    pooled: dict[str, torch.Tensor],
    ridge: float,
    mapridge: float,
) -> ExpertBank:
    """Fit one centered low-rank ligand-response ridge expert per target."""

    if ridge <= 0 or mapridge <= 0:
        raise ValueError("expert and protein-map ridge values must be positive")
    targets: list[str] = []
    coefficients: list[torch.Tensor] = []
    proteins: list[torch.Tensor] = []
    noises: list[torch.Tensor] = []
    dimension = len(mapping.columns)
    eye = torch.eye(dimension, device="cuda")
    for target, group in frame.groupby("target", sort=True):
        if target not in pooled or len(group) < 6:
            continue
        rows = group.source.to_numpy(dtype=np.int64)
        ligand = mapping.transform(feature, rows)
        response = torch.as_tensor(
            group.affinity.to_numpy(dtype=np.float32) - base[rows], device="cuda"
        )
        ligand = ligand - ligand.mean(dim=0)
        response = response - response.mean()
        if len(group) < dimension:
            dual = torch.linalg.solve(
                ligand @ ligand.T + ridge * torch.eye(len(group), device="cuda"),
                response,
            )
            coefficient = ligand.T @ dual
        else:
            coefficient = torch.linalg.solve(
                ligand.T @ ligand + ridge * eye, ligand.T @ response
            )
        error = response - ligand @ coefficient
        noises.append(error.square().mean())
        targets.append(str(target))
        coefficients.append(coefficient)
        proteins.append(pooled[str(target)])
    if len(targets) < 2:
        raise RuntimeError("the fitting role produced fewer than two experts")
    coefficient = torch.stack(coefficients)
    rawprotein = torch.stack(proteins).float()
    proteincenter = rawprotein.mean(dim=0)
    protein = F.normalize(rawprotein - proteincenter, dim=1, eps=1e-8)
    kernel = protein @ protein.T
    mapalpha = torch.linalg.solve(
        kernel + mapridge * torch.eye(len(targets), device="cuda"), coefficient
    )
    fitnoise = float(torch.stack(noises).median().sqrt().clamp_min(0.1).cpu())
    return ExpertBank(
        tuple(targets), coefficient, protein, proteincenter, mapalpha, fitnoise
    )


@dataclass(frozen=True)
class MemoryConfig:
    noise: float = 0.8
    proteinweight: float = 6.0
    functionweight: float = 0.5
    topk: int = 48
    offsetscale: float = 4.0
    dynamicbias: float = 0.5
    nullbias: float = 0.0


class FunctionMemory:
    """Bayesian expert retrieval with an exact five-support evidence update."""

    def __init__(self, bank: ExpertBank, config: MemoryConfig) -> None:
        if config.noise <= 0 or config.offsetscale <= 0:
            raise ValueError("noise and offset scale must be positive")
        if not 1 <= config.topk <= len(bank.targets):
            raise ValueError("top-k retrieval is outside the expert bank")
        self.bank = bank
        self.config = config

    def _evidence(
        self, score: torch.Tensor, residual: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Integrate a Gaussian target offset instead of fitting it by query loss."""

        difference = residual[:, None] - score
        count = len(residual)
        noise2 = self.config.noise**2
        offset2 = self.config.offsetscale**2
        total = difference.sum(dim=0)
        quadratic = difference.square().sum(dim=0) / noise2
        quadratic -= offset2 * total.square() / (
            noise2 * (noise2 + count * offset2)
        )
        logdet = (count - 1) * math.log(noise2) + math.log(
            noise2 + count * offset2
        )
        loglikelihood = -0.5 * (
            count * math.log(2.0 * math.pi) + logdet + quadratic
        )
        offset = offset2 * total / (noise2 + count * offset2)
        return loglikelihood, offset

    @torch.no_grad()
    def predict(
        self,
        support: torch.Tensor,
        supportlabel: torch.Tensor,
        supportbase: torch.Tensor,
        query: torch.Tensor,
        querybase: torch.Tensor,
        protein: torch.Tensor | None,
    ) -> dict[str, torch.Tensor]:
        if len(supportlabel) != 5:
            raise ValueError("function memory requires exactly five support labels")
        if support.device.type != "cuda" or query.device.type != "cuda":
            raise RuntimeError("function memory numerical inference requires CUDA")
        residual = supportlabel - supportbase
        static = self.bank.coefficient
        if protein is None:
            coefficient = torch.cat((static, torch.zeros_like(static[:1])), dim=0)
            prior = torch.zeros(len(coefficient), device="cuda")
            prior[-1] = self.config.nullbias
            dynamicindex = None
        else:
            proteinvalue = self.bank.proteinvalue(protein)
            mapped = self.bank.mappedfunction(protein)
            rawsimilarity = proteinvalue @ self.bank.protein.T
            functionsimilarity = F.cosine_similarity(
                mapped[None, :], static, dim=1, eps=1e-8
            )
            staticprior = self.config.proteinweight * (
                rawsimilarity + self.config.functionweight * functionsimilarity
            )
            if self.config.topk < len(staticprior):
                keep = torch.topk(staticprior, self.config.topk).indices
                mask = torch.full_like(staticprior, -torch.inf)
                mask[keep] = staticprior[keep]
                staticprior = mask
            reference = torch.logsumexp(staticprior, dim=0) - math.log(
                self.config.topk
            )
            coefficient = torch.cat(
                (static, mapped[None, :], torch.zeros_like(static[:1])), dim=0
            )
            prior = torch.cat(
                (
                    staticprior,
                    reference.new_tensor(
                        [
                            float(reference + self.config.dynamicbias),
                            float(reference + self.config.nullbias),
                        ]
                    ),
                )
            )
            dynamicindex = len(static)
        supportscore = support @ coefficient.T
        queryscore = query @ coefficient.T
        loglikelihood, offset = self._evidence(supportscore, residual)
        weight = torch.softmax(prior + loglikelihood, dim=0)
        experts = querybase[:, None] + queryscore + offset[None, :]
        prediction = experts @ weight
        entropy = -(weight * weight.clamp_min(1e-12).log()).sum()
        return {
            "prediction": prediction,
            "weight": weight,
            "entropy": entropy,
            "nullprobability": weight[-1],
            "dynamicprobability": (
                weight[dynamicindex]
                if dynamicindex is not None
                else weight.new_zeros(())
            ),
        }


def pooledproteins() -> dict[str, torch.Tensor]:
    archive = np.load(PROTEINS, allow_pickle=False)
    return {
        str(key): torch.as_tensor(value, device="cuda", dtype=torch.float32)
        for key, value in zip(archive["keys"], archive["pooled"])
    }


def attachcomponents(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    document = json.loads(SEQUENCES.read_text(encoding="utf-8"))
    components, stats = proteincomponents(document["sequences"])
    output = frame[frame.target.isin(components)].copy()
    output["component"] = output.target.map(components)
    if output.empty or output.component.isna().any():
        raise RuntimeError("strict protein mapping did not cover TRAIN rows")
    return output, stats


def episodes(
    frame: pd.DataFrame,
    pooled: dict[str, torch.Tensor],
    feature: np.ndarray,
    targets: int,
    queries: int,
) -> tuple[pd.DataFrame, list]:
    local = frame.reset_index(drop=True).copy()
    local["hcluster"] = local.component
    roster = maketrainroster(
        local,
        pooled,
        feature,
        targets=targets,
        querycap=queries,
        support=5,
    )
    return local, roster


def supportkernel(
    support: torch.Tensor,
    residual: torch.Tensor,
    query: torch.Tensor,
    querybase: torch.Tensor,
) -> torch.Tensor:
    """Five-point ligand-only kernel control in the same frozen coordinates."""

    width = support.shape[1]
    supportsquare = torch.cdist(support, support).square() / width
    querysquare = torch.cdist(query, support).square() / width
    kernel = torch.exp(-0.5 * supportsquare)
    cross = torch.exp(-0.5 * querysquare)
    centered = residual - residual.mean()
    weight = torch.linalg.solve(
        kernel + 0.5 * torch.eye(len(support), device="cuda"), centered
    )
    return querybase + residual.mean() + cross @ weight


@torch.no_grad()
def evaluate(
    memory: FunctionMemory,
    frame: pd.DataFrame,
    roster: list,
    feature: np.ndarray,
    base: np.ndarray,
    mapping: FeatureMap,
    pooled: dict[str, torch.Tensor],
    meanlabel: float,
    seed: int,
    bootstrap: bool,
) -> dict[str, object]:
    names = (
        "correct",
        "wrongprotein",
        "supportonly",
        "supportkernel",
        "intercept",
        "ligand",
        "global",
    )
    prediction: dict[str, list[float]] = {name: [] for name in names}
    labels: list[float] = []
    indices: list[int] = []
    entropy: list[float] = []
    nullprobability: list[float] = []
    dynamicprobability: list[float] = []
    wrong = wrongtargets(frame, seed + 19)
    for item in roster:
        supportrows = frame.source.iloc[list(item.support_indices)].to_numpy(
            dtype=np.int64
        )
        queryrows = frame.source.iloc[list(item.query_indices)].to_numpy(
            dtype=np.int64
        )
        sx = mapping.transform(feature, supportrows)
        qx = mapping.transform(feature, queryrows)
        sy = torch.as_tensor(
            frame.affinity.iloc[list(item.support_indices)].to_numpy(dtype=np.float32),
            device="cuda",
        )
        qy = torch.as_tensor(
            frame.affinity.iloc[list(item.query_indices)].to_numpy(dtype=np.float32),
            device="cuda",
        )
        sb = torch.as_tensor(base[supportrows], device="cuda")
        qb = torch.as_tensor(base[queryrows], device="cuda")
        correct = memory.predict(sx, sy, sb, qx, qb, pooled[item.target_key])
        wrongvalue = memory.predict(
            sx, sy, sb, qx, qb, pooled[wrong[item.target_key]]
        )
        supportvalue = memory.predict(sx, sy, sb, qx, qb, None)
        intercept = qb + (sy - sb).mean()
        arms = {
            "correct": correct["prediction"],
            "wrongprotein": wrongvalue["prediction"],
            "supportonly": supportvalue["prediction"],
            "supportkernel": supportkernel(sx, sy - sb, qx, qb),
            "intercept": intercept,
            "ligand": qb,
            "global": torch.full_like(qb, meanlabel),
        }
        for name, value in arms.items():
            prediction[name].extend(value.cpu().tolist())
        labels.extend(qy.cpu().tolist())
        indices.extend(item.query_indices)
        entropy.append(float(correct["entropy"].cpu()))
        nullprobability.append(float(correct["nullprobability"].cpu()))
        dynamicprobability.append(float(correct["dynamicprobability"].cpu()))
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
    result: dict[str, object] = {
        "metrics": metrics,
        "episodes": len(roster),
        "components": len(set(components.values())),
        "queries": len(labels),
        "posterior": {
            "mean_entropy": float(np.mean(entropy)),
            "mean_null_probability": float(np.mean(nullprobability)),
            "mean_dynamic_probability": float(np.mean(dynamicprobability)),
        },
    }
    if not bootstrap:
        return result
    paired = pairedcomponents(
        predictions=prediction,
        labels=labels,
        episodes=roster,
        prediction_indices=indices,
        component_by_target=components,
        reference="correct",
        seed=seed,
    )
    strictcontrols = ("wrongprotein", "supportonly")
    baselines = ("supportkernel", "intercept", "ligand", "global")
    specificity = all(
        paired[name][metric]["ci95"][0] > 0.0
        for name in strictcontrols
        for metric in ("rmse_gain", "spearman_gain", "pairwise_gain")
    )
    superiority = all(
        paired[name][metric]["ci95"][0] > 0.0
        for name in baselines
        for metric in ("rmse_gain", "spearman_gain", "pairwise_gain")
    )
    result |= {
        "gate": "PASS" if specificity and superiority else "STOP",
        "criteria": {
            "protein_specificity": specificity,
            "absolute_and_ranking_superiority": superiority,
        },
        "paired_component_bootstrap": paired,
    }
    return result


def selectionscore(result: dict[str, object]) -> float:
    metrics = result["metrics"]
    correct = metrics["correct"]
    score = (
        correct["target_macro_rmse"]
        + 0.25 * (1.0 - correct["within_target_spearman"])
        + 0.15 * (1.0 - correct["pairwise_accuracy"])
    )
    for control in ("wrongprotein", "supportonly"):
        score += 2.0 * max(
            0.0,
            correct["target_macro_rmse"] - metrics[control]["target_macro_rmse"],
        )
        score += 0.5 * max(
            0.0,
            metrics[control]["within_target_spearman"]
            - correct["within_target_spearman"],
        )
    return float(score)


def candidateconfigs(experts: int) -> list[MemoryConfig]:
    output = []
    for noise in (0.55, 0.8, 1.1):
        for proteinweight in (3.0, 7.0, 12.0):
            for topk in (24, 64):
                output.append(
                    MemoryConfig(
                        noise=noise,
                        proteinweight=proteinweight,
                        topk=min(topk, experts),
                    )
                )
    return output


def tune(
    bank: ExpertBank,
    frame: pd.DataFrame,
    roster: list,
    feature: np.ndarray,
    base: np.ndarray,
    mapping: FeatureMap,
    pooled: dict[str, torch.Tensor],
    meanlabel: float,
    seed: int,
) -> tuple[MemoryConfig, list[dict[str, object]]]:
    records = []
    for config in candidateconfigs(len(bank.targets)):
        result = evaluate(
            FunctionMemory(bank, config),
            frame,
            roster,
            feature,
            base,
            mapping,
            pooled,
            meanlabel,
            seed,
            bootstrap=False,
        )
        score = selectionscore(result)
        records.append(
            {
                "config": config.__dict__,
                "score": score,
                "correct": result["metrics"]["correct"],
                "wrongprotein": result["metrics"]["wrongprotein"],
                "supportonly": result["metrics"]["supportonly"],
            }
        )
    chosen = min(records, key=lambda value: value["score"])
    return MemoryConfig(**chosen["config"]), records


def jsonsafe(value):
    if isinstance(value, dict):
        return {key: jsonsafe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonsafe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def writemarkdown(result: dict[str, object], path: Path) -> None:
    metrics = result["metrics"]
    lines = [
        "# Target-function memory TRAIN-only gate",
        "",
        f"Decision: **{result['gate']}**",
        "",
        "This isolated route stores one low-dimensional ligand-response function per fitting target. "
        "A frozen protein embedding supplies the prior over those functions; exactly five support "
        "labels update it through an analytic Gaussian marginal likelihood with the target offset "
        "integrated out. A zero-function ligand-only expert is always present.",
        "",
        "## Results",
        "",
        "| Arm | RMSE | Spearman | Pairwise |",
        "|---|---:|---:|---:|",
    ]
    for name in (
        "correct",
        "wrongprotein",
        "supportonly",
        "supportkernel",
        "intercept",
        "ligand",
        "global",
    ):
        item = metrics[name]
        lines.append(
            f"| {name} | {item['target_macro_rmse']:.4f} | "
            f"{item['within_target_spearman']:.4f} | {item['pairwise_accuracy']:.4f} |"
        )
    lines += [
        "",
        "The decision is fail-closed: PASS requires component-bootstrap lower confidence bounds "
        "above zero for RMSE, Spearman, and pairwise gains over wrong-protein and support-only "
        "retrieval, and the same gains over every absolute baseline. A STOP result is not promoted "
        "to the core model.",
        "",
        "## Protocol",
        "",
        f"Only `{', '.join(result['roles_read'])}` affinity labels were read. The outer gate contains "
        f"{result['gate_episodes']} episodes from {result['gate_components']} globally closed "
        "Smith-Waterman homology components. Numerical fitting, posterior inference, metrics, and "
        "bootstrap evaluation ran on CUDA; utilization, power, memory, and wall time are in the JSON.",
        "",
        "## Literature basis",
        "",
        "- HyperPCM motivates protein-conditioned generation of target-specific QSAR functions: "
        "https://doi.org/10.1021/acs.jcim.3c01417",
        "- Graph neural processes motivate inferring a molecular response function from a small "
        "context set: https://doi.org/10.1186/s13321-024-00904-2",
        "- In-context few-shot molecular property prediction motivates retaining a non-parametric "
        "function memory at very small support size: https://arxiv.org/abs/2310.08863",
        "- Bayesian model averaging supplies the posterior expert-combination interpretation: "
        "https://doi.org/10.1214/ss/1009212519",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--latent", type=int, default=96)
    parser.add_argument("--expert-ridge", type=float, default=8.0)
    parser.add_argument("--map-ridge", type=float, default=2.0)
    parser.add_argument("--holdout", type=float, default=0.2)
    parser.add_argument("--inner-holdout", type=float, default=0.2)
    parser.add_argument("--queries", type=int, default=64)
    parser.add_argument("--fit-targets", type=int, default=0)
    parser.add_argument("--gate-targets", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--out", type=Path, default=REPORT)
    args = parser.parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError("the target-function memory gate requires CUDA")
    if args.queries < 1:
        raise ValueError("query cap must be positive")
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    torch.set_float32_matmul_precision("high")
    torch.cuda.reset_peak_memory_stats()
    started = time.perf_counter()
    samples: list[dict[str, float]] = []
    stop = threading.Event()
    monitor = threading.Thread(target=telemetry, args=(stop, samples), daemon=True)
    monitor.start()
    try:
        full = loadframe("pKi")
        feature, _, _ = loadarrays(full)
        full, componentstats = attachcomponents(full)
        pooled = pooledproteins()
        full = full[full.target.isin(pooled)].copy()
        outerfit, outergate = strictsplit(full, args.holdout, args.seed)
        outerfit = targetlimit(outerfit, args.fit_targets, args.seed + 1)
        outergate = targetlimit(outergate, args.gate_targets, args.seed + 2)
        innerfit, innergate = strictsplit(
            outerfit, args.inner_holdout, args.seed + 3
        )

        innerbase, innerbaseinfo = ridgebase(
            innerfit, innergate, feature, ridge=10.0, batchsize=512
        )
        innermap = fitfeatures(
            innerfit, feature, innerbase, args.latent
        )
        innerbank = fitexperts(
            innerfit,
            feature,
            innerbase,
            innermap,
            pooled,
            args.expert_ridge,
            args.map_ridge,
        )
        innergate, innerroster = episodes(
            innergate, pooled, feature, 0, args.queries
        )
        chosen, tuning = tune(
            innerbank,
            innergate,
            innerroster,
            feature,
            innerbase,
            innermap,
            pooled,
            float(innerbaseinfo["label_center"]),
            args.seed + 5,
        )

        base, baseinfo = ridgebase(
            outerfit, outergate, feature, ridge=10.0, batchsize=512
        )
        mapping = fitfeatures(outerfit, feature, base, args.latent)
        bank = fitexperts(
            outerfit,
            feature,
            base,
            mapping,
            pooled,
            args.expert_ridge,
            args.map_ridge,
        )
        outergate, gateroster = episodes(
            outergate, pooled, feature, args.gate_targets, args.queries
        )
        result = evaluate(
            FunctionMemory(bank, replace(chosen, topk=min(chosen.topk, len(bank.targets)))),
            outergate,
            gateroster,
            feature,
            base,
            mapping,
            pooled,
            float(baseinfo["label_center"]),
            args.seed,
            bootstrap=True,
        )
    finally:
        stop.set()
        monitor.join(timeout=2.0)
    result |= {
        "protocol": "TRAIN-only strict-homology-held-out k=5 target-function memory",
        "roles_read": ["train"],
        "endpoint": "pKi",
        "support": 5,
        "seed": args.seed,
        "fit_targets": int(outerfit.target.nunique()),
        "fit_components": int(outerfit.component.nunique()),
        "gate_targets": int(outergate.target.nunique()),
        "gate_components": int(outergate.component.nunique()),
        "gate_episodes": len(gateroster),
        "queries_per_episode": args.queries,
        "expert_count": len(bank.targets),
        "expert_dimension": int(bank.coefficient.shape[1]),
        "expert_fit_noise": bank.fitnoise,
        "selected_config": chosen.__dict__,
        "inner_selection": {
            "fit_targets": int(innerfit.target.nunique()),
            "gate_targets": int(innergate.target.nunique()),
            "gate_episodes": len(innerroster),
            "candidates": tuning,
        },
        "base": {
            "kind": baseinfo["kind"],
            "ridge": baseinfo["ridge"],
            "fit_residual_variance": baseinfo["fit_residual_variance"],
            "label_center": baseinfo["label_center"],
        },
        "strict_homology": componentstats,
        "feature_columns": [int(value) for value in mapping.columns.cpu()],
        "gpu": torch.cuda.get_device_name(),
        "wall_seconds": time.perf_counter() - started,
        "gpu_telemetry": telemetrysummary(
            samples, torch.cuda.max_memory_allocated() / (1024**2)
        ),
        "sources": {
            "registry_sha256": filesha(REGISTRY),
            "ligand_sha256": filesha(LIGANDS),
            "protein_sha256": filesha(PROTEINS),
            "sequence_sha256": filesha(SEQUENCES),
        },
    }
    serialized = jsonsafe(result)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(serialized, indent=2, allow_nan=False), encoding="utf-8"
    )
    writemarkdown(serialized, args.out.with_suffix(".md"))
    print(
        json.dumps(
            {
                "gate": serialized["gate"],
                "metrics": serialized["metrics"],
                "criteria": serialized["criteria"],
                "out": str(args.out),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
