"""Synthetic positive control for the support-identifiable Meta-Section."""
from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import torch

from research.meta_fewshot.meta_section import IdentifiableMetaSection

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "report" / "meta_fewshot" / "SYNTHETIC_META_SECTION.json"


@dataclass(frozen=True)
class SyntheticConfig:
    seed: int = 20260810
    input_dim: int = 12
    true_dim: int = 3
    observations_per_task: int = 12
    source_tasks: int = 96
    evaluation_tasks: int = 48
    steps: int = 700
    tasks_per_step: int = 6
    ridge: float = 0.08
    noise_sd: float = 0.05
    learning_rate: float = 0.025


def make_tasks(generator: torch.Generator, count: int, observations: int,
               input_dim: int, basis: torch.Tensor, noise_sd: float):
    tasks = []
    for _ in range(count):
        phi = torch.randn(observations, input_dim, generator=generator, dtype=torch.float64)
        coefficient = torch.randn(basis.shape[1], generator=generator, dtype=torch.float64)
        y = phi @ basis @ coefficient
        y += noise_sd * torch.randn(observations, generator=generator, dtype=torch.float64)
        tasks.append((phi, y))
    return tasks


def train(config: SyntheticConfig):
    torch.manual_seed(config.seed)
    generator = torch.Generator().manual_seed(config.seed)
    true_basis = torch.linalg.qr(torch.randn(
        config.input_dim, config.true_dim, generator=generator, dtype=torch.float64
    ), mode="reduced").Q
    source = make_tasks(generator, config.source_tasks, config.observations_per_task,
                        config.input_dim, true_basis, config.noise_sd)
    evaluation = make_tasks(generator, config.evaluation_tasks, config.observations_per_task,
                            config.input_dim, true_basis, config.noise_sd)
    model = IdentifiableMetaSection(config.input_dim, config.true_dim, config.ridge)
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate)
    legal_k = (1, 2, 3, 5)
    losses = []
    for step in range(config.steps):
        optimizer.zero_grad()
        task_indices = torch.randint(config.source_tasks, (config.tasks_per_step,),
                                     generator=generator)
        batch_loss = torch.zeros((), dtype=torch.float64)
        for task_index in task_indices:
            phi, y = source[int(task_index)]
            order = torch.randperm(len(phi), generator=generator)
            k = legal_k[int(torch.randint(len(legal_k), (), generator=generator))]
            support, query = order[:k], order[k:]
            state = model.adapt(phi[support], y[support])
            prediction, _, _ = model.query(phi[query], state)
            batch_loss = batch_loss + (prediction - y[query]).square().mean()
        batch_loss = batch_loss / config.tasks_per_step
        batch_loss.backward()
        optimizer.step()
        losses.append(float(batch_loss.detach()))
    return model, true_basis, evaluation, losses


def evaluate(model, true_basis, tasks, config: SyntheticConfig) -> dict:
    correct, zero, foreign, permuted, coverages = [], [], [], [], []
    ranks = []
    for index, (phi, y) in enumerate(tasks):
        support = torch.arange(5)
        query = torch.arange(5, len(phi))
        state = model.adapt(phi[support], y[support])
        prediction, coverage, _ = model.query(phi[query], state)
        correct.append(float((prediction - y[query]).square().mean()))
        zero.append(float(y[query].square().mean()))
        permuted_state = model.adapt(phi[support], y[support].roll(1))
        permuted_prediction, _, _ = model.query(phi[query], permuted_state)
        permuted.append(float((permuted_prediction - y[query]).square().mean()))

        foreign_phi, foreign_y = tasks[(index + 1) % len(tasks)]
        foreign_state = model.adapt(foreign_phi[support], foreign_y[support])
        foreign_prediction, _, _ = model.query(phi[query], foreign_state)
        foreign.append(float((foreign_prediction - y[query]).square().mean()))
        coverages.extend(float(value) for value in coverage)
        ranks.append(state.rank)

    learned = model.basis().detach()
    overlap = float((true_basis.T @ learned).square().sum() / config.true_dim)
    metrics = {
        "correct_mse": float(np.mean(correct)),
        "zero_mse": float(np.mean(zero)),
        "foreign_mse": float(np.mean(foreign)),
        "permuted_mse": float(np.mean(permuted)),
        "subspace_overlap": overlap,
        "mean_coverage": float(np.mean(coverages)),
        "max_section_rank": max(ranks),
    }
    checks = {
        "planted_family_recovered": overlap >= 0.90,
        "d_true_beats_d0": metrics["correct_mse"] < 0.50 * metrics["zero_mse"],
        "correct_beats_foreign": metrics["correct_mse"] < 0.50 * metrics["foreign_mse"],
        "correct_beats_permuted": metrics["correct_mse"] < 0.50 * metrics["permuted_mse"],
        "freedom_bounded_by_k": metrics["max_section_rank"] <= 5,
    }
    return metrics, checks


def run(config: SyntheticConfig = SyntheticConfig(), output: Path = DEFAULT_OUTPUT) -> dict:
    model, true_basis, tasks, losses = train(config)
    metrics, checks = evaluate(model, true_basis, tasks, config)
    result = {
        "schema": "MetaSieve.SyntheticMetaSection.v1",
        "config": asdict(config),
        "environment": {
            "python": subprocess.check_output(["python", "--version"], text=True).strip(),
            "torch": torch.__version__,
            "cuda_available": torch.cuda.is_available(),
            "git_head": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
            ).strip(),
        },
        "query_labels_exposed_to_model": False,
        "training_loss": {"initial_50": float(np.mean(losses[:50])),
                          "final_50": float(np.mean(losses[-50:]))},
        "metrics": metrics,
        "checks": checks,
        "TERMINAL_VERDICT": (
            "SYNTHETIC_META_SECTION_POSITIVE_CONTROL_PASS"
            if all(checks.values()) else "META_SECTION_IMPLEMENTATION_NOT_VALIDATED"
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(output=args.output)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["TERMINAL_VERDICT"].endswith("PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())
