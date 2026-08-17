"""Stage T2 — train the matched arms on the crossed double difference.

Six arms, one code path, one seed, one fixed budget, **no checkpoint selection
of any kind** (the reporting population is the withheld internal-validation
components, so selecting on them would leak).

    A  A_zero                    constant response; D_hat is identically 0
    B  B_transformation_only     R = f(tau); D_hat is identically 0
    C  C_protein                 correct protein
    D  D_protein_shuffled        stable cross-component protein permutation
    E  E_protein_matched_wrong   similarity-matched wrong protein, consistently
    F  F_label_shuffled          arm C with D permuted inside each key

**No "correct beats wrong" counterfactual loss is trained.** Stage S measured
that such a loss passes by inflating the wrong branch. Wrong proteins appear
here as a consistently substituted training input (arm E) and as an
evaluation-time substitution on arm C -- never as an optimization target.

Run:
    python -m tools.research.stageT_mmp.t2_train --arm C_protein
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys
import time

import numpy as np
import torch

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.internal_validation import partition_components
from scripts.qpsmp_data import stable_seed
from tools.research.stageT_mmp.observations import build_observations, load_governed
from tools.research.stageT_mmp.t2_dataset import (
    DoubleDifference, double_differences, effective_independent_units,
    shuffle_within_key, split_by_key_overlap, target_effects,
)
from tools.research.stageT_mmp.t2_model import (
    DiscriminatorConfig, DoubleDifferenceModel, parameter_report,
)

HERE = Path(__file__).resolve().parent

# Frozen before any arm trained.
SEED = 20260820
STEPS = 3000
BATCH = 256
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 1e-4
HUBER_BETA = 1.0
GRADIENT_CLIP = 5.0


@dataclass(frozen=True)
class ArmConfig:
    name: str
    mode: str
    protein_source: str      # "correct" | "shuffled" | "matched_wrong"
    shuffle_labels: bool


ARMS = {
    "A_zero": ArmConfig("A_zero", "zero", "correct", False),
    "B_transformation_only": ArmConfig(
        "B_transformation_only", "transformation", "correct", False),
    "C_protein": ArmConfig("C_protein", "protein", "correct", False),
    "D_protein_shuffled": ArmConfig(
        "D_protein_shuffled", "protein", "shuffled", False),
    "E_protein_matched_wrong": ArmConfig(
        "E_protein_matched_wrong", "protein", "matched_wrong", False),
    "F_label_shuffled": ArmConfig("F_label_shuffled", "protein", "correct", True),
}


def verify_gpu() -> dict:
    available = torch.cuda.is_available()
    return {
        "cuda_available": bool(available),
        "device": torch.cuda.get_device_name(0) if available else "cpu",
        "torch": torch.__version__,
        "note": ("this discriminator is small enough to be CPU-feasible; the "
                 "device is recorded, not required"),
    }


class ProteinTable:
    """Frozen ESM-2 150M pooled + masked-mean residue summary, per target."""

    def __init__(self, data) -> None:
        self.data = data
        self._cache: dict[str, np.ndarray] = {}

    def get(self, target: str) -> np.ndarray:
        if target not in self._cache:
            pooled, residues, mask = self.data.protein_for_target(target)
            pooled = pooled.to(torch.float32).numpy()
            residues = residues.to(torch.float32).numpy()
            weights = mask.to(torch.float32).numpy()[:, None]
            summary = (residues * weights).sum(0) / max(float(weights.sum()), 1.0)
            self._cache[target] = np.concatenate([pooled, summary]).astype(np.float32)
        return self._cache[target]

    @property
    def width(self) -> int:
        any_target = next(iter(self.data.tasks["meta_train"]))
        return int(self.get(any_target).shape[0])


def matched_wrong_map(data, table: ProteinTable, recipients: list[str],
                      candidates: list[str]) -> dict[str, str]:
    """Similarity-matched wrong protein: hardest admissible negative.

    Different CD-HIT40 component; drawn from `meta_train` and chosen on
    `meta_train`-only frozen protein similarity; most similar admissible protein
    by cosine; no shared document programme. Only the protein input changes.
    """
    component = {cell["target_id"]: cell["protein_group_40"] for cell in data.cells}
    documents: dict[str, set[str]] = {}
    for cell in data.cells:
        documents.setdefault(cell["target_id"], set()).update(
            str(panel).split("|")[0] for panel in cell["panel_ids"])
    pool = sorted(set(candidates))
    matrix = np.stack([table.get(target) for target in pool])
    matrix = matrix / np.maximum(np.linalg.norm(matrix, axis=1, keepdims=True), 1e-12)
    out: dict[str, str] = {}
    for recipient in sorted(set(recipients)):
        vector = table.get(recipient)
        vector = vector / max(float(np.linalg.norm(vector)), 1e-12)
        order = np.argsort(-(matrix @ vector))
        chosen = fallback = None
        for position in order:
            candidate = pool[int(position)]
            if candidate == recipient or component[candidate] == component[recipient]:
                continue
            if fallback is None:
                fallback = candidate
            if documents.get(candidate, set()) & documents.get(recipient, set()):
                continue
            chosen = candidate
            break
        out[recipient] = chosen or fallback
    return out


def shuffled_map(data, targets: list[str], seed: int) -> dict[str, str]:
    """Stable cross-component permutation of protein identity."""
    component = {cell["target_id"]: cell["protein_group_40"] for cell in data.cells}
    ordered = sorted(set(targets))
    rng = np.random.default_rng(stable_seed("stageT-protein-shuffle", seed))
    permuted = list(rng.permutation(ordered))
    for offset in range(1, len(ordered)):
        mapping = {ordered[i]: permuted[(i + offset) % len(permuted)]
                   for i in range(len(ordered))}
        if all(component[a] != component[b] for a, b in mapping.items()):
            return mapping
    raise ValueError("no cross-component permutation found")


@dataclass
class Banks:
    train: list[DoubleDifference]
    internal_repeated: list[DoubleDifference]
    internal_disjoint: list[DoubleDifference]
    internal_all: list[DoubleDifference]
    fit_unsampled: list[DoubleDifference]
    fit_targets: list[str]
    internal_targets: list[str]
    fit_components: tuple[str, ...]
    internal_components: tuple[str, ...]


def build_banks(data, seed: int) -> Banks:
    fit, internal = partition_components(data)
    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    fit_targets = sorted(t for t, c in component_of.items() if c in set(fit))
    internal_targets = sorted(t for t, c in component_of.items()
                              if c in set(internal))
    built = build_observations(data, fit_targets + internal_targets)
    observations = built["observations"]
    fit_set, internal_set = set(fit_targets), set(internal_targets)

    fit_effects = target_effects([o for o in observations if o.target in fit_set])
    internal_effects = target_effects(
        [o for o in observations if o.target in internal_set])
    fit_rows = double_differences(fit_effects)
    internal_rows = double_differences(internal_effects)

    fit_keys = {row.key for row in fit_rows}
    repeated, disjoint = split_by_key_overlap(internal_rows, fit_keys)

    # A held-out slice of FIT rows, for the target-key shortcut check (gate 9).
    rng = np.random.default_rng(stable_seed("stageT-fit-holdout", seed))
    order = rng.permutation(len(fit_rows))
    holdout = {int(i) for i in order[:max(1, len(fit_rows) // 10)]}
    train = [row for position, row in enumerate(fit_rows) if position not in holdout]
    fit_unsampled = [row for position, row in enumerate(fit_rows)
                     if position in holdout]
    return Banks(train=train, internal_repeated=repeated,
                 internal_disjoint=disjoint, internal_all=internal_rows,
                 fit_unsampled=fit_unsampled, fit_targets=fit_targets,
                 internal_targets=internal_targets, fit_components=fit,
                 internal_components=internal)


class Trainer:
    def __init__(self, data, table: ProteinTable, banks: Banks, arm: ArmConfig,
                 device: torch.device, seed: int) -> None:
        self.data = data
        self.table = table
        self.banks = banks
        self.arm = arm
        self.device = device
        self.seed = seed
        all_targets = banks.fit_targets + banks.internal_targets
        self.shuffle = shuffled_map(data, all_targets, seed)
        self.matched_wrong = matched_wrong_map(
            data, table, all_targets, banks.fit_targets)
        torch.manual_seed(stable_seed("stageT-init", seed, arm.name) % (2 ** 31))
        self.config = DiscriminatorConfig(
            descriptor_dim=len(banks.train[0].descriptor),
            protein_dim=table.width, mode=arm.mode)
        self.model = DoubleDifferenceModel(self.config).to(device)
        self.optimizer = torch.optim.AdamW(
            self.model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
        self.schedule = torch.optim.lr_scheduler.CosineAnnealingLR(
            self.optimizer, T_max=STEPS)
        rows = banks.train
        if arm.shuffle_labels:
            rows = shuffle_within_key(rows, seed)
        self.train_rows = rows

    def protein_for(self, target: str, override: str | None = None) -> str:
        source = override or self.arm.protein_source
        if source == "correct":
            return target
        if source == "shuffled":
            return self.shuffle[target]
        if source == "matched_wrong":
            return self.matched_wrong[target]
        raise ValueError(f"unknown protein source: {source}")

    def tensors(self, rows: list[DoubleDifference], override: str | None = None):
        descriptor = torch.tensor(
            np.asarray([row.descriptor for row in rows], dtype=np.float32),
            device=self.device)
        left = torch.tensor(np.stack([
            self.table.get(self.protein_for(row.target_left, override))
            for row in rows]), device=self.device)
        right = torch.tensor(np.stack([
            self.table.get(self.protein_for(row.target_right, override))
            for row in rows]), device=self.device)
        truth = torch.tensor(
            np.asarray([row.value for row in rows], dtype=np.float32),
            device=self.device)
        return descriptor, left, right, truth

    def train(self, steps: int, log_every: int = 250) -> list[dict]:
        history: list[dict] = []
        self.model.train()
        rows = self.train_rows
        for step in range(steps):
            # The arm name is deliberately NOT in the batch seed: every arm
            # must see the identical sequence of training rows so the only
            # difference between arms is the arm's own configuration. Arm F
            # then sees the same rows with permuted labels, which is exactly
            # the control it is meant to be.
            rng = np.random.default_rng(stable_seed(
                "stageT-batch", self.seed, step))
            picked = rng.choice(len(rows), size=min(BATCH, len(rows)),
                                replace=False)
            batch = [rows[int(i)] for i in picked]
            descriptor, left, right, truth = self.tensors(batch)
            prediction = self.model(descriptor, left, right)
            loss = torch.nn.functional.huber_loss(prediction, truth,
                                                  delta=HUBER_BETA)
            self.optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), GRADIENT_CLIP)
            self.optimizer.step()
            self.schedule.step()
            if step % log_every == 0 or step == steps - 1:
                history.append({"step": step, "loss": float(loss.detach()),
                                "lr": float(self.schedule.get_last_lr()[0])})
        return history

    @torch.no_grad()
    def evaluate(self, rows: list[DoubleDifference],
                 override: str | None = None) -> dict:
        self.model.eval()
        if not rows:
            return {"row_id": [], "truth": [], "prediction": [], "key": [],
                    "component_left": [], "component_right": [],
                    "activity_cliff": [], "cross_component": []}
        predictions: list[float] = []
        for start in range(0, len(rows), 512):
            chunk = rows[start:start + 512]
            descriptor, left, right, _ = self.tensors(chunk, override)
            predictions.extend(
                self.model(descriptor, left, right).cpu().numpy().tolist())
        return {
            "row_id": [row.row_id for row in rows],
            "truth": [row.value for row in rows],
            "prediction": predictions,
            "key": [row.key for row in rows],
            "component_left": [row.component_left for row in rows],
            "component_right": [row.component_right for row in rows],
            "activity_cliff": [row.activity_cliff for row in rows],
            "cross_component": [row.cross_component for row in rows],
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", required=True, choices=sorted(ARMS))
    parser.add_argument("--steps", type=int, default=STEPS)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--output", type=Path, default=HERE / "runs")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    arm = ARMS[args.arm]
    destination = args.output / arm.name
    if destination.exists() and not args.force:
        raise SystemExit(f"{destination} exists; pass --force to overwrite")
    destination.mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    data, seal = load_governed()
    if not seal["isolation"]["physically_isolated"]:
        raise SystemExit("refusing to train without the physical split view")
    table = ProteinTable(data)
    banks = build_banks(data, args.seed)
    trainer = Trainer(data, table, banks, arm, device, args.seed)

    started = time.time()
    history = trainer.train(args.steps)
    elapsed = time.time() - started

    rows = {
        "internal_repeated_correct": trainer.evaluate(banks.internal_repeated),
        "internal_disjoint_correct": trainer.evaluate(banks.internal_disjoint),
        "internal_all_correct": trainer.evaluate(banks.internal_all),
        "fit_unsampled_correct": trainer.evaluate(banks.fit_unsampled),
    }
    if arm.mode == "protein" and arm.protein_source == "correct":
        # Paired counterfactuals: identical rows, only the protein input changes.
        rows["internal_all_matched_wrong"] = trainer.evaluate(
            banks.internal_all, "matched_wrong")
        rows["internal_all_shuffled"] = trainer.evaluate(
            banks.internal_all, "shuffled")
        rows["internal_repeated_matched_wrong"] = trainer.evaluate(
            banks.internal_repeated, "matched_wrong")
        rows["internal_disjoint_matched_wrong"] = trainer.evaluate(
            banks.internal_disjoint, "matched_wrong")

    record = {
        "schema": "MetaSieve.StageT.T2Arm.v1",
        "arm": asdict(arm),
        "seed": args.seed,
        "steps": args.steps,
        "elapsed_seconds": elapsed,
        "device": verify_gpu(),
        "meta_test": seal,
        "selection": {
            "rule": "fixed_budget_final_checkpoint",
            "leak_free": True,
            "internal_validation_read_during_training": False,
            "disclosure": ("no checkpoint selection of any kind; every arm "
                           "trains the same fixed number of steps"),
        },
        "banks": {
            "train": effective_independent_units(banks.train),
            "internal_repeated": effective_independent_units(banks.internal_repeated),
            "internal_disjoint": effective_independent_units(banks.internal_disjoint),
            "internal_all": effective_independent_units(banks.internal_all),
            "fit_unsampled": effective_independent_units(banks.fit_unsampled),
            "fit_components": len(banks.fit_components),
            "internal_components": len(banks.internal_components),
        },
        "hyperparameters": {
            "learning_rate": LEARNING_RATE, "weight_decay": WEIGHT_DECAY,
            "batch": BATCH, "huber_beta": HUBER_BETA,
            "gradient_clip": GRADIENT_CLIP,
        },
        "parameters": parameter_report(trainer.model),
        "history": history,
    }
    (destination / "RUN.json").write_text(
        json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    for name, value in rows.items():
        (destination / f"{name}.rows.json").write_text(
            json.dumps(value) + "\n", encoding="utf-8")
    print(json.dumps({"arm": arm.name, "elapsed_seconds": round(elapsed, 1),
                      "final_loss": history[-1]["loss"],
                      "parameters": record["parameters"]["total"],
                      "train_rows": len(banks.train)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
