"""Phase 2: the single-seed BindingDB discriminator for the SAR-field family.

Five matched arms, one code path, one seed, one fixed step budget, no
checkpoint selection of any kind:

    A  ligand_only        protein-free potential (learned constant response)
    B  protein            protein-conditioned potential, no counterfactual loss
    C  protein_cf         B plus a correct-vs-hard-wrong-protein relational loss
    D  protein_shuffled   B trained under a stable cross-component permutation
                          of protein identity
    E  label_shuffled     B trained on within-target permuted labels

Training reads only the **fit** components of `meta_train`.  Reporting reads
only the withheld **internal-validation** components of `meta_train`.  `meta_val`
is not opened anywhere in this stage, and the `meta_test` label artifact is not
present on the mounted surface at all.

There is deliberately no early stopping and no best-checkpoint rule: the
reporting population is the internal-validation components, so selecting on them
would leak, and selecting on `meta_val` is the defect this repository already
measured at ~0.62 pK^2.  Every arm trains for exactly `--steps` updates and the
final parameters are evaluated.

Ordinary gradient training only.  No ridge, no closed-form solver, no
pseudoinverse, no test-time query-label gradient, no second dataset.
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
from tools.research.stageS_sar_field.features import (
    LigandFeatureStore, ProteinFeatureStore, cross_component_permutation,
    hard_wrong_protein_map, relabel, within_target_label_shuffle,
)
from tools.research.stageS_sar_field.field import (
    FieldConfig, build_field, parameter_report,
)
from tools.research.stageS_sar_field.pairs import (
    PairSpec, build_target_pairs, component_of_target, load_data,
    target_balanced_bank,
)

HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Frozen before any arm trained (see PREREGISTRATION.md).
SEED = 20260819
STEPS = 4000
TARGETS_PER_BATCH = 16
PAIRS_PER_TARGET_IN_BATCH = 8
TRAIN_PAIRS_PER_TARGET = 512
EVAL_PAIRS_PER_TARGET = 512
LEARNING_RATE = 3e-4
WEIGHT_DECAY = 1e-4
HUBER_BETA = 1.0
SIGN_WEIGHT = 0.2
SIGN_TEMPERATURE = 1.0
COUNTERFACTUAL_WEIGHT = 0.5
COUNTERFACTUAL_MARGIN = 0.25
GRADIENT_CLIP = 5.0
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ArmConfig:
    name: str
    protein_conditioned: bool
    counterfactual_weight: float
    permute_proteins: bool
    shuffle_labels: bool


ARMS = {
    "A_ligand_only": ArmConfig("A_ligand_only", False, 0.0, False, False),
    "B_protein": ArmConfig("B_protein", True, 0.0, False, False),
    "C_protein_cf": ArmConfig("C_protein_cf", True, COUNTERFACTUAL_WEIGHT,
                              False, False),
    "D_protein_shuffled": ArmConfig("D_protein_shuffled", True, 0.0, True, False),
    "E_label_shuffled": ArmConfig("E_label_shuffled", True, 0.0, False, True),
}


def verify_gpu() -> dict:
    if not torch.cuda.is_available():
        raise SystemExit("refusing to train without CUDA (repository rule)")
    device = torch.device("cuda")
    probe = torch.zeros(8, device=device)
    return {
        "cuda_available": True,
        "device": torch.cuda.get_device_name(0),
        "probe_device": str(probe.device),
        "torch": torch.__version__,
    }


@dataclass
class Banks:
    train_by_target: dict[str, list[PairSpec]]
    train_targets: list[str]
    fit_unsampled: list[PairSpec]
    internal_same_panel: list[PairSpec]
    internal_cross_panel: list[PairSpec]
    fit_components: tuple[str, ...]
    internal_components: tuple[str, ...]
    donor_train: dict[str, str]
    donor_internal: dict[str, str]
    protein_permutation: dict[str, str]
    label_shuffle: dict[int, float]


def build_banks(data, proteins: ProteinFeatureStore, seed: int) -> Banks:
    fit, internal = partition_components(data)
    component = component_of_target(data)
    fit_targets = sorted(t for t, c in component.items() if c in set(fit))
    internal_targets = sorted(t for t, c in component.items() if c in set(internal))

    fit_pairs = build_target_pairs(data, "meta_train", fit_targets)
    # Primary supervision: the highest-confidence comparable pairs only.
    # Cross-panel differences are never pooled into the training labels; they
    # are reported as their own evaluation stratum.
    train_by_target = {
        target: [spec for spec in specs if spec.same_panel]
        for target, specs in fit_pairs.items()
    }
    train_by_target = {target: specs for target, specs in train_by_target.items()
                       if specs}
    available = train_by_target
    train_by_target = _cap(train_by_target, seed, TRAIN_PAIRS_PER_TARGET,
                           "sar-field-train")
    # Same-panel fit pairs the target-balanced draw did not take. Same targets,
    # different pairs: a pure training-health monitor that shows whether an arm
    # has overfitted its own bank. Preregistered as DIAGNOSTIC ONLY -- it never
    # selects a checkpoint, an arm or a threshold.
    taken = {(spec.a, spec.b) for specs in train_by_target.values()
             for spec in specs}
    fit_unsampled = _cap({target: [spec for spec in specs
                                   if (spec.a, spec.b) not in taken]
                          for target, specs in available.items()},
                         seed, 64, "sar-field-fit-unsampled")

    internal_pairs = build_target_pairs(data, "meta_train", internal_targets)
    same_panel = _cap({t: [s for s in v if s.same_panel]
                       for t, v in internal_pairs.items()},
                      seed, EVAL_PAIRS_PER_TARGET, "sar-field-eval-same-panel")
    cross_panel = _cap({t: [s for s in v if not s.same_panel]
                        for t, v in internal_pairs.items()},
                       seed, EVAL_PAIRS_PER_TARGET, "sar-field-eval-cross-panel")

    donor_train = hard_wrong_protein_map(
        data, proteins, sorted(train_by_target), fit_targets)
    donor_internal = hard_wrong_protein_map(
        data, proteins, internal_targets, fit_targets)
    permutation = cross_component_permutation(
        data, fit_targets + internal_targets, seed)
    labels = within_target_label_shuffle(data, fit_targets, seed)
    return Banks(
        train_by_target=train_by_target,
        train_targets=sorted(train_by_target),
        fit_unsampled=[s for t in sorted(fit_unsampled) for s in fit_unsampled[t]],
        internal_same_panel=[s for t in sorted(same_panel) for s in same_panel[t]],
        internal_cross_panel=[s for t in sorted(cross_panel)
                              for s in cross_panel[t]],
        fit_components=fit, internal_components=internal,
        donor_train=donor_train, donor_internal=donor_internal,
        protein_permutation=permutation, label_shuffle=labels)


def _cap(by_target: dict[str, list[PairSpec]], seed: int, per_target: int,
         namespace: str) -> dict[str, list[PairSpec]]:
    capped = target_balanced_bank(
        {t: v for t, v in by_target.items() if v}, seed, per_target, namespace)
    out: dict[str, list[PairSpec]] = {}
    for spec in capped:
        out.setdefault(spec.target, []).append(spec)
    return out


def pair_id(spec: PairSpec) -> str:
    return f"{spec.target}:{spec.a}:{spec.b}"


class Trainer:
    def __init__(self, data, ligands: LigandFeatureStore,
                 proteins: ProteinFeatureStore, banks: Banks,
                 arm: ArmConfig, device: torch.device, seed: int) -> None:
        self.data = data
        self.ligands = ligands
        self.proteins = proteins
        self.banks = banks
        self.arm = arm
        self.device = device
        self.seed = seed
        torch.manual_seed(stable_seed("sar-field-init", seed, arm.name) % (2 ** 31))
        self.config = FieldConfig(protein_conditioned=arm.protein_conditioned)
        self.field = build_field(self.config).to(device)
        self.optimizer = torch.optim.AdamW(
            self.field.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
        self.schedule = torch.optim.lr_scheduler.CosineAnnealingLR(
            self.optimizer, T_max=STEPS)
        self.train_bank = self._training_bank()

    def _training_bank(self) -> dict[str, list[PairSpec]]:
        bank = self.banks.train_by_target
        if self.arm.shuffle_labels:
            return {target: relabel(specs, self.banks.label_shuffle)
                    for target, specs in bank.items()}
        return bank

    def protein_for(self, target: str) -> str:
        """The protein this arm associates with a target."""
        if self.arm.permute_proteins:
            return self.banks.protein_permutation[target]
        return target

    # -- forward -----------------------------------------------------------

    def coordinates(self, specs: list[PairSpec]) -> tuple[torch.Tensor, torch.Tensor]:
        """phi for both sides of every pair, encoding each ligand exactly once.

        Pairs inside one target reuse the same ligands heavily, so deduplicating
        before the graph encoder is a pure cost saving: gathering the same
        ligand twice would return the same tensor.
        """
        cells: list[int] = []
        index: dict[int, int] = {}
        for spec in specs:
            for cell in (spec.a, spec.b):
                if cell not in index:
                    index[cell] = len(cells)
                    cells.append(cell)
        atoms, bonds, mask, prints = self.ligands.gather(cells, self.device)
        coordinate = self.field.phi(atoms, bonds, mask, prints)
        rows_a = torch.tensor([index[spec.a] for spec in specs],
                              dtype=torch.long, device=self.device)
        rows_b = torch.tensor([index[spec.b] for spec in specs],
                              dtype=torch.long, device=self.device)
        return coordinate.index_select(0, rows_a), coordinate.index_select(0, rows_b)

    def responses(self, protein_targets: list[str]) -> torch.Tensor:
        unique = sorted(set(protein_targets))
        index = {key: position for position, key in enumerate(unique)}
        pooled, residues, protein_mask = self.proteins.gather(unique, self.device)
        response = self.field.alpha(pooled, residues, protein_mask)
        rows = torch.tensor([index[key] for key in protein_targets],
                            dtype=torch.long, device=self.device)
        return response.index_select(0, rows)

    def predict(self, specs: list[PairSpec], protein_targets: list[str]
                ) -> torch.Tensor:
        phi_a, phi_b = self.coordinates(specs)
        return self.field(phi_a, phi_b, self.responses(protein_targets))

    # -- losses ------------------------------------------------------------

    @staticmethod
    def _huber(prediction: torch.Tensor, truth: torch.Tensor) -> torch.Tensor:
        return torch.nn.functional.huber_loss(
            prediction, truth, delta=HUBER_BETA)

    @staticmethod
    def _sign(prediction: torch.Tensor, truth: torch.Tensor) -> torch.Tensor:
        return torch.nn.functional.softplus(
            -torch.sign(truth) * prediction / SIGN_TEMPERATURE).mean()

    def _counterfactual(self, specs: list[PairSpec], truth: torch.Tensor,
                        correct: torch.Tensor, phi_a: torch.Tensor,
                        phi_b: torch.Tensor) -> torch.Tensor:
        """Relational loss: the correct protein must fit better than a hard wrong one.

        **Only the protein input is replaced.**  `phi_a` and `phi_b` are the
        very tensors used for the correct-protein prediction, so the ligand
        path is bitwise identical on both sides and the contrast is attributable
        to the protein alone.  The hinge saturates once the wrong protein is
        `COUNTERFACTUAL_MARGIN` pK^2 worse, which bounds the reward for simply
        making the wrong-protein prediction explode.
        """
        donors = [self.banks.donor_train[spec.target] for spec in specs]
        wrong = self.field(phi_a, phi_b, self.responses(donors))
        gap = (wrong - truth) ** 2 - (correct - truth) ** 2
        return torch.nn.functional.softplus(COUNTERFACTUAL_MARGIN - gap).mean()

    # -- loop --------------------------------------------------------------

    def sample(self, step: int) -> list[PairSpec]:
        rng = np.random.default_rng(stable_seed(
            "sar-field-batch", self.seed, self.arm.name, step))
        targets = self.banks.train_targets
        picked = rng.choice(len(targets), size=min(TARGETS_PER_BATCH,
                                                   len(targets)), replace=False)
        out: list[PairSpec] = []
        for position in picked:
            specs = self.train_bank[targets[int(position)]]
            take = min(PAIRS_PER_TARGET_IN_BATCH, len(specs))
            order = rng.choice(len(specs), size=take, replace=False)
            out.extend(specs[int(i)] for i in order)
        return out

    def train(self, steps: int, log_every: int = 200) -> list[dict]:
        history: list[dict] = []
        self.field.train()
        for step in range(steps):
            specs = self.sample(step)
            truth = torch.tensor([spec.delta_y for spec in specs],
                                 dtype=torch.float32, device=self.device)
            phi_a, phi_b = self.coordinates(specs)
            response = self.responses([self.protein_for(spec.target)
                                       for spec in specs])
            prediction = self.field(phi_a, phi_b, response)
            regression = self._huber(prediction, truth)
            sign = self._sign(prediction, truth)
            loss = regression + SIGN_WEIGHT * sign
            counterfactual = torch.zeros((), device=self.device)
            if self.arm.counterfactual_weight > 0.0:
                counterfactual = self._counterfactual(
                    specs, truth, prediction, phi_a, phi_b)
                loss = loss + self.arm.counterfactual_weight * counterfactual
            self.optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.field.parameters(), GRADIENT_CLIP)
            self.optimizer.step()
            self.schedule.step()
            if step % log_every == 0 or step == steps - 1:
                history.append({
                    "step": step,
                    "loss": float(loss.detach()),
                    "regression": float(regression.detach()),
                    "sign": float(sign.detach()),
                    "counterfactual": float(counterfactual.detach()),
                    "lr": float(self.schedule.get_last_lr()[0]),
                })
        return history

    # -- evaluation --------------------------------------------------------

    @torch.no_grad()
    def evaluate(self, specs: list[PairSpec], protein: str = "correct",
                 batch: int = 256) -> dict:
        """Prediction rows on an evaluation bank.

        `protein` selects which protein input is used and **nothing else**:
        the recipient's ligands, its `delta_y` and every other input are fixed.
        """
        self.field.eval()
        predictions: list[float] = []
        for start in range(0, len(specs), batch):
            chunk = specs[start:start + batch]
            if protein == "correct":
                targets = [self.protein_for(spec.target) for spec in chunk]
            elif protein == "hard_wrong":
                targets = [self.banks.donor_internal[spec.target] for spec in chunk]
            elif protein == "true_identity":
                targets = [spec.target for spec in chunk]
            else:
                raise ValueError(f"unknown protein condition: {protein}")
            predictions.extend(self.predict(chunk, targets).cpu().numpy().tolist())
        return {
            "pair_id": [pair_id(spec) for spec in specs],
            "delta_y": [spec.delta_y for spec in specs],
            "delta_hat": predictions,
            "target": [spec.target for spec in specs],
            "component": [spec.component for spec in specs],
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

    gpu = verify_gpu()
    device = torch.device("cuda")
    data = load_data()
    seal = data.seal_record()
    if not seal["isolation"]["physically_isolated"]:
        raise SystemExit("refusing to train without the physical split view")

    ligands = LigandFeatureStore(data)
    proteins = ProteinFeatureStore(data)
    banks = build_banks(data, proteins, args.seed)
    trainer = Trainer(data, ligands, proteins, banks, arm, device, args.seed)

    started = time.time()
    history = trainer.train(args.steps)
    elapsed = time.time() - started

    rows = {
        "internal_same_panel_correct": trainer.evaluate(
            banks.internal_same_panel, "correct"),
        "internal_cross_panel_correct": trainer.evaluate(
            banks.internal_cross_panel, "correct"),
        "fit_unsampled_correct": trainer.evaluate(
            banks.fit_unsampled, "correct"),
    }
    if arm.protein_conditioned:
        rows["internal_same_panel_hard_wrong"] = trainer.evaluate(
            banks.internal_same_panel, "hard_wrong")
        rows["internal_cross_panel_hard_wrong"] = trainer.evaluate(
            banks.internal_cross_panel, "hard_wrong")

    record = {
        "schema": "MetaSieve.StageS.Arm.v1",
        "arm": asdict(arm),
        "seed": args.seed,
        "steps": args.steps,
        "elapsed_seconds": elapsed,
        "gpu": gpu,
        "meta_test": seal,
        "selection": {
            "rule": "fixed_budget_final_checkpoint",
            "leak_free": True,
            "meta_val_read": False,
            "internal_validation_read_during_training": False,
            "disclosure": ("no checkpoint selection of any kind; every arm "
                           "trains for the same fixed number of steps and the "
                           "final parameters are evaluated"),
        },
        "partition": {
            "fit_components": len(banks.fit_components),
            "internal_validation_components": len(banks.internal_components),
            "train_targets": len(banks.train_targets),
            "train_pairs": sum(len(v) for v in banks.train_by_target.values()),
            "fit_unsampled_diagnostic_pairs": len(banks.fit_unsampled),
            "internal_same_panel_pairs": len(banks.internal_same_panel),
            "internal_cross_panel_pairs": len(banks.internal_cross_panel),
        },
        "hyperparameters": {
            "learning_rate": LEARNING_RATE, "weight_decay": WEIGHT_DECAY,
            "targets_per_batch": TARGETS_PER_BATCH,
            "pairs_per_target_in_batch": PAIRS_PER_TARGET_IN_BATCH,
            "train_pairs_per_target": TRAIN_PAIRS_PER_TARGET,
            "eval_pairs_per_target": EVAL_PAIRS_PER_TARGET,
            "huber_beta": HUBER_BETA, "sign_weight": SIGN_WEIGHT,
            "sign_temperature": SIGN_TEMPERATURE,
            "counterfactual_weight": arm.counterfactual_weight,
            "counterfactual_margin": COUNTERFACTUAL_MARGIN,
            "gradient_clip": GRADIENT_CLIP,
        },
        "parameters": parameter_report(trainer.field),
        "history": history,
    }
    (destination / "RUN.json").write_text(
        json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    for name, value in rows.items():
        (destination / f"{name}.rows.json").write_text(
            json.dumps(value) + "\n", encoding="utf-8")
    torch.save({"state_dict": trainer.field.state_dict(),
                "config": asdict(trainer.config) | {"dtype": "float32"}},
               destination / "field.pt")
    print(json.dumps({"arm": arm.name, "elapsed_seconds": round(elapsed, 1),
                      "final_loss": history[-1]["loss"],
                      "parameters": record["parameters"]["total"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
