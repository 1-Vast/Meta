"""Attribution diagnostics for a trained arm: is the protein path inert or wrong?

A null gate result has two very different causes, and the repository has already
met both:

* **inert** -- the protein moves the prediction by almost nothing (the A2 exact
  operator: 0.0028 pK of query-specific content against a 0.884 pK label spread);
* **consistently wrong** -- the protein moves the prediction a lot, reproducibly,
  in a direction that carries no affinity information (Stage P: seed-to-seed
  cosine +0.316, alignment with truth +0.022).

Reporting a failure without saying which one it is has no diagnostic value, so
this module measures both directly on the trained field:

1. `protein_induced_spread` -- the sd, across proteins, of `dy_hat` for a fixed
   ligand pair.  Zero means the FiLM path collapsed;
2. `response_spread` -- the sd of `alpha(P)` across targets, which separates "the
   protein encoder produced nothing" from "the potential ignored what it
   produced";
3. `shift_alignment` -- the correlation between the protein-induced shift
   `dy_hat(P_correct) - dy_hat(P_wrong)` and the part of the truth the
   ligand-only field leaves unexplained.  This is the number that decides
   whether protein conditioning is pointing anywhere useful.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from tools.research.stageS_sar_field.features import (
    LigandFeatureStore, ProteinFeatureStore,
)
from tools.research.stageS_sar_field.field import FieldConfig, build_field
from tools.research.stageS_sar_field.pairs import load_data
from tools.research.stageS_sar_field.train import ARMS, Trainer, build_banks, pair_id

HERE = Path(__file__).resolve().parent
SEED = 20260819
DONOR_COUNT = 8


def load_arm(path: Path, arm_name: str, data, ligands, proteins, banks, device):
    payload = torch.load(path / "field.pt", map_location=device,
                         weights_only=False)
    arm = ARMS[arm_name]
    trainer = Trainer(data, ligands, proteins, banks, arm, device, SEED)
    trainer.field.load_state_dict(payload["state_dict"])
    trainer.field.eval()
    return trainer


@torch.no_grad()
def diagnose(trainer: Trainer, banks, specs, donors: list[str]) -> dict:
    device = trainer.device
    correct = []
    alternates = []
    for start in range(0, len(specs), 256):
        chunk = specs[start:start + 256]
        phi_a, phi_b = trainer.coordinates(chunk)
        base = trainer.field(phi_a, phi_b,
                             trainer.responses([s.target for s in chunk]))
        correct.append(base.cpu().numpy())
        rows = []
        for donor in donors:
            response = trainer.responses([donor] * len(chunk))
            rows.append(trainer.field(phi_a, phi_b, response).cpu().numpy())
        alternates.append(np.stack(rows, axis=0))
    correct = np.concatenate(correct)
    alternates = np.concatenate(alternates, axis=1)
    stacked = np.concatenate([correct[None, :], alternates], axis=0)

    targets = sorted({spec.target for spec in specs})
    pooled, residues, mask = trainer.proteins.gather(targets, device)
    response = trainer.field.alpha(pooled, residues, mask).cpu().numpy()

    truth = np.asarray([spec.delta_y for spec in specs], dtype=np.float64)
    hard_wrong = []
    for start in range(0, len(specs), 256):
        chunk = specs[start:start + 256]
        phi_a, phi_b = trainer.coordinates(chunk)
        donor_targets = [banks.donor_internal[s.target] for s in chunk]
        hard_wrong.append(trainer.field(
            phi_a, phi_b, trainer.responses(donor_targets)).cpu().numpy())
    hard_wrong = np.concatenate(hard_wrong)
    shift = correct - hard_wrong
    return {
        "pairs": int(correct.size),
        "label_spread_pK": float(truth.std()),
        "prediction_spread_pK": float(correct.std()),
        "protein_induced_spread_pK": float(stacked.std(axis=0).mean()),
        "protein_induced_spread_over_label_spread": float(
            stacked.std(axis=0).mean() / max(truth.std(), 1e-12)),
        "hard_wrong_shift_sd_pK": float(shift.std()),
        "response_vector_sd": float(response.std(axis=0).mean()),
        "response_vector_pairwise_cosine": float(_mean_cosine(response)),
        "correct": correct.tolist(),
        "hard_wrong": hard_wrong.tolist(),
        "shift": shift.tolist(),
    }


def _mean_cosine(matrix: np.ndarray) -> float:
    normalized = matrix / np.maximum(
        np.linalg.norm(matrix, axis=1, keepdims=True), 1e-12)
    gram = normalized @ normalized.T
    upper = gram[np.triu_indices_from(gram, k=1)]
    return float(upper.mean()) if upper.size else float("nan")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, default=HERE / "runs")
    parser.add_argument("--output", type=Path, default=HERE / "ATTRIBUTION.json")
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    data = load_data()
    ligands = LigandFeatureStore(data)
    proteins = ProteinFeatureStore(data)
    banks = build_banks(data, proteins, SEED)
    specs = banks.internal_same_panel
    truth = np.asarray([spec.delta_y for spec in specs], dtype=np.float64)

    donor_pool = sorted({banks.donor_internal[spec.target] for spec in specs})
    donors = donor_pool[:DONOR_COUNT]

    out: dict = {
        "schema": "MetaSieve.StageS.Attribution.v1",
        "donor_proteins_used": donors,
        "arms": {},
    }
    ligand_only_error = None
    for arm_name in sorted(ARMS):
        path = args.runs / arm_name
        if not (path / "field.pt").exists():
            continue
        trainer = load_arm(path, arm_name, data, ligands, proteins, banks, device)
        record = diagnose(trainer, banks, specs, donors)
        if arm_name == "A_ligand_only":
            ligand_only_error = truth - np.asarray(record["correct"])
        out["arms"][arm_name] = record

    for arm_name, record in out["arms"].items():
        shift = np.asarray(record.pop("shift"))
        record.pop("correct")
        record.pop("hard_wrong")
        if ligand_only_error is not None and shift.std() > 1e-12:
            record["shift_alignment_with_ligand_only_residual"] = float(
                np.corrcoef(shift, ligand_only_error)[0, 1])
        else:
            record["shift_alignment_with_ligand_only_residual"] = float("nan")
        record["shift_alignment_with_truth"] = (
            float(np.corrcoef(shift, truth)[0, 1]) if shift.std() > 1e-12
            else float("nan"))

    args.output.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(out, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
