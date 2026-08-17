"""Why one inner step helps `A1` and destroys `A0`: the step-size conditioning.

At k=1 the inner loop has a closed-form effect that makes the whole pattern
predictable. With one support point and squared error,

    L = (p - y)^2,   dL/dw = 2 (p - y) h,   p' = p - lr * 2 (p - y) |h|^2

so writing `alpha = 2 * lr * |h|^2`, the post-step residual is

    (p' - y) = (p - y) (1 - alpha)

`alpha < 1` converges, `alpha = 2` flips the residual exactly (oscillation
without progress), and `alpha > 2` diverges. `|h|` is the norm of the readout's
hidden activation, which training controls.

This script measures `alpha` for each arm on real meta_val episodes. It decides
nothing — the gates are already settled — but it names the cause of the
inner-step sweep's shape, which the preregistration requires whenever
adaptation misbehaves.

No training, no meta_test.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData                              # noqa: E402
from scripts.stageR6_compare_arms import load_arm                     # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, normalized_episode, training_label_scale,
)
from tools.research.stageA_innerloop.inner_loop import AdaptationConfig  # noqa: E402
from tools.research.stageA_innerloop.train_meta import (              # noqa: E402
    align_atoms, encode_parts, episode_tensors,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
ARMS = ("A0", "A1", "A2")
OUT = Path(__file__).resolve().parent / "CONDITIONING.json"


def hidden_activations(model, task) -> torch.Tensor:
    """The readout's hidden vector `h`, whose norm sets the inner step scale."""
    body = model.interaction_head
    first = F.linear(task.support_hidden, body[0].weight, body[0].bias)
    return body[1](first)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--inner-lr", type=float, default=0.1)
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    arguments = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    specs = data.fixed_nested_episode_banks(
        "meta_val", (1,), 16, 2, 73101, None)[1]
    episodes = [compact_episode(normalized_episode(
        data.materialize(spec), label_scale)) for spec in specs]

    payload = {"schema": "MetaSieve.StageA.Conditioning.v1",
               "date": "2026-08-17", "inner_lr": arguments.inner_lr,
               "support_size": 1, "episodes": len(episodes),
               "definition": "alpha = 2 * inner_lr * |h|^2; "
                             "residual multiplier = |1 - alpha|",
               "arms": {}, "meta_test": data.seal_record()}

    for arm in ARMS:
        checkpoint = arguments.stage / arm / "checkpoint.pt"
        model, _, _ = load_arm(checkpoint, data, arguments.device)
        model.eval()
        alphas = []
        for episode in episodes:
            parts = align_atoms(episode_tensors(
                model, episode, arguments.device, torch.float32))
            with torch.no_grad():
                task = encode_parts(model, parts)
                h = hidden_activations(model, task)
                alphas.append(float(2.0 * arguments.inner_lr
                                    * h.square().sum(-1).mean()))
        values = np.asarray(alphas)
        payload["arms"][arm] = {
            "alpha_mean": float(values.mean()),
            "alpha_median": float(np.median(values)),
            "alpha_p10": float(np.quantile(values, 0.10)),
            "alpha_p90": float(np.quantile(values, 0.90)),
            "residual_multiplier_at_mean": float(abs(1.0 - values.mean())),
            "fraction_overshooting": float((values > 1.0).mean()),
            "fraction_oscillating": float((values > 2.0).mean()),
        }
        print(f"{arm}: alpha mean {values.mean():.3f} "
              f"median {np.median(values):.3f} "
              f"[p10 {np.quantile(values, 0.10):.3f}, "
              f"p90 {np.quantile(values, 0.90):.3f}]  "
              f"|1-alpha| {abs(1.0 - values.mean()):.3f}  "
              f"overshoot {(values > 1.0).mean():.2f}  "
              f"oscillate {(values > 2.0).mean():.2f}")
        del model
        if arguments.device.startswith("cuda"):
            torch.cuda.empty_cache()

    OUT.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
