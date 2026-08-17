"""Corrected inner-step conditioning for the Stage B arms (preregistered gate 8).

Uses `alpha = 2 * lr * (||h_support||^2 + bias)` with `bias = 1` when the bias is
adapted — the Stage A formula omitted that term. `alpha > 1` means one inner
step overshoots the support residual and flips its sign.

Stated precisely: this governs the **support** residual. It does not predict the
query error, which moves by `-2*lr*r_s*(h_q . h_s + bias)` and therefore varies
across the panel. `query_step_delta` reports that directly.

Written as its own script because the Stage B evaluator attached its alpha to
the wrong row (the last appended condition rather than `correct`), so the gate
was not measurable from `STAGE2_meta_val.json`.
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
from tools.research.stageA_innerloop.train_meta import (              # noqa: E402
    align_atoms, encode_parts, episode_tensors,
)
from tools.research.stageB_complementary.arms import (                # noqa: E402
    ADAPTED_BY_MODE, InnerStepSizes, StageBAdaptation,
)
from tools.research.stageB_complementary.residual import (            # noqa: E402
    conditioning_alpha, query_step_delta,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
OUT = Path(__file__).resolve().parent / "CONDITIONING_STAGEB.json"


def hidden(model, task, side: str) -> torch.Tensor:
    body = model.interaction_head
    source = task.support_hidden if side == "support" else task.query_hidden
    return body[1](F.linear(source, body[0].weight, body[0].bias))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--arms", nargs="+", required=True)
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    arguments = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    specs = data.fixed_nested_episode_banks("meta_val", (1,), 16, 2, 73101, None)[1]
    episodes = [compact_episode(normalized_episode(
        data.materialize(spec), label_scale)) for spec in specs]

    payload = {"schema": "MetaSieve.StageB.Conditioning.v1", "date": "2026-08-17",
               "definition": "alpha = 2*lr*(||h_support||^2 + bias_term)",
               "note": ("governs the support residual only; the query effect is "
                        "-2*lr*r_s*(h_q . h_s + bias) and varies per query"),
               "arms": {}, "meta_test": data.seal_record()}

    for arm in arguments.arms:
        checkpoint = arguments.stage / arm / "checkpoint.pt"
        blob = torch.load(checkpoint, map_location="cpu", weights_only=False)
        adaptation = StageBAdaptation.from_dict(blob["adaptation"])
        if adaptation.mode == "T":
            payload["arms"][arm] = {"mode": "T", "adapts": [],
                                    "note": "no inner loop; alpha undefined"}
            print(f"{arm}: no inner loop")
            continue
        scope = ADAPTED_BY_MODE[adaptation.mode]
        adapts_bias = any(n.endswith("bias") for n in scope)
        lr = adaptation.inner_lr
        if blob.get("inner_step_state") is not None:
            steps = InnerStepSizes(adaptation.inner_lr, adaptation.max_step)
            steps.load_state_dict(blob["inner_step_state"])
            lr = float(steps.weight_step())
        model, _, _ = load_arm(checkpoint, data, arguments.device)
        model.eval()
        alphas, deltas = [], []
        for episode in episodes:
            parts = align_atoms(episode_tensors(model, episode,
                                                arguments.device, torch.float32))
            with torch.no_grad():
                task = encode_parts(model, parts)
                hs = hidden(model, task, "support")
                hq = hidden(model, task, "query")
                alphas.append(float(conditioning_alpha(hs, lr, adapts_bias).mean()))
                residual = torch.ones(hs.shape[0], hs.shape[1],
                                      device=hs.device, dtype=hs.dtype)
                deltas.append(float(query_step_delta(
                    hq, hs, residual, lr, adapts_bias).abs().mean()))
        alpha = np.asarray(alphas)
        payload["arms"][arm] = {
            "mode": adaptation.mode, "adapts": list(scope),
            "effective_lr": lr, "bias_adapted": adapts_bias,
            "alpha_mean": float(alpha.mean()), "alpha_max": float(alpha.max()),
            "fraction_overshooting": float((alpha > 1.0).mean()),
            "fraction_oscillating": float((alpha > 2.0).mean()),
            "query_delta_per_unit_residual": float(np.mean(deltas)),
        }
        print(f"{arm} ({adaptation.mode}): alpha {alpha.mean():.4f} "
              f"max {alpha.max():.4f} overshoot {(alpha > 1).mean():.3f} "
              f"| query move per unit support residual "
              f"{np.mean(deltas):.5f}")
        del model
        if arguments.device.startswith("cuda"):
            torch.cuda.empty_cache()

    OUT.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
