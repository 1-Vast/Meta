"""The four Stage 2 arms, and the learned inner step size.

| arm | inner loop | inner target | transport |
|---|---|---|---|
| `T` | none | — | incumbent Tanimoto |
| `M` | weight+bias | raw support label | **disabled** |
| `H` | weight+bias | raw support label | incumbent Tanimoto |
| `C` | weight only | zero-shot + **complementary** residual | incumbent Tanimoto |

`H` reproduces Stage A's `A1`. `C` is the candidate: it differs from `T` by
exactly one additive term, so the contrast isolates the meta-correction.

The identity that makes `C` a clean comparison: the incumbent transport is
`shrink * sum_k w_qk r_k`, and since the softmax weights sum to one that equals
`shrink * (level_raw + sum_k w_qk (r_k - level_raw))`. So splitting the support
residual into a level and a centered part and recombining under one `shrink`
reproduces the incumbent transport **exactly**, and `C = T + meta_correction`.
A test pins that identity.

`C` adapts the weight only. The level already has an explicit closed term inside
the transport, so letting the bias adapt as well would fit the same level twice —
which the governing contract forbids.
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor

from tools.research.stageA_innerloop.inner_loop import base_weights, readout
from tools.research.stageB_complementary.residual import (
    centered_shape, decompose, level_term, query_transport,
)

MODES = ("T", "M", "H", "C")

# Which parameters each arm adapts. `C` holds the bias fixed because the level
# is already carried by an explicit term.
ADAPTED_BY_MODE = {
    "T": (),
    "M": ("interaction_head.2.weight", "interaction_head.2.bias"),
    "H": ("interaction_head.2.weight", "interaction_head.2.bias"),
    "C": ("interaction_head.2.weight",),
}


@dataclass
class StageBAdaptation:
    mode: str = "T"
    inner_steps: int = 1
    inner_lr: float = 0.1
    learned_step: bool = False
    max_step: float = 0.5
    first_order: bool = True

    def to_dict(self) -> dict:
        return {"mode": self.mode, "inner_steps": int(self.inner_steps),
                "inner_lr": float(self.inner_lr),
                "learned_step": bool(self.learned_step),
                "max_step": float(self.max_step),
                "first_order": bool(self.first_order),
                "adapted": list(ADAPTED_BY_MODE[self.mode])}

    @classmethod
    def from_dict(cls, payload: dict) -> "StageBAdaptation":
        return cls(mode=payload["mode"], inner_steps=int(payload["inner_steps"]),
                   inner_lr=float(payload["inner_lr"]),
                   learned_step=bool(payload["learned_step"]),
                   max_step=float(payload["max_step"]),
                   first_order=bool(payload["first_order"]))


class InnerStepSizes(nn.Module):
    """Two bounded positive scalars: one for shape, one for level.

    Meta-SGD-**inspired** and deliberately much smaller than Meta-SGD: that
    method learns a full per-parameter learning-rate vector, which here would be
    97 extra degrees of freedom fitted on as few as two support points. Two
    bounded scalars keep the mechanism identifiable and keep the failure mode
    (an unbounded step that diverges) structurally impossible.

    `step = max_step * sigmoid(raw)`, initialised so `step == inner_lr` exactly.
    Kept outside the model's own `state_dict` so previously recorded checkpoints
    still load strictly.
    """

    def __init__(self, initial: float, max_step: float,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if not 0.0 < initial < max_step:
            raise ValueError(
                f"initial step {initial} must lie strictly inside (0, {max_step})")
        self.max_step = float(max_step)
        ratio = initial / max_step
        raw = float(torch.logit(torch.tensor(ratio, dtype=torch.float64)))
        self.raw_weight = nn.Parameter(torch.tensor(raw, dtype=dtype))
        self.raw_bias = nn.Parameter(torch.tensor(raw, dtype=dtype))

    def weight_step(self) -> Tensor:
        return self.max_step * torch.sigmoid(self.raw_weight)

    def bias_step(self) -> Tensor:
        return self.max_step * torch.sigmoid(self.raw_bias)

    def for_parameter(self, name: str) -> Tensor:
        return self.bias_step() if name.endswith("bias") else self.weight_step()


def adapt_scoped(model, task, target: Tensor, config: StageBAdaptation,
                 steps: InnerStepSizes | None = None
                 ) -> tuple[dict[str, Tensor], list[float]]:
    """Inner loop over this arm's adapted scope, fitting `target` on support.

    `target` is an endpoint-space quantity, not a label: `H`/`M` pass the raw
    support label, `C` passes the zero-shot endpoint plus the complementary
    residual. Adaptation is a gradient computation, so it runs inside an
    explicit `enable_grad` and detaches on the way out when called under
    inference mode.
    """
    scope = ADAPTED_BY_MODE[config.mode]
    weights = {name: value for name, value in base_weights(model).items()}
    if not scope or config.inner_steps <= 0 or target.shape[-1] == 0:
        return weights, []
    fixed = {name: value for name, value in weights.items() if name not in scope}
    adapted = {name: weights[name] for name in scope}
    inference = not torch.is_grad_enabled()
    trace: list[float] = []
    with torch.enable_grad():
        for _ in range(int(config.inner_steps)):
            prediction = readout(model, task.support_hidden,
                                 task.support_additive,
                                 task.support_occupancy, {**fixed, **adapted})
            loss = F.mse_loss(prediction, target)
            trace.append(float(loss.detach()))
            gradients = torch.autograd.grad(
                loss, list(adapted.values()),
                create_graph=(not config.first_order) and not inference)
            if config.first_order:
                gradients = [gradient.detach() for gradient in gradients]
            updated = {}
            for (name, value), gradient in zip(adapted.items(), gradients):
                step = (steps.for_parameter(name) if steps is not None
                        else config.inner_lr)
                updated[name] = value - step * gradient
            adapted = updated
    if inference:
        adapted = {name: value.detach() for name, value in adapted.items()}
    return {**fixed, **adapted}, trace


def restrict(model, adapted: dict[str, Tensor], keep: str) -> dict[str, Tensor]:
    """Retain only part of the inner update, for the parameter-level ablation.

    Note what this does and does not mean. `keep="bias"` really is a pure level
    shift, because the bias moves every query identically. `keep="weight"` is
    **not** pure shape: it moves query `q` by `-2*lr*r_s*(h_q . h_s)`, whose
    within-episode mean is generally non-zero. Shape is measured separately by
    centering the query correction (`residual.centered_shape`), and the two
    ablations answer different questions.
    """
    base = base_weights(model)
    if keep == "both":
        return dict(adapted)
    if keep == "none":
        return dict(base)
    if keep not in {"bias", "weight"}:
        raise ValueError(f"unknown restriction {keep!r}")
    out = {}
    for name, value in adapted.items():
        is_bias = name.endswith("bias")
        out[name] = value if (keep == "bias") == is_bias else base[name]
    return out


def predict(model, parts: dict, task, config: StageBAdaptation,
            steps: InnerStepSizes | None = None, *,
            support_y_override: Tensor | None = None,
            keep: str = "both",
            disable_meta: bool = False,
            disable_transport: bool = False) -> dict:
    """One episode through the arm named by `config.mode`.

    Every additive term is returned separately so the report can ablate them
    without a second forward pass.
    """
    support_y = (parts["support_y"] if support_y_override is None
                 else support_y_override)
    base = base_weights(model)
    support_zero = readout(model, task.support_hidden, task.support_additive,
                           task.support_occupancy, base)
    query_zero = readout(model, task.query_hidden, task.query_additive,
                         task.query_occupancy, base)
    count = support_y.shape[-1]
    zeros = torch.zeros_like(query_zero)

    if count == 0:
        # k = 0 is the ordinary zero-shot path. No support-derived value is
        # constructed at all — not even a zero placeholder that the graph could
        # depend on.
        return {"prediction": query_zero, "zero_shot": query_zero,
                "level": zeros, "transport": zeros, "meta": zeros,
                "inner_trace": [], "support_zero": support_zero,
                "complementary": support_y.new_zeros(support_y.shape)}

    shrink = model.transport.shrinkage(count, support_y)
    scale = model.transport.similarity_scale
    parts_r = decompose(support_y, support_zero, parts["support_fingerprint"],
                        shrink, scale)
    residual = parts_r["residual"]

    level = zeros
    transport = zeros
    if not disable_transport and config.mode != "M":
        # Identical to the incumbent `shrink * sum_k w_qk r_k`, but split so the
        # level and neighbourhood shares are separately reportable.
        raw_level = residual.mean(-1, keepdim=True)
        level = (shrink * raw_level).expand_as(query_zero)
        transport = shrink * query_transport(
            parts["query_fingerprint"], parts["support_fingerprint"],
            residual - raw_level, scale)

    meta = zeros
    trace: list[float] = []
    if config.mode != "T" and not disable_meta:
        from tools.research.stageB_complementary.residual import inner_target
        target = inner_target(config.mode, support_y, support_zero, parts_r)
        fast, trace = adapt_scoped(model, task, target.detach(), config, steps)
        if keep != "both":
            fast = restrict(model, fast, keep)
        meta = readout(model, task.query_hidden, task.query_additive,
                       task.query_occupancy, fast) - query_zero

    prediction = query_zero + level + transport + meta
    return {"prediction": prediction, "zero_shot": query_zero,
            "level": level, "transport": transport, "meta": meta,
            "inner_trace": trace, "support_zero": support_zero,
            "complementary": parts_r["complementary"],
            "meta_shape": centered_shape(meta),
            "meta_level": meta.mean(-1, keepdim=True).expand_as(meta)}
