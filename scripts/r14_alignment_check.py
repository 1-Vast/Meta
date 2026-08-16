"""R14 diagnostic 2: does the ranking term vanish at the regression optimum?

This is the derivation check that can kill the R14 design before any GPU
time. The RCR alignment proof ([arXiv:2211.01494]) is stated for a sigmoid
link with binary relevance. Ours is squared error on continuous pK, so the
identity has to be re-derived and then verified numerically.

The claim, for a within-target panel with labels `y` and predictions `s`,
a fixed shift `m` below the label range, `T(x) = x - m > 0`, and weights
`w_i = y_i - m`:

    ListCE(s, y) = -(1/Σ_j w_j) Σ_i w_i · log[ T(s_i) / Σ_j T(s_j) ]

    ∂ListCE/∂s_k = -(1/Σ_j w_j) [ w_k/T(s_k) - Σ_j w_j / Σ_j T(s_j) ]

which is zero for every k exactly when `T(s) ∝ w`, i.e. when the centered
prediction is proportional to the centered label. The squared-error optimum
`s = y` is one such point — so **at the regression optimum the ranking term
contributes exactly zero gradient.** The two objectives share a minimum
instead of competing for it.

RankNet and hinge-margin losses do *not* have this property: at `s = y` they
still push for wider margins, and that residual gradient is the mechanism
Phase 2 measured as a loss of within-target correlation `r`.

Run::

    python -m scripts.r14_alignment_check
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch

SHIFT = 2.0             # fixed, below the pK label range; not a swept knob


def list_ce(scores: torch.Tensor, labels: torch.Tensor,
            shift: float = SHIFT) -> torch.Tensor:
    """Regression-compatible listwise cross-entropy over one target's panel."""
    weights = labels - shift
    transformed = (scores - shift).clamp_min(1e-6)
    normalized = transformed / transformed.sum()
    return -(weights * normalized.log()).sum() / weights.sum()


def ranknet(scores: torch.Tensor, labels: torch.Tensor) -> torch.Tensor:
    """Pairwise softplus on score differences — the incumbent shape loss."""
    rows, cols = torch.triu_indices(len(labels), len(labels), 1)
    delta_y = labels[rows] - labels[cols]
    delta_s = scores[rows] - scores[cols]
    comparable = delta_y != 0
    signed = torch.sign(delta_y[comparable]) * delta_s[comparable]
    return torch.nn.functional.softplus(-signed).mean()


def hinge(scores: torch.Tensor, labels: torch.Tensor,
          margin: float = 0.1) -> torch.Tensor:
    """The R12 margin-ranking form."""
    rows, cols = torch.triu_indices(len(labels), len(labels), 1)
    delta_y = labels[rows] - labels[cols]
    delta_s = scores[rows] - scores[cols]
    comparable = delta_y != 0
    signed = torch.sign(delta_y[comparable]) * delta_s[comparable]
    return (margin - signed).clamp_min(0.0).mean()


def gradient_norm(loss_fn, scores: torch.Tensor, labels: torch.Tensor) -> float:
    scores = scores.clone().requires_grad_(True)
    loss_fn(scores, labels).backward()
    return float(scores.grad.norm())


def operating_point(labels: torch.Tensor, r: float, dispersion: float,
                    generator: torch.Generator) -> torch.Tensor:
    """A prediction with the measured within-target statistics of the incumbent.

    Alignment at `s = y` says what a term does when the model is already
    perfect. What decides whether a term is *useful* is how much gradient it
    supplies at the model's actual operating point, which Phase 2 measured as
    `r ~ 0.2` and `sd_p/sd_y ~ 0.2`.
    """
    centered = labels - labels.mean()
    spread = centered.std()
    noise = torch.randn(len(labels), generator=generator, dtype=labels.dtype)
    noise = (noise - noise.mean()) / noise.std()
    signal = r * centered / spread + (1.0 - r ** 2) ** 0.5 * noise
    return labels.mean() + dispersion * spread * signal / signal.std()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panels", type=int, default=200)
    parser.add_argument("--panel-size", type=int, default=16)
    parser.add_argument("--operating-r", type=float, default=0.2)
    parser.add_argument("--operating-dispersion", type=float, default=0.2)
    parser.add_argument("--output", type=Path, default=None)
    arguments = parser.parse_args()

    generator = torch.Generator().manual_seed(20260816)
    at_optimum: dict[str, list[float]] = {"list_ce": [], "ranknet": [], "hinge": []}
    at_operating: dict[str, list[float]] = {"list_ce": [], "ranknet": [],
                                            "hinge": [], "squared_error": []}
    proportional: list[float] = []
    for _ in range(arguments.panels):
        # pK labels on a realistic double-cold support: mean ~7.5, sd ~0.85
        labels = 7.5 + 0.85 * torch.randn(arguments.panel_size, generator=generator,
                                          dtype=torch.float64)
        labels = labels.clamp(4.0, 11.0)
        at_optimum["list_ce"].append(gradient_norm(list_ce, labels, labels))
        at_optimum["ranknet"].append(gradient_norm(ranknet, labels, labels))
        at_optimum["hinge"].append(gradient_norm(hinge, labels, labels))
        # ListCE is scale-free: any positive multiple of the centered label
        # must also be a stationary point. This is what lets MSE pin the scale.
        scaled = SHIFT + 1.7 * (labels - SHIFT)
        proportional.append(gradient_norm(list_ce, scaled, labels))

        # The decisive question for usefulness: how much gradient does each
        # term supply where the model actually sits?
        current = operating_point(labels, arguments.operating_r,
                                  arguments.operating_dispersion, generator)
        at_operating["list_ce"].append(gradient_norm(
            lambda s, y: list_ce(s, y), current, labels))
        at_operating["ranknet"].append(gradient_norm(ranknet, current, labels))
        at_operating["hinge"].append(gradient_norm(hinge, current, labels))
        at_operating["squared_error"].append(gradient_norm(
            lambda s, y: (s - y).square().mean(), current, labels))

    report = {
        "schema": "MetaSieve.R14AlignmentCheck.v1",
        "shift": SHIFT,
        "panels": arguments.panels,
        "panel_size": arguments.panel_size,
        "gradient_norm_at_the_regression_optimum": {
            name: {"mean": sum(v) / len(v), "max": max(v)}
            for name, v in at_optimum.items()},
        "list_ce_gradient_norm_at_a_proportional_point": {
            "scale": 1.7,
            "mean": sum(proportional) / len(proportional),
            "max": max(proportional)},
        "gradient_norm_at_the_measured_operating_point": {
            "r": arguments.operating_r,
            "dispersion": arguments.operating_dispersion,
            **{name: {"mean": sum(v) / len(v)}
               for name, v in at_operating.items()}},
    }
    operating = report["gradient_norm_at_the_measured_operating_point"]
    operating["listce_over_ranknet"] = (
        operating["list_ce"]["mean"] / operating["ranknet"]["mean"])
    operating["listce_over_squared_error"] = (
        operating["list_ce"]["mean"] / operating["squared_error"]["mean"])

    print(json.dumps(report, indent=1))
    aligned = report["gradient_norm_at_the_regression_optimum"]["list_ce"]["max"] < 1e-9
    print("\nALIGNED: the regression-compatible ListCE contributes exactly zero "
          "gradient at s = y" if aligned else
          "\nNOT ALIGNED: the design premise fails; do not implement")
    print(f"  ListCE  max |grad| at s=y : "
          f"{report['gradient_norm_at_the_regression_optimum']['list_ce']['max']:.3e}")
    print(f"  RankNet max |grad| at s=y : "
          f"{report['gradient_norm_at_the_regression_optimum']['ranknet']['max']:.3e}")
    print(f"  hinge   max |grad| at s=y : "
          f"{report['gradient_norm_at_the_regression_optimum']['hinge']['max']:.3e}")

    print(f"\nmean |grad| at the measured operating point "
          f"(r={arguments.operating_r}, sd_p/sd_y={arguments.operating_dispersion}):")
    for name in ("squared_error", "ranknet", "hinge", "list_ce"):
        print(f"  {name:<14} {operating[name]['mean']:.4e}")
    print(f"  ListCE / RankNet       : {operating['listce_over_ranknet']:.3f}")
    print(f"  ListCE / squared error : {operating['listce_over_squared_error']:.3f}")

    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps(report, indent=1), encoding="utf-8")
        print(f"\nwrote {arguments.output}")
    return 0 if aligned else 1


if __name__ == "__main__":
    raise SystemExit(main())
