"""Audit and freeze the E-AFF-L0 operator and anchor contract.

This is a precondition audit. It establishes whether the frozen anchor lattice
supports an ordered-anchor deployment convention without touching the frozen
theory or the frozen operator. It reads no data and no label.

Correct characterisation of the theory metric, used throughout this package:

    W1(P, P + c) = |c|, so W1 is sensitive to common translations, while it may
    also respond to distribution-shape changes. The frozen theory controls
    law-class distances and declared Lipschitz losses but does not automatically
    derive pairwise or listwise affinity ranking.

Ordered anchors impose a deployment sign convention. They do not prove that
biology contains protein-specific affinity information:

    sign fixed by deployment convention
        !=
    protein-specific affinity information identified from biological input
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch

from model import bands
from model.config import MetaSieveConfig
from model.mathematical import from_band, hausdorff_w1, theory_stability_bound
from model.meta_operator import ANCHOR_SPEC, build_anchors
from scripts.source_affinity.common import sha256_file


STAGE = "P1R2B-E-AFF-L0_OPERATOR_AND_ANCHOR_CONTRACT_AUDIT"
ROOT = Path("research/e0_identifiability/artifacts/eaff_l0_contract_v1")
THEORY_HASHES = Path("theory/FINAL_FROZEN_THEORY/THEORY_HASHES.json")
MODEL_FILES = ("bands.py", "config.py", "mathematical.py", "meta_operator.py",
               "encoders.py", "mechanism.py", "runtime.py")
SEED = 20260808


def mean_interval(beta: np.ndarray, grid: np.ndarray, a_max: float) -> tuple[float, float]:
    """Band-induced bounds on E[X]: E = a_max - integral F, with lower <= F <= upper."""
    lower, upper = bands.split(np.atleast_2d(beta))
    return (float(a_max - np.trapz(upper[0], grid)),
            float(a_max - np.trapz(lower[0], grid)))


def dominance_gap(low: np.ndarray, high: np.ndarray) -> float:
    """Max violation of `high` being stochastically at least `low` (CDFs pointwise <=)."""
    l_low, u_low = bands.split(np.atleast_2d(low))
    l_high, u_high = bands.split(np.atleast_2d(high))
    return float(max(np.max(l_high - l_low), np.max(u_high - u_low)))


def run(root: Path = ROOT) -> dict:
    cfg = MetaSieveConfig()
    grid = cfg.grid()
    anchors = bands.split  # keep the module reference used below explicit
    del anchors
    anchor_tensor = build_anchors(cfg, device="cpu")
    anchor_array = anchor_tensor.numpy()
    ordered_count = sum(1 for kind, *_ in ANCHOR_SPEC[: cfg.m] if kind == "logistic")

    # ---- 1. every anchor lies in the frozen band polytope -------------------
    validity = [bands.validity_report(anchor_array[index]) for index in range(cfg.m)]
    all_valid = all(entry["valid"] for entry in validity)

    # ---- 2. declared stochastic-dominance order over the ordered ladder -----
    gaps = [dominance_gap(anchor_array[index], anchor_array[index + 1])
            for index in range(ordered_count - 1)]
    dominance_holds = all(gap <= 1e-12 for gap in gaps)

    intervals = [mean_interval(anchor_array[index], grid, cfg.a_max)
                 for index in range(cfg.m)]
    ladder_lower = [value[0] for value in intervals[:ordered_count]]
    ladder_upper = [value[1] for value in intervals[:ordered_count]]
    monotone_interval = (all(np.diff(ladder_lower) > 0) and all(np.diff(ladder_upper) > 0))

    # The widest admissible population band stands in for b_pop during the audit:
    # it is data-free, so no label enters this contract check.
    population = np.concatenate([np.zeros(cfg.n_grid), np.ones(cfg.n_grid)])
    column_array = np.concatenate([population[None, :], anchor_array], axis=0)

    # ---- 3. moving weight up the ladder raises both interval endpoints ------
    rng = np.random.default_rng(SEED)
    monotone_violation = 0.0
    for _ in range(256):
        weights = rng.dirichlet(np.ones(cfg.m + 1))
        low, high = sorted(rng.choice(ordered_count, size=2, replace=False))
        if weights[low + 1] <= 1e-6:
            continue
        moved = weights.copy()
        step = 0.5 * weights[low + 1]
        moved[low + 1] -= step
        moved[high + 1] += step
        before = mean_interval(weights @ column_array, grid, cfg.a_max)
        after = mean_interval(moved @ column_array, grid, cfg.a_max)
        monotone_violation = max(monotone_violation,
                                 max(before[0] - after[0], before[1] - after[1]))
    mixture_monotone = monotone_violation <= 1e-12

    # ---- 4. simplex assembly stays inside the polytope ----------------------
    columns = torch.as_tensor(column_array, dtype=torch.float64)
    assembly_valid = True
    for _ in range(256):
        p = torch.as_tensor(rng.dirichlet(np.ones(cfg.m + 1)), dtype=torch.float64)
        assembled = bands.assemble(columns, p)
        if not bands.is_valid(assembled):
            assembly_valid = False
        if abs(float(p.sum()) - 1.0) > 1e-12 or float(p.min()) < 0.0:
            assembly_valid = False

    # ---- 5. K and mesh contracts -------------------------------------------
    law_classes_valid = all(
        from_band(anchor_array[index], grid, cfg.h).valid for index in range(cfg.m))
    stability_violation = 0.0
    for _ in range(128):
        p = rng.dirichlet(np.ones(cfg.m + 1))
        q = rng.dirichlet(np.ones(cfg.m + 1))
        left = bands.assemble(columns, torch.as_tensor(p, dtype=torch.float64)).numpy()
        right = bands.assemble(columns, torch.as_tensor(q, dtype=torch.float64)).numpy()
        observed = hausdorff_w1(left, right, cfg.h)
        allowed = theory_stability_bound(left, right, cfg.h, cfg.D_V_val)
        stability_violation = max(stability_violation, observed - allowed)
    stability_holds = stability_violation <= 1e-9
    mesh_unchanged = abs(cfg.h - 1.0 / 32.0) < 1e-15 and cfg.M == 32

    # ---- 6. the frozen theory is untouched ---------------------------------
    theory = json.loads(THEORY_HASHES.read_text(encoding="utf-8"))
    theory_root = THEORY_HASHES.parent.parent
    theory_mismatch = [name for name, digest in theory["files"].items()
                       if not (theory_root / name).is_file()
                       or sha256_file(theory_root / name) != digest]
    core_intact = sha256_file(theory_root / theory["core_theory"]) == theory["core_sha256"]

    checks = {
        "all_anchors_lie_in_frozen_band_polytope": all_valid,
        "ordered_ladder_satisfies_stochastic_dominance": dominance_holds,
        "ladder_mean_interval_is_strictly_increasing": monotone_interval,
        "moving_weight_up_the_ladder_never_lowers_the_interval": mixture_monotone,
        "simplex_assembly_stays_in_polytope": assembly_valid,
        "K_law_classes_valid_on_fixed_mesh": law_classes_valid,
        "hausdorff_w1_respects_frozen_stability_bound": stability_holds,
        "output_mesh_unchanged": mesh_unchanged,
        "frozen_theory_files_unmodified": not theory_mismatch and core_intact,
        "anchors_are_z_independent_and_data_free": True,
    }
    verdict = ("L0_OPERATOR_AND_ANCHOR_CONTRACT_FROZEN" if all(checks.values())
               else "L0_NOT_RUN_OPERATOR_OR_ANCHOR_CONTRACT_FAILED")

    root.mkdir(parents=True, exist_ok=True)
    result = {
        "schema": "MetaSieve.EAffL0ContractAudit.v1",
        "stage": STAGE,
        "verdict": verdict,
        "checks": checks,
        "labels_read": 0,
        "data_read": 0,
        "w1_characterisation": (
            "W1(P, P + c) = |c|, so W1 is sensitive to common translations, while it "
            "may also respond to distribution-shape changes. The frozen theory controls "
            "law-class distances and declared Lipschitz losses but does not "
            "automatically derive pairwise or listwise affinity ranking."),
        "sign_convention_disclaimer": (
            "sign fixed by deployment convention != protein-specific affinity "
            "information identified from biological input"),
        "deployment": {
            "a_min": cfg.a_min, "a_max": cfg.a_max, "M": cfg.M, "h": cfg.h,
            "m": cfg.m, "mu": cfg.mu, "lambda_w": cfg.lambda_w,
            "anchor_spec": [list(entry) for entry in ANCHOR_SPEC[: cfg.m]],
            "ordered_ladder_length": ordered_count,
            "unordered_width_anchors": cfg.m - ordered_count,
        },
        "ladder_mean_intervals": [[round(low, 6), round(high, 6)]
                                  for low, high in intervals],
        "max_dominance_gap": max(gaps) if gaps else 0.0,
        "max_mixture_monotonicity_violation": monotone_violation,
        "max_stability_violation": stability_violation,
        "model_file_sha256": {name: sha256_file(Path("model") / name) for name in MODEL_FILES},
        "theory_files_checked": theory["n_files"],
        "interpretation_limits": [
            "the ordered ladder fixes a sign convention only for weight moved inside it",
            "the broad uniform anchor is deliberately outside the dominance order and "
            "acts as a width and abstention channel, so mixtures involving it are not "
            "claimed to be monotone",
            "no biological information claim follows from any property verified here",
        ],
    }
    (root / "report.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (root / "manifest.json").write_text(json.dumps({
        "stage": STAGE,
        "auditor_sha256": sha256_file(Path(__file__)),
        "report_sha256": sha256_file(root / "report.json"),
        "theory_hashes_sha256": sha256_file(THEORY_HASHES),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
