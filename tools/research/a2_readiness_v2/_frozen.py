"""Everything Phase 1 freezes before it observes a single result.

Nothing in this module may be edited in response to an outcome. It exists so
that "we froze the design first" is a checkable property of the repository
rather than a claim in a report: `frozen_manifest()` hashes the checkpoints and
the split assignment, and every probe writes that manifest into its output.

Design and threshold choices belong to `meta_train` component folds
(`MODEL_FOLD_SEED`). `meta_val` is read once, after freezing. `meta_test` is
never constructible from here — the probes use the fail-closed default.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

SPLIT_DIRECTORY = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"

# --- arms ----------------------------------------------------------------
# The three retained A0 seeds, frozen by path and verified by sha256.
A0_CHECKPOINTS = tuple(
    ROOT / "report/meta_fewshot/stageR3R4_level_shape_20260815"
         / f"A0_incumbent_seed{seed}" / "checkpoint.pt"
    for seed in (20260815, 20260816, 20260817))

# Independently initialised comparison models. Ten is the number that makes the
# across-seed reproducibility statistic (Section 3 of the phase brief) estimable
# at all: with three you cannot separate "reproducible across initialisation"
# from "two draws happened to agree".
RANDOM_INIT_SEEDS = tuple(range(20260901, 20260911))

# --- evaluation bank ------------------------------------------------------
EVALUATION_SEED = 73101      # the frozen R6-R14 bank; unchanged
QUERY_SIZE = 16
SUPPORT_SIZES = (0, 1, 2, 3, 5)

# --- donor rule -----------------------------------------------------------
# Donors are drawn from `meta_val` (both proteins equally unseen) and whitened
# with `meta_train` statistics only. Strata are quantiles of the whitened cosine
# similarity to the cross-component candidate pool, so perturbation magnitude is
# a controlled variable rather than a confound.
DONOR_POOL = "meta_val"
WHITENING_POOL = "meta_train"
DONOR_STRATA = ("nearest", "q25", "median", "q75", "farthest")

# --- controls -------------------------------------------------------------
CONTROL_SEED = 20260817      # permutations, random ligand panels, shuffles
REPEATED_FORWARD_PASSES = 2  # numerical floor of the measurement itself

# --- statistics -----------------------------------------------------------
BOOTSTRAP_DRAWS = 9999
BOOTSTRAP_SEED = 20260816
MODEL_FOLD_SEED = 20260818   # meta_train component folds for any fitted probe

# --- thresholds, frozen before observation --------------------------------
# A contrast is RESOLVED only if the component-paired 95% interval excludes
# zero. A null is DECISIVE only if the interval is also narrower than the
# smallest effect the phase would act on; otherwise it is UNDERPOWERED, which
# is a different verdict and must be reported as such.
RESOLVED_RULE = "component-paired 95% bootstrap interval excludes zero"
SMALLEST_EFFECT_OF_INTEREST_R = 0.05   # in within-target Pearson r
DECISIVE_NULL_HALF_WIDTH_R = 0.05      # |hi-lo|/2 must be below this


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def split_assignment_sha256() -> str:
    manifest = json.loads(
        (SPLIT_DIRECTORY / "manifest.json").read_text(encoding="utf-8"))
    return manifest["assignment_sha256"]


def frozen_manifest() -> dict:
    """The pre-registration fingerprint every probe output must carry."""
    return {
        "schema": "MetaSieve.A2ReadinessV2.FrozenDesign.v1",
        "split_directory": str(SPLIT_DIRECTORY.relative_to(ROOT)),
        "split_assignment_sha256": split_assignment_sha256(),
        "a0_checkpoints": {
            str(path.relative_to(ROOT)): sha256(path) for path in A0_CHECKPOINTS},
        "random_init_seeds": list(RANDOM_INIT_SEEDS),
        "bank": {"evaluation_seed": EVALUATION_SEED, "query_size": QUERY_SIZE,
                 "support_sizes": list(SUPPORT_SIZES)},
        "donors": {"pool": DONOR_POOL, "whitening_pool": WHITENING_POOL,
                   "strata": list(DONOR_STRATA)},
        "controls": {"seed": CONTROL_SEED,
                     "repeated_forward_passes": REPEATED_FORWARD_PASSES},
        "statistics": {"bootstrap_draws": BOOTSTRAP_DRAWS,
                       "bootstrap_seed": BOOTSTRAP_SEED,
                       "model_fold_seed": MODEL_FOLD_SEED},
        "thresholds": {
            "resolved_rule": RESOLVED_RULE,
            "smallest_effect_of_interest_r": SMALLEST_EFFECT_OF_INTEREST_R,
            "decisive_null_half_width_r": DECISIVE_NULL_HALF_WIDTH_R},
    }


def verdict(interval: dict) -> str:
    """Four states. 'No effect' is never one of them.

    ``RESOLVED``
        the interval excludes zero **and** the effect is at least as large as
        the smallest one this phase would act on.
    ``RESOLVED_NEGLIGIBLE``
        the interval excludes zero but the effect is smaller than
        `SMALLEST_EFFECT_OF_INTEREST_R`. Statistically real, scientifically
        inert — a distinction that a bare significance test hides, and one that
        matters here because 41 targets x 16 queries can resolve a 0.002
        correlation difference that no mechanism could exploit.
    ``DECISIVE_NULL``
        the interval contains zero and is narrow enough that an effect worth
        acting on would have been detected.
    ``UNDERPOWERED``
        the interval contains zero and is too wide to conclude anything.
    """
    lo, hi = interval["lo"], interval["hi"]
    if lo > 0 or hi < 0:
        return ("RESOLVED" if abs(interval["mean"]) >= SMALLEST_EFFECT_OF_INTEREST_R
                else "RESOLVED_NEGLIGIBLE")
    if (hi - lo) / 2.0 < DECISIVE_NULL_HALF_WIDTH_R:
        return "DECISIVE_NULL"
    return "UNDERPOWERED"
