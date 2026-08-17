"""Frozen design for the method-ladder cycle. Written before any result.

Nothing in this module may be edited in response to an outcome. Thresholds,
populations, panel rules, control definitions and probe budgets are fixed here
so that "we froze the design first" is a checkable property of the repository.

The ladder's shared question, asked once and reused by every family that needs
it, is deliberately *not* "is this representation sensitive to the protein".
Stage P already answered that: the incumbent's protein response can be made
reproducible across seeds (+0.316 cosine) while carrying no affinity
information (+0.022 alignment). Sensitivity is cheap. The question here is
whether a representation carries within-target ordering information that the
protein-blind ligand encoder does not.

Selection happens on `meta_train` component folds. `meta_val` is read once per
frozen probe. `meta_test` is never constructible from here.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]

SPLIT_DIRECTORY = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
VIEWS_DIRECTORY = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views"

# --- frozen reference model ------------------------------------------------
# The three retained A0 seeds. Every frozen-representation probe reports the
# three seeds separately and takes their median; a checkpoint average is a
# different model and is never called "the incumbent".
A0_CHECKPOINTS = tuple(
    ROOT / "report/meta_fewshot/stageR3R4_level_shape_20260815"
         / f"A0_incumbent_seed{seed}" / "checkpoint.pt"
    for seed in (20260815, 20260816, 20260817))

# --- panels ----------------------------------------------------------------
# One panel per target: unique ligands in a deterministic stable-seed order.
# A within-target Pearson r on fewer than six ligands is dominated by its own
# sampling noise, so smaller targets are excluded from the ordering statistic
# rather than being allowed to widen every interval.
PANEL_SEED = 20260816
PANEL_MAX = 32
PANEL_MIN = 6

# --- probe -----------------------------------------------------------------
# A tiny head trained by ordinary SGD. No ridge, no pseudoinverse, no closed
# form, no matrix solve — the prohibition applies to diagnostics too, not only
# to candidate models.
PROBE_STEPS = 400
PROBE_LR = 3e-3
PROBE_BATCH_TARGETS = 32
PROBE_WEIGHT_DECAYS = (1e-4, 1e-3, 1e-2, 1e-1)   # the only selected knob
PROBE_SEEDS = (20260816, 20260817, 20260818)
PROBE_FOLDS = 5
FOLD_SEED = 20260818

# --- controls --------------------------------------------------------------
CONTROL_SEED = 20260817
RANDOM_FEATURE_SEED = 20260819

# --- statistics ------------------------------------------------------------
BOOTSTRAP_DRAWS = 9999
BOOTSTRAP_SEED = 20260816

# --- thresholds, frozen before observation ---------------------------------
# `r` here is the within-target Pearson correlation between a probe's output
# and the centered label. Centering removes the per-target level exactly, so
# every number below speaks to ordering and nothing else.
SMALLEST_EFFECT_OF_INTEREST_R = 0.05
DECISIVE_NULL_HALF_WIDTH_R = 0.05

# The shared information gate. A representation passes only if all three hold.
INFORMATION_GATE = {
    "beats_ligand_only": "lo > 0 and mean >= SMALLEST_EFFECT_OF_INTEREST_R",
    "beats_matched_wrong_protein": "lo > 0",
    "destroyed_by_label_permutation": ">= 80% of the effect removed",
}
LABEL_PERMUTATION_DESTRUCTION = 0.80

# Verdict vocabulary. Every family terminates in exactly one of these.
VERDICTS = (
    "REJECTED_BY_INPUT_CONTRACT",
    "REJECTED_BY_STRUCTURAL_GATE",
    "REJECTED_BY_FROZEN_DISCRIMINATOR",
    "REJECTED_BY_TRAINING_SCREEN",
    "ADMITTED_TO_FULL_EVALUATION",
    "ADMITTED_COMPONENT",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def frozen_manifest() -> dict:
    """The design fingerprint every ladder artifact must carry."""
    manifest = json.loads(
        (SPLIT_DIRECTORY / "manifest.json").read_text(encoding="utf-8"))
    return {
        "schema": "MetaSieve.MethodLadder.FrozenDesign.v1",
        "split_directory": str(SPLIT_DIRECTORY.relative_to(ROOT)),
        "split_assignment_sha256": manifest["assignment_sha256"],
        "a0_checkpoints": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in A0_CHECKPOINTS if path.exists()},
        "panels": {"seed": PANEL_SEED, "max": PANEL_MAX, "min": PANEL_MIN},
        "probe": {"steps": PROBE_STEPS, "lr": PROBE_LR,
                  "batch_targets": PROBE_BATCH_TARGETS,
                  "weight_decays": list(PROBE_WEIGHT_DECAYS),
                  "seeds": list(PROBE_SEEDS), "folds": PROBE_FOLDS,
                  "fold_seed": FOLD_SEED,
                  "solver": "AdamW; no ridge/pseudoinverse/closed form"},
        "controls": {"seed": CONTROL_SEED,
                     "random_feature_seed": RANDOM_FEATURE_SEED},
        "statistics": {"bootstrap_draws": BOOTSTRAP_DRAWS,
                       "bootstrap_seed": BOOTSTRAP_SEED},
        "thresholds": {
            "smallest_effect_of_interest_r": SMALLEST_EFFECT_OF_INTEREST_R,
            "decisive_null_half_width_r": DECISIVE_NULL_HALF_WIDTH_R,
            "label_permutation_destruction": LABEL_PERMUTATION_DESTRUCTION},
        "information_gate": INFORMATION_GATE,
    }


def verdict(interval: dict) -> str:
    """Four states. "No effect" is never one of them.

    ``RESOLVED`` the interval excludes zero and the effect is at least the
    smallest one this cycle would act on. ``RESOLVED_NEGLIGIBLE`` excludes zero
    but is too small to exploit. ``DECISIVE_NULL`` contains zero and is narrow
    enough that a useful effect would have been seen. ``UNDERPOWERED``
    contains zero and is too wide to conclude anything.
    """
    lo, hi = interval["lo"], interval["hi"]
    if lo > 0 or hi < 0:
        return ("RESOLVED"
                if abs(interval["mean"]) >= SMALLEST_EFFECT_OF_INTEREST_R
                else "RESOLVED_NEGLIGIBLE")
    if (hi - lo) / 2.0 < DECISIVE_NULL_HALF_WIDTH_R:
        return "DECISIVE_NULL"
    return "UNDERPOWERED"
