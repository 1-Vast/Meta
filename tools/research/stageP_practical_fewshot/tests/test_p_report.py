"""P1 comparison report structure tests (CPU, dry artifacts allowed)."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
sys.path.insert(0, str(STAGE))

import p_report as PR  # noqa: E402


def test_self_comparison_is_exactly_zero():
    recs = PR.load_arm_records(PR.OUT / "P1_ARM3_ORDINARYFT.json")
    rng = PR.stable_rng("stageP", "pcompare", PR.BOOT_SEED)
    cmp = PR.compare("arm3", recs, "arm3", recs, rng)
    for split_name in ("p_val", "p_test"):
        for k in (0, 5, 10):
            row = cmp.get(f"{split_name}:k{k}")
            if row is None:
                continue
            assert row["delta_mse_mean"] == 0.0
            assert row["delta_mse_ci95"] == [0.0, 0.0]
            assert row["n_pairs"] > 0


def test_baseline_pairing_counts_match_bank():
    bl = PR.load_baseline_records()
    bl_lig = {k[1:]: v for k, v in bl.items() if k[0] == "ligand_only"}
    recs = PR.load_arm_records(PR.OUT / "P1_ARM3_ORDINARYFT.json")
    bank = json.loads((PR.OUT / "P_BANK.json").read_text(encoding="utf-8"))
    n_seeds = len(set(k[0] for k in recs))
    for split_name in ("p_val", "p_test"):
        for k in (5, 10):
            expected = sum(1 for r in bank["records"]
                           if r["split"] == split_name and r["k"] == k)
            got = len(PR.paired_rows(recs, bl_lig, split_name, k,
                                     ref_is_baseline=True)) // n_seeds
            assert got == expected, (split_name, k, got, expected)
