"""Formal adjudication for frozen CIIP-1A control-arm protocol."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PREREG_SHA = "39d02166f69acf235a34d351b649a4cdbf3b828491a0994901bf2378777463f7"
REQUIRED_ARMS = {
    "oracle_local_esm_correct", "family_preserving_shuffle",
    "random_local_window", "ligand_only", "ligand_invariant_shift",
    "random_protein", "free_pairwise",
}


def load_and_adjudicate(result_path: Path) -> dict:
    prereg_path = HERE / "PREREGISTRATION_STAGE1_CONTROLS.md"
    prereg_sha = hashlib.sha256(prereg_path.read_bytes()).hexdigest()
    if prereg_sha != PREREG_SHA:
        raise ValueError(f"frozen preregistration SHA mismatch: {prereg_sha}")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if result.get("preregistration_controls_sha256") != PREREG_SHA:
        raise ValueError("result does not cite the frozen preregistration")
    if set(result.get("arms", {})) != REQUIRED_ARMS:
        raise ValueError("result arms do not match the frozen control matrix")

    arms = result["arms"]
    effects = result["effects"]
    correct = arms["oracle_local_esm_correct"]["agg"]
    random_window = arms["random_local_window"]["agg"]
    family = effects["v_family"]
    random = effects["v_random_window"]
    random_protein = effects["v_random_protein"]
    ligand = effects["v_ligand_only"]

    coverage_exceeds_random = correct["n_nonconstant"] > random_window["n_nonconstant"]
    floor_r2_supported = ligand["bootstrap_ci"]["lo2.5"] > 0
    floor_sign_supported = (correct["sign_acc"] >= 0.55 and
                            arms["ligand_only"]["agg"]["sign_acc"] < correct["sign_acc"])
    family_supported = family["bootstrap_ci"]["lo2.5"] > 0
    random_protein_supported = random_protein["bootstrap_ci"]["lo2.5"] > 0
    random_supported = (random["bootstrap_ci"]["lo2.5"] > 0 and
                        correct["sign_acc"] - random_window["sign_acc"] >= 0.05)
    lopo_supported = random["leave_one_parent_out_sign_stable"]

    supported = all((coverage_exceeds_random, floor_r2_supported,
                     floor_sign_supported, family_supported,
                     random_protein_supported, random_supported, lopo_supported))
    family_absent = (abs(family["observed_pair_mean_effect"]) < 0.02 or
                     family["bootstrap_ci"]["lo2.5"] <= 0)
    random_absent = (abs(random["observed_pair_mean_effect"]) < 0.02 or
                     random["bootstrap_ci"]["lo2.5"] <= 0)
    not_supported = family_absent and random_absent and not coverage_exceeds_random
    verdict = ("ORACLE_LOCAL_SIGNAL_SUPPORTED" if supported else
               "ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED" if not_supported else
               "ORACLE_LOCAL_SIGNAL_UNRESOLVED")
    return {
        "schema": "MetaSieve.StageCIIP1A.Controls.Adjudication.v1",
        "preregistration_controls_sha256": PREREG_SHA,
        "control_result_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
        "scope": result["scope"],
        "verdict": verdict,
        "checks": {
            "correct_nonconstant_coverage_exceeds_random": coverage_exceeds_random,
            "correct_beats_ligand_floor_r2_bootstrap": floor_r2_supported,
            "correct_beats_ligand_floor_sign": floor_sign_supported,
            "correct_beats_family_r2_bootstrap": family_supported,
            "correct_beats_random_protein_r2_bootstrap": random_protein_supported,
            "correct_beats_random_window_r2_and_sign": random_supported,
            "correct_vs_random_window_lopo_sign_stable": lopo_supported,
            "not_supported_family_gap_absent": family_absent,
            "not_supported_random_window_gap_absent": random_absent,
        },
        "evidence": {
            "correct": correct,
            "random_window": random_window,
            "v_family": family,
            "v_random_window": random,
            "v_random_protein": random_protein,
            "v_ligand_only": ligand,
        },
        "authorization": {
            "deployable_protein_representation": "NOT_VALIDATED",
            "ciip_1a_pass": "NOT_AUTHORIZED",
            "ciip_1b": "NOT_AUTHORIZED",
            "bindingdb_potential_bridge": "NOT_AUTHORIZED",
            "production_integration": "NOT_AUTHORIZED",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, default=HERE / "CONTROL_RESULT.json")
    parser.add_argument("--out", type=Path, default=HERE / "CONTROL_ADJUDICATION.json")
    args = parser.parse_args()
    output = load_and_adjudicate(args.result)
    args.out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
