"""Phase 2A — terminal verdict under the frozen precedence rules.

Registered by PREREG_S7_L2B_PHASE2A.md (sha 4e01401d...), section 10.

Reads only the five Phase 2A machine artifacts, applies the precedence in the
order registered, and emits exactly one verdict plus its mandated next action.
No metric is recomputed here and no threshold is re-decided.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

ROOT = Path(r"D:\MetaSieve")
OUT = ROOT / "report" / "s7_l2b_r0r"
RES = ROOT / "research" / "s7_l2b_r0r"
PREREG_SHA = "4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e"

ACTION = {
    "PHASE2A_CONTRACT_OR_ARTIFACT_FAIL_CLOSED":
        "repair only the missing evidence contract; run nothing else",
    "PHASE2A_DATA_NOT_IDENTIFIABLE":
        "preregister a metadata-only census for a new multi-ligand structural "
        "corpus; do not train",
    "TEACHER_GENERIC_POCKET_ONLY":
        "close the exact-coupling claim on the current MONN teacher; preregister a "
        "dense same-protein multi-ligand structural corpus; do not repair B5 on the "
        "same labels",
    "LABEL_SEMANTICS_AMBIGUOUS":
        "preregister one dense continuous-coordinate or audited soft-teacher "
        "reconstruction; do not change the learner yet",
    "LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING":
        "preregister one ligand-conditioned residue residual head",
    "EDGE_COUPLING_PRESENT_B5_ABSENT":
        "preregister the T0-T3 input-observability ladder first; do not train a pair head",
    "EDGE_COUPLING_ALREADY_IDENTIFIED":
        "do not repair; preregister one sealed independent structural confirmation",
}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    files = {
        "contract": OUT / "PHASE2A_INPUT_MANIFEST.json",
        "census": OUT / "PHASE2A_DATA_IDENTIFIABILITY_CENSUS.json",
        "teacher": OUT / "PHASE2A_TEACHER_CONDITIONALITY.json",
        "coupling": OUT / "PHASE2A_MARGINAL_COUPLING_AUDIT.json",
        "labels": OUT / "PHASE2A_LABEL_SEMANTICS.json",
    }
    missing = [k for k, p in files.items() if not p.is_file()]
    if missing:
        print(f"missing artifacts: {missing}")
        return 1
    d = {k: json.loads(p.read_text(encoding="utf-8")) for k, p in files.items()}

    trail = []
    verdict = None

    # rule 1
    r1 = d["contract"]["verdict"] != "PHASE2A_CONTRACT_PASS"
    trail.append({"rule": 1, "test": "any Phase 0 check fails", "fired": r1})
    if r1:
        verdict = "PHASE2A_CONTRACT_OR_ARTIFACT_FAIL_CLOSED"

    # rule 2
    suff = d["census"]["sufficiency"]
    r2 = suff["verdict"] != "DATA_IDENTIFIABLE"
    trail.append({"rule": 2, "test": "D1, D2 or D3 fails", "fired": r2,
                  "D1": suff["D1_components_ge_30"], "D2": suff["D2_scaffold_distinct_pairs_ge_100"],
                  "D3": suff["D3_all_records_masked"]})
    if verdict is None and r2:
        verdict = "PHASE2A_DATA_NOT_IDENTIFIABLE"

    # rule 3
    tv = d["teacher"]["teacher_verdict"]
    r3 = not (tv["T1_pass"] and tv["T6_pass"])
    trail.append({"rule": 3, "test": "T1 fails OR T6 fails", "fired": r3,
                  "T1_pass": tv["T1_pass"], "T6_pass": tv["T6_pass"]})
    if verdict is None and r3:
        verdict = "TEACHER_GENERIC_POCKET_ONLY"

    # rule 4
    amb = d["labels"]["ambiguity_verdict"]
    r4 = bool(amb["LABEL_SEMANTICS_AMBIGUOUS"])
    trail.append({"rule": 4, "test": "label ambiguity positively demonstrated",
                  "fired": r4, "detail": amb})
    if verdict is None and r4:
        verdict = "LABEL_SEMANTICS_AMBIGUOUS"

    BC = bool(d["coupling"]["BC_b5_edge_coupling"]["BC"])
    TC = bool(d["coupling"]["TC_teacher_edge_coupling"]["TC"])
    trail.append({"rule": "5/6/7", "BC": BC, "TC": TC})
    if verdict is None:
        if BC:
            verdict = "EDGE_COUPLING_ALREADY_IDENTIFIED"
        elif TC:
            verdict = "EDGE_COUPLING_PRESENT_B5_ABSENT"
        else:
            verdict = "LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING"

    dep = d["coupling"]["DEPLOYABLE_prediction_decomposition"]["arms"]
    orc = d["coupling"]["ORACLE_label_fitted_ceilings"]
    t1 = d["teacher"]["T1_jaccard"]

    out = {
        "schema": "MetaSieve.S7L2B.P2A.TerminalVerdict.v1",
        "created_utc": "2026-08-10",
        "preregistration": {
            "file": "research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2A.md",
            "sha256": PREREG_SHA,
            "committed": False,
            "chronology_note": "registered by hash, not by a git commit; commit "
                               "authorization was not granted for this run",
            "amendments": {
                f: sha(RES / f) for f in
                ("PREREG_S7_L2B_PHASE2A_AMENDMENT_01.md",
                 "PREREG_S7_L2B_PHASE2A_AMENDMENT_02.md",
                 "PREREG_S7_L2B_PHASE2A_AMENDMENT_03.md")
                if (RES / f).is_file()},
        },
        "repo_commit": "623602e76b7d4f445af069014782278163183d59",
        "input_artifact_sha256": {k: sha(p) for k, p in files.items()},
        "code_sha256": {f.name: sha(f) for f in sorted(RES.glob("pa*.py"))},
        "software": {"python": platform.python_version(),
                     "platform": platform.platform()},
        "seeds": {"census_none": None, "teacher_bootstrap": 20260819,
                  "teacher_permutation": 20260820, "rewiring": 20260817,
                  "component_bootstrap": 20260818, "within_complex_shuffle": 20260821},
        "inference_unit": "protein closure component; atom-residue rows are never "
                          "inference units",
        "label_field_reads": d["contract"]["checks"]["C7_label_and_affinity_reads"],
        "affinity_value_reads": 0,
        "numerical_tolerances_achieved": d["coupling"]["numerics"],
        "precedence_trail": trail,
        "TERMINAL_VERDICT": verdict,
        "authorized_next_action": ACTION[verdict],
        "headline_evidence": {
            "teacher_replicate_jaccard": t1["replicate_mean"],
            "teacher_alternative_ligand_jaccard": t1["alternative_ligand_mean"],
            "teacher_dJ": t1["dJ_paired"]["delta"],
            "teacher_dJ_lcb95": t1["dJ_paired"]["lcb95_one_sided"],
            "teacher_chemistry_association_rho": d["teacher"]["T5_scaffold_distance_sensitivity"]["component_macro_rho"],
            "B5_full_ap": dep["B5"]["full"],
            "B5_additive_ap": dep["B5"]["add"],
            "B5_coupling_ap": dep["B5"]["coup"],
            "B5_coupling_rewiring_null_ap": d["coupling"]["NULL_evaluation_only"]["degree_preserving_rewiring_coupling_ap"]["B5"],
            "label_fitted_additive_ceiling": orc["least_squares_additive_projection_of_Y"],
            "true_residue_margin_ceiling": orc["true_residue_margin_only"],
            "fraction_of_additive_ceiling_reached_by_B5": (
                dep["B5"]["full"] / orc["least_squares_additive_projection_of_Y"]),
        },
        "frozen_boundaries_still_in_force": [
            "real ChEMBL/BindingDB affinity training",
            "DAVIS, KIBA and recipient labels",
            "independent confirmation scoring",
            "new PLM, attention stack, geometry branch, typed-interaction branch, "
            "affinity head, PU loss, knowledge graph or parallel module",
            "few-shot section adaptation and any k-shot claim",
            "admission of any biological statistic into production z",
            "CSMO, Band, mesh, positive ridge, and A(F,z) = K(B(z)F(z))",
            "P2-P4",
        ],
        "claims_not_made": [
            "EXACT_RESIDUE_ATOM_COUPLING_IDENTIFIED",
            "AFFINITY_ENERGETICS_IDENTIFIED",
            "K_SHOT_SECTION_IDENTIFIED",
            "BIOLOGICAL_STATISTIC_ADMITTED_TO_Z",
            "VALIDATED_END_TO_END_DTA_MODEL",
        ],
    }
    (OUT / "PHASE2A_VERDICT.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"TERMINAL_VERDICT": verdict,
                      "authorized_next_action": out["authorized_next_action"],
                      "precedence_trail": trail,
                      "headline_evidence": out["headline_evidence"]}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
