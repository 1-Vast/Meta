import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
sys.path.insert(0, str(STAGE))

from context_propagation import erase_at, load_pairs, summarize_delta  # noqa: E402


def test_erasure_makes_matched_inputs_identical():
    assert erase_at("ACDE", 2) == "AXDE"
    wt, mt = "ACDE", "AKDE"
    assert erase_at(wt, 2) == erase_at(mt, 2)


def test_delta_summary_uses_residue_coordinates():
    wt = np.zeros((20, 640), dtype=np.float32)
    mt = wt.copy()
    mt[3, 0] = 2.0
    summary, distances, norms = summarize_delta(wt, mt, 3)
    assert summary["site_delta_norm"] == 2.0
    assert distances[2] == 0 and norms[2] == 2.0


def test_covered_pairs_have_verified_sequences_and_coordinates():
    data2, selected, _ = load_pairs()
    assert len(selected) == len(data2["covered_pair_indices"]) == 49
    assert all(x["pair"]["pos"] == x["record"]["pos"] for x in selected)


def test_completed_erasure_gate_is_exact():
    import json
    result = json.loads((STAGE / "CONTEXT_PROPAGATION_RESULT.json").read_text())
    assert len(result["mutation_erasure"]) == 49
    assert all(row["inputs_identical"] and row["within_tolerance"]
               for row in result["mutation_erasure"])
    assert result["summary"]["erasure_max_absolute_embedding_delta"] == 0.0
