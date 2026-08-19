"""Stage P split/bank invariant tests (P1 prereg + AD1)."""
import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
sys.path.insert(0, str(STAGE))

ART = STAGE / "artifacts"


@pytest.fixture(scope="module")
def split_art():
    return json.loads((ART / "P_SPLIT.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def bank_art():
    return json.loads((ART / "P_BANK.json").read_text(encoding="utf-8"))


def test_split_clusters_whole_and_disjoint(split_art):
    cells = split_art["cell_split"]
    by_cluster = {}
    for cell, rec in cells.items():
        by_cluster.setdefault(rec["cluster"], set()).add(rec["split"])
    for cl, splits in by_cluster.items():
        assert len(splits) == 1, f"cluster {cl} split across {splits}"
    target_split = {}
    for rec in cells.values():
        target_split.setdefault(rec["target_id"], set()).add(rec["split"])
    for t, s in target_split.items():
        assert len(s) == 1, f"target {t} across {s}"


def test_split_target_fractions(split_art):
    census = split_art["census"]
    n = sum(census[s]["targets"] for s in census)
    assert n > 0
    for s, frac in zip(("p_train", "p_val", "p_test"), (0.60, 0.20, 0.20)):
        assert abs(census[s]["targets"] / n - frac) < 0.05


def test_bank_per_k_eligibility_and_disjoint(bank_art):
    Q = bank_art["Q"]
    for rec in bank_art["records"]:
        assert len(rec["support_cell_ids"]) == rec["k"]
        assert len(rec["query_cell_ids"]) <= Q
        sup = set(rec["support_cell_ids"])
        q = set(rec["query_cell_ids"])
        assert not (sup & q)
        assert len(sup | q) >= 1
        assert rec["n_ligands"] >= rec["k"] + len(rec["query_cell_ids"])


def test_bank_k0_has_no_support(bank_art):
    k0 = [r for r in bank_art["records"] if r["k"] == 0]
    assert k0, "k=0 records must exist"
    for rec in k0:
        assert rec["support_cell_ids"] == []
        assert len(rec["query_cell_ids"]) == bank_art["Q"]


def test_bank_donor_foreign(bank_art, split_art):
    cells = split_art["cell_split"]
    for rec in bank_art["records"]:
        assert rec["donor_target_id"] != rec["target_id"]
        assert rec["donor_cluster"] != rec["cluster"]
        assert cells[rec["query_cell_ids"][0]]["split"] == rec["split"]


def test_bank_manifest_sha():
    manifest = json.loads((ART / "P_BANK.manifest.json").read_text(encoding="utf-8"))
    import hashlib
    sha = hashlib.sha256((ART / "P_BANK.json").read_bytes()).hexdigest()
    assert sha == manifest["sha256"]
