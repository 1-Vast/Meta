"""Repository contract for the split-isolated development surface.

`QPSMPData(include_meta_test=False)` gives **logical exclusion after parsing**:
`cells.jsonl.gz` is one all-label artifact, so every sealed label is
decompressed and parsed on every construction and then discarded. That is
strictly weaker than isolation, and the 2026-08-16 governance incident showed
why it matters — six analysis scripts took the permissive default and two
artifacts remain permanently `process_unsealed`.

`scripts/build_governed_split_views.py` produces the stronger surface: separate
per-split label artifacts, with the `meta_test` artifact written **outside the
development tree** and only its hash recorded in the development manifest.

The load-bearing test here is `test_development_mount_works_with_the_sealed_
artifact_physically_absent`: it moves the sealed artifact away entirely and
proves that training and `meta_val` evaluation still mount. A guarantee that
only holds while the file happens to be present is not isolation.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path
import shutil

import pytest

from scripts.build_governed_split_views import (
    DEVELOPMENT_SPLITS, GovernedSplitView, SEALED_SPLIT, sha256_file,
)

ROOT = Path(__file__).resolve().parents[2]
VIEWS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views"
SEALED = Path("C:/Users/59964/MetaSieve_sealed/bindingdb_ki_double_cold_v1")

pytestmark = pytest.mark.skipif(
    not (VIEWS / "manifest.json").exists(),
    reason="governed split views not built; run "
           "scripts/build_governed_split_views.py")


@pytest.fixture(scope="module")
def manifest() -> dict:
    return json.loads((VIEWS / "manifest.json").read_text(encoding="utf-8"))


# --- the surface withholds the sealed artifact ----------------------------

def test_the_development_tree_contains_no_meta_test_label_artifact():
    assert not (VIEWS / SEALED_SPLIT).exists()
    for path in VIEWS.rglob("*.jsonl.gz"):
        assert SEALED_SPLIT not in path.parts, path


def test_the_manifest_asserts_the_sealed_artifact_was_withheld(manifest):
    assert manifest["meta_test_label_artifact_emitted"] is False
    assert "path" not in manifest.get("meta_test_label_artifact_location", "")


def test_no_meta_test_label_appears_anywhere_in_the_development_tree():
    """Read every development row and check the split field."""
    for split in DEVELOPMENT_SPLITS:
        with gzip.open(VIEWS / split / "cells.jsonl.gz", "rt",
                       encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    assert json.loads(line)["split"] == split


def test_governance_carries_no_labels():
    """`governance.jsonl` maps identity to split and must not leak a pK."""
    rows = [json.loads(line) for line
            in (VIEWS / "governance.jsonl").read_text(encoding="utf-8").splitlines()
            if line.strip()]
    assert rows
    assert all(set(row) == {"cell_id", "split", "target_id", "protein_group_40"}
               for row in rows)
    assert not any("pK" in row for row in rows)


# --- the assignment is unchanged ------------------------------------------

def test_the_split_assignment_is_bit_identical_to_the_frozen_one(manifest):
    """The builder must not have re-cut the split."""
    source = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
    frozen = json.loads((source / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["split_assignment_sha256"] == frozen["assignment_sha256"]


def test_counts_match_the_recorded_protocol(manifest):
    counts = manifest["counts"]
    assert counts["meta_train"]["rows"] == 5643
    assert counts["meta_train"]["targets"] == 346
    assert counts["meta_val"]["targets"] == 41
    assert counts["meta_val"]["components"] == 19
    assert counts[SEALED_SPLIT]["targets"] == 22
    assert counts[SEALED_SPLIT]["components"] == 10


def test_every_emitted_artifact_matches_its_recorded_hash(manifest):
    for split in DEVELOPMENT_SPLITS:
        path = VIEWS / split / "cells.jsonl.gz"
        assert sha256_file(path) == manifest["artifacts"][f"{split}/cells.jsonl.gz"]


# --- the load-bearing negative test ---------------------------------------

def test_development_mount_works_with_the_sealed_artifact_physically_absent(tmp_path):
    """Move the sealed artifact away; development must be unaffected.

    This is the property `include_meta_test=False` cannot have. The sealed
    labels are not merely skipped — they are not on the development process's
    filesystem path at all, and nothing in the mount depends on them.
    """
    relocated = None
    try:
        if SEALED.exists():
            relocated = tmp_path / "relocated"
            shutil.move(str(SEALED), str(relocated))
            assert not SEALED.exists()

        view = GovernedSplitView(VIEWS, visible=DEVELOPMENT_SPLITS)
        assert view.cells
        assert {cell["split"] for cell in view.cells} == set(DEVELOPMENT_SPLITS)
        record = view.isolation_record()
        assert record["physically_isolated"] is True
        assert record["level"] == "physically_isolated"
        assert SEALED_SPLIT not in record["visible_splits"]

        train_only = GovernedSplitView(VIEWS, visible=("meta_train",))
        assert {cell["split"] for cell in train_only.cells} == {"meta_train"}
    finally:
        if relocated is not None and relocated.exists():
            SEALED.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(relocated), str(SEALED))


# --- mounting the sealed split is gated -----------------------------------

def test_mounting_meta_test_requires_authorization_and_an_out_of_tree_path():
    with pytest.raises(ValueError, match="written authorization"):
        GovernedSplitView(VIEWS, visible=(SEALED_SPLIT,))
    with pytest.raises(ValueError, match="out-of-tree"):
        GovernedSplitView(VIEWS, visible=(SEALED_SPLIT,),
                          authorization="test")


def test_meta_test_cannot_be_mounted_beside_a_development_split():
    with pytest.raises(ValueError, match="alongside a development split"):
        GovernedSplitView(VIEWS, visible=("meta_train", SEALED_SPLIT),
                          authorization="test", sealed_directory=SEALED)


def test_authorization_without_requesting_meta_test_is_refused():
    with pytest.raises(ValueError, match="ambiguous mount"):
        GovernedSplitView(VIEWS, visible=DEVELOPMENT_SPLITS,
                          authorization="test")


def test_a_tampered_artifact_fails_the_mount(tmp_path):
    """Hash bindings must be enforced, not merely recorded."""
    copy = tmp_path / "views"
    shutil.copytree(VIEWS, copy)
    path = copy / "meta_val" / "cells.jsonl.gz"
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        rows = [json.loads(line) for line in handle if line.strip()]
    rows = rows[:-1]
    with gzip.GzipFile(path, "wb", mtime=0) as raw:
        for row in rows:
            raw.write((json.dumps(row, sort_keys=True) + "\n").encode("utf-8"))
    with pytest.raises(ValueError, match="hash mismatch"):
        GovernedSplitView(copy, visible=DEVELOPMENT_SPLITS)


def test_no_meta_test_metric_is_recorded_anywhere_in_the_surface():
    """The builder computed a row count and a hash, and nothing else."""
    text = (VIEWS / "manifest.json").read_text(encoding="utf-8")
    payload = json.loads(text)
    sealed_block = payload["counts"][SEALED_SPLIT]
    assert set(sealed_block) == {"rows", "targets", "components"}
    for banned in ("mse", "rmse", "spearman", "concordance", "pK", "score"):
        assert banned not in text.lower().replace("checkpoint_sha256", "")
