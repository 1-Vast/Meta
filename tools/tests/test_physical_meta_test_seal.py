"""The physical meta_test seal: proved by what the process never opens.

`tools/tests/test_governed_split_views.py` pins the *artifact* contract — the
development surface withholds the `meta_test` label file and the manifest says
so. This file pins the *process* contract, which is the property the 2026-08-16
governance incident actually lacked:

1. **file-access spy** — every path a normal train / internal-validation /
   meta_val construction opens is recorded, and none of them is the sealed
   artifact or the all-label `cells.jsonl.gz`. The same spy is run against the
   default corpus surface, where it *does* see `cells.jsonl.gz`, so the
   instrument is not vacuous;
2. **fail-closed authorization** — missing, blank, non-textual, too-short and
   in-tree mounts of `meta_test` all raise;
3. **equivalence to the governed source** — mounting the view yields the same
   cells, task indices and component maps as the corpus construction, and the
   view's own hash and count bindings are enforced rather than merely recorded.

No test here evaluates `meta_test` or computes any statistic of its labels. The
only `meta_test` quantity read anywhere is the row count the builder recorded in
the development manifest.
"""
from __future__ import annotations

import builtins
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import shutil

import pytest

from scripts.build_governed_split_views import (
    AUTHORIZATION_MIN_CHARS, DEVELOPMENT_SPLITS, GovernedSplitView,
    SEALED_SPLIT, check_authorization,
)
from scripts.qpsmp_data import QPSMPData

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
PROTEIN_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank"
LIGAND_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank"
COMPACT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank_compact"
SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
VIEWS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views"

AUTHORIZED = ("contract test: seal-behaviour assertions only, "
              "no meta_test label read and no metric computed")

pytestmark = pytest.mark.skipif(
    not (VIEWS / "manifest.json").exists(),
    reason="governed split views not built; run "
           "scripts/build_governed_split_views.py")


# --- the spy --------------------------------------------------------------

class _AccessSpy:
    """Record every filesystem path the process opens while it is installed.

    `builtins.open`, `io.open`, `gzip.open` and `os.open` are patched
    separately: `pathlib` calls `io.open` by attribute, `numpy` and `torch`
    reach the descriptor layer, and `gzip` has its own entry point. Patching
    only one of them would leave a hole exactly where a label file would be
    read.
    """

    def __init__(self) -> None:
        self.paths: list[str] = []
        self._patched: list[tuple[object, str, object]] = []

    def _record(self, target) -> None:
        try:
            self.paths.append(str(Path(os.fspath(target))))
        except TypeError:  # a descriptor or a file object, not a path
            pass

    def __enter__(self) -> "_AccessSpy":
        for module, name in ((builtins, "open"), (io, "open"),
                             (gzip, "open"), (os, "open")):
            original = getattr(module, name)
            self._patched.append((module, name, original))

            def wrapper(file, *args, _original=original, **kwargs):
                self._record(file)
                return _original(file, *args, **kwargs)

            setattr(module, name, wrapper)
        return self

    def __exit__(self, *exception) -> None:
        for module, name, original in reversed(self._patched):
            setattr(module, name, original)
        self._patched.clear()

    def touched(self, needle: str) -> list[str]:
        return [path for path in self.paths if needle in path]


def _development_workflow(data: QPSMPData) -> None:
    """The ordinary train / internal-validation / meta_val read path."""
    import numpy as np

    rng = np.random.default_rng(20260818)
    for _ in range(5):
        spec = data.draw_episode("meta_train", 3, 8, rng)
        data.materialize(spec)
    banks = data.fixed_nested_episode_banks(
        "meta_val", (0, 1, 2, 3, 5), 8, 1, 73101, 1)
    for size in sorted(banks):
        data.materialize(banks[size][0])


@pytest.fixture(scope="module")
def isolated() -> QPSMPData:
    return QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                     split_directory=SPLIT, split_view=VIEWS)


@pytest.fixture(scope="module")
def corpus_surface() -> QPSMPData:
    return QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                     split_directory=SPLIT)


# --- 1. the file-access spy ----------------------------------------------

def test_the_spy_sees_the_all_label_corpus_on_the_default_surface():
    """Negative control: an instrument that never fires proves nothing."""
    with _AccessSpy() as spy:
        QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                  split_directory=SPLIT)
    assert spy.touched(str(CORPUS / "cells.jsonl.gz")), (
        "the spy did not observe the all-label corpus read that the default "
        "surface certainly performs; the instrument is broken")


def test_a_development_construction_never_touches_the_sealed_path():
    with _AccessSpy() as spy:
        data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                         split_directory=SPLIT, split_view=VIEWS)
        _development_workflow(data)

    assert spy.paths, "the spy recorded no access at all"
    assert spy.touched(SEALED_SPLIT) == [], (
        "a train/meta_val workflow opened a meta_test path: "
        f"{spy.touched(SEALED_SPLIT)}")
    all_label = str(CORPUS / "cells.jsonl.gz")
    assert spy.touched(all_label) == [], (
        "the isolated surface opened the all-label corpus: "
        f"{spy.touched(all_label)}")
    sealed_root = json.loads(
        (VIEWS / "manifest.json").read_text(encoding="utf-8"))
    assert sealed_root["meta_test_label_artifact_emitted"] is False


def test_the_development_read_path_is_exactly_the_two_visible_artifacts():
    """Whitelist, not filter: name every label file the process opened."""
    with _AccessSpy() as spy:
        data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                         split_directory=SPLIT, split_view=VIEWS)
        _development_workflow(data)
    label_files = {Path(path).parent.name for path in spy.paths
                   if path.endswith("cells.jsonl.gz")}
    assert label_files == set(DEVELOPMENT_SPLITS), label_files


def test_the_sealed_artifact_can_be_absent_entirely(tmp_path):
    """Isolation must not depend on the sealed file happening to exist."""
    record = json.loads(
        (VIEWS / "manifest.json").read_text(encoding="utf-8"))
    sealed_hash = record["artifacts"][f"{SEALED_SPLIT}/cells.jsonl.gz"]
    assert sealed_hash, "the manifest records no sealed artifact hash"

    surface = tmp_path / "views"
    shutil.copytree(VIEWS, surface)
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                     split_directory=SPLIT, split_view=surface)
    assert "meta_test" not in data.tasks
    assert all(cell["split"] != SEALED_SPLIT for cell in data.cells)


# --- 2. fail-closed authorization ----------------------------------------

def test_a_missing_authorization_is_refused():
    with pytest.raises(ValueError, match="none given"):
        check_authorization(None)
    with pytest.raises(ValueError, match="written authorization"):
        GovernedSplitView(VIEWS, visible=(SEALED_SPLIT,))


@pytest.mark.parametrize("blank", ["", "   ", "\t\n"])
def test_a_blank_authorization_is_refused(blank):
    with pytest.raises(ValueError, match="blank"):
        check_authorization(blank)


@pytest.mark.parametrize("malformed", [123, True, ["a reason"], {"why": "x"},
                                       Path("reason.txt")])
def test_a_non_textual_authorization_is_refused(malformed):
    with pytest.raises(ValueError, match="malformed"):
        check_authorization(malformed)


def test_a_placeholder_authorization_is_refused():
    short = "x" * (AUTHORIZATION_MIN_CHARS - 1)
    with pytest.raises(ValueError, match="at least"):
        check_authorization(short)
    assert check_authorization("y" * AUTHORIZATION_MIN_CHARS)


def test_opening_meta_test_requires_an_out_of_tree_directory():
    with pytest.raises(ValueError, match="out-of-tree"):
        GovernedSplitView(VIEWS, visible=(SEALED_SPLIT,),
                          authorization=AUTHORIZED)


def test_a_sealed_directory_inside_the_development_surface_is_refused():
    with pytest.raises(ValueError, match="outside the development surface"):
        GovernedSplitView(VIEWS, visible=(SEALED_SPLIT,),
                          authorization=AUTHORIZED,
                          sealed_directory=VIEWS / "smuggled")


def test_qpsmpdata_refuses_an_isolated_open_without_the_sealed_directory():
    with pytest.raises(ValueError, match="out-of-tree"):
        QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                  split_directory=SPLIT, split_view=VIEWS,
                  include_meta_test=True,
                  meta_test_authorization=AUTHORIZED)


def test_a_sealed_directory_without_a_view_is_refused(tmp_path):
    with pytest.raises(ValueError, match="only meaningful with a split_view"):
        QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                  split_directory=SPLIT,
                  sealed_meta_test_directory=tmp_path)


# --- 3. fail-closed manifests --------------------------------------------

def test_an_absent_manifest_is_refused(tmp_path):
    with pytest.raises(ValueError, match="no governed-split-view manifest"):
        GovernedSplitView(tmp_path, visible=DEVELOPMENT_SPLITS)


def test_a_malformed_manifest_is_refused(tmp_path):
    surface = tmp_path / "views"
    shutil.copytree(VIEWS, surface)
    (surface / "manifest.json").write_text("{not json", encoding="utf-8")
    with pytest.raises(ValueError, match="malformed"):
        GovernedSplitView(surface, visible=DEVELOPMENT_SPLITS)


def test_a_manifest_missing_a_binding_is_refused(tmp_path):
    surface = tmp_path / "views"
    shutil.copytree(VIEWS, surface)
    manifest = json.loads(
        (surface / "manifest.json").read_text(encoding="utf-8"))
    manifest.pop("cell_id_index_sha256")
    (surface / "manifest.json").write_text(json.dumps(manifest),
                                           encoding="utf-8")
    with pytest.raises(ValueError, match="missing binding"):
        GovernedSplitView(surface, visible=DEVELOPMENT_SPLITS)


def test_a_tampered_governance_record_is_refused(tmp_path):
    surface = tmp_path / "views"
    shutil.copytree(VIEWS, surface)
    rows = (surface / "governance.jsonl").read_text(
        encoding="utf-8").splitlines()
    (surface / "governance.jsonl").write_text(
        "\n".join(rows[:-1]) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="governance.jsonl hash mismatch"):
        GovernedSplitView(surface, visible=DEVELOPMENT_SPLITS)


def test_a_view_cut_from_another_corpus_is_refused(tmp_path):
    surface = tmp_path / "views"
    shutil.copytree(VIEWS, surface)
    manifest = json.loads(
        (surface / "manifest.json").read_text(encoding="utf-8"))
    manifest["source_manifest_sha256"] = "0" * 64
    (surface / "manifest.json").write_text(json.dumps(manifest),
                                           encoding="utf-8")
    with pytest.raises(ValueError, match="cut from a different corpus"):
        QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                  split_directory=SPLIT, split_view=surface)


def test_a_view_of_another_split_assignment_is_refused(tmp_path):
    surface = tmp_path / "views"
    shutil.copytree(VIEWS, surface)
    manifest = json.loads(
        (surface / "manifest.json").read_text(encoding="utf-8"))
    manifest["split_assignment_sha256"] = "1" * 64
    (surface / "manifest.json").write_text(json.dumps(manifest),
                                           encoding="utf-8")
    with pytest.raises(ValueError, match="governed split assignment"):
        QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                  split_directory=SPLIT, split_view=surface)


# --- 4. equivalence with the governed source ------------------------------

def test_the_isolated_surface_reproduces_the_corpus_construction_exactly(
        isolated, corpus_surface):
    """Cell membership, ordering and identity, element for element.

    This is what makes the migration inert: every recorded episode index, task
    index and component map stays valid, so no stored number depends on which
    surface a future run mounts.
    """
    assert len(isolated.cells) == len(corpus_surface.cells)
    assert isolated.cells == corpus_surface.cells


def test_task_and_component_identity_is_preserved(isolated, corpus_surface):
    assert set(isolated.tasks) == set(corpus_surface.tasks) == {
        "meta_train", "meta_val"}
    for split, targets in corpus_surface.tasks.items():
        assert set(isolated.tasks[split]) == set(targets)
        for target, indices in targets.items():
            assert list(isolated.tasks[split][target]) == list(indices)
    assert isolated.components == corpus_surface.components


def test_counts_match_the_governed_split_assignment(isolated):
    """Count-equivalence against the frozen assignment, not against a rerun."""
    assignment = json.loads(
        (SPLIT / "assignment.json").read_text(encoding="utf-8"))
    expected: dict[str, int] = {}
    for split in assignment.values():
        expected[split] = expected.get(split, 0) + 1
    counts = isolated.split_view.manifest["counts"]
    for split, rows in expected.items():
        assert counts[split]["rows"] == rows, split
    mounted = {split: 0 for split in DEVELOPMENT_SPLITS}
    for cell in isolated.cells:
        mounted[cell["split"]] += 1
    for split in DEVELOPMENT_SPLITS:
        assert mounted[split] == expected[split]
    assert isolated.sealed_cell_count == expected[SEALED_SPLIT]


def test_the_recorded_hashes_are_enforced_not_merely_stored(isolated):
    manifest = isolated.split_view.manifest
    for split in DEVELOPMENT_SPLITS:
        path = VIEWS / split / "cells.jsonl.gz"
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == manifest["artifacts"][f"{split}/cells.jsonl.gz"]
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            ids = sorted(json.loads(line)["cell_id"]
                         for line in handle if line.strip())
        joined = hashlib.sha256("\n".join(ids).encode("utf-8")).hexdigest()
        assert joined == manifest["cell_id_index_sha256"][split]


def test_the_seal_record_states_the_isolated_surface(isolated, corpus_surface):
    isolated_record = isolated.seal_record()
    assert isolated_record["included"] is False
    assert isolated_record["evaluated"] is False
    assert isolated_record["sealed_cells_withheld"] == 768
    isolation = isolated_record["isolation"]
    assert isolation["physically_isolated"] is True
    assert isolation["level"] == "physically_isolated"
    assert isolation["surface"] == "governed_split_view"
    assert isolation["labels_parsed_in_process"] == list(DEVELOPMENT_SPLITS)

    # The two surfaces must never report the same isolation state.
    default_record = corpus_surface.seal_record()
    assert default_record["isolation"]["physically_isolated"] is False
    assert default_record["isolation"]["level"] == (
        "logical_exclusion_after_parsing")
