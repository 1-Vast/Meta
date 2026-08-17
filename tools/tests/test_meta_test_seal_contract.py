"""Repository contract for the meta_test seal (2026-08-16).

The seal this file pins is **logical exclusion after parsing**, not a physical
label seal, and the tests assert that the artifact says so. `cells.jsonl.gz` is
one all-label file: every sealed label is decompressed and parsed on every
construction and then discarded. `test_the_seal_does_not_claim_physical_
isolation` exists specifically to stop the stronger word coming back.


These tests exist because of a recorded governance incident: `QPSMPData`
defaulted to *admitting* the sealed confirmation split, three analysis scripts
took the default, and one of them wrote a hand-authored seal string into every
artifact it produced that contradicted the arguments it had actually passed.
See `tools/research/a2_readiness_v2/GOVERNANCE_INCIDENT.md`.

The contract has four parts, and each is tested here:

1. the seal is **fail-closed** — the default excludes;
2. exclusion happens **before** task indices and materialization exist, not as
   a filter applied at read time by a well-behaved caller;
3. opening the split requires a **written authorization**, not a bare flag;
4. the seal state written into artifacts comes from the object, so it cannot
   disagree with the constructor arguments.

Part 5 is the reproducibility contract that lets the pre-incident numbers
stand: flipping the seal must not change which cells a `meta_val` episode
draws, only the integer indices used to name them.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from scripts.qpsmp_data import QPSMPData

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
PROTEIN_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank"
LIGAND_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank"
COMPACT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank_compact"
SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"

AUTHORIZED = "contract test: seal-behaviour assertions only, no metric computed"


@pytest.fixture(scope="module")
def sealed() -> QPSMPData:
    """The default construction — no seal keyword passed at all."""
    return QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                     split_directory=SPLIT)


@pytest.fixture(scope="module")
def opened() -> QPSMPData:
    return QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                     split_directory=SPLIT, include_meta_test=True,
                     meta_test_authorization=AUTHORIZED)


# --- 1. fail-closed -------------------------------------------------------

def test_the_default_construction_excludes_meta_test(sealed):
    """No keyword passed: the sealed split must be absent."""
    assert sealed.include_meta_test is False
    assert not any(cell["split"] == "meta_test" for cell in sealed.cells)


def test_the_default_construction_withholds_a_nonzero_population(sealed):
    """A seal that withheld nothing would pass every other test vacuously."""
    assert sealed.sealed_cell_count > 0


# --- 2. exclusion precedes indices and materialization --------------------

def test_no_task_index_covers_a_sealed_cell(sealed):
    assert "meta_test" not in sealed.tasks
    assert "meta_test" not in sealed.components
    for split, targets in sealed.tasks.items():
        for target, indices in targets.items():
            for index in indices:
                assert sealed.cells[int(index)]["split"] == split != "meta_test"


def test_every_index_reachable_by_materialize_is_unsealed(sealed):
    """`materialize` reads `self.cells`, so this is the whole reachable set."""
    assert all(cell["split"] in {"meta_train", "meta_val"}
               for cell in sealed.cells)


def test_sealed_targets_are_absent_from_every_split(sealed, opened):
    sealed_targets = set(opened.tasks["meta_test"])
    assert sealed_targets
    reachable = {target for targets in sealed.tasks.values() for target in targets}
    assert not (sealed_targets & reachable)


def test_an_episode_cannot_be_drawn_from_the_sealed_split(sealed):
    with pytest.raises(KeyError):
        sealed.fixed_episode_bank("meta_test", 2, 4, 1, 17)


# --- 3. opening requires a written authorization --------------------------

def test_a_bare_include_flag_is_refused():
    with pytest.raises(ValueError, match="meta_test_authorization"):
        QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                  split_directory=SPLIT, include_meta_test=True)


@pytest.mark.parametrize("blank", ["", "   ", "\t\n"])
def test_a_whitespace_authorization_is_refused(blank):
    with pytest.raises(ValueError, match="meta_test_authorization"):
        QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                  split_directory=SPLIT, include_meta_test=True,
                  meta_test_authorization=blank)


def test_an_authorization_without_the_flag_is_refused():
    """Refuse the ambiguous state rather than silently picking a meaning."""
    with pytest.raises(ValueError, match="ambiguous seal state"):
        QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT,
                  split_directory=SPLIT, include_meta_test=False,
                  meta_test_authorization="please")


def test_an_authorized_open_actually_exposes_the_split(opened):
    """The refusal above must be the authorization, not a broken code path."""
    assert opened.tasks["meta_test"]
    assert opened.meta_test_authorization == AUTHORIZED


# --- 4. artifacts cannot disagree with the constructor --------------------

def test_the_seal_record_reports_the_sealed_state_truthfully(sealed):
    record = sealed.seal_record()
    assert record["included"] is False
    assert record["evaluated"] is False
    assert record["authorization"] is None
    assert record["sealed_cells_withheld"] == sealed.sealed_cell_count
    assert record["seal"].startswith("logical exclusion after parsing:")


def test_the_seal_does_not_claim_physical_isolation(sealed):
    """The terminology guard.

    Logical exclusion after parsing, metric non-consumption, and physical
    isolation are three different properties. Only the first two hold today,
    and an artifact that claims the third would be making the same class of
    false statement as the 2026-08-16 incident.
    """
    record = sealed.seal_record()
    isolation = record["isolation"]
    assert isolation["level"] == "logical_exclusion_after_parsing"
    assert isolation["labels_parsed_in_process"] is True
    assert isolation["physically_isolated"] is False
    assert "all-label artifact" in isolation["why_not"]
    assert "physical" not in record["seal"].lower().replace(
        "not a physical label seal", "")
    assert (ROOT / "tools/research/a2_readiness_v2/SPLIT_ISOLATION_SPEC.md").exists()


def test_an_open_run_reports_its_isolation_state_too(opened):
    record = opened.seal_record(evaluated=False)
    assert record["isolation"]["level"] == "open"
    assert record["isolation"]["physically_isolated"] is False


def test_the_seal_record_cannot_claim_a_seal_that_is_not_held(opened):
    record = opened.seal_record(evaluated=False)
    assert record["included"] is True
    assert record["seal"].startswith("OPEN:")
    assert AUTHORIZED in record["seal"]
    assert record["sealed_cells_withheld"] == 0


def test_an_open_run_must_state_whether_it_evaluated(opened):
    """The incident was an unstated claim; refuse to infer one."""
    with pytest.raises(ValueError, match="evaluated"):
        opened.seal_record()


# --- 5. the seal does not change what meta_val episodes contain -----------

def _episode_content(data: QPSMPData, spec) -> tuple:
    """Identify an episode by its biology, not by its integer cell indices."""
    return tuple(
        (data.cells[index]["ligand_id"], data.cells[index]["pK"])
        for index in (*spec.support, *spec.query))


def test_sealing_changes_indices_but_not_the_meta_val_episodes(sealed, opened):
    """The reproducibility half of the incident report.

    Dropping sealed cells renumbers `self.cells`, so `EpisodeSpec` integers
    differ between the two constructions. What must not differ is the biology:
    the same target, in the same order, drawing the same ligands with the same
    labels. This is the direct evidence that the pre-incident meta_val results
    were computed from an unsealed *process* but a sealed *population*, and
    therefore stand unchanged.
    """
    banks_sealed = sealed.fixed_nested_episode_banks(
        "meta_val", (0, 1, 2, 3, 5), 16, 1, 73101, None)
    banks_opened = opened.fixed_nested_episode_banks(
        "meta_val", (0, 1, 2, 3, 5), 16, 1, 73101, None)

    assert set(banks_sealed) == set(banks_opened)
    for size in banks_sealed:
        left, right = banks_sealed[size], banks_opened[size]
        assert len(left) == len(right) and left
        for a, b in zip(left, right):
            assert (a.target, a.component, a.donor_target) == (
                b.target, b.component, b.donor_target)
            assert _episode_content(sealed, a) == _episode_content(opened, b)
