import hashlib

import pytest

from research.correspondence_router.prepare_r0c_download_panel import (
    MIN_SCORABLE_COMPONENTS,
    _homologous,
    _sequence_id,
    select_fresh_components,
)


def row(index: int, sequence: str):
    return {
        "source_entry_id": f"entry-{index}",
        "sequence": sequence,
        "pdb_id": f"p{index:03d}",
        "ligand_comp_id": f"L{index:03d}",
    }


def test_homology_requires_identity_and_bidirectional_coverage():
    assert _homologous("ACDEFG" * 10, "ACDEFG" * 10)
    assert not _homologous("A" * 60, "C" * 60)
    assert not _homologous("ACDEFG" * 10, "ACDEFG" * 20)


def test_component_closure_removes_whole_exposed_component():
    records = [row(i, chr(65 + (i % 20)) * 60 + str(i)) for i in range(130)]
    ids = [_sequence_id(value["sequence"]) for value in records]
    edges = {(ids[0], ids[1]), (ids[1], ids[2])}
    selected, audit = select_fresh_components(
        records, {ids[2]}, edges, limit=MIN_SCORABLE_COMPONENTS
    )
    selected_ids = {value["source_entry_id"] for value in selected}
    assert not {"entry-0", "entry-1", "entry-2"} & selected_ids
    assert audit["exposure_contaminated_components"] == 1
    assert audit["selected_components"] == MIN_SCORABLE_COMPONENTS


def test_download_limit_cannot_undercut_power_target():
    with pytest.raises(ValueError, match="scorable-component"):
        select_fresh_components([], set(), set(), limit=MIN_SCORABLE_COMPONENTS - 1)
