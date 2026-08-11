from research.correspondence_router.prepare_r0c_candidates import (
    select_fresh_candidates,
)
from scripts.structure_sources.biolip import BioLiPEntry


def entry(name: str, sequence: str, ligand: str, pdb: str | None = None):
    return BioLiPEntry(
        pdb_id=pdb or name[:4],
        receptor_auth_asym_id="A",
        resolution=2.0,
        binding_site_id="BS01",
        ligand_comp_id=ligand,
        ligand_auth_asym_id="B",
        ligand_serial="1",
        ligand_auth_seq_id="1",
        sequence=sequence * 10,
    )


def test_fresh_selection_excludes_all_old_metadata_axes():
    old = [entry("old1", "AAAAAA", "OLD", pdb="old1")]
    rows = [
        old[0],
        entry("same", "CCCCCC", "NEW", pdb="old1"),
        entry("seq1", "AAAAAA", "N01"),
        entry("lig1", "DDDDDD", "OLD"),
        entry("new1", "EEEEEE", "N02"),
        entry("new2", "FFFFFF", "N03"),
    ]
    selected, audit = select_fresh_candidates(
        rows, old, {"OLD", "NEW", "N01", "N02", "N03"}, limit=2
    )
    assert {value.ligand_comp_id for value in selected} == {"N02", "N03"}
    assert audit["selected_unique_pdb"] == 2
    assert audit["affinity_value_reads"] == 0
    assert audit["coordinate_value_reads"] == 0


def test_fresh_selection_is_order_invariant_and_one_to_one():
    old = [entry("old1", "AAAAAA", "OLD")]
    rows = [
        entry("new1", "CCCCCC", "N01", pdb="p001"),
        entry("new2", "CCCCCC", "N02", pdb="p002"),
        entry("new3", "DDDDDD", "N01", pdb="p003"),
        entry("new4", "EEEEEE", "N03", pdb="p004"),
    ]
    left = select_fresh_candidates(rows, old, {"OLD", "N01", "N02", "N03"}, limit=2)[0]
    right = select_fresh_candidates(
        list(reversed(rows)), old, {"OLD", "N01", "N02", "N03"}, limit=2
    )[0]
    assert left == right
    assert len({value.sequence for value in left}) == 2
    assert len({value.ligand_comp_id for value in left}) == 2
    assert len({value.pdb_id for value in left}) == 2
