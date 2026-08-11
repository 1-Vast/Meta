from research.correspondence_router.finalize_r0c_panel import finalize


def row(name, ccd, connectivity, scaffold):
    return {
        "source_entry_id": name,
        "ccd_sha256": ccd,
        "connectivity_sha256": connectivity,
        "murcko_scaffold": scaffold,
    }


def test_finalize_excludes_exposure_and_closes_internal_chemistry():
    exposure = [row("old", "c0", "x0", "s0")]
    records = [
        row("a", "c0", "x1", "s1"),
        row("b", "c2", "x2", "shared"),
        row("c", "c3", "x3", "shared"),
        row("d", "c4", "x4", ""),
        row("e", "c5", "x5", ""),
    ]
    selected, exclusions, audit = finalize(records, exposure)
    assert "a" not in {value["source_entry_id"] for value in selected}
    assert len({value["source_entry_id"] for value in selected} & {"b", "c"}) == 1
    assert {"d", "e"} <= {value["source_entry_id"] for value in selected}
    assert audit["chemistry_clean_records"] == 4
    assert audit["selected_components"] == 3
    assert any(value["reason"] == "ccd" for value in exclusions)
