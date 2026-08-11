from research.correspondence_router.prepare_r0b_panel import select_panel


def row(name, split, group, smiles, scaffold, coverage=1.0):
    return {
        "source_entry_id": name,
        "source_split": split,
        "homology_group_id": group,
        "pdb_id": name[:4],
        "sequence": "ACDEFG",
        "sequence_sha256": f"seq-{name}",
        "protein_mapping_coverage": coverage,
        "canonical_smiles": smiles,
        "murcko_scaffold": scaffold,
        "ccd_sha256": f"ccd-{name}",
        "connectivity_sha256": smiles,
    }


def test_panel_selects_one_train_record_per_family_and_filters_chemistry():
    records = [
        row("tr-a1", "train", "a", "C1", "S1"),
        row("tr-a2", "train", "a", "C2", "S2"),
        row("tr-b1", "train", "b", "C3", "S3"),
        row("tr-z1", "train", "z", "T0", "TZ", coverage=0.9),
        row("va-c1", "val", "c", "V1", "VS"),
        row("va-d1", "val", "d", "C3", "NEW"),
        row("te-e1", "test", "e", "T1", "TS"),
        row("te-f1", "test", "f", "T2", "VS"),
        row("te-g1", "test", "g", "T3", "TG", coverage=0.9),
        row("te-h1", "test", "h", "T0", "TH"),
    ]
    selected, exclusions, audit = select_panel(records)
    train = [value for value in selected if value["r0_split"] == "train"]
    val = [value for value in selected if value["r0_split"] == "val"]
    heldout = [value for value in selected if value["r0_split"] == "heldout_a"]
    assert len(train) == 2
    assert {value["homology_group_id"] for value in train} == {"a", "b"}
    assert [value["source_entry_id"] for value in val] == ["va-c1"]
    assert [value["source_entry_id"] for value in heldout] == ["te-e1"]
    assert any("exact_connectivity" in value["reason"] for value in exclusions)
    assert any("scaffold" in value["reason"] for value in exclusions)
    assert audit["affinity_value_reads"] == 0


def test_panel_selection_is_order_invariant():
    records = [
        row("tr-a1", "train", "a", "C1", "S1"),
        row("tr-a2", "train", "a", "C2", "S2"),
        row("va-b1", "val", "b", "V1", "VS"),
        row("te-c1", "test", "c", "T1", "TS"),
    ]
    left = select_panel(records)[0]
    right = select_panel(list(reversed(records)))[0]
    assert left == right
