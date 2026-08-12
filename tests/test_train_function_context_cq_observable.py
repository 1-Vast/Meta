import gzip
import json

import numpy as np

from research.crossed_interaction.train_function_context_cq_observable import (
    FUNCTION_KEYWORDS,
    materialize_features,
    read_projection_names,
    target_function_descriptor,
)


def test_target_function_descriptor_maps_known_classes():
    descriptor = target_function_descriptor([
        "Tyrosine-protein kinase receptor",
        "Human immunodeficiency virus protease",
    ])
    names = [name for name, _ in FUNCTION_KEYWORDS]
    assert descriptor[names.index("kinase")] == 1.0
    assert descriptor[names.index("protease")] == 1.0
    assert descriptor[names.index("viral")] == 1.0
    assert descriptor[names.index("unknown_or_other")] == 0.0


def test_target_function_descriptor_sets_unknown_when_no_keyword_matches():
    descriptor = target_function_descriptor(["Uncharacterized binding protein"])
    assert descriptor[-1] == 1.0
    assert descriptor[:-1].sum() == 0.0


def test_read_projection_names_filters_targets(tmp_path):
    path = tmp_path / "projection.jsonl.gz"
    rows = [
        {"target_sequence_sha256": "t1", "target_name": "Protein kinase A"},
        {"target_sequence_sha256": "t2", "target_name": "Ignored"},
        {"target_sequence_sha256": "t1", "target_name": "Protein kinase A"},
    ]
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row) + "\n")
    names, metadata = read_projection_names(path, {"t1", "missing"})
    assert names["t1"] == ["Protein kinase A"]
    assert names["missing"] == []
    assert metadata["projection_rows_scanned"] == 3
    assert metadata["projection_rows_matched"] == 2


def test_materialize_features_unknown_control_dimensions(tmp_path):
    corpus = tmp_path / "corpus"
    corpus.mkdir()
    (corpus / "manifest.json").write_text("{}", encoding="utf-8")
    (corpus / "ligands.jsonl").write_text(
        '{"drug_key":"l1","smiles":"CCO"}\n'
        '{"drug_key":"l2","smiles":"CCN"}\n',
        encoding="utf-8")
    cells = [
        {
            "cell_id": "c1", "target_id": "t1", "ligand_id": "l1",
            "protein_group_40": "g1", "scaffold": "s1",
        },
        {
            "cell_id": "c2", "target_id": "t2", "ligand_id": "l2",
            "protein_group_40": "g2", "scaffold": "s2",
        },
    ]
    with gzip.open(corpus / "cells.jsonl.gz", "wt", encoding="utf-8") as handle:
        for cell in cells:
            handle.write(json.dumps(cell) + "\n")
    projection = tmp_path / "projection.jsonl.gz"
    with gzip.open(projection, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps({
            "target_sequence_sha256": "t1",
            "target_name": "Protein kinase A",
        }) + "\n")
    features, metadata = materialize_features(
        corpus, projection, function_mode="unknown_only")
    assert metadata["target_descriptor_dim"] == len(FUNCTION_KEYWORDS)
    assert metadata["ligand_descriptor_dim"] == 40
    assert features["c1"]["correct"].shape == (len(FUNCTION_KEYWORDS) * 40,)
    assert np.isfinite(features["c1"]["correct"]).all()
