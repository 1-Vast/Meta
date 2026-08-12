import gzip
import json

import numpy as np

from research.crossed_interaction.train_physchem_cq_observable import (
    ligand_pharmacophore_descriptor,
    materialize_features,
    protein_physchem_descriptor,
)


def test_protein_physchem_descriptor_has_expected_lag_channels():
    descriptor = protein_physchem_descriptor("ACDEFGHIKLMNPQRSTVWY", max_lag=4)
    assert descriptor.shape == (20,)
    assert np.isfinite(descriptor).all()


def test_protein_physchem_descriptor_is_order_sensitive():
    first = protein_physchem_descriptor("AAAAACCCCCDDDDD", max_lag=4)
    second = protein_physchem_descriptor("ACDACDACDACDACD", max_lag=4)
    assert not np.allclose(first, second)


def test_ligand_pharmacophore_descriptor_is_finite():
    descriptor = ligand_pharmacophore_descriptor("CCO")
    assert descriptor.shape == (40,)
    assert np.isfinite(descriptor).all()
    assert descriptor[35] > 0.0


def test_materialize_features_reports_expected_dimensions(tmp_path):
    corpus = tmp_path
    (corpus / "manifest.json").write_text("{}", encoding="utf-8")
    (corpus / "proteins.jsonl").write_text(
        '{"sequence_sha256":"t1","sequence":"ACDEFGHIKLMNPQRSTVWY"}\n'
        '{"sequence_sha256":"t2","sequence":"WYWYVVVVKKKKDDDD"}\n',
        encoding="utf-8")
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

    features, metadata = materialize_features(corpus, max_lag=3)
    assert metadata["protein_descriptor_dim"] == 15
    assert metadata["ligand_descriptor_dim"] == 40
    assert metadata["feature_dim"] == 600
    assert features["c1"]["correct"].shape == (600,)
