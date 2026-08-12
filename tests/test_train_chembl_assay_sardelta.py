import json
from pathlib import Path

import numpy as np

from research.source_affinity.train_chembl_assay_sardelta import (
    assay_rows,
    build_matched_pairs,
    fit_ridge,
    pair_feature,
    predict,
)


def test_pair_feature_delta_and_concat_shapes():
    left = np.asarray([1.0, 2.0, 3.0], dtype=np.float64)
    right = np.asarray([0.5, 1.0, 1.5], dtype=np.float64)
    delta = pair_feature(left, right, target=None, mode="delta")
    concat = pair_feature(left, right, target=None, mode="concat")
    target = np.asarray([9.0, 8.0], dtype=np.float64)
    delta_target = pair_feature(left, right, target=target, mode="delta_target")
    assert delta.shape == (3,)
    assert concat.shape == (6,)
    assert delta_target.shape == (5,)
    assert np.allclose(delta, [0.5, 1.0, 1.5])


def test_assay_rows_merge_replicate_ligands_by_median(tmp_path: Path):
    path = tmp_path / "CHEMBL1.jsonl"
    rows = [
        {
            "assay_chembl_id": "CHEMBL1",
            "document_chembl_id": "DOC1",
            "target_chembl_id": "T1",
            "target_accession": "P1",
            "ligand_connectivity_key": "AAAA",
            "standard_inchi_key": "AAAA-UHFFFAOYSA-N",
            "canonical_smiles": "CCO",
            "p_value": 7.0,
            "endpoint_family": "Ki",
            "activity_id": 1,
        },
        {
            "assay_chembl_id": "CHEMBL1",
            "document_chembl_id": "DOC1",
            "target_chembl_id": "T1",
            "target_accession": "P1",
            "ligand_connectivity_key": "AAAA",
            "standard_inchi_key": "AAAA-UHFFFAOYSA-N",
            "canonical_smiles": "CCO",
            "p_value": 8.0,
            "endpoint_family": "Ki",
            "activity_id": 2,
        },
    ]
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")
    merged = assay_rows(path)
    assert len(merged) == 1
    assert merged[0]["p_value"] == 7.5
    assert merged[0]["scaffold"] == ""


def test_build_matched_pairs_requires_same_scaffold_and_similarity(tmp_path: Path):
    rows = [
        {
            "assay_chembl_id": "CHEMBL1",
            "document_chembl_id": "DOC1",
            "target_chembl_id": "T1",
            "target_accession": "P1",
            "ligand_connectivity_key": "A",
            "standard_inchi_key": "A-UHFFFAOYSA-N",
            "canonical_smiles": "CCO",
            "p_value": 7.0,
            "endpoint_family": "Ki",
            "activity_id": 1,
            "scaffold": "",
        },
        {
            "assay_chembl_id": "CHEMBL1",
            "document_chembl_id": "DOC1",
            "target_chembl_id": "T1",
            "target_accession": "P1",
            "ligand_connectivity_key": "B",
            "standard_inchi_key": "B-UHFFFAOYSA-N",
            "canonical_smiles": "CCN",
            "p_value": 8.0,
            "endpoint_family": "Ki",
            "activity_id": 2,
            "scaffold": "",
        },
    ]
    pairs = build_matched_pairs(rows, feature_mode="delta_target")
    assert pairs == []


def test_build_matched_pairs_includes_target_context(tmp_path: Path):
    rows = [
        {
            "assay_chembl_id": "CHEMBL1",
            "document_chembl_id": "DOC1",
            "target_chembl_id": "T1",
            "target_accession": "P1",
            "protein_sequence": "ACDEFGHIKLMNPQRSTVWY",
            "ligand_connectivity_key": "A",
            "standard_inchi_key": "A-UHFFFAOYSA-N",
            "canonical_smiles": "CCOc1ccccc1",
            "p_value": 7.0,
            "endpoint_family": "Ki",
            "activity_id": 1,
            "scaffold": "c1ccccc1",
        },
        {
            "assay_chembl_id": "CHEMBL1",
            "document_chembl_id": "DOC1",
            "target_chembl_id": "T1",
            "target_accession": "P1",
            "protein_sequence": "ACDEFGHIKLMNPQRSTVWY",
            "ligand_connectivity_key": "B",
            "standard_inchi_key": "B-UHFFFAOYSA-N",
            "canonical_smiles": "CCOc1ccc(Cl)cc1",
            "p_value": 8.0,
            "endpoint_family": "Ki",
            "activity_id": 2,
            "scaffold": "c1ccccc1",
        },
    ]
    pairs = build_matched_pairs(rows, feature_mode="delta_target")
    assert len(pairs) == 1
    assert pairs[0]["feature"].shape[0] > 3


def test_ridge_fit_and_predict_are_finite():
    x = np.asarray([[0.0, 1.0], [1.0, 0.0], [1.0, 1.0]], dtype=np.float64)
    y = np.asarray([1.0, 2.0, 3.0], dtype=np.float64)
    model = fit_ridge(x, y, ridge=0.1)
    pred = predict(model, x)
    assert pred.shape == (3,)
    assert np.isfinite(pred).all()
