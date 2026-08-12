import json

import numpy as np

from research.crossed_interaction.train_chembl_affinity_teacher_cq_observable import (
    bindingdb_ligand_connectivity_keys,
    fit_affinity_teacher_from_arrays,
    load_external_teacher_training_arrays,
    predict_affinity_teacher,
    teacher_feature,
    two_way_residual_labels,
)


def test_bindingdb_ligand_connectivity_keys_use_inchikey_prefix(tmp_path):
    (tmp_path / "ligands.jsonl").write_text(
        '{"drug_key":"ABCDEFGHIJKLMN-UHFFFAOYSA-N","smiles":"CCO"}\n'
        '{"drug_key":"NOPQRSTUVWXYZA-UHFFFAOYSA-N","smiles":"CCN"}\n',
        encoding="utf-8")

    assert bindingdb_ligand_connectivity_keys(tmp_path) == {
        "ABCDEFGHIJKLMN",
        "NOPQRSTUVWXYZA",
    }


def test_fit_affinity_teacher_requires_positive_ridge_and_predicts_signal():
    x = np.asarray([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
    ], dtype=np.float64)
    y = np.asarray([1.0, 3.0, 2.0, 4.0], dtype=np.float64)

    teacher = fit_affinity_teacher_from_arrays(x, y, ridge=0.01)

    assert teacher["train_mse"] < 0.01
    assert predict_affinity_teacher(teacher, np.asarray([1.0, 1.0])) > 3.8


def test_teacher_feature_modes_are_finite_and_dimensioned():
    x = np.asarray([
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
    ], dtype=np.float64)
    y = np.asarray([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    teacher = fit_affinity_teacher_from_arrays(x, y, ridge=0.1)
    product = np.asarray([2.0, 3.0, 4.0], dtype=np.float64)

    weighted = teacher_feature(product, teacher, mode="weighted_product")
    prediction = teacher_feature(product, teacher, mode="prediction")

    assert weighted.shape == (3,)
    assert prediction.shape == (1,)
    assert np.isfinite(weighted).all()
    assert np.isfinite(prediction).all()


def test_two_way_residual_labels_remove_task_and_ligand_means():
    tasks = ["t1", "t1", "t2", "t2"]
    ligands = ["l1", "l2", "l1", "l2"]
    labels = np.asarray([1.0, 3.0, 4.0, 6.0], dtype=np.float64)

    residual = two_way_residual_labels(tasks, ligands, labels)

    assert abs(residual.mean()) < 1e-9
    for task in sorted(set(tasks)):
        assert abs(residual[[index for index, value in enumerate(tasks) if value == task]].mean()) < 1e-9
    for ligand in sorted(set(ligands)):
        assert abs(residual[[index for index, value in enumerate(ligands) if value == ligand]].mean()) < 1e-9


def test_external_teacher_loader_excludes_bindingdb_targets_and_ligands(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    rows = [
        {
            "task_id": "kept",
            "protein_sequence_sha256": "blocked_target",
            "protein_sequence": "ACDEFGHIKLMNPQRSTVWY",
            "ligand_connectivity_key": "FREEKEY0000001",
            "canonical_smiles": "CCO",
            "p_affinity": 6.0,
        },
        {
            "task_id": "kept",
            "protein_sequence_sha256": "free_target",
            "protein_sequence": "ACDEFGHIKLMNPQRSTVWY",
            "ligand_connectivity_key": "BLOCKED0000001",
            "canonical_smiles": "CCO",
            "p_affinity": 6.0,
        },
    ]
    for index in range(10):
        rows.append({
            "task_id": "kept",
            "protein_sequence_sha256": f"free_target_{index}",
            "protein_sequence": "ACDEFGHIKLMNPQRSTVWY",
            "ligand_connectivity_key": f"FREEKEY{index:07d}",
            "canonical_smiles": "CCO",
            "p_affinity": 6.0 + index / 10.0,
        })
    (source_dir / "canonical_rows.jsonl").write_text(
        "".join(json.dumps(row) + "\n" for row in rows),
        encoding="utf-8")

    x, y, counters = load_external_teacher_training_arrays(
        source_dir=source_dir,
        governed_task_ids={"kept"},
        blocked_targets={"blocked_target"},
        blocked_ligands={"BLOCKED0000001"},
        max_source_rows=10)

    assert x.shape[0] == 10
    assert y.shape == (10,)
    assert counters["skipped_blocked_target"] == 1
    assert counters["skipped_blocked_ligand"] == 1
    assert counters["source_rows_used"] == 10
