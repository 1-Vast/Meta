import numpy as np

from research.crossed_interaction.train_bindingdb_rectangle_descriptor import (
    crossed_feature,
    fit_ridge,
    load_descriptors,
    predict,
    protein_plm_slot_descriptors,
    shuffled_protein_features,
    split_summary,
    wrong_target_features,
)


def test_crossed_feature_changes_sign_under_one_axis_swap():
    protein_delta = np.asarray([1.0, -2.0])
    ligand_delta = np.asarray([3.0, 4.0])

    feature = crossed_feature(protein_delta, ligand_delta)

    assert np.allclose(crossed_feature(-protein_delta, ligand_delta), -feature)
    assert np.allclose(crossed_feature(protein_delta, -ligand_delta), -feature)
    assert np.allclose(crossed_feature(-protein_delta, -ligand_delta), feature)


def test_ridge_quotient_model_has_zero_origin():
    x = np.asarray([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]])
    y = np.asarray([1.0, 2.0, -1.0, -2.0])

    model = fit_ridge(x, y, ridge=0.1)

    assert np.allclose(predict(model, np.zeros((1, 2))), [0.0])


def test_plm_slot_descriptor_metadata_uses_frozen_bank(monkeypatch, tmp_path):
    def fake_read_jsonl(path):
        return [{"sequence_sha256": "p"}]

    def fake_load_protein_bank(path, required):
        assert required == {"p"}
        return {
            "p": {
                "residues": np.asarray([[1.0, 3.0], [5.0, 7.0]]),
                "mask": np.asarray([True, True]),
            }
        }, {
            "model_id": "esm",
            "model_revision": "rev",
            "slot_policy": "slots",
        }

    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(
        "research.crossed_interaction.train_bindingdb_rectangle_descriptor.read_jsonl",
        fake_read_jsonl)
    monkeypatch.setattr(
        "research.crossed_interaction.train_bindingdb_rectangle_descriptor.load_protein_bank",
        fake_load_protein_bank)

    descriptors, metadata = protein_plm_slot_descriptors(
        tmp_path, tmp_path, slot_segments=2, hidden_blocks=1)

    assert np.allclose(descriptors["p"], [2.0, 6.0])
    assert metadata["protein_descriptor"] == "frozen_esm2_slot_region_means"
    assert metadata["protein_descriptor_dim"] == 2
    assert metadata["protein_bank_model_id"] == "esm"


def test_load_descriptors_rejects_unknown_protein_mode(tmp_path):
    try:
        load_descriptors(
            tmp_path,
            protein_descriptor_mode="unknown",
            protein_bank=tmp_path,
            slot_segments=1,
            hidden_blocks=1)
    except ValueError as error:
        assert "unknown protein descriptor mode" in str(error)
    else:
        raise AssertionError("expected unknown protein descriptor mode to fail")


def test_wrong_target_features_use_other_dependency_component():
    examples = [
        {
            "dependency_component": "a",
            "protein_delta": np.asarray([1.0, 0.0]),
            "ligand_delta": np.asarray([1.0]),
        },
        {
            "dependency_component": "b",
            "protein_delta": np.asarray([0.0, 2.0]),
            "ligand_delta": np.asarray([1.0]),
        },
    ]

    features, metadata = wrong_target_features(examples)

    assert features.shape == (2, 2)
    assert metadata["dependency_components"] == 2
    assert metadata["self_donors"] == 0
    assert np.allclose(features[0], [0.0, 2.0])
    assert np.allclose(features[1], [1.0, 0.0])


def test_shuffled_protein_features_do_not_keep_fixed_points():
    examples = [
        {
            "protein_delta": np.asarray([float(index)]),
            "ligand_delta": np.asarray([1.0]),
        }
        for index in range(5)
    ]

    features, metadata = shuffled_protein_features(examples, seed=7)

    assert features.shape == (5, 1)
    assert metadata["fixed_points"] == 0


def test_split_summary_counts_rectangle_entities():
    examples = [
        {
            "panel_id": "p",
            "target_pair": "ta|tb",
            "ligand_a": "la",
            "ligand_b": "lb",
            "dependency_component": "c",
        },
        {
            "panel_id": "p",
            "target_pair": "ta|tc",
            "ligand_a": "la",
            "ligand_b": "lc",
            "dependency_component": "c",
        },
    ]

    summary = split_summary(examples)

    assert summary["rectangles"] == 2
    assert summary["panels"] == 1
    assert summary["target_pairs"] == 2
    assert summary["ligand_pairs"] == 2
    assert summary["dependency_components"] == 1
