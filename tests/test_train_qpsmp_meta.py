from types import SimpleNamespace

import numpy as np
import torch

from model.qpsmp_meta import QPSMPBioModel
from scripts.train_qpsmp import (
    LabelScale, freeze_for_section_training, training_label_scale)


def test_label_scale_uses_meta_train_only():
    data = SimpleNamespace(cells=[
        {"split": "meta_train", "pK": 2.0},
        {"split": "meta_train", "pK": 4.0},
        {"split": "meta_val", "pK": 1000.0},
        {"split": "meta_test", "pK": -1000.0},
    ])

    scale = training_label_scale(data)

    assert np.isclose(scale.mean, 3.0)
    assert np.isclose(scale.scale, 1.0)


def test_label_scale_reports_squared_error_in_raw_pk_units():
    scale = LabelScale(mean=6.0, scale=2.0)
    prediction = torch.tensor([0.0, 1.0])
    truth = torch.tensor([1.0, 0.0])

    assert torch.equal(scale.squared_error_pk(prediction, truth), torch.tensor([4.0, 4.0]))


def test_label_scale_aligns_truth_dtype_to_prediction():
    scale = LabelScale(mean=6.0, scale=2.0)
    prediction = torch.tensor([1.0], dtype=torch.float64)
    truth = torch.tensor([0.0], dtype=torch.float32)

    error = scale.squared_error_pk(prediction, truth)

    assert error.dtype == prediction.dtype
    assert torch.allclose(error, torch.tensor([4.0], dtype=torch.float64))


def test_section_only_freezes_potential_and_trains_only_geometry_and_ridge():
    model = QPSMPBioModel(8, 6, 3, ligand_layers=1, section_mode="ridge")

    freeze_for_section_training(model)

    trainable = {name for name, parameter in model.named_parameters()
                 if parameter.requires_grad}
    assert trainable == {
        "meta.section_head.weight", "meta.support_span.ridge_raw"}


def test_neural_section_only_trains_geometry_and_amortized_adapter():
    model = QPSMPBioModel(8, 6, 3, ligand_layers=1, section_mode="neural")

    freeze_for_section_training(model)

    trainable = {name for name, parameter in model.named_parameters()
                 if parameter.requires_grad}
    assert "meta.section_head.weight" in trainable
    assert any(name.startswith("meta.adapter.") for name in trainable)
    assert "meta.support_span.ridge_raw" not in trainable
    assert all(name.startswith(("meta.section_head.", "meta.adapter."))
               for name in trainable)


def test_qp_ams_section_only_trains_only_geometry_and_set_adapter():
    model = QPSMPBioModel(8, 6, 3, ligand_layers=1, section_mode="qp_ams")

    freeze_for_section_training(model)

    trainable = {name for name, parameter in model.named_parameters()
                 if parameter.requires_grad}
    assert "meta.section_head.weight" in trainable
    assert any(name.startswith("meta.qp_ams.") for name in trainable)
    assert all(name.startswith(("meta.section_head.", "meta.qp_ams."))
               for name in trainable)


def test_section_former_stage_freezes_zero_shot_potential():
    model = QPSMPBioModel(
        8, 8, 4, ligand_layers=1, section_mode="section_former",
        interaction_mode="bpsf", pair_blocks=1, pair_heads=2)

    freeze_for_section_training(model)

    trainable = {name for name, parameter in model.named_parameters()
                 if parameter.requires_grad}
    assert any(name.startswith("pair_section.latent.cross.") for name in trainable)
    assert any(name.startswith("meta.section_former.") for name in trainable)
    assert not any(name.startswith("pair_section.latent.endpoint.") for name in trainable)
    assert not any(name.startswith("meta.zero_shot_head.") for name in trainable)
