from types import SimpleNamespace

import numpy as np
import torch

from model.qpsmp_meta import QPSMPBioModel
from scripts.train_qpsmp import LabelScale, freeze_for_section_training, training_label_scale


def test_label_scale_uses_meta_train_only():
    data = SimpleNamespace(cells=[
        {"split": "meta_train", "pK": 2.0},
        {"split": "meta_train", "pK": 4.0},
        {"split": "meta_val", "pK": 1000.0},
    ])
    scale = training_label_scale(data)
    assert np.isclose(scale.mean, 3.0)
    assert np.isclose(scale.scale, 1.0)


def test_label_scale_reports_raw_pk_squared_error():
    scale = LabelScale(6.0, 2.0)
    assert torch.equal(
        scale.squared_error_pk(torch.tensor([0.0, 1.0]), torch.tensor([1.0, 0.0])),
        torch.tensor([4.0, 4.0]))


def test_section_stage_freezes_scalar_potential_and_trains_learned_operator():
    model = QPSMPBioModel(
        8, hidden_dim=16, task_dim=4, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2,
        support_hidden_dim=16, support_blocks=1)
    freeze_for_section_training(model)
    trainable = {name for name, parameter in model.named_parameters()
                 if parameter.requires_grad}
    assert any(name.startswith("pair_section.latent.section.") for name in trainable)
    assert any(name.startswith("meta.section_operator.") for name in trainable)
    assert not any(name.startswith("pair_section.latent.endpoint.") for name in trainable)
    assert all(name.startswith("pair_section.latent.section.") or
               name.startswith("meta.section_operator.") for name in trainable)
