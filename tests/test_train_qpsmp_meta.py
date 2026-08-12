from types import SimpleNamespace

import numpy as np
import torch

from scripts.train_qpsmp import LabelScale, training_label_scale


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
