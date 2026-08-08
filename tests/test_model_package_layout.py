import model
import torch

from model import (DEFAULT, LigandEncoder, MechanisticInteractionBridge,
                   ProteinEncoder)
from model.meta_operator import BandOperator, CSMO, build_band_operator


def test_model_import_surface_contains_only_verified_components():
    assert not hasattr(model, "MetaSieveDTA")
    assert not hasattr(model, "MetaSieveFrontend")
    assert ProteinEncoder.__module__ == "model.encoders"
    assert LigandEncoder.__module__ == "model.encoders"
    assert MechanisticInteractionBridge.__module__ == "model.mechanism"
    assert BandOperator.__module__ == "model.meta_operator"
    assert CSMO.__module__ == "model.meta_operator"
    assert build_band_operator.__module__ == "model.meta_operator"


def test_band_context_retains_the_frozen_operator_shape():
    operator = BandOperator.__new__(BandOperator)
    torch.nn.Module.__init__(operator)
    operator.cfg = DEFAULT
    context = operator.context(torch.zeros(2, DEFAULT.d_z, dtype=torch.float64))
    assert context.shape == (2,)
