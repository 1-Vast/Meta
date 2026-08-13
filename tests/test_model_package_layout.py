import model

from model import (BipartitePairSectionFormer, LigandEncoder, ProteinEncoder,
                   QPSMPBioModel, QuotientSupportSetOperator)


def test_active_model_surface_excludes_retired_paths():
    assert QPSMPBioModel.__module__ == "model.qpsmp_meta"
    assert BipartitePairSectionFormer.__module__ == "model.bpsf"
    assert QuotientSupportSetOperator.__module__ == "model.bpsf"
    assert ProteinEncoder.__module__ == "model.encoders"
    assert LigandEncoder.__module__ == "model.encoders"
    assert not hasattr(model, "MechanisticInteractionBridge")
    assert not hasattr(model, "CenteredRidgeSection")
