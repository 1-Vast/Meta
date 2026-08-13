import model

from model import (BipartitePairSectionFormer, LigandEncoder,
                   MechanismEvidenceMetaTransformer, ProteinEncoder,
                   QPSMPBioModel)


def test_active_model_surface_excludes_retired_paths():
    assert QPSMPBioModel.__module__ == "model.qpsmp_meta"
    assert BipartitePairSectionFormer.__module__ == "model.bpsf"
    assert MechanismEvidenceMetaTransformer.__module__ == "model.qpsmp_meta"
    assert ProteinEncoder.__module__ == "model.encoders"
    assert LigandEncoder.__module__ == "model.encoders"
    assert not hasattr(model, "MechanisticInteractionBridge")
    assert not hasattr(model, "CenteredRidgeSection")
    assert not hasattr(model, "DifferentiableRidgeSection")
    assert not hasattr(model, "QuotientSupportSetOperator")
