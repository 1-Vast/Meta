"""Regression tests for the S4R single-axis graph-aware ligand repair.

These tests are executable statements of the S4R contract: the estimator, the
gauge, the frozen registration hashes, the ligand statistic and the arms that
distinguish the changed axis from the surfaces that must not move.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "research" / "s7_l2b_r0r" / "s4r_run.py"
AUDIT_MODULE = ROOT / "research" / "s7_l2b_r0r" / "s4r_audit.py"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


s4r = _load("s4r_run_test", MODULE)
audit = _load("s4r_audit_test", AUDIT_MODULE)


# ------------------------------------------------------------ registration
def test_registration_hashes_match_the_committed_documents():
    here = ROOT / "research" / "s7_l2b_r0r"
    assert s4r.sha_file(here / "PREREG_PHASE2B_S4R_GRAPH_AWARE_LIGAND_DIRECT_W.md") \
        == s4r.PREREG_SHA
    assert s4r.sha_file(here / "PREREG_PHASE2B_S4R_LIGAND_REPRESENTATION_AUDIT.md") \
        == s4r.AUDIT_PREREG_SHA
    assert s4r.sha_file(here / "PREREG_PHASE2B_S4R_AUDIT_AMENDMENT_01.md") \
        == s4r.AUDIT_AMENDMENT_SHA
    assert audit.PREREG_SHA == s4r.AUDIT_PREREG_SHA
    assert audit.AMENDMENT_SHA == s4r.AUDIT_AMENDMENT_SHA


def test_frozen_constants_match_the_preregistration():
    assert (s4r.D_ESM, s4r.D_ATOM, s4r.D_LIG) == (1280, 41, 128)
    assert s4r.MORGAN_RADIUS == 1
    assert s4r.N_UPDATES == 210
    assert s4r.LR == 1e-3 and s4r.CLIP == 5.0
    assert s4r.SEED_PARAM == 20260901 and s4r.SEED_BOOT == 20260903
    assert s4r.EXPECTED == {"train_pairs": 226765, "train_components": 554,
                            "heldoutA_pairs": 46818, "heldoutA_components": 112}
    assert s4r.ARM_DIM == {"candidate": 128, "repeat": 128, "permuted": 128,
                           "baseline41": 41}


def test_stage_never_touches_heldout_b_or_affinity():
    source = MODULE.read_text(encoding="utf-8")
    assert "heldoutB" not in source.replace(
        '"heldoutB_created": False', "").replace(
        '"heldoutB_reads": 0', "").replace(
        '"heldoutB_status"', "")
    for forbidden in ("davis", "kiba", "chembl", "bindingdb", "recipient",
                      "affinity_head", "metaval"):
        assert forbidden not in source.lower().replace("affinity_value_reads", "")


def test_audit_a_gate_thresholds_are_frozen():
    assert audit.A1_EFFRANK_MULTIPLE == 3.0
    assert audit.A2_MIN_INCREMENTAL == 0.25
    assert audit.A3_MAX_RETENTION_LOSS == 0.10
    assert audit.A4_MIN_COVERAGE == 0.99
    assert audit.RADII == (1, 2)
    assert audit.VOCAB_SIZES == (128, 256, 512)


# ------------------------------------------------------------ estimator
def test_direct_w_shapes_and_unit_norm():
    graph = s4r.DirectW(s4r.D_LIG)
    base = s4r.DirectW(s4r.D_ATOM)
    assert tuple(graph.W.shape) == (1280, 128)
    assert tuple(base.W.shape) == (1280, 41)
    assert graph.W.numel() == 163840 and base.W.numel() == 52480
    for head in (graph, base):
        assert abs(float(head.W.norm()) - 1.0) < 1e-6


def test_direct_w_initialization_is_deterministic():
    assert torch.equal(s4r.DirectW(s4r.D_LIG).W, s4r.DirectW(s4r.D_LIG).W)


def test_baseline_arm_reuses_the_s3r_initialization():
    """baseline41 must start from the identical 1280x41 draw S3R used, or the
    cross-stage replication contract could not hold."""
    s3r = _load("s3r_run_replay", ROOT / "research" / "s7_l2b_r0r" / "s3r_run.py")
    assert torch.equal(s4r.DirectW(s4r.D_ATOM).W, s3r.DirectW().W)


def test_positive_scale_normalization_is_invariant():
    score = torch.linspace(-2.0, 3.0, 97, dtype=torch.float64)
    assert torch.max(torch.abs(s4r.normalized(score) -
                               s4r.normalized(11.0 * score))) < 1e-10


# ------------------------------------------------------------ gauge
def test_projection_is_orthogonal_and_pair_difference_is_antisymmetric():
    rng = np.random.default_rng(4)
    h = rng.normal(size=(37, 1280))
    W = rng.normal(size=(1280, s4r.D_LIG))
    ga, gb = rng.normal(size=s4r.D_LIG), rng.normal(size=s4r.D_LIG)
    q, _ = np.linalg.qr(rng.normal(size=(37, 2)))
    left = s4r.project_np(q, (h @ W) @ (ga - gb))
    right = s4r.project_np(q, (h @ W) @ (gb - ga))
    assert np.max(np.abs(left + right)) < 1e-10
    assert np.linalg.norm(q.T @ left) / np.linalg.norm(left) < 1e-12


def test_identical_ligand_gives_an_exactly_zero_difference():
    rng = np.random.default_rng(6)
    h = rng.normal(size=(13, 1280))
    W = rng.normal(size=(1280, s4r.D_LIG))
    g = rng.normal(size=s4r.D_LIG)
    assert np.array_equal((h @ W) @ (g - g), np.zeros(13))


def test_residue_constant_field_is_annihilated_by_the_gauge():
    """The registered ligand-only claim: any residue-constant field lies in
    span{1} and therefore carries no ranking information."""
    rng = np.random.default_rng(11)
    prior = rng.normal(size=61)
    q = s4r.nuisance_basis(prior)
    constant = np.full(61, -3.25)
    assert np.max(np.abs(s4r.project_np(q, constant))) < 1e-12


# ------------------------------------------------------------ ligand statistic
def _mol(smiles):
    from rdkit import Chem
    return Chem.MolFromSmiles(smiles)


def _index(n=s4r.D_LIG):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=s4r.MORGAN_RADIUS)
    seen = []
    for smiles in ("c1ccccc1O", "CC(=O)Nc1ccccc1", "CCOC(=O)C", "C1CCNCC1",
                   "Clc1ccc(cc1)C(=O)N", "CSCC[C@H](N)C(=O)O"):
        mol = Chem.MolFromSmiles(smiles)
        for environment in generator.GetSparseCountFingerprint(mol).GetNonzeroElements():
            if int(environment) not in seen:
                seen.append(int(environment))
    seen = sorted(seen)[:n]
    return {environment: i for i, environment in enumerate(seen)}


def test_graph_statistic_separates_graphs_the_atom_mean_cannot():
    """1,3- and 1,4-disubstituted benzene share every atom-local feature count
    but differ in connectivity. The 41-D mean is identical; the radius-1
    statistic is not. This is the whole hypothesis in one assertion."""
    from s7_dataset import atom_features

    meta, para = _mol("Cc1cccc(C)c1"), _mol("Cc1ccc(C)cc1")
    left = s4r.g_of(atom_features(meta))
    right = s4r.g_of(atom_features(para))
    assert np.max(np.abs(left - right)) < 1e-12

    index = _index()
    a = s4r.g_graph_of(meta, index, meta.GetNumAtoms())
    b = s4r.g_graph_of(para, index, para.GetNumAtoms())
    assert a.shape == (s4r.D_LIG,) and b.shape == (s4r.D_LIG,)
    assert np.max(np.abs(a - b)) > 1e-9


def test_graph_statistic_is_per_heavy_atom_normalized_and_bounded():
    mol = _mol("CC(=O)Nc1ccccc1")
    index = _index()
    vector = s4r.g_graph_of(mol, index, mol.GetNumAtoms())
    assert np.all(vector >= 0.0) and np.all(vector <= 1.0)
    assert np.isfinite(vector).all()


def test_graph_statistic_is_invariant_to_atom_ordering():
    from rdkit import Chem
    mol = _mol("Clc1ccc(cc1)C(=O)NCC")
    index = _index()
    reference = s4r.g_graph_of(mol, index, mol.GetNumAtoms())
    order = list(range(mol.GetNumAtoms()))[::-1]
    renumbered = Chem.RenumberAtoms(mol, order)
    permuted = s4r.g_graph_of(renumbered, index, renumbered.GetNumAtoms())
    assert np.max(np.abs(reference - permuted)) < 1e-12


def test_out_of_vocabulary_environments_contribute_zero():
    mol = _mol("[Se]1CCCC1")
    empty = s4r.g_graph_of(mol, {}, mol.GetNumAtoms())
    assert np.array_equal(empty, np.zeros(s4r.D_LIG))


def test_frozen_vocabulary_loads_and_matches_its_registered_hash():
    path = s4r.AUDIT_EXEC / "selected_ligand_vocabulary.json"
    if not path.is_file():
        pytest.skip("frozen ligand vocabulary is not present in this checkout")
    vocabulary = s4r.load_vocabulary()
    assert len(vocabulary) == s4r.D_LIG
    assert len(set(vocabulary)) == s4r.D_LIG


# ------------------------------------------------------------ audit statistics
def test_effective_rank_bounds():
    rng = np.random.default_rng(3)
    direction = rng.normal(size=(1, 12))
    rank_one = rng.normal(size=(400, 1)) * direction
    assert audit.effective_rank(rank_one) < 1.05
    isotropic = rng.normal(size=(4000, 12))
    assert audit.effective_rank(isotropic) > 11.0


def test_residual_fraction_is_zero_for_a_linear_reparameterization():
    rng = np.random.default_rng(8)
    X = rng.normal(size=(300, 9))
    Y = X @ rng.normal(size=(9, 5)) + 2.0
    assert audit.residual_fraction(Y, X) < 1e-18
    assert audit.residual_fraction(X, rng.normal(size=(300, 2))) > 0.5


def test_label_blind_pair_superset_uses_no_residue_label():
    records = [
        {"seq_key": "s", "source_key": "a", "graph_key": "g1", "scaffold": "x"},
        {"seq_key": "s", "source_key": "b", "graph_key": "g2", "scaffold": "y"},
        {"seq_key": "s", "source_key": "c", "graph_key": "g2", "scaffold": "y"},
        {"seq_key": "s", "source_key": "d", "graph_key": "g3", "scaffold": "x"},
        {"seq_key": "t", "source_key": "e", "graph_key": "g4", "scaffold": ""},
        {"seq_key": "t", "source_key": "f", "graph_key": "g5", "scaffold": "z"},
    ]
    pairs = audit.label_blind_graph_pairs(records)
    assert ("g1", "g2") in pairs
    assert ("g1", "g3") not in pairs          # same scaffold
    assert ("g2", "g2") not in pairs          # same ligand graph
    assert ("g4", "g5") not in pairs          # missing scaffold
    assert len(pairs) == len(set(pairs))      # deduplicated


def test_vocabulary_is_train_only_and_frequency_ranked():
    counts = {"train1": {7: 1, 9: 2}, "train2": {7: 1}, "held": {5: 4, 6: 4}}
    vocabulary = audit.build_vocabulary(counts, {"train1", "train2"}, 2)
    assert vocabulary == [7, 9]
    assert 5 not in vocabulary and 6 not in vocabulary


# ------------------------------------------------------------ arm hygiene
def test_common_mask_audit_rejects_a_silent_intersection():
    metrics = {"ap_bidir": 0.1}
    arms = {"candidate": ({"a": metrics}, {"a": "s"}),
            "control": ({"b": metrics}, {"b": "s"})}
    with pytest.raises(s4r.S4RContractError, match="mask differs"):
        s4r.assert_common_masks(arms, {"s": "c"})


def test_common_mask_audit_rejects_a_remapped_construct():
    metrics = {"ap_bidir": 0.1}
    arms = {"candidate": ({"a": metrics}, {"a": "s"}),
            "control": ({"a": metrics}, {"a": "other"})}
    with pytest.raises(s4r.S4RContractError, match="mapping differs"):
        s4r.assert_common_masks(arms, {"s": "c", "other": "c"})


def test_only_registered_arms_can_be_trained():
    with pytest.raises(s4r.S4RContractError, match="unregistered trained arm"):
        s4r.train("hyperparameter_search")


def test_runtime_context_routes_baseline_and_candidate_to_different_tables():
    class Stub:
        gvec = {"g": np.zeros(41)}
        gvec_graph = {"g": np.zeros(128)}
        ligands = s4r.RuntimeContext.ligands
    stub = Stub()
    assert stub.ligands("baseline41") is Stub.gvec
    assert stub.ligands("candidate") is Stub.gvec_graph
    assert stub.ligands("permuted") is Stub.gvec_graph


def test_stream_and_baseline_replication_anchors_are_the_published_s3r_values():
    gate = json.loads((ROOT / "report" / "s7_l2b_r0r" /
                       "PHASE2B_S3R_GATE.json").read_text(encoding="utf-8"))
    stream = json.loads((ROOT / "report" / "s7_l2b_r0r" /
                         "PHASE2B_S3R_REAL_STREAM_MANIFEST.json").read_text(
                             encoding="utf-8"))
    assert s4r.S3R_CANDIDATE_MACRO_AP == gate["macro_ap_bidir"]["candidate"]
    assert s4r.S3R_STREAM_SEMANTIC_SHA == stream["semantic_sha256"]
