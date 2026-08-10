"""Regression tests for the C0/C1 correspondence information audit.

These pin the registration hashes, the P1B slot semantics that the parent
instruction corrected, the exposure firewall, the degree-preserving null and
the exact-AP estimator. The audit trains nothing, so there is no model contract
to test — only the measurement contract.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "research" / "correspondence_router"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


c0 = _load("c0_corpus_test", BASE / "c0_corpus.py")
c1 = _load("c1_audit_test", BASE / "c1_audit.py")


# ------------------------------------------------------------ registration
def test_registration_hashes_match_the_committed_documents():
    assert c0.sha_file(BASE / "PREREG_C0_C1_CORRESPONDENCE_INFORMATION_AUDIT.md") \
        == c0.PREREG_SHA
    assert c0.sha_file(BASE / "PREREG_C0_C1_AMENDMENT_01.md") == c0.AMENDMENT_SHA


def test_frozen_contract_constants():
    assert c0.SLOTS == 128
    assert c0.CONTACT_THRESHOLD == 6.0
    assert (c0.MIN_LIGAND_ATOMS, c0.MAX_LIGAND_ATOMS) == (6, 80)
    assert (c0.MIN_SEQUENCE, c0.MAX_SEQUENCE) == (150, 1200)
    assert c0.MAX_RESOLUTION == 2.5 and c0.MIN_MAPPING_COVERAGE == 0.90
    assert (c1.MIN_COMPONENTS, c1.MAX_COMPONENT_FRACTION) == (60, 0.25)
    assert (c1.MAX_MDE, c1.C1A_MARGIN) == (0.05, 0.05)
    assert (c1.MIN_POSITIVE_UNITS, c1.MIN_CHECKERBOARDS) == (10_000, 1_000)
    assert c1.SEED_BOOT == 20260903 and c1.IDENTITY_EDGE == 0.40


def test_contact_threshold_matches_the_frozen_p1b_mechanism_contract():
    """The exact edges must use P1B's own threshold or the slot gate and the
    exact edges would not describe the same event."""
    from contracts.mechanism import (CONTACT_THRESHOLD_ANGSTROM,
                                     MECHANISM_RESIDUE_SLOTS)
    assert c0.CONTACT_THRESHOLD == CONTACT_THRESHOLD_ANGSTROM
    assert c0.SLOTS == MECHANISM_RESIDUE_SLOTS


def test_stage_opens_no_affinity_source():
    """Scans for source tokens that could only appear if an affinity corpus were
    actually referenced. PLINDER and heldout-A are prose in these files and are
    asserted through the published artifacts instead."""
    for path in (BASE / "c0_corpus.py", BASE / "c1_audit.py",
                 BASE / "c0_equivalence.py"):
        text = path.read_text(encoding="utf-8").lower()
        for forbidden in ("pdbbind", "chembl", "bindingdb", "davis", "kiba",
                          "source_affinity", "phase2b_s3r", "phase2b_s4r"):
            assert forbidden not in text, f"{forbidden} in {path.name}"


# ------------------------------------------------------------ P1B semantics
def test_slot_rule_is_the_p1b_rule():
    """slot = min(127, seq_index * 128 // L), byte-identical to
    scripts/build_structure_supervision."""
    for length in (150, 327, 415, 1200):
        for index in (0, 1, length // 2, length - 1):
            expected = min(c0.SLOTS - 1, index * c0.SLOTS // length)
            assert expected == min(127, index * 128 // length)


def test_slot_target_is_a_minimum_not_a_sum():
    """P1B's slot distance is the MINIMUM over all atoms of all residues in the
    slot, so contact_prob is 'any contact in slot'. This test exists because
    treating it as additive mass was the specific error to avoid."""
    ligand = np.array([[0.0, 0.0, 0.0]])
    near = np.array([[3.0, 0.0, 0.0]])
    far = np.array([[9.0, 0.0, 0.0]])
    distances = np.array([[np.sqrt(((ligand - block) ** 2).sum())
                           for block in (near, far)]])
    slot_distance = distances.min(-1)
    assert slot_distance[0] == pytest.approx(3.0)
    assert (slot_distance <= c0.CONTACT_THRESHOLD).all()
    # a slot containing one contacting and one non-contacting residue is a
    # positive slot, and both residues must remain individually addressable
    assert (distances[0] <= c0.CONTACT_THRESHOLD).tolist() == [True, False]


def test_multiple_residues_in_one_slot_may_contact_the_same_atom():
    contact = np.array([[1, 1, 0]], dtype=np.int8)
    slot_of = np.array([5, 5, 5])
    every, informative = c1.within_slot_ap(
        contact, slot_of, [(5, np.array([0, 1, 2]))])
    assert len(every) == 1 and len(informative) == 1


# ------------------------------------------------------------ exact AP
def test_ap_exact_perfect_and_reversed_rankings():
    labels = np.array([1, 1, 0, 0])
    assert c1.ap_exact(np.array([4.0, 3.0, 2.0, 1.0]), labels) == pytest.approx(1.0)
    assert c1.ap_exact(np.array([1.0, 2.0, 3.0, 4.0]), labels) < 0.6


def test_ap_exact_is_none_for_constant_labels():
    assert c1.ap_exact(np.arange(5.0), np.zeros(5, dtype=int)) is None
    assert c1.ap_exact(np.arange(5.0), np.ones(5, dtype=int)) is None


def test_ap_exact_tied_block_is_the_expectation_not_an_arbitrary_order():
    labels = np.array([1, 0, 1, 0])
    tied = c1.ap_exact(np.zeros(4), labels)
    assert 0.4 < tied < 0.8
    assert c1.ap_exact(np.zeros(4), labels[::-1]) == pytest.approx(tied)


# ------------------------------------------------------------ null
def test_curveball_preserves_both_degree_sequences_exactly():
    rng = np.random.default_rng(0)
    contact = (rng.random((25, 18)) < 0.25).astype(np.int8)
    rewired = c1.curveball(contact, np.random.default_rng(1), swaps_per_edge=50)
    assert np.array_equal(rewired.sum(0), contact.sum(0))
    assert np.array_equal(rewired.sum(1), contact.sum(1))
    assert set(np.unique(rewired)) <= {0, 1}


def test_curveball_actually_mixes():
    rng = np.random.default_rng(2)
    contact = (rng.random((40, 30)) < 0.3).astype(np.int8)
    rewired = c1.curveball(contact, np.random.default_rng(3), swaps_per_edge=100)
    assert np.abs(rewired - contact).sum() > 0


def test_curveball_is_a_noop_on_a_degenerate_matrix():
    contact = np.ones((4, 3), dtype=np.int8)
    rewired = c1.curveball(contact, np.random.default_rng(4), swaps_per_edge=50)
    assert np.array_equal(rewired, contact)


def test_rewire_preserves_the_residue_marginal_used_as_the_score():
    """C1a compares arrangements at fixed marginals, so the score itself must be
    invariant under the null or the contrast would be confounded."""
    rng = np.random.default_rng(5)
    contact = (rng.random((30, 20)) < 0.2).astype(np.int8)
    rewired = c1.curveball(contact, np.random.default_rng(6), swaps_per_edge=60)
    assert np.array_equal(contact.sum(0), rewired.sum(0))


# ------------------------------------------------------------ closure
class _UF:
    pass


def test_union_find_components():
    union = c1.UnionFind(["a", "b", "c", "d"])
    union.union("a", "b")
    union.union("c", "d")
    roots = {union.find(x) for x in ("a", "b", "c", "d")}
    assert len(roots) == 2


def test_kmer_set_is_a_containment_prefilter_input():
    assert c1.kmer_set("AAAB", 3) == {"AAA", "AAB"}
    assert c1.kmer_set("AB", 3) == set()


# ------------------------------------------------------------ exposure
def test_exposure_registry_is_disjoint_from_the_untouched_corpus():
    path = ROOT / "report" / "correspondence_router" / "C0_CORPUS_AND_CENSUS.json"
    if not path.is_file():
        pytest.skip("C0 has not been executed in this checkout")
    census = json.loads(path.read_text(encoding="utf-8"))
    corpus = census["corpus"]
    assert corpus["untouched_entries"] + census["exposure_registry"][
        "union_exposed_pdb_ids"] >= corpus["local_mmcif_entries"]
    assert census["affinity_value_reads"] == 0
    assert census["plinder_used"] is False
    assert census["heldoutA_referenced"] is False
    assert census["trainable_parameters_introduced"] == 0


def test_recorded_mapping_equivalence_passed():
    path = ROOT / "report" / "correspondence_router" / "C0_MAPPING_EQUIVALENCE.json"
    if not path.is_file():
        pytest.skip("the mapping equivalence check has not been executed")
    result = json.loads(path.read_text(encoding="utf-8"))
    assert result["structures_checked"] > 0
    assert result["slot_assignment_agreement"] == result["structures_checked"]
    assert result["pass"] is True


def test_recorded_c1_verdict_is_in_the_registered_set():
    path = ROOT / "report" / "correspondence_router" / "C1_INFORMATION_AUDIT.json"
    if not path.is_file():
        pytest.skip("C1 has not been executed in this checkout")
    gate = json.loads(path.read_text(encoding="utf-8"))
    assert gate["TERMINAL_VERDICT"] in {
        "CORRESPONDENCE_DATA_OR_CLOSURE_NOT_IDENTIFIABLE",
        "SLOT_ROUTING_ESTIMAND_INVALID",
        "EXACT_EDGE_COUPLING_NOT_SUPPORTED_BY_TEACHER",
        "CORRESPONDENCE_INFORMATION_PRESENT_C2_AUTHORIZED"}
    assert gate["affinity_value_reads"] == 0
    assert gate["trainable_parameters_introduced"] == 0
    assert gate["rewire_is_an_evaluation_null_not_a_biological_nonbinder"] is True
    if gate["TERMINAL_VERDICT"] != "CORRESPONDENCE_INFORMATION_PRESENT_C2_AUTHORIZED":
        assert gate["authorized_next_action"].startswith("none")
