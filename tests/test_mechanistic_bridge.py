import gzip
import json

import numpy as np
import pytest
import torch

from model.geometry_supervision import PairGeometryTeacher
from scripts.build_holo_complex_index import (_protein_sequence_mapping,
    build_holo_complex_index)
from scripts.build_protein_bank import residue_slot_mapping
from scripts.cache_structure_proteins import _record_sequence_key, _sequence_chunks
from scripts.build_structure_supervision import (audit_plinder_mlsb,
    build_structure_supervision)
from scripts.data_contract import read_jsonl
from scripts.evaluate_pair_geometry import _metrics, _paired_bootstrap


def test_residue_slot_mapping_covers_sequence_without_overlap():
    mapping = residue_slot_mapping(10, 4)
    assert mapping["residue_slot_start"].tolist() == [0, 2, 5, 7]
    assert mapping["residue_slot_end"].tolist() == [2, 5, 7, 10]
    assert mapping["residue_mask"].tolist() == [True, True, True, True]


def test_esm_long_sequences_are_chunked_without_truncation():
    sequence = "A" * 2050
    chunks = _sequence_chunks(sequence)
    assert list(map(len, chunks)) == [1022, 1022, 6]
    assert "".join(chunks) == sequence
    assert _record_sequence_key({"target_key": "target"}) == "target"


def test_mmcif_residues_align_to_canonical_sequence_indices():
    names = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS"]
    rows = [{"label_seq_id": str(index), "label_comp_id": name}
            for index, name in enumerate(names, start=1)]
    labels, indices, coverage = _protein_sequence_mapping(rows, "ACDMEFGHIK")
    assert labels == [str(index) for index in range(1, 10)]
    assert indices == [0, 1, 2, 4, 5, 6, 7, 8, 9]
    assert coverage == pytest.approx(0.9)


def test_structure_sidecar_compiles_canonical_mmcif_and_ccd(tmp_path):
    pytest.importorskip("gemmi")
    pytest.importorskip("rdkit")
    raw = tmp_path / "raw"
    (raw / "mmcif").mkdir(parents=True)
    (raw / "ccd").mkdir()
    (raw / "biolip2").mkdir()
    ccd = raw / "ccd" / "CYC.cif"
    ccd.write_text("""data_CYC
_chem_comp.type 'NON-POLYMER'
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
C1 C 0
C2 C 0
C3 C 0
C4 C 0
C5 C 0
C6 C 0
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
C1 C2 SING
C2 C3 SING
C3 C4 SING
C4 C5 SING
C5 C6 SING
C6 C1 SING
""", encoding="utf-8")
    tags = ["group_PDB", "label_atom_id", "label_comp_id", "label_asym_id",
            "label_seq_id", "auth_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y",
            "Cartn_z", "type_symbol", "label_alt_id", "occupancy", "pdbx_PDB_model_num"]
    rows = []
    for index in range(1, 51):
        rows.append(f"ATOM CA ALA A {index} A {index} {index * 3:.1f} 0 0 C . 1.0 1")
    for index in range(1, 7):
        rows.append(f"HETATM C{index} CYC B . L 101 {3 + index * .2:.1f} 1 0 C . 1.0 1")
    cif_text = "data_test\n_exptl.method 'X-RAY DIFFRACTION'\n_refine.ls_d_res_high 2.0\nloop_\n"
    cif_text += "\n".join(f"_atom_site.{tag}" for tag in tags) + "\n"
    cif_text += "\n".join(rows) + "\n"
    mmcif = raw / "mmcif" / "test.cif.gz"
    with gzip.open(mmcif, "wt", encoding="utf-8") as handle:
        handle.write(cif_text)
    columns = ["test", "A", "2.0", "BS01", "CYC", "L", "1", "", "", "",
               "", "", "", "", "", "", "", "", "", "101", "A" * 50]
    annotation = raw / "biolip2" / "BioLiP.txt.gz"
    with gzip.open(annotation, "wt", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
    with gzip.open(raw / "biolip2" / "ligand.tsv.gz", "wt", encoding="utf-8") as handle:
        handle.write("CYC\tC6H12\t-\t-\tC1CCCCC1\tcyclohexane\t-\t-\t-\n")
    index_dir = tmp_path / "index"
    index = build_holo_complex_index(annotation, raw, index_dir, candidate_limit=1)
    assert index["valid_holo_complexes"] == 1
    output = tmp_path / "sidecar"
    manifest = build_structure_supervision(index_dir / "complexes.jsonl", output)
    pair = read_jsonl(output / "pairs.jsonl")[0]
    with np.load(output / pair["shard"], allow_pickle=False) as shard:
        assert shard["contact"].shape == (1, 128, 128)
        assert shard["contact"][0, 0, 0] == 1.0
        assert shard["distance_bin"][0, 0, 0] == 0
        assert shard["contact"].dtype == np.uint8
        assert shard["distance_bin"].dtype == np.uint8
    assert manifest["pairs"] == 1
    assert manifest["residue_slots"] == 128
    assert pair["atom_mapping_hash"] and pair["residue_mapping_hash"]


def test_local_plinder_mlsb_is_rejected_as_holo_supervision():
    root = "dataset/raw/plinder/plinder_2024-06_v2_metadata/splits/mlsb/inputs"
    result = audit_plinder_mlsb(root)
    assert result["systems"] == 346
    assert result["ligand_coordinate_files"] == 0
    assert result["admissible_for_geometry_supervision"] is False


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA bridge validation")
def test_pair_geometry_teacher_shapes_masks_and_gradients_on_cuda():
    bridge = PairGeometryTeacher(
        hidden_dim=8, section_dim=4, pair_dim=8, blocks=1,
        latents=4, heads=2, chunk_size=8).cuda()
    atoms = torch.randn(2, 5, 8, device="cuda", requires_grad=True)
    residues = torch.randn(2, 7, 8, device="cuda", requires_grad=True)
    atom_mask = torch.tensor([[1, 1, 1, 0, 0], [1, 1, 1, 1, 1]], device="cuda")
    residue_mask = torch.tensor([[1, 1, 0, 0, 0, 0, 0], [1] * 7], device="cuda")
    adjacency = torch.eye(5, device="cuda").expand(2, -1, -1)
    output = bridge(atoms, atom_mask, residues, residue_mask, adjacency)
    assert output.contact_logits.shape == (2, 5, 7)
    assert output.distance_logits.shape == (2, 5, 7, 5)
    assert output.contact_prob[0, 3:].count_nonzero() == 0
    (output.contact_logits.mean() + output.distance_logits.mean()).backward()
    assert atoms.grad is not None and residues.grad is not None


def test_mechanism_gate_metrics_use_frozen_distance_centers_and_top_l():
    labels = np.array([1, 0, 1, 0], dtype=np.uint8)
    contact_prob = np.array([0.9, 0.8, 0.7, 0.1])
    distance_labels = np.array([0, 1, 2, 4], dtype=np.uint8)
    distance_prob = np.eye(5)[distance_labels]
    result = _metrics(labels, contact_prob, distance_labels, distance_prob,
                      [(labels, contact_prob, 2)])
    assert result["contact_precision_at_top_l"] == pytest.approx(0.5)
    assert result["expected_distance_mae_angstrom"] == pytest.approx(0.0)
    assert result["distance_bin_accuracy"] == pytest.approx(1.0)
    correct = [{"source_entry_id": "a", "contact_auprc": 0.8},
               {"source_entry_id": "b", "contact_auprc": 0.7}]
    control = [{"source_entry_id": "a", "contact_auprc": 0.4},
               {"source_entry_id": "b", "contact_auprc": 0.3}]
    interval = _paired_bootstrap(correct, control, "contact_auprc",
                                 higher_is_better=True, seed=17)
    assert interval["point"] == pytest.approx(0.4)
    assert interval["lower_95"] > 0
