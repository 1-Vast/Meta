import numpy as np

from research.crossed_interaction.train_kmer_cq_observable import (
    materialize_features,
    protein_kmer_descriptor,
)


def test_protein_kmer_descriptor_is_normalized():
    descriptor = protein_kmer_descriptor("ACDACD", k=3, bins=16)
    assert descriptor.shape == (16,)
    assert np.isfinite(descriptor).all()
    assert np.isclose(descriptor.sum(), 1.0)


def test_protein_kmer_descriptor_captures_local_motif_order():
    first = protein_kmer_descriptor("ACDACD", k=3, bins=64)
    second = protein_kmer_descriptor("AAACCC", k=3, bins=64)
    assert not np.array_equal(first, second)


def test_short_sequence_uses_single_stable_bin():
    descriptor = protein_kmer_descriptor("AC", k=3, bins=16)
    assert np.isclose(descriptor.sum(), 1.0)
    assert np.count_nonzero(descriptor) == 1


def test_materialize_features_reports_expected_dimensions(tmp_path):
    corpus = tmp_path
    (corpus / "manifest.json").write_text("{}", encoding="utf-8")
    (corpus / "proteins.jsonl").write_text(
        '{"sequence_sha256":"t1","sequence":"ACDACD"}\n'
        '{"sequence_sha256":"t2","sequence":"WYWYWY"}\n',
        encoding="utf-8")
    (corpus / "ligands.jsonl").write_text(
        '{"drug_key":"l1","smiles":"CCO"}\n'
        '{"drug_key":"l2","smiles":"CCN"}\n',
        encoding="utf-8")

    import gzip
    import json

    cells = [
        {
            "cell_id": "c1", "target_id": "t1", "ligand_id": "l1",
            "protein_group_40": "g1", "scaffold": "s1",
        },
        {
            "cell_id": "c2", "target_id": "t2", "ligand_id": "l2",
            "protein_group_40": "g2", "scaffold": "s2",
        },
    ]
    with gzip.open(corpus / "cells.jsonl.gz", "wt", encoding="utf-8") as handle:
        for cell in cells:
            handle.write(json.dumps(cell) + "\n")

    features, metadata = materialize_features(
        corpus, kmer_size=3, kmer_bins=8, ligand_fp_bits=16)
    assert metadata["protein_descriptor_dim"] == 8
    assert metadata["ligand_descriptor_dim"] == 24
    assert metadata["feature_dim"] == 192
    assert features["c1"]["correct"].shape == (192,)
