import json
from dataclasses import replace

import numpy as np
import pytest
import torch

from model.config import DEFAULT
from scripts.data import CompiledEpisodes
from scripts.data_contract import digest, read_jsonl, write_jsonl
from test.audit_sealed_dataset import audit_sealed_dataset
from scripts.build_ligand_bank import build_ligand_bank
from scripts.build_protein_bank import migrate_protein_cache
from scripts.seal_compiled_dataset import SealedCompiledDataset, seal_compiled_dataset


def _compiled_rows(tmp_path):
    compiled = tmp_path / "compiled"
    rows = []
    for split, target in (("train", "source-a"), ("train", "source-b"),
                          ("val", "metaval-a"), ("test", "recipient-sentinel")):
        for drug in range(7):
            rows.append({
                "row_id": f"{target}-{drug}", "drug_key": f"{target}-drug-{drug}",
                "target_key": target, "pair_key": f"{target}-pair-{drug}",
                "task_key": f"{target}-Ki", "smiles": "CCO", "sequence": "ACDEFGHIK",
                "gamma_key": "Ki", "gamma": 0.5, "split": split, "y": drug / 10,
                "provenance": [str(drug)],
            })
    write_jsonl(compiled / "rows.jsonl", rows)
    (compiled / "manifest.json").write_text(json.dumps({
        "dataset": "synthetic", "normalization": {"low": 0.0, "high": 10.0},
    }), encoding="utf-8")
    return compiled


def _cache(path, keys, dim=4):
    torch.save({key: (torch.ones(dim), torch.ones(2, dim)) for key in keys}, path)


def test_compiled_seal_physically_excludes_recipient_labels(tmp_path):
    sealed = tmp_path / "sealed"
    manifest = seal_compiled_dataset(_compiled_rows(tmp_path), sealed)

    source = (sealed / "source" / "rows.jsonl").read_text(encoding="utf-8")
    metaval = (sealed / "metaval" / "rows.jsonl").read_text(encoding="utf-8")
    assert "recipient-sentinel" not in source
    assert "recipient-sentinel" not in metaval
    assert not (sealed / "recipient").exists()
    assert manifest["recipient_label_artifact_emitted"] is False

    view = SealedCompiledDataset(sealed, "source", "synthetic")
    assert {row["split"] for row in view.rows} == {"source"}
    assert view.audit_snapshot()["recipient_label_reads"] == 0
    assert audit_sealed_dataset(sealed)["valid"]


def test_seal_permits_one_pair_in_distinct_assay_tasks(tmp_path):
    compiled = _compiled_rows(tmp_path)
    rows = read_jsonl(compiled / "rows.jsonl")
    duplicate = dict(rows[0])
    duplicate["row_id"] = "source-a-kd-repeat"
    duplicate["task_key"] = "source-a-Kd"
    duplicate["gamma_key"] = "Kd"
    duplicate["gamma"] = 0.25
    rows.append(duplicate)
    write_jsonl(compiled / "rows.jsonl", rows)

    manifest = seal_compiled_dataset(compiled, tmp_path / "sealed")
    assert manifest["counts"]["source"]["rows"] == 15


def test_legacy_cache_migration_emits_only_the_requested_view_targets(tmp_path):
    rows = tmp_path / "source.jsonl"
    write_jsonl(rows, [{"target_key": digest("ACD"), "sequence": "ACD"},
                       {"target_key": digest("EFG"), "sequence": "EFG"}])
    index = tmp_path / "index.tsv"
    index.write_text("target\tsequence\nlegacy-a\tACD\nlegacy-b\tEFG\n", encoding="utf-8")
    legacy = tmp_path / "legacy.pt"
    torch.save({"legacy-a": (torch.ones(2), torch.ones(1, 2)),
                "legacy-b": (torch.ones(2), torch.ones(1, 2))}, legacy)
    output = tmp_path / "view.pt"

    result = migrate_protein_cache(rows, legacy, index, output,
                                   target_column="target", sequence_column="sequence")
    assert result["targets"] == 2
    assert set(torch.load(output, weights_only=False)) == {digest("ACD"), digest("EFG")}


@pytest.mark.skipif(not torch.cuda.is_available(), reason="sealed runtime requires CUDA")
def test_runtime_mounts_one_view_and_rejects_cross_view_protein_cache(tmp_path):
    sealed = tmp_path / "sealed"
    seal_compiled_dataset(_compiled_rows(tmp_path), sealed)
    build_ligand_bank(sealed / "source" / "rows.jsonl", sealed / "source_ligand_bank")
    build_ligand_bank(sealed / "metaval" / "rows.jsonl", sealed / "metaval_ligand_bank")
    source_cache = tmp_path / "source.pt"
    _cache(source_cache, ["source-a", "source-b"])

    source = CompiledEpisodes(
        "synthetic", sealed_dir=sealed, protein_cache_path=source_cache, device="cuda",
    )
    episode, target = source.sample(8, np.random.default_rng(3))
    assert not hasattr(episode, "query_y") and not hasattr(episode, "Y")
    assert tuple(target.query_y.shape) == (8,)
    for row in range(len(episode)):
        valid = episode.support_mask[row].bool()
        assert not bool((episode.support_pair_idx[row, valid]
                         == episode.query_pair_idx[row]).any())
    with pytest.raises(PermissionError, match="not mounted"):
        source.sample(1, np.random.default_rng(4), "metaval")

    leaked_cache = tmp_path / "leaked.pt"
    _cache(leaked_cache, ["source-a", "source-b", "metaval-a"])
    with pytest.raises(ValueError, match="exactly match"):
        CompiledEpisodes(
            "synthetic", sealed_dir=sealed, protein_cache_path=leaked_cache, device="cuda",
        )

    metaval_cache = tmp_path / "metaval.pt"
    _cache(metaval_cache, ["metaval-a"])
    metaval = CompiledEpisodes(
        "synthetic", sealed_dir=sealed, protein_cache_path=metaval_cache, device="cuda",
        visible_splits=("metaval",), normalization_bounds=(0.0, 10.0),
    )
    assert metaval.fixed_support_tasks("metaval-a-Ki", 2, np.random.default_rng(7)) is not None
    with pytest.raises(PermissionError, match="not mounted"):
        metaval.sample(1, np.random.default_rng(8), "source")
