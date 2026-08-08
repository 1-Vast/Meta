import csv
import json

from scripts.preprocess_dataset import audit_dataset, compile_dataset
from test.compiled_dataset_tools import (batch_from_episode_indices, build_episodes,
                                         freeze_population_bands)
from scripts.data_contract import digest, read_jsonl


def test_dataset_contract_compiles_without_dataset_specific_code(tmp_path):
    raw = tmp_path / "raw.csv"
    with raw.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["smiles", "sequence", "affinity", "relation", "endpoint"])
        writer.writeheader()
        for target in range(30):
            sequence = "M" + "ACDEFGHIK" + "A" * target
            for drug in range(6):
                writer.writerow({"smiles": "C" * (drug + 1), "sequence": sequence,
                                 "affinity": 10 + drug, "relation": "=",
                                 "endpoint": "ki" if drug % 2 else "kd"})
        writer.writerow({"smiles": "CC", "sequence": "MACDEFGHIK", "affinity": 11,
                         "relation": "=", "endpoint": "kd"})
        writer.writerow({"smiles": "CCC", "sequence": "MACDEFGHIK", "affinity": 12,
                         "relation": "<", "endpoint": "ki"})
    spec = {"name": "synthetic", "input": {"path": str(raw), "columns": {
        "smiles": "smiles", "sequence": "sequence", "label": "affinity",
        "relation": "relation", "gamma": "endpoint"}},
        "label": {"transform": "identity", "allowed_relations": ["="],
                  "normalization": {"mode": "explicit_affine", "min": 0, "max": 20}},
        "split": {"fractions": {"train": 0.6, "val": 0.15, "test": 0.25}},
        "duplicates": {"label_tolerance": 0.0, "aggregate": "median", "on_conflict": "fail"}}
    spec_path = tmp_path / "spec.json"
    spec_path.write_text(json.dumps(spec), encoding="utf-8")
    compiled = tmp_path / "compiled"
    manifest = compile_dataset(spec_path, compiled)
    assert manifest["schema"] == "MetaSieve.CanonicalDTA.v2"
    assert manifest["biological_context"]["endpoint"]["mapping"]["Ki"] == 1
    assert all("gamma" not in row for row in read_jsonl(compiled / "rows.jsonl"))
    assert manifest["dropped"]["nonpoint_relation"] == 1
    assert audit_dataset(compiled)["valid"]

    episode_path = compiled / "episodes.npz"
    result = build_episodes(compiled / "rows.jsonl", episode_path, k=2, draws_per_task=1)
    assert result["episodes"] > 0
    assert result["contains_query_label"] is False
    batch, query_y = batch_from_episode_indices(compiled / "rows.jsonl", episode_path)
    assert not hasattr(batch, "Y")
    assert len(query_y) == len(batch.query_pair_idx)

    bands = freeze_population_bands(compiled / "rows.jsonl", compiled / "bands.json", minimum_count=100)
    assert {value["status"] for value in bands["context_bands"].values()} == {"vacuous"}


def test_label_blind_governance_fits_normalization_from_declared_source_only(tmp_path):
    raw = tmp_path / "raw.csv"
    raw.write_text(
        "smiles,sequence,affinity\nCCO,ACD,10\nCCN,EFG,20\nCCC,HIK,100\nCCCC,LMN,200\n",
        encoding="utf-8",
    )
    spec = {
        "name": "governed", "input": {"path": str(raw), "columns": {
            "smiles": "smiles", "sequence": "sequence", "label": "affinity"}},
        "label": {"transform": "identity", "normalization": {"mode": "train_minmax"}},
        "split": {"fractions": {"train": 0.6, "val": 0.15, "test": 0.25}},
        "duplicates": {"label_tolerance": 0.0, "aggregate": "median", "on_conflict": "fail"},
    }
    spec_path = tmp_path / "spec.json"
    spec_path.write_text(json.dumps(spec), encoding="utf-8")
    governance = tmp_path / "governance.tsv"
    governance.write_text(
        "sequence_sha256\tsplit\n"
        f"{digest('ACD')}\tsource\n{digest('EFG')}\tsource\n"
        f"{digest('HIK')}\tmetaval\n{digest('LMN')}\trecipient\n",
        encoding="utf-8",
    )
    compiled = tmp_path / "compiled"
    manifest = compile_dataset(spec_path, compiled, governance)
    rows = {row["sequence"]: row for row in read_jsonl(compiled / "rows.jsonl")}

    assert manifest["split_governance"]["mode"] == "external_label_blind"
    assert manifest["normalization"] == {"mode": "train_minmax", "low": 10.0, "high": 20.0}
    assert rows["ACD"]["split"] == "train" and rows["ACD"]["y"] == 0.0
    assert rows["EFG"]["split"] == "train" and rows["EFG"]["y"] == 1.0
    assert rows["HIK"]["split"] == "val"
