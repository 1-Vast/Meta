import csv
import importlib.util
import json
from pathlib import Path

import pytest

from dataset.splitseal import (AccessLedger, SealedDataset, scan_raw_metadata,
                               seal_source_metaval)


def _raw_table(path):
    rows = [
        ["drug_id", "smiles", "target_id", "sequence", "affinity"],
        ["d1", "CC", "source-target", "AAAA", "1.25"],
        ["d2", "CCC", "metaval-target", "BBBB", "2.50"],
        ["SECRET_RECIPIENT_DRUG", "SECRET_RECIPIENT_SMILES", "recipient-target", "CCCC", "SECRET_RECIPIENT_LABEL"],
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)


def _sealed(tmp_path):
    raw = tmp_path / "raw.tsv"
    out = tmp_path / "sealed"
    _raw_table(raw)
    mapping = {
        "source-target": "source",
        "metaval-target": "metaval",
        "recipient-target": "recipient",
    }
    seal_source_metaval(raw, out, "FIXTURE", mapping)
    return out, mapping


def test_seal_never_emits_recipient_labels(tmp_path):
    out, mapping = _sealed(tmp_path)
    source_labels = (out / "source.tsv").read_text(encoding="utf-8")
    metaval_labels = (out / "metaval.tsv").read_text(encoding="utf-8")
    governance = (out / "governance.tsv").read_text(encoding="utf-8")
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))

    assert "SECRET_RECIPIENT_LABEL" not in source_labels
    assert "SECRET_RECIPIENT_LABEL" not in metaval_labels
    assert "SECRET_RECIPIENT_LABEL" not in governance
    assert "SECRET_RECIPIENT_DRUG" not in governance
    assert "SECRET_RECIPIENT_SMILES" not in governance
    assert "recipient-target" not in source_labels
    assert "recipient-target" not in metaval_labels
    assert "recipient-target" in governance
    assert governance.splitlines()[0].split("\t") == [
        "target_id", "sequence", "sequence_sha256", "split",
    ]
    assert "\ty\t" not in f"\t{governance.splitlines()[0]}\t"
    assert manifest["recipient_label_columns_emitted"] == 0
    assert manifest["rows"]["recipient"] == 1

    sealed = SealedDataset(out, "FIXTURE", mapping)
    assert len(sealed.frame("source", "training")) == 1
    assert len(sealed.frame("source", "hyperparameter")) == 1
    metaval = SealedDataset(out, "FIXTURE", mapping, visible_splits=("metaval",))
    assert len(metaval.frame("metaval", "validation")) == 1
    assert sealed.audit_snapshot() == {
        "recipient_label_reads": 0,
        "recipient_training_rows": 0,
        "recipient_validation_rows": 0,
        "recipient_hyperparameter_access": 0,
        "permitted_row_ids_sha256": manifest["label_row_ids_sha256"]["source"],
        "episode_inputs_sha256": AccessLedger("").snapshot()["episode_inputs_sha256"],
        "episode_count": 0,
    }


def test_split_metadata_scan_never_loads_label_column(tmp_path):
    raw = tmp_path / "raw.tsv"
    _raw_table(raw)
    metadata = scan_raw_metadata(raw)
    assert list(metadata.columns) == ["drug_id", "smiles", "target_id", "sequence"]
    assert "SECRET_RECIPIENT_LABEL" not in metadata.to_string()


@pytest.mark.parametrize("purpose", ["training", "validation", "hyperparameter"])
def test_recipient_access_fails_closed(tmp_path, purpose):
    out, mapping = _sealed(tmp_path)
    sealed = SealedDataset(out, "FIXTURE", mapping)
    with pytest.raises(PermissionError, match="recipient labels"):
        sealed.frame("recipient", purpose)


def test_metaval_cannot_be_used_for_training_or_tuning(tmp_path):
    out, mapping = _sealed(tmp_path)
    sealed = SealedDataset(out, "FIXTURE", mapping, visible_splits=("metaval",))
    with pytest.raises(PermissionError, match="source-only"):
        sealed.frame("metaval", "training")
    with pytest.raises(PermissionError, match="source-only"):
        sealed.frame("metaval", "hyperparameter")


def test_split_map_and_artifact_tampering_are_detected(tmp_path):
    out, mapping = _sealed(tmp_path)
    wrong = dict(mapping)
    wrong["metaval-target"] = "source"
    with pytest.raises(ValueError, match="different split map"):
        SealedDataset(out, "FIXTURE", wrong)

    with (out / "source.tsv").open("a", encoding="utf-8") as handle:
        handle.write("tampered\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        SealedDataset(out, "FIXTURE", mapping)


def test_source_mount_does_not_open_or_validate_metaval_labels(tmp_path):
    out, mapping = _sealed(tmp_path)
    with (out / "metaval.tsv").open("a", encoding="utf-8") as handle:
        handle.write("tampered-metaval-label\n")

    source = SealedDataset(out, "FIXTURE", mapping, visible_splits=("source",))
    assert set(source.rows["split"]) == {"source"}
    with pytest.raises(ValueError, match="hash mismatch"):
        SealedDataset(out, "FIXTURE", mapping, visible_splits=("metaval",))


def test_runtime_cannot_mount_source_and_metaval_labels_together(tmp_path):
    out, mapping = _sealed(tmp_path)
    with pytest.raises(ValueError, match="exactly one"):
        SealedDataset(out, "FIXTURE", mapping, visible_splits=("source", "metaval"))


def test_episode_input_hash_is_deterministic_and_order_sensitive():
    first = AccessLedger("rows")
    second = AccessLedger("rows")
    for ledger in (first, second):
        ledger.record_episode_inputs(
            "training", "t1", ["q1", "q2"], [["s1"], ["s2"]])
    assert (first.snapshot()["episode_inputs_sha256"]
            == second.snapshot()["episode_inputs_sha256"])

    second.record_episode_inputs("training", "t1", ["q3"], [["s1"]])
    assert (first.snapshot()["episode_inputs_sha256"]
            != second.snapshot()["episode_inputs_sha256"])


def test_external_split_map_must_close_every_homology_component(tmp_path, monkeypatch):
    path = Path(__file__).resolve().parents[1] / "test" / "sealdata.py"
    spec = importlib.util.spec_from_file_location("phase_z_sealdata", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    raw = tmp_path / "raw.tsv"
    _raw_table(raw)
    mapping = {
        "source-target": "source",
        "metaval-target": "metaval",
        "recipient-target": "recipient",
    }

    def straddling(_sequences, identity, verbose):
        assert identity > 0 and verbose is False
        return {
            "source-target": 0,
            "metaval-target": 0,
            "recipient-target": 1,
        }, {"n_clusters": 2}

    monkeypatch.setattr(module, "cluster_targets", straddling)
    with pytest.raises(ValueError, match="straddles 1 homology components"):
        module.audit_external_split(raw, mapping)
