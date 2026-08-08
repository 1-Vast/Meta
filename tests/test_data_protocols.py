import json
from pathlib import Path

import pytest

from scripts.data_contract import read_jsonl
from scripts.data_protocols import PROTOCOLS, build_protocols


def _write_jsonl(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
        encoding="utf-8",
    )


def _inputs(tmp_path: Path):
    identities = []
    labels = []
    assignments = []
    for protein, component in (("p1", "h1"), ("p2", "h1"), ("p3", "h2")):
        assignments.append({
            "protein_sequence_sha256": protein,
            "homology_component_id": component,
        })
        for ligand_index in range(18):
            row_id = f"{protein}-{ligand_index}"
            identities.append({
                "row_id": row_id,
                "protein_sequence_sha256": protein,
                "ligand_connectivity_key": f"ligand-{ligand_index}",
                "task_keys": {"protein_assay_context": f"{protein}-Ki-assay"},
                "canonical_smiles": "CCO",
                "murcko_scaffold": f"scaffold-{ligand_index % 3}",
                "pocket_id": f"pocket-{protein}",
            })
            labels.append({"row_id": row_id, "p_value": 6.0 + ligand_index / 10})

    # Exact duplicates may have distinct administrative row identifiers.
    duplicate = dict(identities[0], row_id="duplicate-row")
    identities.append(duplicate)
    labels.append({"row_id": "duplicate-row", "p_value": labels[0]["p_value"]})

    paths = {
        "identity": tmp_path / "identity.jsonl",
        "labels": tmp_path / "labels.jsonl",
        "assignments": tmp_path / "assignments.jsonl",
        "recipient": tmp_path / "recipient.jsonl",
        "annotations": tmp_path / "protein_annotations.jsonl",
    }
    _write_jsonl(paths["identity"], identities)
    _write_jsonl(paths["labels"], labels)
    _write_jsonl(paths["assignments"], assignments)
    _write_jsonl(paths["recipient"], [{"sequence_sha256": "p3"}])
    _write_jsonl(paths["annotations"], [
        {"protein_sequence_sha256": "p1", "split": "obsolete", "go_terms": ["GO:1"]},
        {"protein_sequence_sha256": "p2", "split": "obsolete", "go_terms": ["GO:2"]},
        {"protein_sequence_sha256": "p3", "split": "obsolete", "go_terms": ["GO:3"]},
    ])
    return paths


def _build(tmp_path: Path, name: str = "out"):
    paths = _inputs(tmp_path)
    output = tmp_path / name
    summaries = build_protocols(
        paths["identity"],
        paths["labels"],
        paths["assignments"],
        output,
        recipient_identity_paths=[paths["recipient"]],
        protein_annotations_path=paths["annotations"],
    )
    return output, summaries


def test_dual_protocol_preserves_estimator_boundaries(tmp_path):
    output, summaries = _build(tmp_path)

    for protocol in PROTOCOLS:
        directory = output / protocol
        identities = read_jsonl(directory / "identity.jsonl")
        labels = read_jsonl(directory / "labels.jsonl")
        assignments = read_jsonl(directory / "assignments.jsonl")
        tasks = read_jsonl(directory / "tasks.jsonl")
        annotations = read_jsonl(directory / "annotations.jsonl")
        protein_annotations = read_jsonl(directory / "protein_annotations.jsonl")

        assert len(identities) == len(labels) == 36
        assert {row["protein_sequence_sha256"] for row in identities} == {"p1", "p2"}
        assert all("p_value" not in row for row in identities)
        assert all(set(row) == {"row_id", "p_value"} for row in labels)
        assert all(row["eligible_k5_q12"] for row in tasks)
        assert all(row["support_compounds"] >= 5 for row in tasks)
        assert all(row["query_compounds"] >= 12 for row in tasks)

        for task in tasks:
            rows = [row for row in identities if row["task_id"] == task["task_id"]]
            support = {row["ligand_connectivity_key"] for row in rows if row["episode_role"] == "support_pool"}
            query = {row["ligand_connectivity_key"] for row in rows if row["episode_role"] == "query_pool"}
            assert support.isdisjoint(query)
            assert {row["split"] for row in rows} == {task["split"]}

        assert {row["row_id"] for row in annotations} == {row["row_id"] for row in identities}
        assert all("murcko_scaffold" in row and "pocket_id" in row for row in annotations)
        rebound = {row["protein_sequence_sha256"]: row for row in protein_annotations}
        assignment = {row["protein_sequence_sha256"]: row for row in assignments}
        assert set(rebound) == {"p1", "p2"}
        assert all(rebound[p]["split"] == assignment[p]["split"] for p in rebound)
        assert all(rebound[p]["protocol"] == protocol for p in rebound)

        manifest = summaries[protocol]
        assert manifest["recipient_rows_excluded"] == 18
        assert manifest["recipient_registry_exact_targets"] == 1
        assert manifest["recipient_corpus_exact_targets_excluded"] == 1
        assert manifest["exact_duplicate_rows_removed"] == 1
        assert manifest["recipient_label_reads"] == 0
        assert manifest["homology_components_deleted"] == 0
        assert manifest["scaffolds_deleted"] == 0

    stress = read_jsonl(output / "novel_family_stress_v1" / "assignments.jsonl")
    assert len({row["split"] for row in stress}) == 1


def test_output_is_immutable_and_build_is_deterministic(tmp_path):
    output, first = _build(tmp_path, "first")
    paths = _inputs(tmp_path)
    with pytest.raises(FileExistsError):
        build_protocols(paths["identity"], paths["labels"], paths["assignments"], output)

    second_output = tmp_path / "second"
    second = build_protocols(
        paths["identity"], paths["labels"], paths["assignments"], second_output,
        [paths["recipient"]], paths["annotations"],
    )
    for protocol in PROTOCOLS:
        assert first[protocol]["artifacts"] == second[protocol]["artifacts"]


def test_conflicting_row_identity_fails_closed(tmp_path):
    paths = _inputs(tmp_path)
    identities = read_jsonl(paths["identity"])
    identities.append(dict(identities[0], canonical_smiles="CCC"))
    _write_jsonl(paths["identity"], identities)
    with pytest.raises(ValueError, match="conflicting identities"):
        build_protocols(
            paths["identity"], paths["labels"], paths["assignments"], tmp_path / "out"
        )


def test_conflicting_replicate_labels_are_registered_and_retained(tmp_path):
    paths = _inputs(tmp_path)
    identities = read_jsonl(paths["identity"])
    labels = read_jsonl(paths["labels"])
    identities.append(dict(identities[0], row_id="conflicting-replicate"))
    labels.append({"row_id": "conflicting-replicate", "p_value": 9.0})
    _write_jsonl(paths["identity"], identities)
    _write_jsonl(paths["labels"], labels)
    output = tmp_path / "out"
    summaries = build_protocols(
        paths["identity"], paths["labels"], paths["assignments"], output,
        [paths["recipient"]], paths["annotations"],
    )
    for protocol in PROTOCOLS:
        conflicts = read_jsonl(output / protocol / "replicate_conflicts.jsonl")
        assert len(conflicts) == 1
        assert conflicts[0]["p_values"] == [6.0, 9.0]
        assert summaries[protocol]["conflicting_replicate_groups_retained"] == 1
        assert summaries[protocol]["rows"] == 37


def test_label_blind_governance_tsv_can_seal_recipient_targets(tmp_path):
    paths = _inputs(tmp_path)
    governance = tmp_path / "governance.tsv"
    governance.write_text(
        "target_id\tsequence\tsequence_sha256\tsplit\n"
        "recipient\tSEQ\tp3\trecipient\n",
        encoding="utf-8",
    )
    summaries = build_protocols(
        paths["identity"], paths["labels"], paths["assignments"], tmp_path / "out",
        [governance], paths["annotations"],
    )
    assert all(value["recipient_rows_excluded"] == 18 for value in summaries.values())
