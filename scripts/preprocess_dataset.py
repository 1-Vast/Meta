"""Compile a declared raw affinity table into the canonical DTA row contract."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

from scripts.data_contract import (canonical_sequence, canonical_smiles, digest,
    biological_context_registry, duplicate_governance, encode_biological_context,
    fit_normalizer, load_spec, normalize, read_jsonl, target_split,
    transform_label, write_jsonl)
from scripts.seal_compiled_dataset import seal_compiled_dataset


def _raw_rows(spec: dict) -> list[dict[str, str]]:
    source = Path(spec["input"]["path"])
    delimiter = spec["input"].get("delimiter", ",")
    with source.open(encoding=spec["input"].get("encoding", "utf-8"), newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _governance_splits(path: str | Path) -> tuple[dict[str, str], str]:
    """Read only sequence hashes and declared splits from a label-free TSV."""
    governance = Path(path)
    with governance.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sequence_sha256", "split"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError("governance TSV must contain sequence_sha256 and split columns")
        mapping = {}
        for row in reader:
            key, split = str(row["sequence_sha256"]), str(row["split"])
            if split not in {"source", "metaval", "recipient"}:
                raise ValueError(f"invalid governance split {split!r}")
            if key in mapping and mapping[key] != split:
                raise ValueError("governance TSV assigns one sequence to multiple splits")
            mapping[key] = split
    return mapping, sha256(governance.read_bytes()).hexdigest()


def compile_dataset(spec_path: str | Path, output_dir: str | Path,
                    governance_path: str | Path | None = None) -> dict:
    destination = Path(output_dir)
    if destination.exists():
        raise FileExistsError(f"compiled dataset output already exists: {destination}")
    spec = load_spec(spec_path)
    columns = spec["input"]["columns"]
    allowed = set(spec["label"].get("allowed_relations", ["="]))
    raw = _raw_rows(spec)
    context_registry = biological_context_registry(spec)
    has_external_governance = governance_path is not None
    governance, governance_sha256 = ({}, None) if governance_path is None \
        else _governance_splits(governance_path)
    governance_to_canonical = {"source": "train", "metaval": "val", "recipient": "test"}
    dropped: dict[str, int] = {}
    staged = []
    for index, source in enumerate(raw):
        try:
            relation = source.get(columns.get("relation", ""), "=").strip() or "="
            if relation not in allowed:
                raise ValueError("nonpoint_relation")
            smiles = canonical_smiles(source[columns["smiles"]])
            sequence = canonical_sequence(source[columns["sequence"]])
            context = encode_biological_context(source, context_registry)
            target_key, drug_key = digest(sequence), digest(smiles)
            task_key = digest(f"{target_key}|{context['context_key']}")
            pair_key = digest(f"{target_key}|{drug_key}")
            split = target_split(target_key, spec["split"])
            if governance:
                try:
                    split = governance_to_canonical[governance[target_key]]
                except KeyError as error:
                    raise ValueError("governance TSV does not cover every canonical target") from error
            staged.append({"row_id": str(index), "smiles": smiles, "sequence": sequence,
                           "drug_key": drug_key, "target_key": target_key, "pair_key": pair_key,
                           "task_key": task_key, **context,
                           "split": split,
                           "label_value": transform_label(source[columns["label"]], spec["label"]),
                           "provenance": (str(index),)})
        except (KeyError, TypeError, ValueError) as error:
            if has_external_governance and "governance TSV" in str(error):
                raise
            key = str(error) or "invalid_row"
            dropped[key] = dropped.get(key, 0) + 1
    governed = duplicate_governance(staged, spec)
    if has_external_governance:
        present = {row["target_key"] for row in governed}
        unexpected = set(governance) - present
        if unexpected:
            raise ValueError("governance TSV contains targets absent from the raw dataset")
    fitted = fit_normalizer([row["label_value"] for row in governed if row["split"] == "train"],
                            spec["label"].get("normalization", {}))
    canonical = [{key: value for key, value in row.items() if key != "label_value"} |
                 {"y": normalize(row["label_value"], fitted)}
                 for row in governed]
    write_jsonl(destination / "rows.jsonl", canonical)
    manifest = {"schema": "MetaSieve.CanonicalDTA.v2", "dataset": spec["name"],
                "spec_sha256": sha256(Path(spec_path).read_bytes()).hexdigest(),
                "counts": {"raw_rows": len(raw), "canonical_rows": len(canonical),
                           "drugs": len({row["drug_key"] for row in canonical}),
                           "targets": len({row["target_key"] for row in canonical}),
                           "tasks": len({row["task_key"] for row in canonical})},
                "dropped": dropped, "normalization": fitted,
                "biological_context": context_registry,
                "contract": {"value_range": [0.0, 1.0], "target_split_unit": "exact_sequence_sha256"},
                "split_governance": {"mode": "external_label_blind" if governance else "hash",
                                     "sha256": governance_sha256}}
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def audit_dataset(directory: str | Path) -> dict:
    """Fail closed when canonical staging rows violate the row contract."""
    rows = read_jsonl(Path(directory) / "rows.jsonl")
    errors = []
    if not rows:
        errors.append("no canonical rows")
    target_splits: dict[str, set[str]] = {}
    duplicates = set()
    for row in rows:
        if not 0.0 <= float(row["y"]) <= 1.0:
            errors.append(f"row {row['row_id']} has y outside [0,1]")
        target_splits.setdefault(row["target_key"], set()).add(row["split"])
        key = (row["task_key"], row["pair_key"])
        if key in duplicates:
            errors.append(f"unresolved duplicate {key}")
        duplicates.add(key)
    crossed = [key for key, values in target_splits.items() if len(values) != 1]
    if crossed:
        errors.append(f"{len(crossed)} exact targets cross split boundaries")
    return {"valid": not errors, "rows": len(rows), "targets": len(target_splits), "errors": errors}


def build_dataset(spec_path: str | Path, compiled_dir: str | Path, *,
                  governance_path: str | Path | None = None,
                  sealed_dir: str | Path | None = None) -> dict:
    """Compile, audit, and optionally seal one declarative DTA dataset."""
    manifest = compile_dataset(spec_path, compiled_dir, governance_path)
    audit = audit_dataset(compiled_dir)
    if not audit["valid"]:
        raise RuntimeError(f"compiled dataset audit failed: {audit['errors']}")
    result = {"compiled": manifest, "audit": audit}
    if sealed_dir is not None:
        result["sealed"] = seal_compiled_dataset(compiled_dir, sealed_dir)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("spec")
    parser.add_argument("output")
    parser.add_argument("--governance", help="label-free target split TSV")
    parser.add_argument("--sealed-dir", help="optional source/metaval runtime seal output")
    args = parser.parse_args()
    allowed_output = (Path(__file__).resolve().parents[1] / "dataset/processed").resolve()
    outputs = [Path(args.output)]
    if args.sealed_dir is not None:
        outputs.append(Path(args.sealed_dir))
    for raw_output in outputs:
        output = raw_output if raw_output.is_absolute() else Path.cwd() / raw_output
        output = output.resolve()
        if output == allowed_output or not output.is_relative_to(allowed_output):
            raise ValueError(
                f"data output must be a new child of {allowed_output}: {output}")
        if output.exists():
            raise FileExistsError(f"data output already exists: {output}")
    print(json.dumps(build_dataset(args.spec, args.output, governance_path=args.governance,
                                   sealed_dir=args.sealed_dir), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
