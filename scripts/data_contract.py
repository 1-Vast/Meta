"""Dataset-agnostic canonical contract for few-shot DTA corpora."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable

from contracts.biological_context import (CONTEXT_SCHEMA,
    DEFAULT_ENDPOINT_VOCABULARY, UNKNOWN_ENDPOINT)


POINT_RELATION = "="
VALID_AA = frozenset("ACDEFGHIKLMNPQRSTVWYXBZUO")


@dataclass(frozen=True)
class CanonicalRow:
    row_id: str
    drug_key: str
    target_key: str
    pair_key: str
    task_key: str
    smiles: str
    sequence: str
    endpoint_key: str
    context_key: str
    context_id: int
    context_cont: tuple[float, ...]
    context_mask: tuple[int, ...]
    split: str
    y: float
    provenance: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["provenance"] = list(self.provenance)
        value["context_cont"] = list(self.context_cont)
        value["context_mask"] = list(self.context_mask)
        return value


def load_spec(path: str | Path) -> dict[str, Any]:
    """Load a JSON data specification; dataset names never enter the code path."""
    with Path(path).open(encoding="utf-8") as handle:
        value = json.load(handle)
    required = {"name", "input", "label", "split"}
    missing = required - set(value)
    if missing:
        raise ValueError(f"dataset spec missing keys: {sorted(missing)}")
    return value


def canonical_sequence(value: object) -> str:
    sequence = "".join(str(value).upper().split())
    if not sequence or set(sequence) - VALID_AA:
        raise ValueError("invalid protein sequence")
    return sequence


def canonical_smiles(value: object) -> str:
    """Canonicalize with RDKit when available; otherwise preserve a validated string."""
    smiles = str(value).strip()
    if not smiles:
        raise ValueError("empty SMILES")
    try:
        from rdkit import Chem
    except ImportError:
        return smiles
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("invalid SMILES")
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def digest(value: str) -> str:
    return sha256(value.encode("utf-8")).hexdigest()


def target_split(target_key: str, split: dict[str, Any]) -> str:
    fractions = split.get("fractions", {})
    names = ("train", "val", "test")
    if set(fractions) != set(names):
        raise ValueError("split.fractions must declare train, val, and test")
    if not math.isclose(sum(float(fractions[name]) for name in names), 1.0, abs_tol=1e-9):
        raise ValueError("split fractions must sum to 1")
    salt = str(split.get("salt", "MetaSieve-DTA/v1"))
    value = int(digest(f"{salt}:{target_key}")[:16], 16) / 2**64
    boundary = 0.0
    for name in names:
        boundary += float(fractions[name])
        if value < boundary:
            return name
    return "test"


def transform_label(raw: object, label_spec: dict[str, Any]) -> float:
    value = float(raw)
    transform = label_spec.get("transform", "identity")
    if transform == "identity":
        return value
    if transform == "negate":
        return -value
    if transform == "log10":
        if value <= 0:
            raise ValueError("log10 label must be positive")
        return math.log10(value)
    if transform == "paffinity":
        molar = value * float(label_spec.get("unit_to_molar", 1.0))
        if molar <= 0:
            raise ValueError("paffinity label must be positive")
        return -math.log10(molar)
    raise ValueError(f"unsupported label transform: {transform}")


def biological_context_registry(specification: dict[str, Any]) -> dict[str, Any]:
    declared = specification.get("biological_context", {})
    endpoint = declared.get("categorical", {}).get("endpoint", {})
    legacy_endpoint_column = specification.get("input", {}).get("columns", {}).get("gamma")
    vocabulary = tuple(str(value) for value in
                       endpoint.get("vocabulary", DEFAULT_ENDPOINT_VOCABULARY))
    unknown = str(endpoint.get("unknown", UNKNOWN_ENDPOINT))
    if not vocabulary or len(set(vocabulary)) != len(vocabulary):
        raise ValueError("endpoint vocabulary must contain unique values")
    if unknown not in vocabulary:
        raise ValueError("endpoint vocabulary must contain its unknown value")
    continuous = declared.get("continuous", {})
    features = []
    for name, item in continuous.items():
        low, high = float(item["min"]), float(item["max"])
        if not high > low:
            raise ValueError(f"continuous context {name!r} requires max > min")
        features.append({"name": str(name), "column": str(item["column"]),
                         "min": low, "max": high})
    return {"schema": CONTEXT_SCHEMA, "endpoint": {
        "column": endpoint.get("column", legacy_endpoint_column), "vocabulary": list(vocabulary),
        "unknown": unknown, "mapping": {key: i for i, key in enumerate(vocabulary)},
    }, "continuous": features}


def encode_biological_context(source: dict[str, Any], registry: dict[str, Any]) -> dict[str, Any]:
    endpoint = registry["endpoint"]
    raw = source.get(endpoint.get("column"), "") if endpoint.get("column") else ""
    lookup = {key.casefold(): key for key in endpoint["vocabulary"]}
    endpoint_key = lookup.get(str(raw).strip().casefold(), endpoint["unknown"])
    context_cont, context_mask = [], []
    for feature in registry["continuous"]:
        value = source.get(feature["column"], "")
        if value is None or not str(value).strip():
            context_cont.append(0.0)
            context_mask.append(0)
            continue
        numeric = float(value)
        scaled = (numeric - feature["min"]) / (feature["max"] - feature["min"])
        context_cont.append(min(1.0, max(0.0, scaled)))
        context_mask.append(1)
    signature = json.dumps({"endpoint": endpoint_key, "continuous": context_cont,
                            "mask": context_mask}, sort_keys=True, separators=(",", ":"))
    return {"endpoint_key": endpoint_key, "context_key": digest(signature),
            "context_id": int(endpoint["mapping"][endpoint_key]),
            "context_cont": tuple(context_cont), "context_mask": tuple(context_mask)}


def fit_normalizer(train_values: list[float], specification: dict[str, Any]) -> dict[str, float]:
    if not train_values:
        raise ValueError("normalization requires at least one training label")
    mode = specification.get("mode", "train_minmax")
    if mode == "explicit_affine":
        low, high = float(specification["min"]), float(specification["max"])
    elif mode == "train_minmax":
        low, high = min(train_values), max(train_values)
    elif mode == "train_quantile_affine":
        ordered = sorted(train_values)
        low_q, high_q = float(specification.get("low_quantile", 0.01)), float(specification.get("high_quantile", 0.99))
        low = ordered[round((len(ordered) - 1) * low_q)]
        high = ordered[round((len(ordered) - 1) * high_q)]
    else:
        raise ValueError(f"unsupported normalization mode: {mode}")
    if not high > low:
        raise ValueError("normalization bounds must satisfy high > low")
    return {"mode": mode, "low": low, "high": high}


def normalize(value: float, fitted: dict[str, float]) -> float:
    return min(1.0, max(0.0, (value - fitted["low"]) / (fitted["high"] - fitted["low"])))


def duplicate_governance(rows: list[dict[str, Any]], specification: dict[str, Any]) -> list[dict[str, Any]]:
    """Aggregate exact task/pair repeats before model indices or splits are emitted."""
    config = specification.get("duplicates", {})
    tolerance = float(config.get("label_tolerance", 0.0))
    aggregate = config.get("aggregate", "median")
    conflict = config.get("on_conflict", "fail")
    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault((row["task_key"], row["pair_key"]), []).append(row)
    resolved = []
    for key, values in groups.items():
        labels = sorted(row["label_value"] for row in values)
        if labels[-1] - labels[0] > tolerance:
            if conflict == "drop":
                continue
            raise ValueError(f"conflicting duplicate measurement for {key}")
        middle = len(labels) // 2
        label = labels[middle] if len(labels) % 2 else (labels[middle - 1] + labels[middle]) / 2
        if aggregate != "median":
            raise ValueError(f"unsupported duplicate aggregation: {aggregate}")
        chosen = dict(values[0])
        chosen["label_value"] = label
        chosen["provenance"] = tuple(row["row_id"] for row in values)
        resolved.append(chosen)
    return resolved


def write_jsonl(path: str | Path, rows: Iterable[dict[str, Any]]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")


def read_jsonl(path: str | Path) -> list[dict[str, Any]]:
    with Path(path).open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]
