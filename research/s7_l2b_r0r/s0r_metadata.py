"""Build the label-blind pair universe for the Phase 2B synthetic replay.

The source JSONL also contains structural interaction labels. This module
copies an explicit metadata whitelist and never exposes those label fields to
the synthetic runner. Pair eligibility depends only on exact construct,
ligand graph, Murcko scaffold, the frozen protein closure, and the frozen split.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import pickle
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(r"D:\MetaSieve")
RAW = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r1_raw_corpus"
MONN = ROOT / "dataset" / "raw" / "monn" / "MONN" / "data"
ESM_INDEX = (ROOT / "dataset" / "processed" / "s7_l2b_r0r" /
             "esm2_650M" / "esm2_650M_index.json")
OUT = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "phase2b_s0r"
REPORT = ROOT / "report" / "s7_l2b_r0r" / "S0R_METADATA_ONLY_CENSUS.json"

SEED_SPLIT = 20260810
TEST_TARGET_FRACTION = 0.20
METADATA_FIELDS = (
    "source_key", "pdb_id", "uniprot_id", "uniprot_sequence", "ligand_ccd",
)


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_jsonl(path: Path, rows) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for row in rows:
            f.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def read_metadata(path: Path, cohort: str):
    """Return only the declared metadata whitelist from the mixed source."""
    rows = []
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            source = json.loads(line)
            row = {key: source[key] for key in METADATA_FIELDS}
            row["cohort"] = cohort
            rows.append(row)
    return rows


def load_molecules():
    return [
        pickle.load((MONN / "mol_dict").open("rb"), encoding="bytes"),
        pickle.load((MONN / "independent_dataset_mol_dict").open("rb"),
                    encoding="bytes"),
    ]


def ligand_identity(ccd: str, dictionaries, cache):
    if ccd in cache:
        return cache[ccd]
    mol = None
    key = ccd.encode("ascii", "ignore")
    for dictionary in dictionaries:
        mol = dictionary.get(key)
        if mol is not None:
            break
    if mol is None:
        cache[ccd] = None
        return None
    work = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(work)
        smiles = Chem.MolToSmiles(work, isomericSmiles=False, canonical=True)
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=work)
    except Exception:
        cache[ccd] = None
        return None
    value = {
        "graph_key": sha_bytes(smiles.encode()),
        "scaffold": scaffold,
        "n_atoms": work.GetNumAtoms(),
    }
    cache[ccd] = value
    return value


def protein_components(rows):
    # Importing this function is safe: it uses only the metadata fields below
    # and the already frozen homology alignment table.
    from s7_dataset import protein_components as frozen_components
    return frozen_components(rows)


def split_rows(rows, component_of):
    dev = [row for row in rows if row["cohort"] == "development"]
    by_component = defaultdict(list)
    for row in dev:
        by_component[component_of[row["source_key"]]].append(row)

    import numpy as np
    rng = np.random.default_rng(SEED_SPLIT)
    order = sorted(by_component, key=lambda c: (-len(by_component[c]), c))
    rng.shuffle(order)
    target = int(TEST_TARGET_FRACTION * len(dev))
    test_components, count = set(), 0
    for component in order:
        if count >= target:
            break
        test_components.add(component)
        count += len(by_component[component])

    train = [row for row in dev if component_of[row["source_key"]] not in test_components]
    held = [row for row in dev if component_of[row["source_key"]] in test_components]
    train_graphs = {row["graph_key"] for row in train}
    held_a = [row for row in held if row["graph_key"] not in train_graphs]
    return train, held_a


def build_pairs(rows):
    by_construct = defaultdict(list)
    for row in rows:
        by_construct[row["seq_key"]].append(row)
    pairs, excluded = [], Counter()
    for seq_key in sorted(by_construct):
        records = sorted(by_construct[seq_key], key=lambda row: row["source_key"])
        for i, left in enumerate(records):
            for right in records[i + 1:]:
                if left["graph_key"] == right["graph_key"]:
                    excluded["same_graph"] += 1
                    continue
                if not (left["scaffold"] and right["scaffold"] and
                        left["scaffold"] != right["scaffold"]):
                    excluded["scaffold_not_distinct"] += 1
                    continue
                pairs.append({
                    "seq_key": seq_key,
                    "a": left["source_key"],
                    "b": right["source_key"],
                })
    return pairs, dict(excluded)


def in_train_panel(pair) -> bool:
    token = f"{pair['a']}|{pair['b']}".encode()
    return int(hashlib.sha256(token).hexdigest()[:16], 16) % 16 == 0


def main() -> int:
    RDLogger.DisableLog("rdApp.*")
    rows = read_metadata(RAW / "monn_development_edge_corpus.jsonl.gz", "development")
    rows += read_metadata(RAW / "monn_additional_pdb_edge_corpus.jsonl.gz", "additional_pdb")
    dictionaries = load_molecules()
    esm_index = json.loads(ESM_INDEX.read_text(encoding="utf-8"))
    identity_cache, kept, exclusions = {}, [], Counter()

    for row in rows:
        sequence = row["uniprot_sequence"]
        row["seq_key"] = sha_bytes(sequence.encode())
        if row["seq_key"] not in esm_index:
            exclusions["sequence_not_in_frozen_esm_cache"] += 1
            continue
        identity = ligand_identity(row["ligand_ccd"], dictionaries, identity_cache)
        if identity is None:
            exclusions["ligand_graph_unavailable"] += 1
            continue
        row.update(identity)
        row["n_res"] = len(sequence)
        kept.append(row)

    component_of = protein_components(kept)
    for row in kept:
        row["component"] = component_of[row["source_key"]]
    train, held_a = split_rows(kept, component_of)
    train_pairs, train_pair_exclusions = build_pairs(train)
    held_pairs, held_pair_exclusions = build_pairs(held_a)
    component_by_construct = {}
    for row in kept:
        component_by_construct.setdefault(row["seq_key"], row["component"])
    train_pair_constructs = {pair["seq_key"] for pair in train_pairs}
    held_pair_constructs = {pair["seq_key"] for pair in held_pairs}
    train_pair_components = {component_by_construct[key] for key in train_pair_constructs}
    held_pair_components = {component_by_construct[key] for key in held_pair_constructs}
    train_panel = [pair for pair in train_pairs if in_train_panel(pair)]
    train_panel_constructs = {pair["seq_key"] for pair in train_panel}
    train_panel_components = {component_by_construct[key] for key in train_panel_constructs}

    OUT.mkdir(parents=True, exist_ok=True)
    metadata_path = OUT / "metadata_only_records.jsonl"
    train_pairs_path = OUT / "train_pairs.jsonl"
    held_pairs_path = OUT / "heldoutA_pairs.jsonl"
    write_jsonl(metadata_path, sorted(kept, key=lambda row: row["source_key"]))
    write_jsonl(train_pairs_path, train_pairs)
    write_jsonl(held_pairs_path, held_pairs)

    train_components = {row["component"] for row in train}
    held_components = {row["component"] for row in held_a}
    report = {
        "schema": "MetaSieve.S7L2B.P2B.S0R.MetadataOnlyCensus.v1",
        "source_policy": {
            "whitelist": list(METADATA_FIELDS),
            "structural_label_fields_exposed_to_runner": 0,
            "affinity_value_reads": 0,
            "pair_eligibility": [
                "exact_sequence", "different_exact_ligand_graph",
                "different_nonempty_murcko_scaffold",
            ],
        },
        "records": {
            "source": len(rows), "kept": len(kept), "exclusions": dict(exclusions),
            "train": len(train), "heldoutA": len(held_a),
        },
        "pairs": {
            "train": len(train_pairs), "heldoutA": len(held_pairs),
            "train_constructs": len(train_pair_constructs),
            "heldoutA_constructs": len(held_pair_constructs),
            "train_components": len(train_pair_components),
            "heldoutA_components": len(held_pair_components),
            "hash_stratified_train_panel": len(train_panel),
            "hash_stratified_train_constructs": len(train_panel_constructs),
            "hash_stratified_train_components": len(train_panel_components),
            "train_exclusions": train_pair_exclusions,
            "heldoutA_exclusions": held_pair_exclusions,
        },
        "components": {
            "train": len(train_components), "heldoutA": len(held_components),
            "overlap": len(train_components & held_components),
        },
        "artifacts": {
            str(metadata_path.relative_to(ROOT)).replace("\\", "/"): sha_file(metadata_path),
            str(train_pairs_path.relative_to(ROOT)).replace("\\", "/"): sha_file(train_pairs_path),
            str(held_pairs_path.relative_to(ROOT)).replace("\\", "/"): sha_file(held_pairs_path),
        },
    }
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
