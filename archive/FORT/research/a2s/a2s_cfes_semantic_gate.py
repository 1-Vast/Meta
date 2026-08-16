"""Outcome-free PLINDER semantic gate for A2S-CFES.

The gate tests whether a compact ligand-by-pocket term predicts held-cluster
structural contact profiles beyond matched additive controls. It never reads an
affinity column or PLINDER test rows, and it excludes every A2S source target.
Passing this gate authorizes only the fit-only CFES representation experiment.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import time
from typing import Callable, Iterable, Sequence
from urllib.parse import quote
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
import pyarrow.dataset as pads
import pyarrow.parquet as pq
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.decomposition import PCA
import torch
from torch import nn

from research.a2s.a2s_information_gate import (
    canonical,
    load_metadata,
    sha256_file,
    verify_lock,
)
from research.a2s.a2s_trace_stratum import DEFAULT_LOCK


ROOT = Path(__file__).resolve().parents[2]
PLINDER = ROOT / "dataset" / "public" / "plinder_2024_06_v2"
ANNOTATIONS = PLINDER / "raw" / "annotation_table.parquet"
SPLITS = PLINDER / "raw" / "split.parquet"
LOCAL_SEQUENCES = PLINDER / "processed" / "dualcold" / "uniprot_sequences.json"
SEQUENCE_CACHE = (
    PLINDER / "processed" / "dualcold" / "cfes_uniprot_sequences_2026-08-02.json"
)
TARGET_UNIPROT = ROOT / "dataset" / "public" / "chembl_37" / "processed" / "target_uniprot.json"
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_cfes_semantic_gate_2026-08-02.json"
DEFAULT_RECORDS = (
    ROOT / "reports" / "active" / "a2s_cfes_semantic_gate_records_2026-08-02.parquet"
)
DEFAULT_WEIGHTS = (
    ROOT / "reports" / "active" / "a2s_cfes_semantic_gate_weights_2026-08-02.pt"
)

SAFE_ANNOTATION_COLUMNS = (
    "system_id",
    "entry_pdb_id",
    "system_id_no_biounit",
    "system_pocket_UniProt",
    "ligand_rdkit_canonical_smiles",
    "ligand_neighboring_residues",
    "ligand_interactions",
    "ligand_is_proper",
)
SAFE_SPLIT_COLUMNS = (
    "system_id",
    "uniqueness",
    "split",
    "cluster",
    "cluster_for_val_split",
    "system_pass_validation_criteria",
    "system_pass_statistics_criteria",
    "system_proper_num_ligand_chains",
    "system_proper_pocket_num_residues",
    "system_proper_num_interactions",
)
CONTACT_TYPES = (
    "hydrogen_bonds",
    "hydrophobic_contacts",
    "water_bridges",
    "salt_bridges",
    "pi_stacks",
    "pi_cation",
    "metal_complexes",
    "halogen_bonds",
)
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: index for index, aa in enumerate(AMINO_ACIDS)}
AA_GROUPS = (
    frozenset("AILMFWVY"),
    frozenset("STNQCGP"),
    frozenset("KRH"),
    frozenset("DE"),
    frozenset("FWYH"),
    frozenset("CGP"),
)

SEEDS = (1729, 1731, 1733, 1741, 1753)
PCA_DIM = 32
MORGAN_BITS = 256
POSITION_BINS = 8
HIDDEN = 64
CROSS_RANK = 16
BATCH_SIZE = 1024
EPOCHS = 140
PATIENCE = 18
LEARNING_RATE = 2e-3
WEIGHT_DECAY = 1e-4
BOOTSTRAP_DRAWS = 10_000
INNER_PERCENT = 15
MIN_MAPPED_RESIDUES = 4
MIN_MAPPED_FRACTION = 0.80
REMOVAL_THRESHOLD = 0.70
USER_AGENT = "A2S-CFES-C0B/1.0 (outcome-free structural gate)"
UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query={query}"


def stable_seed(*parts: object) -> int:
    return int(sha256(canonical(parts).encode()).hexdigest()[:8], 16)


def assert_safe_columns(columns: Iterable[str], allowed: Sequence[str]) -> None:
    selected = tuple(columns)
    forbidden = [column for column in selected if "affinity" in column.lower()]
    if forbidden:
        raise AssertionError(f"affinity-bearing columns are prohibited: {forbidden}")
    unexpected = sorted(set(selected) - set(allowed))
    if unexpected:
        raise AssertionError(f"unregistered structural columns requested: {unexpected}")


def assert_allowed_splits(frame: pd.DataFrame) -> None:
    values = set(frame["split"].astype(str).unique())
    if values - {"train", "val"}:
        raise AssertionError(f"forbidden PLINDER split entered C0B: {sorted(values)}")


def parse_fasta(text: str) -> dict[str, str]:
    sequences: dict[str, str] = {}
    key: str | None = None
    parts: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if key is not None:
                sequences[key] = "".join(parts).upper()
            fields = line[1:].split("|", maxsplit=2)
            key = fields[1] if len(fields) >= 2 else fields[0].split()[0]
            parts = []
        elif key is not None:
            parts.append(line)
    if key is not None:
        sequences[key] = "".join(parts).upper()
    return sequences


def fetch_fasta(url: str) -> str:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=45) as response:
        return response.read().decode("utf-8")


def load_sequences(
    accessions: Iterable[str],
    fetcher: Callable[[str], str] = fetch_fasta,
) -> tuple[dict[str, str], dict[str, object]]:
    required = sorted(set(str(value) for value in accessions))
    sequences: dict[str, str] = {}
    for path in (LOCAL_SEQUENCES, SEQUENCE_CACHE):
        if path.exists():
            payload = json.loads(path.read_text(encoding="utf-8"))
            sequences.update(
                {str(key): str(value) for key, value in payload.items() if value}
            )
    missing = [accession for accession in required if accession not in sequences]
    fetched = 0
    failures: list[str] = []
    for start in range(0, len(missing), 40):
        chunk = missing[start : start + 40]
        expression = "(" + " OR ".join(f"accession:{item}" for item in chunk) + ")"
        url = UNIPROT_STREAM.format(query=quote(expression))
        error: Exception | None = None
        for attempt in range(3):
            try:
                parsed = parse_fasta(fetcher(url))
                sequences.update(parsed)
                fetched += len(parsed)
                error = None
                break
            except Exception as exc:  # network errors are recorded, not hidden
                error = exc
                if attempt < 2:
                    time.sleep(0.5 * (attempt + 1))
        if error is not None:
            failures.extend(chunk)

    retained = {key: sequences[key] for key in required if key in sequences}
    SEQUENCE_CACHE.parent.mkdir(parents=True, exist_ok=True)
    SEQUENCE_CACHE.write_text(
        json.dumps(retained, indent=2, sort_keys=True), encoding="utf-8"
    )
    metadata = {
        "requested": len(required),
        "available": len(retained),
        "fetched_entries": fetched,
        "failed_requests_accessions": len(failures),
        "missing": sorted(set(required) - set(retained)),
    }
    return retained, metadata


def a2s_accessions() -> tuple[set[str], dict[str, int]]:
    metadata = load_metadata(verify_lock(DEFAULT_LOCK))
    if "affinity" in metadata.columns:
        raise AssertionError("source metadata unexpectedly contains affinity")
    mapping = json.loads(TARGET_UNIPROT.read_text(encoding="utf-8"))
    targets = set(metadata["target"].astype(str))
    accessions = {str(mapping[target]) for target in targets if mapping.get(target)}
    return accessions, {
        "source_targets": len(targets),
        "mapped_source_accessions": len(accessions),
    }


def populated(value: object) -> bool:
    try:
        return value is not None and len(value) > 0  # type: ignore[arg-type]
    except TypeError:
        return False


def load_raw_rows() -> tuple[pd.DataFrame, dict[str, object]]:
    assert_safe_columns(SAFE_SPLIT_COLUMNS, SAFE_SPLIT_COLUMNS)
    split = pq.read_table(SPLITS, columns=list(SAFE_SPLIT_COLUMNS)).to_pandas()
    quality = (
        split["split"].isin(["train", "val"])
        & split["system_pass_validation_criteria"].astype(bool)
        & split["system_pass_statistics_criteria"].astype(bool)
        & split["system_proper_num_ligand_chains"].eq(1)
        & split["system_proper_pocket_num_residues"].gt(0)
        & split["system_proper_num_interactions"].gt(0)
    )
    split = split.loc[quality].copy()
    assert_allowed_splits(split)

    assert_safe_columns(SAFE_ANNOTATION_COLUMNS, SAFE_ANNOTATION_COLUMNS)
    identifiers = set(split["system_id"].astype(str))
    annotations = pads.dataset(ANNOTATIONS, format="parquet").to_table(
        columns=list(SAFE_ANNOTATION_COLUMNS),
        filter=pads.field("system_id").isin(identifiers),
    ).to_pandas()
    frame = annotations.merge(split, on="system_id", how="inner", validate="many_to_one")
    frame = frame.loc[
        frame["ligand_is_proper"].astype(bool)
        & frame["system_pocket_UniProt"].notna()
        & frame["ligand_rdkit_canonical_smiles"].notna()
        & frame["ligand_neighboring_residues"].map(populated)
        & frame["ligand_interactions"].map(populated)
    ].copy()
    frame["system_pocket_UniProt"] = frame["system_pocket_UniProt"].astype(str)
    frame["ligand_rdkit_canonical_smiles"] = frame[
        "ligand_rdkit_canonical_smiles"
    ].astype(str)
    excluded, source_summary = a2s_accessions()
    frame = frame.loc[~frame["system_pocket_UniProt"].isin(excluded)].copy()
    frame = frame.drop_duplicates(
        ["system_id", "system_pocket_UniProt", "ligand_rdkit_canonical_smiles"]
    ).reset_index(drop=True)
    assert_allowed_splits(frame)
    return frame, {
        "quality_split_rows": split["split"].value_counts().sort_index().to_dict(),
        "candidate_ligand_rows": frame["split"].value_counts().sort_index().to_dict(),
        "excluded_a2s": source_summary,
    }


def molecule_features(smiles: str) -> tuple[np.ndarray, np.ndarray, str, int] | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(
        molecule, radius=2, nBits=MORGAN_BITS
    )
    bits = np.zeros(MORGAN_BITS, dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fingerprint, bits)
    descriptors = np.asarray(
        [
            Descriptors.MolWt(molecule),
            Crippen.MolLogP(molecule),
            Lipinski.NumHDonors(molecule),
            Lipinski.NumHAcceptors(molecule),
            rdMolDescriptors.CalcTPSA(molecule),
            Lipinski.NumRotatableBonds(molecule),
            Lipinski.RingCount(molecule),
            molecule.GetNumHeavyAtoms(),
            rdMolDescriptors.CalcFractionCSP3(molecule),
            Chem.GetFormalCharge(molecule),
        ],
        dtype=np.float32,
    )
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(
        mol=molecule, includeChirality=False
    )
    return bits, descriptors, scaffold or "<ACYCLIC>", molecule.GetNumHeavyAtoms()


def residue_index(token: str) -> int | None:
    match = re.search(r"_(-?\d+)$", str(token))
    return int(match.group(1)) if match else None


def aa_group(aa: str) -> int:
    for index, group in enumerate(AA_GROUPS):
        if aa in group:
            return index
    return 1


def pocket_features(
    neighboring_residues: Sequence[object],
    sequence: str,
    randomize_seed: int | None = None,
) -> tuple[np.ndarray, int, float] | None:
    mapped: list[tuple[str, float]] = []
    unique_tokens = list(dict.fromkeys(str(token) for token in neighboring_residues))
    denominator = max(1, len(sequence) - 1)
    for token in unique_tokens:
        index = residue_index(token)
        if index is None or index < 0 or index >= len(sequence):
            continue
        aa = sequence[index]
        if aa in AA_INDEX:
            mapped.append((aa, index / denominator))
    fraction = len(mapped) / max(1, len(unique_tokens))
    if len(mapped) < MIN_MAPPED_RESIDUES or fraction < MIN_MAPPED_FRACTION:
        return None

    amino_acids = [item[0] for item in mapped]
    positions = [item[1] for item in mapped]
    if randomize_seed is not None and len(amino_acids) > 1:
        rng = np.random.default_rng(randomize_seed)
        amino_acids = [amino_acids[index] for index in rng.permutation(len(amino_acids))]

    aa_composition = np.zeros(len(AMINO_ACIDS), dtype=np.float32)
    group_composition = np.zeros(len(AA_GROUPS), dtype=np.float32)
    position_occupancy = np.zeros(POSITION_BINS, dtype=np.float32)
    group_positions = np.zeros((POSITION_BINS, len(AA_GROUPS)), dtype=np.float32)
    for aa, position in zip(amino_acids, positions, strict=True):
        group = aa_group(aa)
        position_bin = min(POSITION_BINS - 1, int(position * POSITION_BINS))
        aa_composition[AA_INDEX[aa]] += 1.0
        group_composition[group] += 1.0
        position_occupancy[position_bin] += 1.0
        group_positions[position_bin, group] += 1.0
    count = float(len(mapped))
    vector = np.concatenate(
        (
            aa_composition / count,
            group_composition / count,
            position_occupancy / count,
            group_positions.reshape(-1) / count,
            np.asarray([math.log1p(count)], dtype=np.float32),
        )
    ).astype(np.float32)
    return vector, len(mapped), fraction


def contact_profile(interactions: Sequence[object]) -> np.ndarray:
    counts = {name: 0 for name in CONTACT_TYPES}
    for interaction in interactions:
        text = str(interaction)
        if "_type:" not in text:
            continue
        kind = text.split("_type:", maxsplit=1)[1].split("__", maxsplit=1)[0]
        if kind in counts:
            counts[kind] += 1
    return np.log1p(np.asarray([counts[name] for name in CONTACT_TYPES], dtype=np.float32))


@dataclass(frozen=True)
class RawFeatures:
    frame: pd.DataFrame
    bits: np.ndarray
    descriptors: np.ndarray
    pocket: np.ndarray
    pocket_randomized: np.ndarray
    labels: np.ndarray


def materialize_raw_features(
    frame: pd.DataFrame, sequences: dict[str, str]
) -> tuple[RawFeatures, dict[str, object]]:
    molecule_cache: dict[str, tuple[np.ndarray, np.ndarray, str, int] | None] = {}
    rows: list[int] = []
    bits: list[np.ndarray] = []
    descriptors: list[np.ndarray] = []
    pockets: list[np.ndarray] = []
    randomized: list[np.ndarray] = []
    labels: list[np.ndarray] = []
    scaffolds: list[str] = []
    heavy_atoms: list[int] = []
    pocket_sizes: list[int] = []
    mapped_fractions: list[float] = []

    for index, row in frame.iterrows():
        accession = str(row["system_pocket_UniProt"])
        sequence = sequences.get(accession)
        if sequence is None:
            continue
        smiles = str(row["ligand_rdkit_canonical_smiles"])
        if smiles not in molecule_cache:
            molecule_cache[smiles] = molecule_features(smiles)
        ligand = molecule_cache[smiles]
        if ligand is None:
            continue
        pocket = pocket_features(row["ligand_neighboring_residues"], sequence)
        pocket_random = pocket_features(
            row["ligand_neighboring_residues"],
            sequence,
            stable_seed("residue-randomization", row["system_id"], smiles),
        )
        if pocket is None or pocket_random is None:
            continue
        rows.append(index)
        bits.append(ligand[0])
        descriptors.append(ligand[1])
        scaffolds.append(ligand[2])
        heavy_atoms.append(ligand[3])
        pockets.append(pocket[0])
        randomized.append(pocket_random[0])
        pocket_sizes.append(pocket[1])
        mapped_fractions.append(pocket[2])
        labels.append(contact_profile(row["ligand_interactions"]))

    selected = frame.loc[rows].copy().reset_index(drop=True)
    selected["scaffold"] = scaffolds
    selected["ligand_heavy_atoms"] = heavy_atoms
    selected["mapped_pocket_residues"] = pocket_sizes
    selected["mapped_residue_fraction"] = mapped_fractions
    features = RawFeatures(
        selected,
        np.stack(bits),
        np.stack(descriptors),
        np.stack(pockets),
        np.stack(randomized),
        np.stack(labels),
    )
    return features, {
        "valid_rows": selected["split"].value_counts().sort_index().to_dict(),
        "unique_ligands": len(molecule_cache),
        "median_mapped_residue_fraction": float(selected["mapped_residue_fraction"].median()),
    }


PURGE_COLUMNS = (
    "system_pocket_UniProt",
    "ligand_rdkit_canonical_smiles",
    "scaffold",
    "entry_pdb_id",
    "system_id_no_biounit",
    "cluster",
    "cluster_for_val_split",
    "uniqueness",
)


def dual_cold_indices(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    train = frame.index[frame["split"].eq("train")].to_numpy(dtype=np.int64)
    audit = frame.index[frame["split"].eq("val")].to_numpy(dtype=np.int64)
    if len(train) == 0 or len(audit) == 0:
        raise RuntimeError("C0B requires non-empty train and validation rows")
    keep = np.ones(len(train), dtype=bool)
    overlap_before: dict[str, int] = {}
    removed_by_axis: dict[str, int] = {}
    audit_frame = frame.loc[audit]
    for column in PURGE_COLUMNS:
        audit_values = set(audit_frame[column].dropna())
        train_values = frame.loc[train, column]
        overlap_before[column] = len(set(train_values.dropna()) & audit_values)
        axis_remove = train_values.isin(audit_values).to_numpy()
        removed_by_axis[column] = int(axis_remove.sum())
        keep &= ~axis_remove
    retained = train[keep]
    overlap_after = {
        column: len(
            set(frame.loc[retained, column].dropna())
            & set(audit_frame[column].dropna())
        )
        for column in PURGE_COLUMNS
    }
    if any(overlap_after.values()):
        raise AssertionError(f"dual-cold purge left overlap: {overlap_after}")
    return retained, audit, {
        "train_before": int(len(train)),
        "train_after": int(len(retained)),
        "audit_rows": int(len(audit)),
        "overlap_before": overlap_before,
        "removed_by_axis": removed_by_axis,
        "overlap_after": overlap_after,
    }


def hash_bucket(value: object, buckets: int, namespace: str) -> int:
    return stable_seed(namespace, str(value)) % buckets


@dataclass(frozen=True)
class Transform:
    mean: np.ndarray
    scale: np.ndarray

    def apply(self, values: np.ndarray) -> np.ndarray:
        return ((values - self.mean) / self.scale).astype(np.float32)


def fit_transform(values: np.ndarray) -> Transform:
    mean = values.mean(axis=0, keepdims=True).astype(np.float32)
    scale = values.std(axis=0, keepdims=True).astype(np.float32)
    scale[scale < 1e-6] = 1.0
    return Transform(mean, scale)


@dataclass(frozen=True)
class PreparedData:
    train_frame: pd.DataFrame
    audit_frame: pd.DataFrame
    train_ligand: np.ndarray
    audit_ligand: np.ndarray
    train_pocket: np.ndarray
    audit_pocket: np.ndarray
    audit_pocket_randomized: np.ndarray
    train_label: np.ndarray
    audit_label: np.ndarray
    fit_mask: np.ndarray
    inner_mask: np.ndarray


def prepare_data(raw: RawFeatures) -> tuple[PreparedData, dict[str, object]]:
    train_index, audit_index, purge = dual_cold_indices(raw.frame)
    train_frame = raw.frame.loc[train_index].reset_index(drop=True)
    audit_frame = raw.frame.loc[audit_index].reset_index(drop=True)
    train_clusters = train_frame["cluster_for_val_split"].astype(str)
    inner_mask = np.asarray(
        [hash_bucket(value, 100, "cfes-inner") < INNER_PERCENT for value in train_clusters],
        dtype=bool,
    )
    fit_mask = ~inner_mask
    if int(inner_mask.sum()) == 0 or int(fit_mask.sum()) == 0:
        raise RuntimeError("deterministic inner split is empty")

    pca = PCA(n_components=PCA_DIM, svd_solver="randomized", random_state=1729)
    pca.fit(raw.bits[train_index][fit_mask])
    train_ligand_raw = np.concatenate(
        (pca.transform(raw.bits[train_index]), raw.descriptors[train_index]), axis=1
    ).astype(np.float32)
    audit_ligand_raw = np.concatenate(
        (pca.transform(raw.bits[audit_index]), raw.descriptors[audit_index]), axis=1
    ).astype(np.float32)
    ligand_transform = fit_transform(train_ligand_raw[fit_mask])
    pocket_transform = fit_transform(raw.pocket[train_index][fit_mask])
    label_transform = fit_transform(raw.labels[train_index][fit_mask])
    prepared = PreparedData(
        train_frame=train_frame,
        audit_frame=audit_frame,
        train_ligand=ligand_transform.apply(train_ligand_raw),
        audit_ligand=ligand_transform.apply(audit_ligand_raw),
        train_pocket=pocket_transform.apply(raw.pocket[train_index]),
        audit_pocket=pocket_transform.apply(raw.pocket[audit_index]),
        audit_pocket_randomized=pocket_transform.apply(raw.pocket_randomized[audit_index]),
        train_label=label_transform.apply(raw.labels[train_index]),
        audit_label=label_transform.apply(raw.labels[audit_index]),
        fit_mask=fit_mask,
        inner_mask=inner_mask,
    )
    summary = {
        "purge": purge,
        "fit_rows": int(fit_mask.sum()),
        "inner_validation_rows": int(inner_mask.sum()),
        "fit_clusters": int(train_frame.loc[fit_mask, "cluster_for_val_split"].nunique()),
        "inner_clusters": int(train_frame.loc[inner_mask, "cluster_for_val_split"].nunique()),
        "audit_clusters": int(audit_frame["cluster_for_val_split"].nunique()),
        "ligand_dimension": int(prepared.train_ligand.shape[1]),
        "pocket_dimension": int(prepared.train_pocket.shape[1]),
        "contact_dimension": int(prepared.train_label.shape[1]),
        "pca_explained_variance": float(pca.explained_variance_ratio_.sum()),
    }
    return prepared, summary


class Branch(nn.Module):
    def __init__(self, input_dim: int, output_dim: int) -> None:
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(input_dim, HIDDEN),
            nn.SiLU(),
            nn.Linear(HIDDEN, output_dim),
        )

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return self.network(values)


class LigandOnly(nn.Module):
    def __init__(self, ligand_dim: int, output_dim: int) -> None:
        super().__init__()
        self.ligand = Branch(ligand_dim, output_dim)

    def forward(self, ligand: torch.Tensor) -> torch.Tensor:
        return self.ligand(ligand)


class AdditiveContactModel(nn.Module):
    def __init__(self, ligand_dim: int, pocket_dim: int, output_dim: int) -> None:
        super().__init__()
        self.ligand = Branch(ligand_dim, output_dim)
        self.pocket = Branch(pocket_dim, output_dim)
        self.bias = nn.Parameter(torch.zeros(output_dim))

    def forward(self, ligand: torch.Tensor, pocket: torch.Tensor) -> torch.Tensor:
        return self.ligand(ligand) + self.pocket(pocket) + self.bias


class BilinearResidual(nn.Module):
    def __init__(
        self,
        ligand_dim: int,
        pocket_dim: int,
        output_dim: int,
        frozen_random: bool = False,
    ) -> None:
        super().__init__()
        self.ligand = nn.Linear(ligand_dim, CROSS_RANK, bias=False)
        self.pocket = nn.Linear(pocket_dim, CROSS_RANK, bias=False)
        self.output = nn.Linear(CROSS_RANK, output_dim, bias=False)
        nn.init.zeros_(self.output.weight)
        if frozen_random:
            for parameter in (*self.ligand.parameters(), *self.pocket.parameters()):
                parameter.requires_grad_(False)

    def forward(self, ligand: torch.Tensor, pocket: torch.Tensor) -> torch.Tensor:
        return self.output(self.ligand(ligand) * self.pocket(pocket))


class NoCrossResidual(nn.Module):
    def __init__(self, ligand_dim: int, pocket_dim: int, output_dim: int) -> None:
        super().__init__()
        self.ligand = nn.Linear(ligand_dim, CROSS_RANK, bias=False)
        self.pocket = nn.Linear(pocket_dim, CROSS_RANK, bias=False)
        self.output = nn.Linear(CROSS_RANK, output_dim, bias=False)
        nn.init.zeros_(self.output.weight)

    def forward(self, ligand: torch.Tensor, pocket: torch.Tensor) -> torch.Tensor:
        return self.output(self.ligand(ligand) + self.pocket(pocket))


def parameter_count(model: nn.Module) -> int:
    return sum(parameter.numel() for parameter in model.parameters())


def device() -> torch.device:
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def fit_network(
    model: nn.Module,
    inputs: tuple[np.ndarray, ...],
    labels: np.ndarray,
    fit_mask: np.ndarray,
    inner_mask: np.ndarray,
    seed: int,
    base: np.ndarray | None = None,
    epochs: int = EPOCHS,
) -> tuple[nn.Module, dict[str, object]]:
    torch.manual_seed(seed)
    run_device = device()
    model = model.to(run_device)
    tensors = tuple(torch.as_tensor(value, device=run_device) for value in inputs)
    target = torch.as_tensor(labels, device=run_device)
    offset = torch.zeros_like(target) if base is None else torch.as_tensor(base, device=run_device)
    fit_rows = torch.as_tensor(np.flatnonzero(fit_mask), device=run_device)
    inner_rows = torch.as_tensor(np.flatnonzero(inner_mask), device=run_device)
    trainable = [parameter for parameter in model.parameters() if parameter.requires_grad]
    optimizer = torch.optim.AdamW(
        trainable, lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY
    )
    best_state = deepcopy(model.state_dict())
    best_loss = float("inf")
    best_epoch = -1
    stale = 0
    history: list[dict[str, float]] = []

    for epoch in range(epochs):
        model.train()
        generator = torch.Generator(device="cpu").manual_seed(
            stable_seed("fit-order", seed, epoch)
        )
        order = fit_rows[torch.randperm(len(fit_rows), generator=generator).to(run_device)]
        losses: list[float] = []
        for start in range(0, len(order), BATCH_SIZE):
            rows = order[start : start + BATCH_SIZE]
            prediction = offset[rows] + model(*(values[rows] for values in tensors))
            loss = torch.mean((prediction - target[rows]) ** 2)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(trainable, 5.0)
            optimizer.step()
            losses.append(float(loss.detach()))
        model.eval()
        with torch.no_grad():
            prediction = offset[inner_rows] + model(
                *(values[inner_rows] for values in tensors)
            )
            validation = float(torch.mean((prediction - target[inner_rows]) ** 2))
        history.append(
            {"epoch": epoch, "train_loss": float(np.mean(losses)), "inner_loss": validation}
        )
        if validation < best_loss - 1e-6:
            best_loss = validation
            best_epoch = epoch
            best_state = deepcopy(model.state_dict())
            stale = 0
        else:
            stale += 1
        if stale >= PATIENCE:
            break
    model.load_state_dict(best_state)
    model.eval()
    return model, {
        "best_epoch": int(best_epoch),
        "best_inner_loss": float(best_loss),
        "epochs_run": len(history),
        "parameter_count": parameter_count(model),
        "trainable_parameter_count": sum(
            parameter.numel() for parameter in model.parameters() if parameter.requires_grad
        ),
        "history": history,
    }


@torch.no_grad()
def predict(model: nn.Module, *inputs: np.ndarray) -> np.ndarray:
    run_device = device()
    model.eval()
    tensors = tuple(torch.as_tensor(value, device=run_device) for value in inputs)
    chunks: list[np.ndarray] = []
    for start in range(0, len(inputs[0]), 4096):
        stop = start + 4096
        chunks.append(
            model(*(value[start:stop] for value in tensors)).detach().cpu().numpy()
        )
    return np.concatenate(chunks)


def matched_donors(
    values: np.ndarray,
    forbidden: Sequence[object],
    seed: int,
    bins: int = 10,
) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    forbidden_values = np.asarray([str(value) for value in forbidden], dtype=object)
    if len(values) < 2:
        raise ValueError("a destruction requires at least two audit rows")
    quantiles = np.unique(np.quantile(values, np.linspace(0.0, 1.0, bins + 1)))
    if len(quantiles) <= 2:
        bucket = np.zeros(len(values), dtype=int)
    else:
        bucket = np.digitize(values, quantiles[1:-1], right=True)
    rng = np.random.default_rng(seed)
    donors = np.empty(len(values), dtype=np.int64)
    all_rows = np.arange(len(values))
    for row in all_rows:
        candidates = all_rows[
            (bucket == bucket[row]) & (forbidden_values != forbidden_values[row])
        ]
        if len(candidates) == 0:
            candidates = all_rows[forbidden_values != forbidden_values[row]]
        if len(candidates) == 0:
            raise ValueError("no physically distinct donor is available")
        distance = np.abs(values[candidates] - values[row])
        nearest = candidates[distance == distance.min()]
        donors[row] = int(rng.choice(nearest))
    return donors


def normalized_state_pool(states: np.ndarray) -> np.ndarray:
    values = np.asarray(states)
    if values.ndim < 2:
        raise ValueError("state pooling requires a state axis and feature axis")
    unique = np.unique(values, axis=0)
    return unique.mean(axis=0)


@dataclass(frozen=True)
class FittedModels:
    ligand_only: LigandOnly
    additive: AdditiveContactModel
    cross: BilinearResidual
    no_cross: NoCrossResidual
    random_cross: BilinearResidual
    training: dict[str, object]


def fit_all_models(data: PreparedData, seed: int, epochs: int = EPOCHS) -> FittedModels:
    ligand_dim = data.train_ligand.shape[1]
    pocket_dim = data.train_pocket.shape[1]
    output_dim = data.train_label.shape[1]
    torch.manual_seed(seed)
    ligand, ligand_log = fit_network(
        LigandOnly(ligand_dim, output_dim),
        (data.train_ligand,),
        data.train_label,
        data.fit_mask,
        data.inner_mask,
        stable_seed(seed, "ligand"),
        epochs=epochs,
    )
    torch.manual_seed(seed)
    additive, additive_log = fit_network(
        AdditiveContactModel(ligand_dim, pocket_dim, output_dim),
        (data.train_ligand, data.train_pocket),
        data.train_label,
        data.fit_mask,
        data.inner_mask,
        stable_seed(seed, "additive"),
        epochs=epochs,
    )
    additive_train = predict(additive, data.train_ligand, data.train_pocket)

    residual_logs: dict[str, object] = {}
    residual_models: dict[str, nn.Module] = {}
    constructors: dict[str, Callable[[], nn.Module]] = {
        "cross": lambda: BilinearResidual(ligand_dim, pocket_dim, output_dim),
        "no_cross": lambda: NoCrossResidual(ligand_dim, pocket_dim, output_dim),
        "random_cross": lambda: BilinearResidual(
            ligand_dim, pocket_dim, output_dim, frozen_random=True
        ),
    }
    for name, constructor in constructors.items():
        torch.manual_seed(seed)
        fitted, log = fit_network(
            constructor(),
            (data.train_ligand, data.train_pocket),
            data.train_label,
            data.fit_mask,
            data.inner_mask,
            stable_seed(seed, name),
            base=additive_train,
            epochs=epochs,
        )
        residual_models[name] = fitted
        residual_logs[name] = log
    return FittedModels(
        ligand_only=ligand,  # type: ignore[arg-type]
        additive=additive,  # type: ignore[arg-type]
        cross=residual_models["cross"],  # type: ignore[arg-type]
        no_cross=residual_models["no_cross"],  # type: ignore[arg-type]
        random_cross=residual_models["random_cross"],  # type: ignore[arg-type]
        training={
            "ligand_only": ligand_log,
            "additive": additive_log,
            **residual_logs,
        },
    )


def squared_error(label: np.ndarray, prediction: np.ndarray) -> np.ndarray:
    return (prediction - label) ** 2


def evaluate_seed(data: PreparedData, models: FittedModels, seed: int) -> pd.DataFrame:
    frame = data.audit_frame[
        [
            "system_id",
            "entry_pdb_id",
            "system_pocket_UniProt",
            "ligand_rdkit_canonical_smiles",
            "scaffold",
            "cluster_for_val_split",
            "mapped_pocket_residues",
            "ligand_heavy_atoms",
        ]
    ].copy()
    frame["model_seed"] = seed
    frame["audit_fold"] = [
        hash_bucket(value, 5, "cfes-audit-fold")
        for value in frame["cluster_for_val_split"]
    ]

    ligand_donor = matched_donors(
        frame["ligand_heavy_atoms"].to_numpy(),
        frame["scaffold"],
        stable_seed(seed, "ligand-shuffle"),
    )
    pocket_donor = matched_donors(
        frame["mapped_pocket_residues"].to_numpy(),
        frame["system_pocket_UniProt"],
        stable_seed(seed, "pocket-shuffle"),
    )
    transplant_donor = matched_donors(
        frame["mapped_pocket_residues"].to_numpy(),
        frame["system_pocket_UniProt"],
        stable_seed(seed, "structure-transplant"),
        bins=20,
    )

    ligand_prediction = predict(models.ligand_only, data.audit_ligand)
    additive = predict(models.additive, data.audit_ligand, data.audit_pocket)
    cross_residual = predict(models.cross, data.audit_ligand, data.audit_pocket)
    no_cross = additive + predict(models.no_cross, data.audit_ligand, data.audit_pocket)
    random_cross = additive + predict(
        models.random_cross, data.audit_ligand, data.audit_pocket
    )
    cross = additive + cross_residual
    pocket_shuffle = additive + predict(
        models.cross, data.audit_ligand, data.audit_pocket[pocket_donor]
    )
    ligand_shuffle = additive + predict(
        models.cross, data.audit_ligand[ligand_donor], data.audit_pocket
    )
    transplant_additive = predict(
        models.additive, data.audit_ligand, data.audit_pocket[transplant_donor]
    )
    transplant_cross = transplant_additive + predict(
        models.cross, data.audit_ligand, data.audit_pocket[transplant_donor]
    )
    residue_additive = predict(
        models.additive, data.audit_ligand, data.audit_pocket_randomized
    )
    residue_cross = residue_additive + predict(
        models.cross, data.audit_ligand, data.audit_pocket_randomized
    )

    predictions = {
        "ligand_only": ligand_prediction,
        "additive": additive,
        "cross": cross,
        "no_cross_capacity": no_cross,
        "frozen_random_cross": random_cross,
        "pocket_shuffle_cross": pocket_shuffle,
        "ligand_shuffle_cross": ligand_shuffle,
        "transplant_additive": transplant_additive,
        "transplant_cross": transplant_cross,
        "residue_additive": residue_additive,
        "residue_cross": residue_cross,
    }
    for name, prediction_value in predictions.items():
        errors = squared_error(data.audit_label, prediction_value)
        frame[f"{name}__loss"] = errors.mean(axis=1)
        for index, contact in enumerate(CONTACT_TYPES):
            frame[f"{name}__{contact}"] = errors[:, index]
    return frame


def cluster_bootstrap(
    records: pd.DataFrame,
    values: pd.Series,
    seed: int,
    draws: int = BOOTSTRAP_DRAWS,
) -> dict[str, float]:
    working = records[["cluster_for_val_split", "model_seed"]].copy()
    working["value"] = np.asarray(values, dtype=float)
    by_seed = working.groupby(
        ["cluster_for_val_split", "model_seed"], sort=True
    )["value"].mean()
    by_cluster = by_seed.groupby(level=0).mean()
    array = by_cluster.to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    sampled = rng.integers(0, len(array), size=(draws, len(array)))
    means = array[sampled].mean(axis=1)
    return {
        "mean": float(array.mean()),
        "lower95": float(np.quantile(means, 0.025)),
        "upper95": float(np.quantile(means, 0.975)),
        "clusters": int(len(array)),
    }


def effect(records: pd.DataFrame, left: str, right: str, name: str) -> dict[str, float]:
    values = records[f"{right}__loss"] - records[f"{left}__loss"]
    return cluster_bootstrap(records, values, stable_seed("bootstrap", name))


def mean_cluster_effect(frame: pd.DataFrame, values: pd.Series) -> float:
    copy = frame[["cluster_for_val_split"]].copy()
    copy["value"] = np.asarray(values, dtype=float)
    return float(copy.groupby("cluster_for_val_split")["value"].mean().mean())


def summarise_records(records: pd.DataFrame) -> dict[str, object]:
    effects = {
        "cross_minus_additive": effect(records, "cross", "additive", "cross-additive"),
        "cross_minus_ligand_only": effect(
            records, "cross", "ligand_only", "cross-ligand"
        ),
        "cross_minus_no_cross_capacity": effect(
            records, "cross", "no_cross_capacity", "cross-no-cross"
        ),
        "cross_minus_frozen_random": effect(
            records, "cross", "frozen_random_cross", "cross-random"
        ),
    }
    destroyed = {
        "pocket_shuffle": effect(
            records, "pocket_shuffle_cross", "additive", "pocket-shuffle"
        ),
        "ligand_shuffle": effect(
            records, "ligand_shuffle_cross", "additive", "ligand-shuffle"
        ),
        "structure_transplant": effect(
            records, "transplant_cross", "transplant_additive", "transplant"
        ),
        "residue_randomization": effect(
            records, "residue_cross", "residue_additive", "residue-random"
        ),
    }
    primary = effects["cross_minus_additive"]["mean"]
    removals = {
        name: float(1.0 - value["mean"] / primary) if abs(primary) > 1e-12 else float("nan")
        for name, value in destroyed.items()
    }
    per_seed: dict[str, float] = {}
    for seed in SEEDS:
        cell = records.loc[records["model_seed"] == seed]
        per_seed[str(seed)] = mean_cluster_effect(
            cell, cell["additive__loss"] - cell["cross__loss"]
        )
    per_fold: dict[str, float] = {}
    for fold in range(5):
        cell = records.loc[records["audit_fold"] == fold]
        per_fold[str(fold)] = mean_cluster_effect(
            cell, cell["additive__loss"] - cell["cross__loss"]
        )
    per_contact: dict[str, float] = {}
    for contact in CONTACT_TYPES:
        per_contact[contact] = mean_cluster_effect(
            records,
            records[f"additive__{contact}"] - records[f"cross__{contact}"],
        )
    return {
        "effects": effects,
        "destroyed_effects": destroyed,
        "removal_fraction": removals,
        "per_seed_cross_minus_additive": per_seed,
        "per_fold_cross_minus_additive": per_fold,
        "per_contact_cross_minus_additive": per_contact,
    }


def synthetic_control(
    ligand_dim: int,
    pocket_dim: int,
    output_dim: int,
    train_size: int,
    audit_clusters: Sequence[object],
    epochs: int = 90,
) -> dict[str, object]:
    seed = 99173
    rng = np.random.default_rng(seed)
    audit_clusters = np.asarray([str(value) for value in audit_clusters])
    audit_size = len(audit_clusters)
    train_size = max(3000, train_size)
    train_l = rng.normal(size=(train_size, ligand_dim)).astype(np.float32)
    train_p = rng.normal(size=(train_size, pocket_dim)).astype(np.float32)
    audit_l = rng.normal(size=(audit_size, ligand_dim)).astype(np.float32)
    audit_p = rng.normal(size=(audit_size, pocket_dim)).astype(np.float32)
    left = rng.normal(scale=1 / math.sqrt(ligand_dim), size=(ligand_dim, 4))
    right = rng.normal(scale=1 / math.sqrt(pocket_dim), size=(pocket_dim, 4))
    output = rng.normal(scale=0.8, size=(4, output_dim))
    ligand_main = rng.normal(scale=0.15, size=(ligand_dim, output_dim))
    pocket_main = rng.normal(scale=0.15, size=(pocket_dim, output_dim))

    def make_label(ligand: np.ndarray, pocket: np.ndarray) -> np.ndarray:
        interaction = ((ligand @ left) * (pocket @ right)) @ output
        additive = ligand @ ligand_main + pocket @ pocket_main
        noise = rng.normal(scale=0.45, size=interaction.shape)
        return (additive + interaction + noise).astype(np.float32)

    train_y = make_label(train_l, train_p)
    audit_y = make_label(audit_l, audit_p)
    inner = np.asarray(
        [stable_seed("synthetic-inner", index) % 100 < INNER_PERCENT for index in range(train_size)]
    )
    fit = ~inner
    synthetic_frame = pd.DataFrame(
        {
            "cluster_for_val_split": audit_clusters,
            "model_seed": seed,
            "audit_fold": [hash_bucket(value, 5, "cfes-audit-fold") for value in audit_clusters],
        }
    )
    synthetic_data = PreparedData(
        train_frame=pd.DataFrame(index=np.arange(train_size)),
        audit_frame=pd.DataFrame(index=np.arange(audit_size)),
        train_ligand=train_l,
        audit_ligand=audit_l,
        train_pocket=train_p,
        audit_pocket=audit_p,
        audit_pocket_randomized=audit_p.copy(),
        train_label=train_y,
        audit_label=audit_y,
        fit_mask=fit,
        inner_mask=inner,
    )
    models = fit_all_models(synthetic_data, seed, epochs=epochs)
    additive = predict(models.additive, audit_l, audit_p)
    cross = additive + predict(models.cross, audit_l, audit_p)
    pocket_donor = np.roll(np.arange(audit_size), 1)
    ligand_donor = np.roll(np.arange(audit_size), -1)
    pocket_shuffle = additive + predict(models.cross, audit_l, audit_p[pocket_donor])
    ligand_shuffle = additive + predict(models.cross, audit_l[ligand_donor], audit_p)
    for name, prediction_value in {
        "additive": additive,
        "cross": cross,
        "pocket_shuffle_cross": pocket_shuffle,
        "ligand_shuffle_cross": ligand_shuffle,
    }.items():
        synthetic_frame[f"{name}__loss"] = squared_error(audit_y, prediction_value).mean(axis=1)
    gain = cluster_bootstrap(
        synthetic_frame,
        synthetic_frame["additive__loss"] - synthetic_frame["cross__loss"],
        stable_seed("synthetic-bootstrap"),
        draws=2000,
    )
    destroyed = {}
    removals = {}
    for name in ("pocket_shuffle", "ligand_shuffle"):
        value = mean_cluster_effect(
            synthetic_frame,
            synthetic_frame["additive__loss"]
            - synthetic_frame[f"{name}_cross__loss"],
        )
        destroyed[name] = value
        removals[name] = 1.0 - value / gain["mean"]
    fold_gain = {}
    for fold in range(5):
        cell = synthetic_frame.loc[synthetic_frame.audit_fold == fold]
        fold_gain[str(fold)] = mean_cluster_effect(
            cell, cell["additive__loss"] - cell["cross__loss"]
        )
    passed = (
        gain["lower95"] > 0.0
        and all(value > 0.0 for value in fold_gain.values())
        and all(value >= REMOVAL_THRESHOLD for value in removals.values())
    )
    return {
        "pass": bool(passed),
        "gain": gain,
        "destroyed_gain": destroyed,
        "removal_fraction": removals,
        "fold_gain": fold_gain,
        "train_rows": train_size,
        "audit_rows": audit_size,
        "true_rank": 4,
        "noise_sd": 0.45,
    }


def invariance_checks() -> dict[str, object]:
    rng = np.random.default_rng(1729)
    states = rng.normal(size=(3, len(CONTACT_TYPES))).astype(np.float64)
    original = normalized_state_pool(states)
    permuted = normalized_state_pool(states[[2, 0, 1]])
    duplicated = normalized_state_pool(np.concatenate((states, states[[1]]), axis=0))
    order_error = float(np.max(np.abs(original - permuted)))
    duplication_error = float(np.max(np.abs(original - duplicated)))
    return {
        "state_order_max_abs_error": order_error,
        "state_duplication_max_abs_error": duplication_error,
        "pass": order_error == 0.0 and duplication_error == 0.0,
    }


def decide(
    synthetic: dict[str, object],
    summary: dict[str, object],
    invariance: dict[str, object],
    firewall: dict[str, object],
) -> dict[str, object]:
    effects = summary["effects"]
    removals = summary["removal_fraction"]
    contacts = summary["per_contact_cross_minus_additive"]
    checks = {
        "synthetic_positive_control": bool(synthetic["pass"]),
        "cross_beats_additive": effects["cross_minus_additive"]["lower95"] > 0.0,
        "cross_beats_ligand_only": effects["cross_minus_ligand_only"]["lower95"] > 0.0,
        "cross_beats_no_cross_capacity": effects["cross_minus_no_cross_capacity"]["lower95"] > 0.0,
        "cross_beats_frozen_random": effects["cross_minus_frozen_random"]["lower95"] > 0.0,
        "every_seed_positive": all(
            value > 0.0 for value in summary["per_seed_cross_minus_additive"].values()
        ),
        "every_audit_fold_positive": all(
            value > 0.0 for value in summary["per_fold_cross_minus_additive"].values()
        ),
        "contact_breadth": (
            sum(value > 0.0 for value in contacts.values()) >= 4
            and contacts["hydrogen_bonds"] > 0.0
            and contacts["hydrophobic_contacts"] > 0.0
        ),
        "physical_destruction": all(
            value >= REMOVAL_THRESHOLD for value in removals.values()
        ),
        "state_invariances": bool(invariance["pass"]),
        "firewall": bool(firewall["pass"]),
    }
    passed = all(checks.values())
    if not synthetic["pass"]:
        verdict = "CFES_C0B_HARNESS_INVALID"
    elif passed:
        verdict = "CFES_C0B_SEMANTICS_ADMITTED_PROCEED_C1"
    else:
        verdict = "CFES_C0B_SEMANTICS_NOT_ADMITTED_STOP_CFES"
    return {"checks": checks, "pass": bool(passed), "verdict": verdict}


def run(output: Path, records_path: Path, weights_path: Path) -> dict[str, object]:
    started = time.perf_counter()
    frame, loading_summary = load_raw_rows()
    sequences, sequence_summary = load_sequences(frame["system_pocket_UniProt"])
    raw, feature_summary = materialize_raw_features(frame, sequences)
    data, preparation_summary = prepare_data(raw)
    firewall = {
        "annotation_columns": list(SAFE_ANNOTATION_COLUMNS),
        "split_columns": list(SAFE_SPLIT_COLUMNS),
        "affinity_columns_read": False,
        "plinder_splits_read": ["train", "val"],
        "processed_registry_read": False,
        "source_affinity_read": False,
        "source_probe_labels_read": False,
        "source_locked_labels_read": False,
        "recipient_labels_read": False,
        "a2s_accessions_excluded": loading_summary["excluded_a2s"]["mapped_source_accessions"],
        "overlap_after_purge": preparation_summary["purge"]["overlap_after"],
    }
    firewall["pass"] = (
        not firewall["affinity_columns_read"]
        and not firewall["processed_registry_read"]
        and firewall["a2s_accessions_excluded"] > 0
        and not any(firewall["overlap_after_purge"].values())
    )

    synthetic = synthetic_control(
        data.train_ligand.shape[1],
        data.train_pocket.shape[1],
        data.train_label.shape[1],
        len(data.train_frame),
        data.audit_frame["cluster_for_val_split"],
    )
    if not synthetic["pass"]:
        records = pd.DataFrame()
        training: dict[str, object] = {}
        summary: dict[str, object] = {
            "effects": {},
            "removal_fraction": {},
            "per_seed_cross_minus_additive": {},
            "per_fold_cross_minus_additive": {},
            "per_contact_cross_minus_additive": {},
        }
    else:
        record_parts: list[pd.DataFrame] = []
        training = {}
        checkpoints: dict[str, dict[str, torch.Tensor]] = {}
        for seed in SEEDS:
            models = fit_all_models(data, seed)
            record_parts.append(evaluate_seed(data, models, seed))
            training[str(seed)] = models.training
            checkpoints[str(seed)] = {
                f"additive.{key}": value.detach().cpu()
                for key, value in models.additive.state_dict().items()
            } | {
                f"cross.{key}": value.detach().cpu()
                for key, value in models.cross.state_dict().items()
            }
        records = pd.concat(record_parts, ignore_index=True)
        summary = summarise_records(records)
        weights_path.parent.mkdir(parents=True, exist_ok=True)
        torch.save(checkpoints, weights_path)

    invariance = invariance_checks()
    decision = decide(synthetic, summary, invariance, firewall)
    output.parent.mkdir(parents=True, exist_ok=True)
    records_path.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    payload: dict[str, object] = {
        "schema": "a2s-cfes-c0b-semantic-gate-v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "device": str(device()),
        "cuda_device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "runtime_seconds": float(time.perf_counter() - started),
        "protocol": {
            "seeds": list(SEEDS),
            "contact_types": list(CONTACT_TYPES),
            "cross_rank": CROSS_RANK,
            "hidden": HIDDEN,
            "epochs": EPOCHS,
            "patience": PATIENCE,
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "removal_threshold": REMOVAL_THRESHOLD,
        },
        "firewall": firewall,
        "data": {
            "loading": loading_summary,
            "sequences": sequence_summary,
            "features": feature_summary,
            "preparation": preparation_summary,
        },
        "synthetic_positive_control": synthetic,
        "training": training,
        "summary": summary,
        "invariance": invariance,
        "decision": decision,
        "artifacts": {
            "annotation_sha256": sha256_file(ANNOTATIONS),
            "split_sha256": sha256_file(SPLITS),
            "source_lock_sha256": sha256_file(DEFAULT_LOCK),
            "target_uniprot_sha256": sha256_file(TARGET_UNIPROT),
            "sequence_cache_sha256": sha256_file(SEQUENCE_CACHE),
            "records": str(records_path.relative_to(ROOT)),
            "records_sha256": sha256_file(records_path),
            "weights": str(weights_path.relative_to(ROOT)) if weights_path.exists() else None,
            "weights_sha256": sha256_file(weights_path) if weights_path.exists() else None,
        },
    }
    payload["content_sha256"] = sha256(canonical(payload).encode()).hexdigest()
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--weights", type=Path, default=DEFAULT_WEIGHTS)
    args = parser.parse_args()
    payload = run(args.output, args.records, args.weights)
    print(payload["decision"]["verdict"])
    print(canonical(payload["decision"]["checks"]))


if __name__ == "__main__":
    main()
