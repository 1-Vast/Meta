"""Governed biological banks and target-disjoint episodes for QPSMP."""
from __future__ import annotations

from collections import OrderedDict, defaultdict
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM, GRAPH_SCHEMA, MAX_ATOMS
from scripts.build_ligand_bank import featurize_smiles


@dataclass(frozen=True)
class EpisodeSpec:
    split: str
    component: str
    target: str
    support: tuple[int, ...]
    query: tuple[int, ...]
    donor_target: str


@dataclass(frozen=True)
class EpisodeBatch:
    spec: EpisodeSpec
    protein_pooled: torch.Tensor
    protein_tokens: torch.Tensor
    protein_mask: torch.Tensor
    support_atoms: torch.Tensor
    support_bonds: torch.Tensor
    support_mask: torch.Tensor
    support_y: torch.Tensor
    query_atoms: torch.Tensor
    query_bonds: torch.Tensor
    query_mask: torch.Tensor
    query_y: torch.Tensor


def stable_seed(*parts: object) -> int:
    raw = "|".join(map(str, parts)).encode("utf-8")
    return int(hashlib.sha256(raw).hexdigest()[:16], 16)


class _ShardedBank:
    def __init__(self, directory: Path, fields: tuple[str, ...], cache_size: int = 2):
        self.directory = directory
        self.fields = fields
        self.cache_size = cache_size
        self.index: dict[str, tuple[Path, int]] = {}
        self.cache: OrderedDict[Path, dict[str, np.ndarray]] = OrderedDict()
        manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
        self.manifest = manifest
        shards = manifest.get("shards", [])
        for item in shards:
            name = item["path"] if isinstance(item, dict) else item
            path = directory / name
            with np.load(path, allow_pickle=False) as stored:
                missing = {"keys", *fields}.difference(stored.files)
                if missing:
                    raise ValueError(f"{path} is missing bank fields {sorted(missing)}")
                for local, key in enumerate(stored["keys"].astype(str)):
                    if key in self.index:
                        raise ValueError(f"duplicate bank key: {key}")
                    self.index[key] = path, local

    def __len__(self) -> int:
        return len(self.index)

    def get(self, key: str) -> tuple[np.ndarray, ...]:
        if key not in self.index:
            raise KeyError(f"bank does not contain {key}")
        path, local = self.index[key]
        if path not in self.cache:
            with np.load(path, allow_pickle=False) as stored:
                self.cache[path] = {field: stored[field] for field in self.fields}
            self.cache.move_to_end(path)
            while len(self.cache) > self.cache_size:
                self.cache.popitem(last=False)
        return tuple(self.cache[path][field][local] for field in self.fields)


class QPSMPData:
    """Read-only main-v0 corpus with lazy, schema-checked biological banks."""

    def __init__(self, corpus: Path, protein_bank: Path, ligand_bank: Path):
        self.corpus = corpus
        manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
        if manifest.get("schema") != "MetaSieve.MetaFewshot.MainV0Corpus.v1":
            raise ValueError("unsupported corpus schema")
        if manifest.get("training_authorized") is not True:
            raise ValueError("corpus manifest does not authorize training")
        self.manifest = manifest
        self.cells = self._read_cells(corpus / "cells.jsonl.gz")
        proteins = self._read_jsonl(corpus / "proteins.jsonl")
        ligands = self._read_jsonl(corpus / "ligands.jsonl")
        self.ligand_smiles = {row["drug_key"]: row["smiles"] for row in ligands}
        self.graph_cache: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        protein_keys = {row["sequence_sha256"] for row in proteins}
        self.protein_bank = _ShardedBank(protein_bank, ("pooled", "residues", "mask"))
        self.ligand_bank = _ShardedBank(ligand_bank, ("X", "A", "mask"))
        self._validate_manifests(proteins, ligands)
        self.tasks, self.components = self._build_tasks()
        target_keys = {cell["target_id"] for cell in self.cells}
        if target_keys != protein_keys:
            raise ValueError("target IDs do not equal the governed protein-bank keys")

    @staticmethod
    def _read_jsonl(path: Path) -> list[dict]:
        with path.open(encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]

    @staticmethod
    def _read_cells(path: Path) -> list[dict]:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]

    def _validate_manifests(self, proteins: list[dict], ligands: list[dict]) -> None:
        pm, lm = self.protein_bank.manifest, self.ligand_bank.manifest
        if pm.get("schema") != "MetaSieve.StructureProteinBank.v1":
            raise ValueError("unsupported protein-bank schema")
        if lm.get("schema") != GRAPH_SCHEMA:
            raise ValueError("unsupported ligand-bank schema")
        expected_graph = (ATOM_FEAT_DIM, BOND_FEAT_DIM, MAX_ATOMS)
        actual_graph = (lm.get("atom_feature_dim"), lm.get("bond_feature_dim"), lm.get("max_atoms"))
        if actual_graph != expected_graph:
            raise ValueError(f"ligand graph contract mismatch: {actual_graph}")
        protein_keys = {row["sequence_sha256"] for row in proteins}
        ligand_keys = {row["drug_key"] for row in ligands}
        if len(protein_keys) != self.manifest["targets"] or len(ligand_keys) != self.manifest["ligands"]:
            raise ValueError("corpus key counts disagree with its manifest")
        if not protein_keys.issubset(self.protein_bank.index):
            raise ValueError("protein bank does not cover the corpus")
        if not ligand_keys.issubset(self.ligand_bank.index):
            raise ValueError("ligand bank does not cover the corpus")

    def _build_tasks(self):
        tasks: dict[str, dict[str, np.ndarray]] = defaultdict(dict)
        grouped: dict[tuple[str, str], list[int]] = defaultdict(list)
        target_component: dict[str, str] = {}
        split_components: dict[str, set[str]] = defaultdict(set)
        for index, cell in enumerate(self.cells):
            split, target, component = cell["split"], cell["target_id"], cell["protein_group_40"]
            previous = target_component.setdefault(target, component)
            if previous != component:
                raise ValueError("one target belongs to multiple components")
            grouped[(split, target)].append(index)
            split_components[split].add(component)
        split_names = ("meta_train", "meta_val", "meta_test")
        for left, right in ((0, 1), (0, 2), (1, 2)):
            if split_components[split_names[left]] & split_components[split_names[right]]:
                raise ValueError("homology component crosses a hard split")
        components: dict[str, dict[str, tuple[str, ...]]] = {}
        for (split, target), indices in grouped.items():
            tasks[split][target] = np.asarray(indices, dtype=np.int64)
        for split, split_tasks in tasks.items():
            by_component: dict[str, list[str]] = defaultdict(list)
            for target in split_tasks:
                by_component[target_component[target]].append(target)
            components[split] = {key: tuple(sorted(value)) for key, value in by_component.items()}
        return dict(tasks), components

    def draw_episode(self, split: str, support_size: int, query_size: int,
                     rng: np.random.Generator) -> EpisodeSpec:
        eligible = {
            component: tuple(target for target in targets
                             if len(self.tasks[split][target]) >= support_size + 1)
            for component, targets in self.components[split].items()
        }
        eligible = {key: value for key, value in eligible.items() if value}
        if not eligible:
            raise ValueError("no task can provide the requested episode")
        component = sorted(eligible)[int(rng.integers(len(eligible)))]
        targets = eligible[component]
        target = targets[int(rng.integers(len(targets)))]
        indices = self.tasks[split][target]
        order = rng.permutation(indices)
        support = order[:support_size]
        query = order[support_size:support_size + min(query_size, len(order) - support_size)]
        donor_components = sorted(set(eligible).difference({component}))
        donor_component = donor_components[int(rng.integers(len(donor_components)))]
        donor_targets = eligible[donor_component]
        donor = donor_targets[int(rng.integers(len(donor_targets)))]
        return EpisodeSpec(split, component, target, tuple(map(int, support)),
                           tuple(map(int, query)), donor)

    def fixed_episode_bank(self, split: str, support_size: int, query_size: int,
                           draws: int, seed: int,
                           max_targets_per_component: int | None = None) -> tuple[EpisodeSpec, ...]:
        episodes = []
        eligible = {target: indices for target, indices in self.tasks[split].items()
                    if len(indices) >= support_size + 1}
        by_component = {
            component: tuple(target for target in targets if target in eligible)
            for component, targets in self.components[split].items()
        }
        by_component = {key: value for key, value in by_component.items() if value}
        if len(by_component) < 2:
            raise ValueError("fixed controls require at least two eligible components")
        component_order = sorted(by_component)
        for component_index, component in enumerate(component_order):
            donor_component = component_order[(component_index + 1) % len(component_order)]
            targets = by_component[component]
            if max_targets_per_component is not None:
                if max_targets_per_component < 1:
                    raise ValueError("max_targets_per_component must be positive")
                order = np.random.default_rng(
                    stable_seed("target-bank", seed, split, component)).permutation(len(targets))
                targets = tuple(targets[int(index)] for index in order[:max_targets_per_component])
            for target in targets:
                for draw in range(draws):
                    rng = np.random.default_rng(stable_seed("episode", seed, split, target, draw))
                    indices = eligible[target]
                    order = rng.permutation(indices)
                    support = order[:support_size]
                    query = order[support_size:support_size + min(query_size, len(order) - support_size)]
                    donor_targets = by_component[donor_component]
                    donor = donor_targets[int(rng.integers(len(donor_targets)))]
                    episodes.append(EpisodeSpec(
                        split, component, target, tuple(map(int, support)),
                        tuple(map(int, query)), donor))
        return tuple(episodes)

    def materialize(self, spec: EpisodeSpec) -> EpisodeBatch:
        if set(spec.support) & set(spec.query):
            raise ValueError("support and query cells overlap")
        rows = [self.cells[index] for index in (*spec.support, *spec.query)]
        if any(row["target_id"] != spec.target or row["split"] != spec.split for row in rows):
            raise ValueError("episode cells do not match the declared task")
        pooled, residues, protein_mask = self.protein_bank.get(spec.target)

        def graphs(indices: tuple[int, ...]):
            values = []
            for index in indices:
                ligand = self.cells[index]["ligand_id"]
                if ligand not in self.graph_cache:
                    graph = featurize_smiles(self.ligand_smiles[ligand])
                    self.graph_cache[ligand] = graph["X"], graph["A"], graph["mask"]
                values.append(self.graph_cache[ligand])
            return tuple(torch.from_numpy(np.stack(field).copy()) for field in zip(*values))

        support_atoms, support_bonds, support_mask = graphs(spec.support)
        query_atoms, query_bonds, query_mask = graphs(spec.query)
        return EpisodeBatch(
            spec, torch.from_numpy(pooled.copy()), torch.from_numpy(residues.copy()),
            torch.from_numpy(protein_mask.copy()), support_atoms, support_bonds, support_mask,
            torch.tensor([self.cells[i]["pK"] for i in spec.support], dtype=torch.float32),
            query_atoms, query_bonds, query_mask,
            torch.tensor([self.cells[i]["pK"] for i in spec.query], dtype=torch.float32))

    def protein_for_target(self, target: str) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        index = int(self.tasks[next(split for split in self.tasks if target in self.tasks[split])][target][0])
        values = self.protein_bank.get(target)
        return tuple(torch.from_numpy(value.copy()) for value in values)
