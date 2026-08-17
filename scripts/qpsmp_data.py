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
    protein_chemistry: torch.Tensor | None = None
    # Morgan fingerprints of the support/query ligands. These are an ordinary
    # 2D ligand descriptor, not an extra label source, and they carry the
    # chemical-similarity structure that learned kernel keys failed to find.
    support_fingerprint: torch.Tensor | None = None
    query_fingerprint: torch.Tensor | None = None
    # Optional common-frame Cartesian inputs.  Coordinates contain padded
    # ligand atoms followed by protein residues for each support/query pair.
    # Edges use the packed flattened indexing consumed by model.cartesian.
    support_coordinates: torch.Tensor | None = None
    query_coordinates: torch.Tensor | None = None
    geometry_edge_index: torch.Tensor | None = None
    geometry_available: torch.Tensor | None = None
    geometry_common_frame: torch.Tensor | None = None


def protein_slot_chemistry(sequence: str, slot_count: int) -> torch.Tensor:
    """Aggregate governed amino-acid roles into the protein-bank slot grid."""
    roles = (
        set("KRH"), set("DE"), set("FWY"), set("AVILMFWY"),
    )
    values = torch.zeros(slot_count, len(roles), dtype=torch.float32)
    length = len(sequence)
    for slot in range(slot_count):
        start = slot * length // slot_count
        stop = (slot + 1) * length // slot_count
        if stop <= start:
            continue
        fragment = sequence[start:stop]
        values[slot] = torch.tensor([
            sum(residue in role for residue in fragment) / len(fragment)
            for role in roles
        ])
    return values


def stable_seed(*parts: object) -> int:
    raw = "|".join(map(str, parts)).encode("utf-8")
    return int(hashlib.sha256(raw).hexdigest()[:16], 16)


class _ShardedBank:
    def __init__(self, directory: Path, fields: tuple[str, ...], cache_size: int | None = None):
        self.directory = directory
        self.fields = fields
        self.index: dict[str, tuple[Path, int]] = {}
        self.cache: OrderedDict[Path, dict[str, np.ndarray]] = OrderedDict()
        manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
        self.manifest = manifest
        shards = manifest.get("shards", [])
        self.cache_size = len(shards) if cache_size is None else cache_size
        for item in shards:
            name = item["path"] if isinstance(item, dict) else item
            path = directory / name
            with np.load(path, allow_pickle=False) as stored:  # noqa: SIM117
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


class _CompactLigandBank:
    FIELDS = ("sizes", "atom_offsets", "bond_offsets", "X", "A")

    def __init__(self, directory: Path, cache_size: int | None = None):
        manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
        if manifest.get("schema") != "MetaSieve.QPSMPCompactLigandBank.v1":
            raise ValueError("unsupported compact ligand-bank schema")
        self.manifest = manifest
        self.index: dict[str, tuple[Path, int]] = {}
        self.cache: OrderedDict[Path, dict[str, np.ndarray]] = OrderedDict()
        # An episode draws ligands uniformly, so it touches most shards. A
        # single-shard cache reloads and decompresses ~0.6 MB npz archives
        # dozens of times per episode; resident shards remove that entirely.
        self.cache_size = (len(manifest["shards"]) if cache_size is None
                           else max(1, int(cache_size)))
        for item in manifest["shards"]:
            path = directory / item["path"]
            with np.load(path, allow_pickle=False) as stored:
                for local, key in enumerate(stored["keys"].astype(str)):
                    self.index[key] = path, local

    def __len__(self) -> int:
        return len(self.index)

    def get(self, key: str) -> tuple[np.ndarray, ...]:
        path, local = self.index[key]
        if path not in self.cache:
            with np.load(path, allow_pickle=False) as stored:
                self.cache[path] = {name: stored[name] for name in self.FIELDS}
        self.cache.move_to_end(path)
        while len(self.cache) > self.cache_size:
            self.cache.popitem(last=False)
        values = self.cache[path]
        size = int(values["sizes"][local])
        atom_start, atom_stop = values["atom_offsets"][local:local + 2]
        bond_start, bond_stop = values["bond_offsets"][local:local + 2]
        atoms = values["X"][atom_start:atom_stop]
        bonds = values["A"][bond_start:bond_stop].reshape(size, size, -1)
        return atoms, bonds, np.ones(size, dtype=np.float32)


class QPSMPData:
    """Read-only main-v0 corpus with lazy, schema-checked biological banks.

    The class has **two label surfaces**, and `seal_record()` always reports
    which one was mounted.

    **Default surface — the all-label corpus.** `include_meta_test` provides
    **logical exclusion after parsing**, which is strictly weaker than a
    physical label seal and must never be described as one. `cells.jsonl.gz` is
    a single all-label artifact: this constructor decompresses and parses every
    row, including every `meta_test` label, and then discards the sealed ones.
    The labels are in process memory for the duration of `_governed_cells`.
    Every recorded R0-R14 and Stage A-Q artifact was produced on this surface,
    which is why it is kept unchanged.

    What that contract does guarantee, and what it does not:

    * ✅ **fail-closed** — the default excludes, so a caller cannot admit the
      sealed split by omission;
    * ✅ **unreachable after construction** — `cells` is never *bound* to a list
      containing a sealed row, so no task index, component map or
      `materialize` call can reach one;
    * ✅ **authorized** — admitting the split requires a written reason in
      `meta_test_authorization`; a bare `include_meta_test=True` raises;
    * ❌ **not physically isolated** — the sealed labels are read off disk and
      parsed on every construction. A process that never asks for them still
      holds them transiently.

    **Isolated surface — `split_view`.** Passing the governed split-view
    directory built by `scripts/build_governed_split_views.py` mounts separately
    hashed per-split label artifacts instead. `cells.jsonl.gz` is then never
    opened, and the `meta_test` label artifact lives outside the development
    tree entirely, so a development process has no filesystem path to it. On
    this surface `seal_record()["isolation"]["physically_isolated"]` is `True`.
    The row order is restored from the view's identity-only `governance.jsonl`,
    so `cells` is element-for-element identical to the corpus construction and
    no recorded episode index changes. Opening `meta_test` from this surface
    additionally requires an explicit out-of-tree `sealed_meta_test_directory`.

    `seal_record()` reports the constructed object's *actual* state and is the
    only sanctioned source for the `meta_test` block of a result artifact.
    Hand-written seal strings caused the incident this contract repairs; see
    `tools/research/a2_readiness_v2/GOVERNANCE_INCIDENT.md`.
    """

    def __init__(self, corpus: Path, protein_bank: Path, ligand_bank: Path,
                 compact_ligand_bank: Path | None = None,
                 split_directory: Path | None = None,
                 include_meta_test: bool = False,
                 meta_test_authorization: str | None = None,
                 split_view: Path | None = None,
                 sealed_meta_test_directory: Path | None = None):
        if include_meta_test and not (meta_test_authorization or "").strip():
            raise ValueError(
                "including meta_test requires an explicit written "
                "meta_test_authorization; the sealed confirmation population "
                "cannot be opened by a bare flag (contract 2026-08-16)")
        if meta_test_authorization and not include_meta_test:
            raise ValueError(
                "meta_test_authorization was supplied without "
                "include_meta_test=True; refusing an ambiguous seal state")
        if sealed_meta_test_directory is not None and split_view is None:
            raise ValueError(
                "sealed_meta_test_directory is only meaningful with a "
                "split_view; the all-label corpus has no separate meta_test "
                "artifact")
        self.include_meta_test = bool(include_meta_test)
        self.meta_test_authorization = (
            meta_test_authorization if self.include_meta_test else None)
        self.corpus = corpus
        if split_directory is not None:
            split_directory = Path(split_directory)
        manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
        if manifest.get("schema") != "MetaSieve.MetaFewshot.MainV0Corpus.v1":
            raise ValueError("unsupported corpus schema")
        if manifest.get("training_authorized") is not True:
            raise ValueError("corpus manifest does not authorize training")
        self.manifest = manifest
        # An alternative governed split replaces `split` on every cell and drops
        # the cells it does not assign, so downstream indices stay consistent.
        # The corpus itself is never modified.
        self.split_manifest = None
        assignment: dict[str, str] | None = None
        if split_directory is not None:
            self.split_manifest = json.loads(
                (split_directory / "manifest.json").read_text(encoding="utf-8"))
            if self.split_manifest.get("schema") != "MetaSieve.DoubleColdSplit.v1":
                raise ValueError("unsupported governed-split schema")
            assignment = json.loads(
                (split_directory / "assignment.json").read_text(encoding="utf-8"))
            recorded = self.split_manifest.get("assignment_sha256")
            actual = hashlib.sha256(
                json.dumps(assignment, sort_keys=True).encode("utf-8")).hexdigest()
            if recorded != actual:
                raise ValueError(
                    f"split assignment hash mismatch: {actual} != {recorded}")
        self.split_view = None
        if split_view is not None:
            # Isolated surface: mount exactly the per-split label artifacts
            # named. `cells.jsonl.gz` is never opened, so no meta_test label is
            # decompressed or parsed by a development construction.
            self.split_view = self._mount_split_view(
                Path(split_view), sealed_meta_test_directory, assignment)
            self.cells = list(self.split_view.cells)
            self.sealed_cell_count = self.split_view.sealed_cell_count
        else:
            # The seal is applied in the same pass that assigns the split, so
            # `self.cells` is never bound to a list containing a sealed cell.
            # Every downstream structure — task indices, component maps,
            # `materialize` — is built from this attribute, so none of them can
            # reach meta_test while the seal holds. The sealed rows are still
            # *parsed* out of the gzip stream here; the defensible claim is that
            # they are discarded during construction and are unreachable
            # afterwards, not that the file is never decompressed.
            self.cells = self._governed_cells(
                corpus / "cells.jsonl.gz", assignment, self.include_meta_test)
            self.sealed_cell_count = self._sealed_cell_count
        proteins = self._read_jsonl(corpus / "proteins.jsonl")
        self._protein_sequences = {
            row["sequence_sha256"]: row["sequence"] for row in proteins}
        ligands = self._read_jsonl(corpus / "ligands.jsonl")
        protein_keys = {row["sequence_sha256"] for row in proteins}
        self.protein_bank = _ShardedBank(
            protein_bank, ("pooled", "residues", "mask"), cache_size=4)
        self.ligand_bank = (
            _CompactLigandBank(compact_ligand_bank)
            if compact_ligand_bank is not None and compact_ligand_bank.exists()
            else _ShardedBank(ligand_bank, ("X", "A", "mask"), cache_size=1))
        self._protein_tensors: dict[str, tuple[torch.Tensor, ...]] = {}
        self._protein_chemistry: dict[str, torch.Tensor] = {}
        self._ligand_tensors: OrderedDict[str, tuple[torch.Tensor, ...]] = OrderedDict()
        self._ligand_tensor_cache_size = 128
        self._ligand_smiles = {row["drug_key"]: row.get("smiles")
                               for row in ligands}
        self._fingerprints: dict[str, torch.Tensor] | None = None
        self.fingerprint_bits = 1024
        self.fingerprint_radius = 2
        self._validate_manifests(proteins, ligands)
        self.tasks, self.components = self._build_tasks()
        target_keys = {cell["target_id"] for cell in self.cells}
        if self.split_manifest is None and self.split_view is None:
            if target_keys != protein_keys:
                raise ValueError(
                    "target IDs do not equal the governed protein-bank keys")
        elif not target_keys.issubset(protein_keys):
            raise ValueError("governed split references targets outside the bank")

    @staticmethod
    def _read_jsonl(path: Path) -> list[dict]:
        with path.open(encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]

    @staticmethod
    def _read_cells(path: Path) -> list[dict]:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]

    def _governed_cells(self, path: Path, assignment: dict[str, str] | None,
                        include_meta_test: bool) -> list[dict]:
        """Stream the corpus, apply the governed split, and drop sealed cells.

        Filtering happens here rather than in a second pass so that no attribute
        of this object is ever bound to a list containing a `meta_test` cell.
        """
        kept: list[dict] = []
        sealed = 0
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                cell = json.loads(line)
                if assignment is not None:
                    if cell["cell_id"] not in assignment:
                        continue
                    cell["split"] = assignment[cell["cell_id"]]
                if cell["split"] == "meta_test":
                    sealed += 1
                    if not include_meta_test:
                        continue
                kept.append(cell)
        self._sealed_cell_count = sealed
        return kept

    def _mount_split_view(self, directory: Path,
                          sealed_directory: Path | None,
                          assignment: dict[str, str] | None):
        """Mount the governed split view and bind it to this corpus and split.

        The view carries the source-corpus manifest hash and the frozen split
        assignment hash, so provenance is verified from the manifest rather than
        from a re-read of the all-label artifact.
        """
        from scripts.build_governed_split_views import (  # local: optional path
            DEVELOPMENT_SPLITS, GovernedSplitView, SEALED_SPLIT,
        )

        visible = ((SEALED_SPLIT,) if self.include_meta_test
                   else DEVELOPMENT_SPLITS)
        view = GovernedSplitView(
            directory, visible=visible,
            authorization=self.meta_test_authorization,
            sealed_directory=sealed_directory)
        corpus_manifest_sha256 = hashlib.sha256(
            (self.corpus / "manifest.json").read_bytes()).hexdigest()
        if view.manifest["source_manifest_sha256"] != corpus_manifest_sha256:
            raise ValueError(
                "split view was cut from a different corpus: "
                f"{view.manifest['source_manifest_sha256']} != "
                f"{corpus_manifest_sha256}")
        if assignment is not None:
            recorded = hashlib.sha256(
                json.dumps(assignment, sort_keys=True).encode("utf-8")).hexdigest()
            if view.manifest["split_assignment_sha256"] != recorded:
                raise ValueError(
                    "split view does not match the governed split assignment: "
                    f"{view.manifest['split_assignment_sha256']} != {recorded}")
        return view

    def seal_record(self, evaluated: bool | None = None) -> dict:
        """The truthful `meta_test` block for a result artifact.

        Artifacts must call this instead of writing a seal string by hand: the
        2026-08-16 incident was a hand-written claim that contradicted the
        constructor arguments actually used.

        While the seal holds, `evaluated=False` is a provable consequence of the
        data being absent. Once the split is open the object cannot know what
        the caller did with it, so `evaluated` becomes mandatory.
        """
        if self.include_meta_test and evaluated is None:
            raise ValueError(
                "an open meta_test run must state evaluated=True/False "
                "explicitly; the dataset cannot infer it")
        if self.split_view is not None:
            return self._isolated_seal_record(evaluated)
        return {
            "evaluated": bool(evaluated) if self.include_meta_test else False,
            "included": self.include_meta_test,
            "authorization": self.meta_test_authorization,
            "sealed_cells_withheld": (
                0 if self.include_meta_test else self.sealed_cell_count),
            # Three distinct states, never collapsed into one word:
            #   parsed_in_process   — the labels were decompressed and read
            #   metric_unconsumed   — no sealed value entered any computation
            #   physically_isolated — the labels were never on this process's
            #                         disk path at all
            # Only the first two are true today.
            "isolation": {
                "level": ("open" if self.include_meta_test
                          else "logical_exclusion_after_parsing"),
                "labels_parsed_in_process": True,
                "physically_isolated": False,
                "why_not": ("cells.jsonl.gz is a single all-label artifact; "
                            "every meta_test label is decompressed and parsed "
                            "on construction before the sealed rows are "
                            "discarded"),
                "specification_for_a_real_seal": (
                    "tools/research/a2_readiness_v2/SPLIT_ISOLATION_SPEC.md "
                    "(not authorized, not implemented)"),
            },
            "seal": ("OPEN: QPSMPData(include_meta_test=True) — "
                     f"authorization: {self.meta_test_authorization}"
                     if self.include_meta_test else
                     "logical exclusion after parsing: "
                     "QPSMPData(include_meta_test=False). Sealed rows are "
                     "parsed from the all-label corpus and discarded during "
                     "construction; they are unreachable from cells, tasks, "
                     "components and materialize afterwards. This is NOT a "
                     "physical label seal."),
        }

    def _isolated_seal_record(self, evaluated: bool | None) -> dict:
        """The `meta_test` block for a construction on the isolated surface.

        The difference from the default surface is one provable fact: this
        process never opened the `meta_test` label artifact, because the
        artifact is not in the mounted directory at all. The wording stays
        specific — "physically isolated (governed split view)" — so a reader
        can tell which surface produced any given artifact.
        """
        isolation = self.split_view.isolation_record()
        return {
            "evaluated": bool(evaluated) if self.include_meta_test else False,
            "included": self.include_meta_test,
            "authorization": self.meta_test_authorization,
            "sealed_cells_withheld": self.sealed_cell_count,
            "isolation": {
                "level": isolation["level"],
                "labels_parsed_in_process": isolation["labels_parsed_in_process"],
                "physically_isolated": isolation["physically_isolated"],
                "surface": "governed_split_view",
                "view_directory": str(self.split_view.directory),
                "source_manifest_sha256": isolation["source_manifest_sha256"],
                "split_assignment_sha256": isolation["split_assignment_sha256"],
                "builder": "scripts/build_governed_split_views.py",
                "specification": ("tools/research/a2_readiness_v2/"
                                  "SPLIT_ISOLATION_SPEC.md"),
            },
            "seal": ("OPEN: authorized meta_test mount from an out-of-tree "
                     f"artifact — authorization: {self.meta_test_authorization}"
                     if self.include_meta_test else
                     "physically isolated (governed split view): this process "
                     "mounted only the meta_train/meta_val label artifacts. "
                     "The meta_test label artifact is not present in the "
                     "development surface, so no sealed label was opened, "
                     "decompressed or parsed."),
        }

    def _validate_manifests(self, proteins: list[dict], ligands: list[dict]) -> None:
        pm, lm = self.protein_bank.manifest, self.ligand_bank.manifest
        if pm.get("schema") != "MetaSieve.StructureProteinBank.v1":
            raise ValueError("unsupported protein-bank schema")
        compact = lm.get("schema") == "MetaSieve.QPSMPCompactLigandBank.v1"
        if not compact and lm.get("schema") != GRAPH_SCHEMA:
            raise ValueError("unsupported ligand-bank schema")
        expected_graph = (ATOM_FEAT_DIM, BOND_FEAT_DIM)
        actual_graph = (lm.get("atom_feature_dim"), lm.get("bond_feature_dim"))
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

    def _unique_ligand_order(self, indices: np.ndarray,
                             rng: np.random.Generator) -> np.ndarray:
        """Randomize cells, retaining at most one observation per ligand."""
        ordered = []
        seen = set()
        for index in rng.permutation(indices):
            ligand = self.cells[int(index)]["ligand_id"]
            if ligand in seen:
                continue
            seen.add(ligand)
            ordered.append(int(index))
        return np.asarray(ordered, dtype=np.int64)

    def _unique_ligand_count(self, indices: np.ndarray) -> int:
        return len({self.cells[int(index)]["ligand_id"] for index in indices})

    def draw_episode(self, split: str, support_size: int, query_size: int,
                     rng: np.random.Generator,
                     min_query_size: int = 1) -> EpisodeSpec:
        if support_size < 0:
            raise ValueError("support size cannot be negative")
        if min_query_size < 1:
            raise ValueError("an episode needs at least one query")
        # Within-target centered and ranking objectives are undefined on a
        # one-query episode, so a caller that uses them can require a wider
        # panel. The default reproduces the historical behaviour exactly.
        needed = support_size + min_query_size
        eligible = {
            component: tuple(target for target in targets
                             if self._unique_ligand_count(
                                 self.tasks[split][target]) >= needed)
            for component, targets in self.components[split].items()
        }
        eligible = {key: value for key, value in eligible.items() if value}
        if not eligible:
            raise ValueError("no task can provide the requested episode")
        component = sorted(eligible)[int(rng.integers(len(eligible)))]
        targets = eligible[component]
        target = targets[int(rng.integers(len(targets)))]
        indices = self.tasks[split][target]
        order = self._unique_ligand_order(indices, rng)
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
                    if self._unique_ligand_count(indices) >= support_size + 1}
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
                    order = self._unique_ligand_order(indices, rng)
                    support = order[:support_size]
                    query = order[support_size:support_size + min(query_size, len(order) - support_size)]
                    donor_targets = by_component[donor_component]
                    donor = donor_targets[int(rng.integers(len(donor_targets)))]
                    episodes.append(EpisodeSpec(
                        split, component, target, tuple(map(int, support)),
                        tuple(map(int, query)), donor))
        return tuple(episodes)

    def fixed_nested_episode_banks(
            self, split: str, support_sizes: tuple[int, ...], query_size: int,
            draws: int, seed: int,
            max_targets_per_component: int | None = None,
    ) -> dict[int, tuple[EpisodeSpec, ...]]:
        """Build nested support prefixes with one common query per target/draw."""
        sizes = tuple(sorted(set(map(int, support_sizes))))
        if not sizes or sizes[0] < 0:
            raise ValueError("support sizes cannot be negative")
        max_support = sizes[-1]
        eligible = {
            target: indices for target, indices in self.tasks[split].items()
            if self._unique_ligand_count(indices) >= max_support + 1
        }
        by_component = {
            component: tuple(target for target in targets if target in eligible)
            for component, targets in self.components[split].items()
        }
        by_component = {key: value for key, value in by_component.items() if value}
        if len(by_component) < 2:
            raise ValueError("nested controls require at least two eligible components")
        banks: dict[int, list[EpisodeSpec]] = {size: [] for size in sizes}
        component_order = sorted(by_component)
        for component_index, component in enumerate(component_order):
            donor_component = component_order[(component_index + 1) % len(component_order)]
            targets = by_component[component]
            if max_targets_per_component is not None:
                if max_targets_per_component < 1:
                    raise ValueError("max_targets_per_component must be positive")
                target_order = np.random.default_rng(stable_seed(
                    "nested-target-bank", seed, split, component)).permutation(len(targets))
                targets = tuple(targets[int(index)] for index in
                                target_order[:max_targets_per_component])
            for target in targets:
                for draw in range(draws):
                    rng = np.random.default_rng(stable_seed(
                        "nested-episode", seed, split, target, draw))
                    order = self._unique_ligand_order(eligible[target], rng)
                    query = tuple(map(int, order[
                        max_support:max_support + min(query_size, len(order) - max_support)]))
                    donor_targets = by_component[donor_component]
                    donor = donor_targets[int(rng.integers(len(donor_targets)))]
                    for size in sizes:
                        banks[size].append(EpisodeSpec(
                            split, component, target,
                            tuple(map(int, order[:size])), query, donor))
        return {size: tuple(episodes) for size, episodes in banks.items()}

    def materialize(self, spec: EpisodeSpec) -> EpisodeBatch:
        if set(spec.support) & set(spec.query):
            raise ValueError("support and query cells overlap")
        support_ligands = {self.cells[index]["ligand_id"] for index in spec.support}
        query_ligands = {self.cells[index]["ligand_id"] for index in spec.query}
        if support_ligands & query_ligands:
            raise ValueError("support and query ligand identities overlap")
        rows = [self.cells[index] for index in (*spec.support, *spec.query)]
        if any(row["target_id"] != spec.target or row["split"] != spec.split for row in rows):
            raise ValueError("episode cells do not match the declared task")
        if spec.target not in self._protein_tensors:
            self._protein_tensors[spec.target] = tuple(
                torch.from_numpy(value.copy()) for value in self.protein_bank.get(spec.target))
        pooled, residues, protein_mask = self._protein_tensors[spec.target]

        def graphs(indices: tuple[int, ...]):
            if not indices:
                atom_dim = int(self.ligand_bank.manifest["atom_feature_dim"])
                bond_dim = int(self.ligand_bank.manifest["bond_feature_dim"])
                return (torch.zeros(0, 1, atom_dim),
                        torch.zeros(0, 1, 1, bond_dim),
                        torch.zeros(0, 1, dtype=torch.bool))
            values = []
            for index in indices:
                ligand = self.cells[index]["ligand_id"]
                if ligand not in self._ligand_tensors:
                    self._ligand_tensors[ligand] = tuple(
                        torch.from_numpy(value.copy()) for value in self.ligand_bank.get(ligand))
                    while len(self._ligand_tensors) > self._ligand_tensor_cache_size:
                        self._ligand_tensors.popitem(last=False)
                self._ligand_tensors.move_to_end(ligand)
                values.append(self._ligand_tensors[ligand])
            max_atoms = max(value[0].shape[0] for value in values)
            atoms, bonds, masks = [], [], []
            for atom, bond, mask in values:
                missing = max_atoms - atom.shape[0]
                atoms.append(torch.nn.functional.pad(atom, (0, 0, 0, missing)))
                bonds.append(torch.nn.functional.pad(
                    bond, (0, 0, 0, missing, 0, missing)))
                masks.append(torch.nn.functional.pad(mask, (0, missing)))
            return torch.stack(atoms), torch.stack(bonds), torch.stack(masks)

        support_atoms, support_bonds, support_mask = graphs(spec.support)
        query_atoms, query_bonds, query_mask = graphs(spec.query)
        if spec.target not in self._protein_chemistry:
            self._protein_chemistry[spec.target] = protein_slot_chemistry(
                self._protein_sequences[spec.target], residues.shape[0])
        return EpisodeBatch(
            spec, pooled, residues, protein_mask, support_atoms, support_bonds, support_mask,
            torch.tensor([self.cells[i]["pK"] for i in spec.support], dtype=torch.float32),
            query_atoms, query_bonds, query_mask,
            torch.tensor([self.cells[i]["pK"] for i in spec.query], dtype=torch.float32),
            protein_chemistry=self._protein_chemistry[spec.target],
            support_fingerprint=self.fingerprint_rows(spec.support),
            query_fingerprint=self.fingerprint_rows(spec.query))

    @property
    def fingerprints(self) -> dict[str, torch.Tensor]:
        """Lazily built Morgan fingerprints, keyed by ligand id.

        Ligands whose SMILES cannot be parsed get an all-zero vector. That
        yields Tanimoto 0 against every other ligand, so such a support gets
        the *lowest* similarity logit — it is not excluded, because a softmax
        over the supports still assigns it finite weight. All 9,880 ligands in
        the active corpus parse, so no episode currently exercises this path;
        if a corpus with unparsable ligands is ever used, explicit masking
        should be added rather than relying on the zero vector.
        """
        if self._fingerprints is None:
            from rdkit import Chem, RDLogger
            from rdkit.Chem import rdFingerprintGenerator
            RDLogger.DisableLog("rdApp.*")
            generator = rdFingerprintGenerator.GetMorganGenerator(
                radius=self.fingerprint_radius, fpSize=self.fingerprint_bits)
            built: dict[str, torch.Tensor] = {}
            for key, smiles in self._ligand_smiles.items():
                vector = torch.zeros(self.fingerprint_bits, dtype=torch.float32)
                molecule = Chem.MolFromSmiles(smiles) if smiles else None
                if molecule is not None:
                    bits = list(generator.GetFingerprint(molecule).GetOnBits())
                    if bits:
                        vector[torch.tensor(bits, dtype=torch.long)] = 1.0
                built[key] = vector
            self._fingerprints = built
        return self._fingerprints

    def fingerprint_rows(self, indices: tuple[int, ...]) -> torch.Tensor:
        if not indices:
            return torch.zeros(0, self.fingerprint_bits, dtype=torch.float32)
        table = self.fingerprints
        return torch.stack([table[self.cells[i]["ligand_id"]] for i in indices])

    def protein_for_target(self, target: str) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        index = int(self.tasks[next(split for split in self.tasks if target in self.tasks[split])][target][0])
        if target not in self._protein_tensors:
            self._protein_tensors[target] = tuple(
                torch.from_numpy(value.copy()) for value in self.protein_bank.get(target))
        return self._protein_tensors[target]

    def protein_chemistry_for_target(self, target: str) -> torch.Tensor:
        if target not in self._protein_chemistry:
            residues = self.protein_for_target(target)[1]
            self._protein_chemistry[target] = protein_slot_chemistry(
                self._protein_sequences[target], residues.shape[0])
        return self._protein_chemistry[target]
