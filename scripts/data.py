"""Runtime episode assembly over one physically isolated label view."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch

from contracts.ligand_graph import MAX_ATOMS
from model.runtime import require_cuda
from scripts.build_ligand_bank import load_ligand_bank
from scripts.seal_compiled_dataset import SealedCompiledDataset


DTYPE = torch.float64


@dataclass
class ModelEpisode:
    target_idx: torch.Tensor
    query_pair_idx: torch.Tensor
    support_pair_idx: torch.Tensor
    support_y: torch.Tensor
    support_mask: torch.Tensor
    context_id: torch.Tensor
    context_cont: torch.Tensor
    context_mask: torch.Tensor
    deranged_target_idx: torch.Tensor | None = None

    def __len__(self) -> int:
        return int(self.target_idx.shape[0])


@dataclass
class EpisodeTarget:
    query_y: torch.Tensor


class Banks:
    """Immutable source-or-metaval graph and protein tensors on CUDA."""

    def __init__(self, X, A, m, prot_pooled, prot_bank, device):
        self.X, self.A, self.m = X, A, m
        self.prot_pooled, self.prot_bank = prot_pooled, prot_bank
        self.device = device

    def drug_batch(self, index):
        return self.X[index], self.A[index], self.m[index]

    def protein_batch(self, index):
        return self.prot_pooled[index], self.prot_bank[index]


def _cache_parts(value):
    if isinstance(value, dict):
        return value["pooled"], value["residues"]
    if isinstance(value, tuple) and len(value) == 2:
        return value
    raise TypeError("protein cache entries must be (pooled, residues) or matching dictionaries")


class CompiledEpisodes:
    """Sample disjoint support/query episodes from one sealed view only."""

    def __init__(self, name: str | None = None, *, protein_cache_path: str | Path,
                 device=None, max_atoms: int = MAX_ATOMS, dtype=DTYPE,
                 sealed_dir: str | Path | None = None,
                 ligand_bank_path: str | Path | None = None,
                 visible_splits: tuple[str, ...] = ("source",),
                 normalization_bounds: tuple[float, float] | None = None):
        if sealed_dir is None:
            raise ValueError("sealed_dir is required; raw or v1 seals are not runtime inputs")
        if len(visible_splits) != 1:
            raise ValueError("each runtime must mount exactly one label split")
        self.name = name
        self.device = require_cuda(device)
        self.visible_splits = tuple(visible_splits)
        self.sealed = SealedCompiledDataset(sealed_dir, self.visible_splits[0], name)
        self.label_split = self.sealed.label_split
        self.rows = self.sealed.rows
        self._validate_view_rows()

        fitted = self.sealed.manifest.get("normalization") or {}
        self.eta_lo = float(fitted.get("low", 0.0))
        self.eta_hi = float(fitted.get("high", 1.0))
        if not self.eta_hi > self.eta_lo:
            raise ValueError("seal is missing valid source-fitted normalization bounds")
        if normalization_bounds is not None and not np.allclose(
                normalization_bounds, (self.eta_lo, self.eta_hi), rtol=0, atol=1e-12):
            raise ValueError("checkpoint normalization bounds do not match the sealed dataset")

        if max_atoms != MAX_ATOMS:
            raise ValueError("runtime graph width is fixed by the ligand graph contract")
        self.drug_ids = sorted({row["drug_key"] for row in self.rows})
        self.didx = {key: index for index, key in enumerate(self.drug_ids)}
        ligand_bank_path = Path(ligand_bank_path) if ligand_bank_path is not None else \
            Path(sealed_dir) / f"{self.label_split}_ligand_bank"
        graph_by_key = load_ligand_bank(ligand_bank_path)
        if set(graph_by_key) != set(self.drug_ids):
            raise ValueError("ligand bank keys must exactly match the mounted label view")
        graphs = [graph_by_key[key] for key in self.drug_ids]
        self._load_protein_bank(protein_cache_path, dtype)
        X = torch.as_tensor(np.stack([graph["X"] for graph in graphs]), dtype=dtype, device=self.device)
        A = torch.as_tensor(np.stack([graph["A"] for graph in graphs]), dtype=dtype, device=self.device)
        m = torch.as_tensor(np.stack([graph["mask"] for graph in graphs]), dtype=dtype, device=self.device)
        self.banks = Banks(X, A, m, self.prot_pooled, self.prot_bank, self.device)

        self.targets: dict[str, dict] = {}
        for row in self.rows:
            task = self.targets.setdefault(row["task_key"], {
                "split": self.label_split, "target_key": row["target_key"],
                "drug_rows": [], "row_ids": [], "y": [], "context_id": [],
                "context_cont": [], "context_mask": [],
            })
            task["drug_rows"].append(self.didx[row["drug_key"]])
            task["row_ids"].append(row["row_id"])
            task["y"].append(float(row["y"]))
            task["context_id"].append(int(row["context_id"]))
            task["context_cont"].append(tuple(float(v) for v in row["context_cont"]))
            task["context_mask"].append(tuple(int(v) for v in row["context_mask"]))
        for task in self.targets.values():
            task["drug_rows"] = np.asarray(task["drug_rows"], dtype=np.int64)
            task["row_ids"] = np.asarray(task["row_ids"], dtype=str)
            task["y"] = np.asarray(task["y"], dtype=float)
            for field in ("context_id", "context_cont", "context_mask"):
                values = task.pop(field)
                if any(value != values[0] for value in values[1:]):
                    raise ValueError(f"a task mixes {field} specifications")
                task[field] = values[0]
        self.by_split = {
            split: (sorted(self.targets) if split == self.label_split else [])
            for split in ("source", "metaval", "recipient")
        }

    def _validate_view_rows(self) -> None:
        seen = set()
        for row in self.rows:
            if row["row_id"] in seen:
                raise ValueError("duplicate row_id inside sealed view")
            if not 0.0 <= float(row["y"]) <= 1.0:
                raise ValueError("sealed label is outside [0,1]")
            if len(row["context_cont"]) != len(row["context_mask"]):
                raise ValueError("sealed continuous context and mask widths differ")
            if any(float(value) < 0.0 or float(value) > 1.0 for value in row["context_cont"]):
                raise ValueError("sealed continuous context is outside [0,1]")
            if any(int(value) not in {0, 1} for value in row["context_mask"]):
                raise ValueError("sealed context mask is not binary")
            seen.add(row["row_id"])

    def _load_protein_bank(self, path: str | Path, dtype) -> None:
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(path)
        cache = torch.load(path, map_location="cpu", weights_only=False)
        expected = sorted({row["target_key"] for row in self.rows})
        if set(cache) != set(expected):
            raise ValueError(
                "protein cache keys must exactly match the mounted label view; "
                "build a separate cache for source and metaval"
            )
        self.target_ids = expected
        self.tidx = {key: index for index, key in enumerate(self.target_ids)}
        parts = [_cache_parts(cache[key]) for key in self.target_ids]
        pooled, residues = zip(*parts)
        if any(value.ndim != 1 for value in pooled) or any(value.ndim != 2 for value in residues):
            raise ValueError("protein cache has invalid pooled or residue tensor ranks")
        if len({tuple(value.shape) for value in pooled}) != 1 or \
                len({tuple(value.shape) for value in residues}) != 1:
            raise ValueError("protein cache entries must have uniform shapes")
        self.prot_pooled = torch.stack(list(pooled)).to(device=self.device, dtype=dtype)
        self.prot_bank = torch.stack(list(residues)).to(device=self.device, dtype=dtype)

    def seal_audit(self) -> dict:
        return self.sealed.audit_snapshot()

    def _require_split(self, split: str, purpose: str) -> None:
        if split != self.label_split:
            raise PermissionError(f"{split} labels are not mounted in this runtime")
        if purpose in {"training", "hyperparameter"} and split != "source":
            raise PermissionError(f"{purpose} is source-only")

    def _make(self, t_idx, q_idx, s_idx, sup_y, mask, labels,
              context_id, context_cont, context_mask):
        tensor = lambda value, dtype=DTYPE: torch.as_tensor(
            np.asarray(value), dtype=dtype, device=self.device)
        index = lambda value: torch.as_tensor(
            np.asarray(value), dtype=torch.long, device=self.device)
        episode = ModelEpisode(
            index(t_idx), index(q_idx), index(s_idx), tensor(sup_y), tensor(mask),
            index(context_id), tensor(context_cont), tensor(context_mask),
        )
        return episode, EpisodeTarget(tensor(labels))

    def sample(self, n: int, rng: np.random.Generator, split="source", ks=(1, 3, 5)):
        self._require_split(split, "training")
        if n < 1:
            raise ValueError("sample count must be positive")
        eligible = [task for task in self.by_split[split] if len(self.targets[task]["y"]) >= 2]
        if not eligible:
            raise ValueError("mounted label view has no task with disjoint support/query rows")
        kmax = max(ks)
        t_idx, q_idx, s_idx, sup_y, mask, labels = [], [], [], [], [], []
        context_id, context_cont, context_mask = [], [], []
        for task_key in rng.choice(eligible, size=n, replace=True):
            task = self.targets[str(task_key)]
            valid_ks = [k for k in ks if k < len(task["y"])]
            k = int(rng.choice(valid_ks))
            chosen = rng.choice(len(task["y"]), size=k + 1, replace=False)
            query, support = int(chosen[0]), chosen[1:]
            padded = np.pad(support, (0, kmax - k), constant_values=int(support[0]))
            t_idx.append(self.tidx[task["target_key"]])
            q_idx.append(task["drug_rows"][query])
            s_idx.append(task["drug_rows"][padded])
            sup_y.append(task["y"][padded])
            mask.append([1.0] * k + [0.0] * (kmax - k))
            labels.append(task["y"][query])
            context_id.append(task["context_id"])
            context_cont.append(task["context_cont"])
            context_mask.append(task["context_mask"])
        return self._make(t_idx, q_idx, s_idx, sup_y, mask, labels,
                          context_id, context_cont, context_mask)

    def fixed_support_tasks(self, task_key: str, k: int, rng: np.random.Generator,
                            mode="correct", max_q: int | None = None):
        if hasattr(self, "label_split"):
            self._require_split(self.label_split, "validation")
        task = self.targets[task_key]
        n_rows = len(task["y"])
        if n_rows <= k + 1:
            return None
        support = rng.choice(n_rows, size=k, replace=False)
        query = np.setdiff1d(np.arange(n_rows), support, assume_unique=False)
        if max_q is not None and len(query) > max_q:
            query = rng.choice(query, size=max_q, replace=False)
        support_y = task["y"][support]
        if mode == "permuted" and k > 1:
            support_y = support_y[np.roll(np.arange(k), 1)]
        elif mode != "correct":
            raise ValueError(f"unsupported fixed-support mode {mode!r}")
        args = (
            np.full(len(query), self.tidx[task.get("target_key", task_key)], dtype=np.int64),
            task["drug_rows"][query], np.tile(task["drug_rows"][support], (len(query), 1)),
            np.tile(support_y, (len(query), 1)), np.ones((len(query), k)), task["y"][query],
            np.full(len(query), task["context_id"], dtype=np.int64),
            np.tile(np.asarray(task["context_cont"], dtype=float), (len(query), 1)),
            np.tile(np.asarray(task["context_mask"], dtype=float), (len(query), 1)),
        )
        return self._make(*args)


DTAEpisodes = CompiledEpisodes
MetaSieveBatch = ModelEpisode


def load_episodes(dataset: str | None, *, device=None, sealed_dir: str | Path | None = None,
                  protein_cache_path: str | Path | None = None,
                  ligand_bank_path: str | Path | None = None,
                  visible_splits: tuple[str, ...] = ("source",),
                  normalization_bounds: tuple[float, float] | None = None) -> CompiledEpisodes:
    if protein_cache_path is None:
        raise ValueError("protein_cache_path is required for the mounted label view")
    return CompiledEpisodes(
        dataset, device=device, sealed_dir=sealed_dir, protein_cache_path=protein_cache_path,
        ligand_bank_path=ligand_bank_path,
        visible_splits=visible_splits, normalization_bounds=normalization_bounds,
    )
