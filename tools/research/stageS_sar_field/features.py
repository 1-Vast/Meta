"""Feature gathering, protein donors and the two falsification permutations.

Everything here is deterministic and process-stable: the permutations are
derived from `scripts.qpsmp_data.stable_seed` (SHA-256), never from Python's
`hash()`.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np
import torch

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.qpsmp_data import QPSMPData, stable_seed
from tools.research.stageS_sar_field.pairs import (
    PairSpec, component_of_target, target_panel_documents,
)


class LigandFeatureStore:
    """Padded graph + fingerprint tensors for a set of ligands.

    Padding is done per batch to the batch maximum.  The encoder masks padded
    atoms and their (all-zero) bonds out of every term, so a prediction does not
    depend on which other ligands share its batch; `tests/test_structural.py`
    measures that rather than assuming it.
    """

    def __init__(self, data: QPSMPData) -> None:
        self.data = data
        self._cache: OrderedDict[str, tuple[torch.Tensor, ...]] = OrderedDict()

    def _ligand(self, ligand_id: str) -> tuple[torch.Tensor, ...]:
        if ligand_id not in self._cache:
            atoms, bonds, mask = self.data.ligand_bank.get(ligand_id)
            self._cache[ligand_id] = (
                torch.from_numpy(atoms.copy()).to(torch.float32),
                torch.from_numpy(bonds.copy()).to(torch.float32),
                torch.from_numpy(mask.copy()).to(torch.float32),
                self.data.fingerprints[ligand_id])
        self._cache.move_to_end(ligand_id)
        return self._cache[ligand_id]

    def gather(self, cell_indices: list[int], device: torch.device):
        ligands = [self.data.cells[index]["ligand_id"] for index in cell_indices]
        values = [self._ligand(key) for key in ligands]
        width = max(value[0].shape[0] for value in values)
        atoms, bonds, masks, prints = [], [], [], []
        for atom, bond, mask, fingerprint in values:
            missing = width - atom.shape[0]
            atoms.append(torch.nn.functional.pad(atom, (0, 0, 0, missing)))
            bonds.append(torch.nn.functional.pad(
                bond, (0, 0, 0, missing, 0, missing)))
            masks.append(torch.nn.functional.pad(mask, (0, missing)))
            prints.append(fingerprint)
        return (torch.stack(atoms).to(device), torch.stack(bonds).to(device),
                torch.stack(masks).to(device), torch.stack(prints).to(device))


class ProteinFeatureStore:
    """Frozen ESM-2 150M pooled + residue-slot tensors, cached per target."""

    def __init__(self, data: QPSMPData) -> None:
        self.data = data
        self._cache: dict[str, tuple[torch.Tensor, torch.Tensor, torch.Tensor]] = {}

    def get(self, target: str) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        if target not in self._cache:
            pooled, residues, mask = self.data.protein_for_target(target)
            self._cache[target] = (pooled.to(torch.float32),
                                   residues.to(torch.float32),
                                   mask.to(torch.float32))
        return self._cache[target]

    def gather(self, targets: list[str], device: torch.device):
        rows = [self.get(target) for target in targets]
        return (torch.stack([row[0] for row in rows]).to(device),
                torch.stack([row[1] for row in rows]).to(device),
                torch.stack([row[2] for row in rows]).to(device))

    def pooled_matrix(self, targets: list[str]) -> np.ndarray:
        return np.stack([self.get(target)[0].numpy() for target in targets])


# -- hard wrong proteins -----------------------------------------------------


@dataclass(frozen=True)
class DonorRule:
    """The frozen selection rule for a hard wrong protein.

    A donor must be
      * from a **different CD-HIT40 component** than the recipient;
      * drawn from `meta_train` and selected with `meta_train`-only protein
        features (frozen PLM pooled vectors and the governed component map);
      * **similar in frozen PLM space where possible** -- the most similar
        admissible protein by cosine, which is the hardest available negative,
        not a random one;
      * **free of a shared assay/document programme** -- zero DOI overlap
        between the recipient's and the donor's panel documents, so the swap
        cannot be solved by recognising a shared testing history.

    Only the protein input is replaced.  The recipient's ligands and its
    `delta_y` are untouched.
    """

    namespace: str = "sar-field-hard-wrong-protein"


def hard_wrong_protein_map(data: QPSMPData, proteins: ProteinFeatureStore,
                           recipients: list[str], candidates: list[str],
                           ) -> dict[str, str]:
    component = component_of_target(data)
    documents = target_panel_documents(data)
    pool = sorted(set(candidates))
    matrix = proteins.pooled_matrix(pool)
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    normalized = matrix / np.maximum(norms, 1e-12)
    out: dict[str, str] = {}
    for recipient in sorted(set(recipients)):
        vector = proteins.get(recipient)[0].numpy()
        vector = vector / max(float(np.linalg.norm(vector)), 1e-12)
        similarity = normalized @ vector
        order = np.argsort(-similarity)
        chosen = None
        relaxed = None
        for position in order:
            candidate = pool[int(position)]
            if candidate == recipient:
                continue
            if component[candidate] == component[recipient]:
                continue
            if relaxed is None:
                relaxed = candidate
            if documents.get(candidate, set()) & documents.get(recipient, set()):
                continue
            chosen = candidate
            break
        if chosen is None:  # pragma: no cover - not reached on this corpus
            chosen = relaxed
        if chosen is None:
            raise ValueError(f"no admissible hard wrong protein for {recipient}")
        out[recipient] = chosen
    return out


# -- controls ----------------------------------------------------------------


def cross_component_permutation(data: QPSMPData, targets: list[str],
                                seed: int) -> dict[str, str]:
    """A stable derangement of protein identity that always crosses components.

    Used by arm D.  Built by repeatedly rotating a seeded ordering until every
    target is paired with a protein from a different CD-HIT40 component, so the
    control cannot accidentally hand a target a homologous protein.
    """
    component = component_of_target(data)
    ordered = sorted(set(targets))
    rng = np.random.default_rng(stable_seed("sar-field-protein-shuffle", seed))
    permuted = list(rng.permutation(ordered))
    for offset in range(1, len(ordered)):
        mapping = {ordered[i]: permuted[(i + offset) % len(permuted)]
                   for i in range(len(ordered))}
        if all(component[key] != component[value] for key, value in mapping.items()):
            return mapping
    # Fall back to a greedy repair; recorded in the artifact if it ever runs.
    mapping = {ordered[i]: permuted[i] for i in range(len(ordered))}
    for key in ordered:  # pragma: no cover - not reached on this corpus
        if component[mapping[key]] == component[key]:
            for other in ordered:
                if (component[mapping[other]] != component[key]
                        and component[mapping[key]] != component[other]):
                    mapping[key], mapping[other] = mapping[other], mapping[key]
                    break
    return mapping


def within_target_label_shuffle(data: QPSMPData, targets: list[str],
                                seed: int) -> dict[int, float]:
    """Permute pK **inside** each target: level preserved, ordering destroyed.

    Used by arm E.  Returns a cell-index -> label map; `delta_y` is recomputed
    from it, so the arm trains on a target whose within-target SAR carries no
    information while every marginal (target mean, spread, label histogram)
    is unchanged.
    """
    out: dict[int, float] = {}
    by_target: dict[str, list[int]] = {}
    for index, cell in enumerate(data.cells):
        if cell["target_id"] in set(targets):
            by_target.setdefault(cell["target_id"], []).append(index)
    for target in sorted(by_target):
        indices = by_target[target]
        labels = [float(data.cells[i]["pK"]) for i in indices]
        rng = np.random.default_rng(
            stable_seed("sar-field-label-shuffle", seed, target))
        order = rng.permutation(len(labels))
        for position, index in enumerate(indices):
            out[index] = labels[int(order[position])]
    return out


def relabel(specs: list[PairSpec], labels: dict[int, float]) -> list[PairSpec]:
    """Rebuild `delta_y` from a substitute label map (arm E only)."""
    out: list[PairSpec] = []
    for spec in specs:
        gap = labels[spec.b] - labels[spec.a]
        out.append(PairSpec(spec.split, spec.component, spec.target, spec.a,
                            spec.b, gap, spec.tanimoto, spec.same_panel,
                            spec.same_document, spec.stratum))
    return out
