"""Deterministic within-target ligand-pair banks for the SAR-field stage.

Every construction rule in this module is fixed by SHA-256 derived seeds
(`scripts.qpsmp_data.stable_seed`).  Python's builtin `hash()` is never used:
it is salted per process by `PYTHONHASHSEED`, and the same defect has already
been caught twice in this repository (Stage R, Stage A).

Contract:

* a pair is **always** two distinct ligands measured against **one** target.
  There is no cross-target pair construction anywhere in this stage;
* the pair bank is a function of (split, component set, seed) only, so two
  processes build the identical bank;
* `same_panel` means the two cells share at least one governed `panel_ids`
  entry, i.e. the same document, the same endpoint and the same target.  That
  is the highest-confidence comparable context available in this corpus;
* `same_document` is the weaker DOI-level relation, kept separately because two
  panels of one paper can still be different assay protocols.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np
import torch

if __package__ in {None, ""}:  # pragma: no cover - direct execution helper
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.qpsmp_data import QPSMPData, stable_seed

ROOT = Path(__file__).resolve().parents[3]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
PROTEIN_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank"
LIGAND_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank"
COMPACT_LIGAND_BANK = (
    ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank_compact")
SPLIT_DIRECTORY = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
SPLIT_VIEW = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views"

# Frozen before any Phase 2 arm trained.  Matches the Stage L2 definition so
# the activity-cliff stratum is comparable with the recorded directional-SAR
# measurement.
CLIFF_TANIMOTO = 0.6
CLIFF_GAP = 1.0
# Pair-distance strata.  "local" is chemically close, "distant" is the tail.
LOCAL_TANIMOTO = 0.5
MEDIUM_TANIMOTO = 0.25


def load_data(include_meta_test: bool = False) -> QPSMPData:
    """Mount the physically isolated governed split view.

    `split_view` is mandatory in this stage: on that surface the `meta_test`
    label artifact is not present in the development tree at all, so no sealed
    label is decompressed or parsed by this process.
    """
    if include_meta_test:
        raise ValueError(
            "this stage never opens meta_test; the argument exists only so a "
            "reader can see that the refusal is explicit")
    return QPSMPData(
        CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
        split_directory=SPLIT_DIRECTORY, split_view=SPLIT_VIEW)


@dataclass(frozen=True)
class PairSpec:
    """One directed within-target ligand pair.

    `a` and `b` are indices into `QPSMPData.cells`.  Both cells carry the same
    `target`; this is enforced at construction and re-checked by the structural
    tests.
    """

    split: str
    component: str
    target: str
    a: int
    b: int
    delta_y: float
    tanimoto: float
    same_panel: bool
    same_document: bool
    stratum: str

    @property
    def reversed(self) -> "PairSpec":
        return PairSpec(self.split, self.component, self.target, self.b, self.a,
                        -self.delta_y, self.tanimoto, self.same_panel,
                        self.same_document, self.stratum)


def _documents(panel_ids) -> set[str]:
    return {str(panel).split("|")[0] for panel in panel_ids}


def target_cell_index(data: QPSMPData, split: str) -> dict[str, list[int]]:
    """One cell per (target, ligand), deterministic in corpus order.

    Duplicate ligand rows inside a target would make `delta_y` ambiguous, so
    the first occurrence in corpus order wins and the rest are dropped.  The
    audit reports how many rows this removes.
    """
    index: dict[str, list[int]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    for position, cell in enumerate(data.cells):
        if cell["split"] != split:
            continue
        key = (cell["target_id"], cell["ligand_id"])
        if key in seen:
            continue
        seen.add(key)
        index[cell["target_id"]].append(position)
    return dict(index)


def tanimoto_matrix(data: QPSMPData, indices: list[int]) -> np.ndarray:
    rows = data.fingerprint_rows(tuple(indices)).numpy()
    intersection = rows @ rows.T
    counts = rows.sum(axis=1)
    union = counts[:, None] + counts[None, :] - intersection
    with np.errstate(divide="ignore", invalid="ignore"):
        value = np.where(union > 0, intersection / np.maximum(union, 1e-12), 0.0)
    return value.astype(np.float64)


def classify(tanimoto: float, gap: float) -> str:
    """Assign one pair to exactly one sampling stratum.

    The cliff test comes first: an activity cliff is a *local* pair with a large
    gap, and it must not be absorbed into the ordinary local stratum, otherwise
    target-balanced sampling would drown the pairs the hypothesis is about.
    """
    if tanimoto >= CLIFF_TANIMOTO and abs(gap) >= CLIFF_GAP:
        return "cliff"
    if tanimoto >= LOCAL_TANIMOTO:
        return "local"
    if tanimoto >= MEDIUM_TANIMOTO:
        return "medium"
    return "distant"


def build_target_pairs(data: QPSMPData, split: str, targets: list[str] | None = None,
                       ) -> dict[str, list[PairSpec]]:
    """All unordered within-target pairs, in one canonical orientation.

    The canonical orientation is corpus order (`a < b` as cell indices), which
    is independent of the labels, so `sign(delta_y)` is not biased by it.  The
    model built in this stage is exactly antisymmetric, so a canonical bank and
    a both-orientations bank give identical metrics; `expand_orientations`
    exists to verify that rather than to change it.
    """
    cells_by_target = target_cell_index(data, split)
    component_of = {}
    for cell in data.cells:
        component_of[cell["target_id"]] = cell["protein_group_40"]
    wanted = set(targets) if targets is not None else set(cells_by_target)
    out: dict[str, list[PairSpec]] = {}
    for target, indices in sorted(cells_by_target.items()):
        if target not in wanted or len(indices) < 2:
            continue
        similarity = tanimoto_matrix(data, indices)
        panels = [set(data.cells[i]["panel_ids"]) for i in indices]
        documents = [_documents(data.cells[i]["panel_ids"]) for i in indices]
        labels = [float(data.cells[i]["pK"]) for i in indices]
        specs: list[PairSpec] = []
        for left in range(len(indices)):
            for right in range(left + 1, len(indices)):
                gap = labels[right] - labels[left]
                value = float(similarity[left, right])
                specs.append(PairSpec(
                    split, component_of[target], target,
                    indices[left], indices[right], gap, value,
                    bool(panels[left] & panels[right]),
                    bool(documents[left] & documents[right]),
                    classify(value, gap)))
        out[target] = specs
    return out


def expand_orientations(specs: list[PairSpec]) -> list[PairSpec]:
    """Balanced (a,b)+(b,a) bank.  Used by the orientation-isolation test."""
    out: list[PairSpec] = []
    for spec in specs:
        out.append(spec)
        out.append(spec.reversed)
    return out


# Target-balanced sampling ---------------------------------------------------

STRATA = ("local", "medium", "cliff", "distant")


def target_balanced_bank(pairs_by_target: dict[str, list[PairSpec]], seed: int,
                         per_target: int, namespace: str = "sar-field-train",
                         ) -> list[PairSpec]:
    """Draw at most `per_target` pairs from every target, stratum-balanced.

    Global pair sampling would let one 200-ligand target contribute ~20,000 of
    the pairs and dominate the objective.  Here each target contributes the same
    quota, split as evenly as its own strata allow, so local pairs,
    medium-distance pairs, activity cliffs and distant pairs are all represented
    without letting large targets set the loss.
    """
    if per_target < 1:
        raise ValueError("per_target must be positive")
    drawn: list[PairSpec] = []
    for target in sorted(pairs_by_target):
        specs = pairs_by_target[target]
        if not specs:
            continue
        buckets = {name: [s for s in specs if s.stratum == name] for name in STRATA}
        present = [name for name in STRATA if buckets[name]]
        quota = {name: per_target // len(present) for name in present}
        for extra in range(per_target - sum(quota.values())):
            quota[present[extra % len(present)]] += 1
        for name in present:
            bucket = buckets[name]
            rng = np.random.default_rng(stable_seed(namespace, seed, target, name))
            take = min(quota[name], len(bucket))
            order = rng.permutation(len(bucket))[:take]
            drawn.extend(bucket[int(i)] for i in sorted(order))
    return drawn


def evaluation_bank(pairs_by_target: dict[str, list[PairSpec]], seed: int,
                    per_target: int, namespace: str = "sar-field-eval",
                    ) -> list[PairSpec]:
    """A frozen, target-capped evaluation bank over every stratum.

    Same construction as the training bank but a different namespace, so a
    training draw and an evaluation draw are independent yet both reproducible.
    """
    return target_balanced_bank(pairs_by_target, seed, per_target, namespace)


# Protein-side helpers -------------------------------------------------------


def component_of_target(data: QPSMPData) -> dict[str, str]:
    return {cell["target_id"]: cell["protein_group_40"] for cell in data.cells}


def target_panel_documents(data: QPSMPData) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for cell in data.cells:
        out[cell["target_id"]].update(_documents(cell["panel_ids"]))
    return dict(out)


def protein_features(data: QPSMPData, target: str) -> torch.Tensor:
    """Frozen pooled PLM vector for one target, float32 on CPU."""
    pooled, _residues, _mask = data.protein_for_target(target)
    return pooled.to(torch.float32)
