"""Core-inclusive MMP transformations and double differences.

The one substantive difference from Stage T: `tau` carries the **shared core**,
so two targets are compared only when they realise the identical chemical
transformation. Stage 0 measured what Stage T's core-blind key cost — a
within-target across-core nuisance of median 0.269 pK (p95 1.268) contaminating
40.4% of its training rows — which is the whole reason this module exists.

Chemistry is inherited from Stage U's `mmp.py`, which already implements the
core-inclusive key and a core-consuming descriptor. Nothing here reads a label.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.qpsmp_data import QPSMPData, stable_seed
from tools.research.stageU_mmp_interaction.mmp import (
    DESCRIPTOR_WIDTH, EDIT_WIDTH, Transformation, descriptor, edit_features,
    fragment, transformation,
)

ROOT = Path(__file__).resolve().parents[3]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
PROTEIN_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank"
LIGAND_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank"
COMPACT_LIGAND_BANK = (
    ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank_compact")
SPLIT_DIRECTORY = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
SPLIT_VIEW = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views"

CLIFF_TANIMOTO = 0.6
CLIFF_GAP = 1.0
SAME_PANEL_STRATA = ("S1_same_panel_single", "S2_same_panel_multi")


def load_governed() -> tuple[QPSMPData, dict]:
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT_DIRECTORY, split_view=SPLIT_VIEW)
    return data, data.seal_record()


@dataclass(frozen=True)
class Observation:
    target: str
    component: str
    core: str
    exact_key: str
    coarse_key: str
    cell_a: int
    cell_b: int
    delta_y: float
    same_panel: bool
    stratum: str
    tanimoto: float
    activity_cliff: bool
    stereo_edit: bool
    charge_change: int
    edit: tuple[float, ...]


def _stratum(cell_a: dict, cell_b: dict) -> tuple[bool, str]:
    panels_a, panels_b = set(cell_a["panel_ids"]), set(cell_b["panel_ids"])
    if not (panels_a & panels_b):
        return False, "S3_cross_panel"
    if cell_a["panel_count"] == 1 and cell_b["panel_count"] == 1:
        return True, "S1_same_panel_single"
    return True, "S2_same_panel_multi"


def build_observations(data: QPSMPData, targets: list[str]) -> dict:
    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    wanted = set(targets)
    by_target: dict[str, list[int]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    for index, cell in enumerate(data.cells):
        if cell["split"] != "meta_train" or cell["target_id"] not in wanted:
            continue
        key = (cell["target_id"], cell["ligand_id"])
        if key in seen:
            continue
        seen.add(key)
        by_target[cell["target_id"]].append(index)

    fingerprints = data.fingerprints
    observations: list[Observation] = []
    no_cut = 0
    for target in sorted(by_target):
        indices = sorted(by_target[target])
        pieces: dict[int, tuple] = {}
        for index in indices:
            smiles = data._ligand_smiles.get(data.cells[index]["ligand_id"])
            if not smiles:
                continue
            parts = fragment(smiles)
            if not parts:
                no_cut += 1
                continue
            pieces[index] = parts
        by_core: dict[str, list[tuple[int, object]]] = defaultdict(list)
        for index, parts in pieces.items():
            for part in parts:
                by_core[part.core].append((index, part))
        emitted: set[tuple[str, str, str]] = set()
        for core, entries in sorted(by_core.items()):
            for position, (left_index, left) in enumerate(entries):
                for right_index, right in entries[position + 1:]:
                    built = transformation(left, right)
                    if built is None:
                        continue
                    item, flipped = built
                    index_a, index_b = ((right_index, left_index) if flipped
                                        else (left_index, right_index))
                    signature = (core, item.r_a, item.r_b)
                    if signature in emitted:
                        continue
                    emitted.add(signature)
                    cell_a, cell_b = data.cells[index_a], data.cells[index_b]
                    same_panel, stratum = _stratum(cell_a, cell_b)
                    left_fp = fingerprints[cell_a["ligand_id"]].numpy()
                    right_fp = fingerprints[cell_b["ligand_id"]].numpy()
                    intersection = float(left_fp @ right_fp)
                    union = float(left_fp.sum() + right_fp.sum() - intersection)
                    tanimoto = intersection / union if union > 0 else 0.0
                    delta = float(cell_b["pK"]) - float(cell_a["pK"])
                    observations.append(Observation(
                        target=target, component=component_of[target], core=core,
                        exact_key=item.exact_key, coarse_key=item.coarse_key,
                        cell_a=index_a, cell_b=index_b, delta_y=delta,
                        same_panel=same_panel, stratum=stratum,
                        tanimoto=tanimoto,
                        activity_cliff=bool(tanimoto >= CLIFF_TANIMOTO
                                            and abs(delta) >= CLIFF_GAP),
                        stereo_edit=item.stereo_edit,
                        charge_change=item.charge_change,
                        edit=tuple(edit_features(item))))
    return {
        "observations": observations,
        "construction": {
            "machinery": "rdkit.Chem.rdMMPA.FragmentMol (Hussain-Rea), single cut",
            "exact_key": ("sha256(core_isomeric | repr(attachment environment) "
                          "| R_a >> R_b)"),
            "coarse_key": ("sha256(core_stereo_stripped | element | aromatic | "
                           "R_a_stripped >> R_b_stripped)"),
            "edit_feature_width": EDIT_WIDTH,
            "descriptor_width": DESCRIPTOR_WIDTH,
            "targets_considered": len(by_target),
            "ligands_with_no_admissible_cut": no_cut,
            "canonical_direction": "delta_y = y(R_b ligand) - y(R_a ligand)",
        },
    }


@dataclass(frozen=True)
class TargetEffect:
    key: str
    coarse_key: str
    target: str
    component: str
    delta_y: float
    observations: int
    cores: int
    activity_cliff: bool
    edit: tuple[float, ...]


def target_effects(observations: list[Observation]) -> list[TargetEffect]:
    grouped: dict[tuple[str, str], list[Observation]] = defaultdict(list)
    for item in observations:
        if item.stratum not in SAME_PANEL_STRATA:
            continue
        grouped[(item.exact_key, item.target)].append(item)
    out: list[TargetEffect] = []
    for (key, target), rows in sorted(grouped.items()):
        out.append(TargetEffect(
            key=key, coarse_key=rows[0].coarse_key, target=target,
            component=rows[0].component,
            delta_y=float(np.median([r.delta_y for r in rows])),
            observations=len(rows), cores=len({r.core for r in rows}),
            activity_cliff=any(r.activity_cliff for r in rows),
            edit=rows[0].edit))
    return out


@dataclass(frozen=True)
class DoubleDifference:
    key: str
    coarse_key: str
    target_left: str
    target_right: str
    component_left: str
    component_right: str
    value: float
    cross_component: bool
    activity_cliff: bool
    edit: tuple[float, ...]

    @property
    def row_id(self) -> str:
        return f"{self.key}|{self.target_left}|{self.target_right}"


def double_differences(effects: list[TargetEffect]) -> list[DoubleDifference]:
    by_key: dict[str, list[TargetEffect]] = defaultdict(list)
    for effect in effects:
        by_key[effect.key].append(effect)
    out: list[DoubleDifference] = []
    for key in sorted(by_key):
        rows = sorted(by_key[key], key=lambda value: value.target)
        for left, right in combinations(rows, 2):
            out.append(DoubleDifference(
                key=key, coarse_key=left.coarse_key,
                target_left=left.target, target_right=right.target,
                component_left=left.component, component_right=right.component,
                value=left.delta_y - right.delta_y,
                cross_component=left.component != right.component,
                activity_cliff=left.activity_cliff or right.activity_cliff,
                edit=left.edit))
    return out


def effective_independent_units(rows: list[DoubleDifference]) -> dict:
    components = {c for row in rows
                  for c in (row.component_left, row.component_right)}
    keys = {row.key for row in rows}
    return {
        "rows": len(rows),
        "cross_component_rows": sum(1 for row in rows if row.cross_component),
        "protein_components": len(components),
        "transformation_keys": len(keys),
        "effective_independent_units": min(len(components), len(keys)),
    }
