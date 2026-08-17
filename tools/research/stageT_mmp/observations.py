"""Directed MMP observations over the governed corpus.

An observation is one transformation `tau` applied inside one target `t`,
carrying `delta_y(t, tau)` in the transformation's canonical direction. The
canonical direction is a function of structure only, so the sign of the label
is never chosen by the label.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.qpsmp_data import QPSMPData

ROOT = Path(__file__).resolve().parents[3]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
PROTEIN_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank"
LIGAND_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank"
COMPACT_LIGAND_BANK = (
    ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank_compact")
SPLIT_DIRECTORY = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
SPLIT_VIEW = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views"

# Frozen, matching the Stage L2 / Stage S definition so the cliff stratum stays
# comparable with the recorded measurements.
CLIFF_TANIMOTO = 0.6
CLIFF_GAP = 1.0


def load_governed() -> tuple[QPSMPData, dict]:
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT_DIRECTORY, split_view=SPLIT_VIEW)
    return data, data.seal_record()


@dataclass(frozen=True)
class MMPObservation:
    target: str
    component: str
    core: str
    exact_key: str
    coarse_key: str
    cell_a: int
    cell_b: int
    ligand_a: str
    ligand_b: str
    delta_y: float
    same_panel: bool
    stratum: str
    tanimoto: float
    activity_cliff: bool
    stereo_edit: bool
    charge_change: int

    @property
    def pair_id(self) -> str:
        return f"{self.target}:{self.cell_a}:{self.cell_b}"


def _stratum(cell_a: dict, cell_b: dict) -> tuple[bool, str]:
    """Frozen, label-blind confidence stratum (PREREGISTRATION.md section 2)."""
    panels_a = set(cell_a["panel_ids"])
    panels_b = set(cell_b["panel_ids"])
    shared = bool(panels_a & panels_b)
    if not shared:
        return False, "S3_cross_panel"
    if cell_a["panel_count"] == 1 and cell_b["panel_count"] == 1:
        return True, "S1_same_panel_single"
    return True, "S2_same_panel_multi"


def build_observations(data: QPSMPData, targets: list[str]) -> dict:
    """All directed MMP observations for the given meta_train targets."""
    from tools.research.stageT_mmp.mmp import fragment, transformation

    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    by_target: dict[str, list[int]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    for index, cell in enumerate(data.cells):
        if cell["split"] != "meta_train" or cell["target_id"] not in set(targets):
            continue
        key = (cell["target_id"], cell["ligand_id"])
        if key in seen:
            continue
        seen.add(key)
        by_target[cell["target_id"]].append(index)

    fingerprints = data.fingerprints
    observations: list[MMPObservation] = []
    unparsable = 0
    no_fragmentation = 0
    considered_ligands = 0

    for target in sorted(by_target):
        indices = sorted(by_target[target])
        fragmentations: dict[int, tuple] = {}
        for index in indices:
            smiles = data._ligand_smiles.get(data.cells[index]["ligand_id"])
            considered_ligands += 1
            if not smiles:
                unparsable += 1
                continue
            pieces = fragment(smiles)
            if not pieces:
                no_fragmentation += 1
                continue
            fragmentations[index] = pieces
        # Index by core so only genuine matched pairs are compared.
        by_core: dict[str, list[tuple[int, object]]] = defaultdict(list)
        for index, pieces in fragmentations.items():
            for piece in pieces:
                by_core[piece.core].append((index, piece))
        emitted: set[tuple[str, str, str]] = set()
        for core, entries in sorted(by_core.items()):
            for position, (left_index, left) in enumerate(entries):
                for right_index, right in entries[position + 1:]:
                    if left_index == right_index:
                        continue
                    built = transformation(left, right)
                    if built is None:
                        continue
                    item, flipped = built
                    # delta_y is y(ligand carrying r_b) - y(ligand carrying r_a).
                    if flipped:
                        index_a, index_b = right_index, left_index
                    else:
                        index_a, index_b = left_index, right_index
                    # Deduplicate a repeated (target, core, r_a, r_b) by the
                    # lower cell index -- never by label value.
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
                    observations.append(MMPObservation(
                        target=target, component=component_of[target], core=core,
                        exact_key=item.exact_key, coarse_key=item.coarse_key,
                        cell_a=index_a, cell_b=index_b,
                        ligand_a=cell_a["ligand_id"], ligand_b=cell_b["ligand_id"],
                        delta_y=delta, same_panel=same_panel, stratum=stratum,
                        tanimoto=tanimoto,
                        activity_cliff=bool(tanimoto >= CLIFF_TANIMOTO
                                            and abs(delta) >= CLIFF_GAP),
                        stereo_edit=item.stereo_edit,
                        charge_change=item.charge_change))

    return {
        "observations": observations,
        "construction": {
            "machinery": "rdkit.Chem.rdMMPA.FragmentMol (Hussain-Rea), single cut",
            "targets_considered": len(by_target),
            "ligand_slots_considered": considered_ligands,
            "ligands_without_smiles": unparsable,
            "ligands_with_no_admissible_cut": no_fragmentation,
            "deduplication": ("one row per (target, core, r_a, r_b); the lower "
                              "cell index wins, never the label"),
            "canonical_direction": ("canonical SMILES sort of the two R groups; "
                                    "delta_y = y(r_b ligand) - y(r_a ligand)"),
            "cliff_definition": (f"tanimoto >= {CLIFF_TANIMOTO} and "
                                 f"|delta_y| >= {CLIFF_GAP}"),
        },
    }
