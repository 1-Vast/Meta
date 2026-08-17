"""Soft controlled chemical-change families for Stage W.

The exact-MMP branch was closed on BindingDB-Ki because the strict key does not
recur across cold protein components. This module builds a **soft family**:
single-cut MMP fragments are kept, but families aggregate different exact cores
and R identities through a frozen, structure-only key:

    sha256( murcko_core | attachment_element | attachment_aromatic |
            attachment_in_ring | category(R_a) >> category(R_b) )

Everything here is label-blind; `delta_y` is computed by the census module.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import sys

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMMPA
from rdkit.Chem.Scaffolds import MurckoScaffold

from tools.research.stageU_mmp_interaction.mmp import fragment, transformation

ATTACHMENT = "[*:1]"


def murcko_core(core: str) -> str:
    """Murcko scaffold of the core with the cut dummy replaced by H."""
    molecule = Chem.MolFromSmiles(core.replace(ATTACHMENT, "[H]"))
    if molecule is None:
        return core
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return Chem.MolToSmiles(molecule)
    return Chem.MolToSmiles(scaffold)


def category(r_smiles: str) -> tuple:
    """Quantized pharmacophore/change class of an R fragment."""
    molecule = Chem.MolFromSmiles(r_smiles)
    if molecule is None:
        return (0, 0, 0, 0, 0, 0)
    heavy = sum(1 for atom in molecule.GetAtoms()
                if atom.GetAtomicNum() > 1)
    aromatic = any(atom.GetIsAromatic() for atom in molecule.GetAtoms())
    rings = molecule.GetRingInfo().NumRings() > 0
    hbd = Lipinski.NumHDonors(molecule)
    hba = Lipinski.NumHAcceptors(molecule)
    charge = Chem.GetFormalCharge(molecule)
    def bucket(value):
        return 0 if value == 0 else 1 if value <= 3 else 2 if value <= 7 else 3
    return (bucket(heavy), int(aromatic), int(rings),
            min(hbd, 2), min(hba, 2),
            0 if charge < 0 else 1 if charge == 0 else 2)


@dataclass(frozen=True)
class SoftTransformation:
    exact_key: str
    family_key: str
    murcko_core: str
    category_a: tuple
    category_b: tuple
    r_a: str
    r_b: str
    flipped: bool

    @property
    def family_string(self) -> str:
        return (f"{self.murcko_core}|{self.category_a}>>{self.category_b}")


def soft_transformation(left, right):
    """Canonical soft family for two single-cut MMP fragments."""
    built = transformation(left, right)
    if built is None:
        return None
    item, flipped = built
    # Family direction sorts by (category, SMILES), structure only.
    cat_left, cat_right = category(item.r_a), category(item.r_b)
    if (cat_left, item.r_a) <= (cat_right, item.r_b):
        cat_a, cat_b, r_a, r_b, family_flip = cat_left, cat_right, item.r_a, item.r_b, False
    else:
        cat_a, cat_b, r_a, r_b, family_flip = cat_right, cat_left, item.r_b, item.r_a, True
    core = murcko_core(item.core)
    payload = (f"{core}|{item.context[0]}|{item.context[1]}|"
               f"{item.context[2]}|{cat_a}>>{cat_b}")
    family_key = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    return SoftTransformation(
        exact_key=item.exact_key, family_key=family_key, murcko_core=core,
        category_a=cat_a, category_b=cat_b, r_a=r_a, r_b=r_b,
        flipped=(flipped != family_flip))
