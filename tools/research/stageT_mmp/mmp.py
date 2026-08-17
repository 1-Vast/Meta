"""True matched-molecular-pair transformations via RDKit's MMPA machinery.

Every rule here is **label-blind**: nothing in this module reads a `pK`. That is
the property that lets the T1 census and the T2 splits be built without a label
path, and `tests/` asserts it by AST parse rather than by inspection.

Construction (frozen in PREREGISTRATION.md section 3):

* `rdkit.Chem.rdMMPA.FragmentMol(mol, minCuts=1, maxCuts=1, ...)` -- the
  supported Hussain-Rea implementation. No ad-hoc SMILES string surgery;
* of the two fragments the one with more heavy atoms is the **core**, ties
  broken by canonical SMILES sort;
* a pair is admissible when two ligands of **one target** share an identical
  core (isomeric SMILES, including the `[*:1]` label) and differ in R;
* the **attachment context** is the core atom bearing `[*:1]`, recorded as
  `(element, aromatic, in_ring, degree, formal_charge)`;
* stereochemistry is retained everywhere (isomeric SMILES); a transformation
  whose R groups differ only by stereochemistry is flagged, never collapsed;
* the formal-charge change `q(R_b) - q(R_a)` is recorded on every transformation;
* the **canonical direction** is the canonical-SMILES sort order of the two R
  groups, so direction is a function of structure alone and the inverse maps to
  the negated label by construction.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import re
import sys

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMMPA

RDLogger.DisableLog("rdApp.*")

# Frozen fragmentation settings.
MIN_CUTS = 1
MAX_CUTS = 1
MAX_CUT_BONDS = 30
ATTACHMENT = "[*:1]"
# RDKit's default acyclic single-bond pattern. Passed explicitly because
# `FragmentMol` is overloaded: four positional integers bind to the *first*
# overload (maxCuts, maxCutBonds, pattern) and silently mis-typed `pattern`,
# which raised inside RDKit and -- under a broad `except` -- would have produced
# an empty census that looked like a real negative result.
CUT_PATTERN = "[#6+0;!$(*=,#[!#6])]!@!=!#[*]"


@dataclass(frozen=True)
class Fragmentation:
    """One single-cut decomposition of a molecule into core + R."""

    core: str
    r_group: str
    attachment_element: str
    attachment_aromatic: bool
    attachment_in_ring: bool
    attachment_degree: int
    attachment_charge: int

    @property
    def context(self) -> tuple:
        return (self.attachment_element, self.attachment_aromatic,
                self.attachment_in_ring, self.attachment_degree,
                self.attachment_charge)

    @property
    def coarse_context(self) -> tuple:
        return (self.attachment_element, self.attachment_aromatic)


def _heavy_atoms(smiles: str) -> int:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return 0
    return sum(1 for atom in molecule.GetAtoms()
               if atom.GetAtomicNum() > 1 and atom.GetAtomicNum() != 0)


def _attachment_context(smiles: str) -> tuple | None:
    """Properties of the heavy atom the `[*:1]` dummy is bonded to.

    Returns plain Python data, never an RDKit `Atom`. An `Atom` borrows a
    pointer into its parent `Mol`; returning one lets the `Mol` be collected
    while the caller still holds the atom, and the next attribute access reads
    freed memory. That is a segfault, not an exception, and it took down this
    module until the ownership was cut here.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        neighbours = list(atom.GetNeighbors())
        if len(neighbours) != 1:
            continue
        neighbour = neighbours[0]
        return (neighbour.GetSymbol(), bool(neighbour.GetIsAromatic()),
                bool(neighbour.IsInRing()), int(neighbour.GetDegree()),
                int(neighbour.GetFormalCharge()))
    return None


def formal_charge(smiles: str) -> int:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return 0
    return Chem.GetFormalCharge(molecule)


def strip_stereochemistry(smiles: str) -> str:
    """Canonical SMILES with stereochemistry removed, for the coarse key."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return smiles
    Chem.RemoveStereochemistry(molecule)
    return Chem.MolToSmiles(molecule)


def has_stereocentre(smiles: str) -> bool:
    return smiles != strip_stereochemistry(smiles)


@lru_cache(maxsize=32768)
def fragment(smiles: str) -> tuple[Fragmentation, ...]:
    """All single-cut core/R decompositions of one ligand.

    Returns an empty tuple for a molecule RDKit cannot parse or that has no
    admissible acyclic single bond to cut. Both are ordinary outcomes and the
    census counts them.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return ()
    # Fully positional against the 6-argument overload. Both alternatives are
    # unsafe on this Boost.Python signature: four positional integers bind to
    # the *first* overload and mis-type `pattern`, and naming every argument
    # (including `pattern`) segfaults the interpreter outright. Only this form
    # is verified. No broad `except` either -- an RDKit failure must surface
    # rather than silently shrink the census into a false negative.
    raw = rdMMPA.FragmentMol(molecule, MIN_CUTS, MAX_CUTS, MAX_CUT_BONDS,
                             CUT_PATTERN, False)
    out: list[Fragmentation] = []
    for _core, chains in raw:
        if not chains or "." not in chains:
            continue
        left, right = chains.split(".", 1)
        if ATTACHMENT not in left or ATTACHMENT not in right:
            continue
        # Core = larger fragment; ties broken by canonical SMILES sort so the
        # choice is deterministic and independent of input order.
        sizes = (_heavy_atoms(left), _heavy_atoms(right))
        if sizes[0] > sizes[1]:
            core, r_group = left, right
        elif sizes[1] > sizes[0]:
            core, r_group = right, left
        else:
            core, r_group = tuple(sorted((left, right)))
        context = _attachment_context(core)
        if context is None:
            continue
        element, aromatic, in_ring, degree, charge = context
        out.append(Fragmentation(
            core=core, r_group=r_group, attachment_element=element,
            attachment_aromatic=aromatic, attachment_in_ring=in_ring,
            attachment_degree=degree, attachment_charge=charge))
    # Deduplicate deterministically.
    return tuple(sorted(set(out), key=lambda value: (value.core, value.r_group)))


@dataclass(frozen=True)
class Transformation:
    """A directed MMP transformation, canonicalised without reading any label."""

    core: str
    r_a: str
    r_b: str
    context: tuple
    coarse_context: tuple
    charge_change: int
    stereo_edit: bool

    @property
    def exact_key(self) -> str:
        return f"{self.context}|{self.r_a}>>{self.r_b}"

    @property
    def coarse_key(self) -> str:
        return (f"{self.coarse_context}|{strip_stereochemistry(self.r_a)}>>"
                f"{strip_stereochemistry(self.r_b)}")

    @property
    def inverse(self) -> "Transformation":
        return Transformation(self.core, self.r_b, self.r_a, self.context,
                              self.coarse_context, -self.charge_change,
                              self.stereo_edit)


def transformation(left: Fragmentation, right: Fragmentation
                   ) -> tuple[Transformation, bool] | None:
    """Build the canonically directed transformation for two fragmentations.

    Returns `(transformation, flipped)` where `flipped` is True when the
    canonical direction is `right -> left`, i.e. the caller must negate the
    label it computed as `y(left_ligand) -> y(right_ligand)`. Returns None when
    the two fragmentations are not a matched pair.
    """
    if left.core != right.core or left.r_group == right.r_group:
        return None
    if left.context != right.context:
        return None
    # Canonical direction: sort the two R groups as strings. Structure only.
    if left.r_group <= right.r_group:
        r_a, r_b, flipped = left.r_group, right.r_group, False
    else:
        r_a, r_b, flipped = right.r_group, left.r_group, True
    stereo = strip_stereochemistry(r_a) == strip_stereochemistry(r_b)
    return Transformation(
        core=left.core, r_a=r_a, r_b=r_b, context=left.context,
        coarse_context=left.coarse_context,
        charge_change=formal_charge(r_b) - formal_charge(r_a),
        stereo_edit=bool(stereo)), flipped


def descriptor(item: Transformation) -> list[float]:
    """Structured, label-blind numeric descriptor of a transformation.

    Deliberately small and interpretable: heavy-atom and ring deltas, element
    counts of the removed and added fragments, the attachment context, the
    formal-charge change and a stereo flag. No learned ligand embedding, so the
    T2 discriminator cannot smuggle a whole-molecule model back in.
    """
    elements = ("C", "N", "O", "S", "F", "Cl", "Br", "I", "P")

    def counts(smiles: str) -> list[float]:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return [0.0] * (len(elements) + 4)
        symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()
                   if atom.GetAtomicNum() > 0]
        aromatic = sum(1 for atom in molecule.GetAtoms() if atom.GetIsAromatic())
        rings = molecule.GetRingInfo().NumRings()
        rotatable = sum(1 for bond in molecule.GetBonds()
                        if bond.GetBondType() == Chem.BondType.SINGLE
                        and not bond.IsInRing())
        return ([float(symbols.count(element)) for element in elements]
                + [float(len(symbols)), float(aromatic), float(rings),
                   float(rotatable)])

    left, right = counts(item.r_a), counts(item.r_b)
    delta = [b - a for a, b in zip(left, right)]
    context = [
        float(ord(item.context[0][0])) / 100.0,
        float(item.context[1]), float(item.context[2]),
        float(item.context[3]), float(item.context[4]),
    ]
    return (left + right + delta + context
            + [float(item.charge_change), float(item.stereo_edit)])


DESCRIPTOR_WIDTH = len(descriptor(Transformation(
    "c1ccccc1[*:1]", "C[*:1]", "CC[*:1]", ("C", True, True, 3, 0), ("C", True),
    0, False)))
