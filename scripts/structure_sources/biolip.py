"""BioLiP2 annotation parsing; no coordinate or model concerns live here."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import gzip
from pathlib import Path
from typing import Iterator, TextIO


BIOLIP_ANNOTATION_URL = "https://zhanggroup.org/BioLiP/download/BioLiP.txt.gz"
BIOLIP_PROTEIN_FASTA_URL = "https://zhanggroup.org/BioLiP/data/protein.fasta.gz"
BIOLIP_LIGAND_SUMMARY_URL = "https://zhanggroup.org/BioLiP/data/ligand.tsv.gz"
BIOLIP_DOWNLOAD_PAGE = "https://zhanggroup.org/BioLiP/download.html"
BIOLIP_README_URL = "https://zhanggroup.org/BioLiP/download/readme.txt"


@dataclass(frozen=True)
class BioLiPEntry:
    pdb_id: str
    receptor_auth_asym_id: str
    resolution: float
    binding_site_id: str
    ligand_comp_id: str
    ligand_auth_asym_id: str
    ligand_serial: str
    ligand_auth_seq_id: str
    sequence: str

    @property
    def source_entry_id(self) -> str:
        return ":".join((self.pdb_id, self.receptor_auth_asym_id,
                         self.ligand_comp_id, self.ligand_auth_asym_id,
                         self.ligand_auth_seq_id, self.binding_site_id))

    def to_dict(self) -> dict:
        return asdict(self) | {"source_entry_id": self.source_entry_id}


def _open_text(path: str | Path) -> TextIO:
    source = Path(path)
    if source.suffix == ".gz":
        return gzip.open(source, "rt", encoding="utf-8", errors="strict")
    return source.open(encoding="utf-8")


def parse_biolip_line(line: str) -> BioLiPEntry:
    columns = line.rstrip("\r\n").split("\t")
    if len(columns) != 21:
        raise ValueError(f"BioLiP row must contain 21 tab-separated columns, got {len(columns)}")
    sequence = "".join(columns[20].upper().split())
    if not sequence:
        raise ValueError("BioLiP receptor sequence is empty")
    return BioLiPEntry(
        pdb_id=columns[0].strip().lower(),
        receptor_auth_asym_id=columns[1].strip(),
        resolution=float(columns[2].strip() or "-1"),
        binding_site_id=columns[3].strip(),
        ligand_comp_id=columns[4].strip().upper(),
        ligand_auth_asym_id=columns[5].strip(),
        ligand_serial=columns[6].strip(),
        ligand_auth_seq_id=columns[19].strip(),
        sequence=sequence,
    )


def iter_biolip(path: str | Path) -> Iterator[BioLiPEntry]:
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                yield parse_biolip_line(line)
            except Exception as error:
                raise ValueError(f"invalid BioLiP row {line_number}: {error}") from error


def regular_ligand_ids(path: str | Path) -> set[str]:
    """Use the BioLiP ligand summary for a cheap heavy-atom prefilter only."""
    from rdkit import Chem, RDLogger
    accepted = set()
    RDLogger.DisableLog("rdApp.error")
    try:
        with _open_text(path) as handle:
            for line in handle:
                columns = line.rstrip("\r\n").split("\t")
                if len(columns) < 5 or columns[0].startswith("#"):
                    continue
                for smiles in columns[4].split(";"):
                    molecule = Chem.MolFromSmiles(smiles.strip())
                    if molecule is not None and 6 <= molecule.GetNumHeavyAtoms() <= 96:
                        accepted.add(columns[0].strip().upper())
                        break
    finally:
        RDLogger.EnableLog("rdApp.error")
    return accepted


def pilot_candidates(path: str | Path, *, limit: int | None = None,
                     allowed_ligands: set[str] | None = None,
                     ligand_diversity_target: int = 4000) -> list[BioLiPEntry]:
    """Apply coordinate-free QC and choose a deterministic diversity-first pilot."""
    unique: dict[str, BioLiPEntry] = {}
    for entry in iter_biolip(path):
        if not (0 < entry.resolution <= 3.0 and 50 <= len(entry.sequence) <= 1022):
            continue
        if not entry.ligand_auth_seq_id.lstrip("-").isdigit():
            continue
        if allowed_ligands is not None and entry.ligand_comp_id not in allowed_ligands:
            continue
        unique.setdefault(entry.source_entry_id, entry)
    entries = list(unique.values())
    if limit is None or len(entries) <= limit:
        return entries

    selected: list[BioLiPEntry] = []
    sequences: set[str] = set()
    ligands: set[str] = set()
    for entry in entries:
        if (len(ligands) < min(limit, ligand_diversity_target) and
                entry.sequence not in sequences and entry.ligand_comp_id not in ligands):
            selected.append(entry)
            sequences.add(entry.sequence)
            ligands.add(entry.ligand_comp_id)
        if len(selected) == limit:
            return selected
    selected_ids = {entry.source_entry_id for entry in selected}
    for entry in entries:
        if (entry.source_entry_id not in selected_ids and entry.sequence not in sequences and
                entry.ligand_comp_id in ligands):
            selected.append(entry)
            sequences.add(entry.sequence)
            selected_ids.add(entry.source_entry_id)
        if len(selected) == limit:
            return selected
    for entry in entries:
        if entry.source_entry_id not in selected_ids and entry.sequence not in sequences:
            selected.append(entry)
            sequences.add(entry.sequence)
            selected_ids.add(entry.source_entry_id)
        if len(selected) == limit:
            return selected
    selected.extend(entry for entry in entries if entry.source_entry_id not in selected_ids)
    return selected[:limit]
