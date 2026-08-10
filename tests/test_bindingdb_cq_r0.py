import csv
import gzip
import json
import zipfile
from pathlib import Path

from research.crossed_interaction.bindingdb_cq_r0 import (
    extract_projection,
    run_census,
)


HEADER = [
    "BindingDB Reactant_set_id",
    "Ligand SMILES",
    "Ligand InChI Key",
    "Target Name",
    "Target Source Organism According to Curator or DataSource",
    "Ki (nM)",
    "IC50 (nM)",
    "Kd (nM)",
    "EC50 (nM)",
    "kon (M-1-s-1)",
    "koff (s-1)",
    "pH",
    "Temp (C)",
    "Curation/DataSource",
    "Article DOI",
    "BindingDB Entry DOI",
    "PMID",
    "PubChem AID",
    "Patent Number",
    "Authors",
    "Date of publication",
    "Date in BindingDB",
    "Institution",
    "Link to Ligand in BindingDB",
    "Link to Target in BindingDB",
    "Link to Ligand-Target Pair in BindingDB",
    "Ligand HET ID in PDB",
    "PDB ID(s) for Ligand-Target Complex",
    "PubChem CID",
    "PubChem SID",
    "ChEBI ID of Ligand",
    "ChEMBL ID of Ligand",
    "DrugBank ID of Ligand",
    "IUPHAR_GRAC ID of Ligand",
    "KEGG ID of Ligand",
    "ZINC ID of Ligand",
    "Number of Protein Chains in Target (>1 implies a multichain complex)",
    "BindingDB Target Chain Sequence 1",
]


def write_zip(path: Path, name: str, rows):
    text = "\n".join("\t".join(row) for row in rows) + "\n"
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(name, text)


def test_projection_does_not_expose_affinity(tmp_path):
    articles = tmp_path / "articles.zip"
    assays = tmp_path / "assays.zip"
    projection = tmp_path / "projection.jsonl.gz"
    row = [""] * len(HEADER)
    values = {
        "BindingDB Reactant_set_id": "1",
        "Ligand SMILES": "CC",
        "Ligand InChI Key": "L1",
        "Target Name": "T1",
        "Target Source Organism According to Curator or DataSource": "Human",
        "Ki (nM)": "SECRET_12.5",
        "Curation/DataSource": "BindingDB",
        "Article DOI": "10.test/a",
        "Date of publication": "2020",
        "Number of Protein Chains in Target (>1 implies a multichain complex)": "1",
        "BindingDB Target Chain Sequence 1": "AAAA",
    }
    for key, value in values.items():
        row[HEADER.index(key)] = value
    write_zip(articles, "articles.tsv", [HEADER, row])
    write_zip(
        assays,
        "assays.tsv",
        [["ENTRYID", "ASSAYID", "ASSAY_NAME", "DESCRIPTION"], ["1", "7", "binding", "same protocol"]],
    )
    manifest = extract_projection(articles, assays, projection)
    payload = gzip.open(projection, "rt", encoding="utf-8").read()
    assert "SECRET_12.5" not in payload
    assert manifest["numeric_affinity_values_parsed"] == 0


def test_document_protocol_recovers_cycle_but_entry_assay_does_not(tmp_path):
    projection = tmp_path / "projection.jsonl.gz"
    rows = []
    for target in ("AAAA", "BBBB"):
        for ligand in ("L1", "L2"):
            rows.append(
                {
                    "source_row_id": f"{target}-{ligand}",
                    "document_id": "doi:10.test/panel",
                    "endpoint_available": ["Ki"],
                    "chain_count": 1,
                    "target_sequence": target,
                    "target_sequence_sha256": target,
                    "ligand_inchikey": ligand,
                    "ligand_smiles": "",
                    "assays": [{"assay_id": "1", "assay_name_norm": "binding", "protocol_sha256": "p"}],
                }
            )
    with gzip.GzipFile(filename="", mode="wb", fileobj=projection.open("wb"), mtime=0) as raw:
        with __import__("io").TextIOWrapper(raw, encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(row) + "\n")
    result = run_census(projection)
    assert result["panel_definitions"]["document_protocol"]["total_cycle_rank"] == 1
    assert result["panel_definitions"]["entry_assay_negative_control"]["total_cycle_rank"] == 0


def test_zero_cycle_census_is_defined(tmp_path):
    projection = tmp_path / "projection.jsonl.gz"
    row = {
        "source_row_id": "1",
        "document_id": "doi:x",
        "endpoint_available": ["Kd"],
        "chain_count": 1,
        "target_sequence": "AAAA",
        "target_sequence_sha256": "t",
        "ligand_inchikey": "l",
        "ligand_smiles": "",
        "assays": [],
    }
    with gzip.open(projection, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(row) + "\n")
    result = run_census(projection)
    assert result["panel_definitions"]["document_protocol"]["median_positive_cycle_rank"] == 0
    assert not result["development_training_ready_preclosure"]
