import hashlib
import json
from pathlib import Path
import sqlite3

from scripts.source_affinity.canonicalize import (
    build_corpus, canonicalize_row, p_affinity, task_manifest,
)


def _raw(activity_id=1, ligand="AAAAAAAAAAAAAA-BBBBBBBBBB-C", value=10.0):
    return {
        "activity_id": activity_id, "assay_chembl_id": "CHEMBL_A",
        "document_chembl_id": "CHEMBL_D", "document_doi": "10/x",
        "document_pmid": 1, "document_year": 2020, "parent_molregno": activity_id,
        "molecule_chembl_id": f"CHEMBL_M{activity_id}", "canonical_smiles": "CCO",
        "standard_inchi_key": ligand, "target_chembl_id": "CHEMBL_T",
        "target_component_id": 7, "target_accession": "P1", "protein_sequence": "ACD",
        "target_variant_id": None, "target_variant_mutation": None,
        "endpoint_family": "Kd", "standard_relation": "=", "standard_value": value,
        "standard_units": "nM", "published_type": "Kd", "published_relation": "=",
        "published_value": value, "published_units": "nM", "pchembl_value_reported": 8.0,
        "bao_endpoint": None, "activity_comment": None, "src_id": 1,
        "assay_type": "B", "assay_confidence": 9, "assay_description": "binding",
        "assay_organism": "Homo sapiens", "relationship_type": "D", "bao_format": "BAO",
        "cell_id": None, "tissue_id": None, "subcellular_fraction": None,
        "variant_id": None,
    }


def test_p_affinity_unit_conversion():
    assert p_affinity(10.0, "nM") == 8.0
    assert p_affinity(0.01, "uM") == 8.0


def test_canonicalize_row_is_pair_local_and_deterministic():
    row = canonicalize_row(_raw(), [])
    assert row["protein_sequence_sha256"] == hashlib.sha256(b"ACD").hexdigest()
    assert row["canonical_smiles"] == "CCO"
    assert row["ligand_connectivity_key"] == "AAAAAAAAAAAAAA"
    assert row["p_affinity"] == 8.0
    assert len(row["task_id"]) == 64


def test_task_manifest_requires_compounds_and_non_tied_comparisons():
    rows = []
    for index in range(20):
        key = f"{index:014d}-BBBBBBBBBB-C"
        rows.append(canonicalize_row(_raw(index + 1, key, value=index + 1), []))
    task = task_manifest(rows)[0]
    assert task["exact_compound_count"] == 20
    assert task["non_tied_pair_comparisons"] == 190
    assert task["eligible_e0_core"] is True


def test_frozen_sql_builds_canonical_rows_from_static_schema(tmp_path):
    database = tmp_path / "chembl.db"
    connection = sqlite3.connect(database)
    connection.executescript("""
        CREATE TABLE activities (
          activity_id INTEGER, assay_id INTEGER, doc_id INTEGER, molregno INTEGER,
          standard_type TEXT, standard_relation TEXT, standard_value REAL,
          standard_units TEXT, type TEXT, relation TEXT, value REAL, units TEXT,
          pchembl_value REAL, bao_endpoint TEXT, activity_comment TEXT,
          standard_flag INTEGER, data_validity_comment TEXT,
          potential_duplicate INTEGER, src_id INTEGER);
        CREATE TABLE assays (
          assay_id INTEGER, chembl_id TEXT, doc_id INTEGER, tid INTEGER,
          assay_type TEXT, confidence_score INTEGER, description TEXT,
          assay_organism TEXT, relationship_type TEXT, bao_format TEXT,
          cell_id INTEGER, tissue_id INTEGER, assay_subcellular_fraction TEXT,
          variant_id INTEGER);
        CREATE TABLE target_dictionary (tid INTEGER, chembl_id TEXT, target_type TEXT);
        CREATE TABLE target_components (tid INTEGER, component_id INTEGER);
        CREATE TABLE component_sequences (
          component_id INTEGER, accession TEXT, sequence TEXT);
        CREATE TABLE variant_sequences (
          variant_id INTEGER, sequence TEXT, mutation TEXT);
        CREATE TABLE molecule_hierarchy (molregno INTEGER, parent_molregno INTEGER);
        CREATE TABLE molecule_dictionary (molregno INTEGER, chembl_id TEXT);
        CREATE TABLE compound_structures (
          molregno INTEGER, canonical_smiles TEXT, standard_inchi_key TEXT);
        CREATE TABLE docs (
          doc_id INTEGER, chembl_id TEXT, doi TEXT, pubmed_id INTEGER, year INTEGER);
        CREATE TABLE assay_parameters (
          assay_id INTEGER, standard_type TEXT, standard_relation TEXT,
          standard_value REAL, standard_units TEXT, standard_text_value TEXT);
        INSERT INTO activities VALUES
          (1, 2, 3, 4, 'Kd', '=', 10.0, 'nM', 'Kd', '=', 10.0, 'nM',
           8.0, NULL, NULL, 1, NULL, 0, 1);
        INSERT INTO assays VALUES
          (2, 'CHEMBL_A', 3, 5, 'B', 9, 'binding', 'Homo sapiens', 'D',
           'BAO', NULL, NULL, NULL, NULL);
        INSERT INTO target_dictionary VALUES (5, 'CHEMBL_T', 'SINGLE PROTEIN');
        INSERT INTO target_components VALUES (5, 6);
        INSERT INTO component_sequences VALUES (6, 'P1', 'ACD');
        INSERT INTO molecule_hierarchy VALUES (4, 7);
        INSERT INTO molecule_dictionary VALUES (7, 'CHEMBL_M');
        INSERT INTO compound_structures VALUES
          (7, 'CCO', 'AAAAAAAAAAAAAA-BBBBBBBBBB-C');
        INSERT INTO docs VALUES (3, 'CHEMBL_D', '10/x', 1, 2020);
        INSERT INTO assay_parameters VALUES (2, 'pH', '=', 7.4, NULL, NULL);
    """)
    connection.commit()
    connection.close()
    sql = Path(__file__).parents[2] / "contracts/source_affinity/chembl37_e0_core.sql"
    result = build_corpus(database, sql, tmp_path / "output")
    assert result["canonical_rows"]["rows"] == 1
    row = json.loads((tmp_path / "output/canonical_rows.jsonl").read_text())
    assert row["activity_id"] == 1
    assert row["assay_context"]["parameters"][0]["type"] == "pH"
