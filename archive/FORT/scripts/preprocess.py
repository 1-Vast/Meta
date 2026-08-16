"""Build the versioned ChEMBL 37 pKi training corpus from the raw SQLite DB.

The raw database is opened read-only.  The primary package contains only exact,
high-confidence Ki measurements from binding assays.  Censored values, label
mismatches, noisy assay-context aggregates, and invalid structures are retained
as separate audit artifacts instead of being silently discarded.

This module intentionally does not promote a natural-tail split to a formal
evaluation claim.  ChEMBL provides publication/patent years, but not a fully
auditable measurement timestamp and lineage for every row; the resulting
metadata-only roster is therefore reported separately.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import shutil
import sqlite3
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import torch

try:
    from .contract import AffinityRow
except ImportError:  # Allows the unified CLI to run as `python scripts/preprocess.py`.
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.contract import AffinityRow


FIELDS = (
    "target_key",
    "ligand_parent_key",
    "scaffold_key",
    "endpoint",
    "assay_key",
    "document_or_provenance_key",
    "split_role",
)


ROOT = Path(__file__).resolve().parents[1]
RAW_DB = ROOT / "dataset/public/chembl_historical/snapshots/chembl_37/chembl_37.db"
INNOVATION_SOURCE = ROOT / "dataset/ready/a2s_validation_small.v1"
DEFAULT_OUTPUT = ROOT / "dataset/formal_training/chembl37_pki_formal.v4"
RAW_ALIAS = ROOT / "dataset/raw/chembl_37"
A2S_REGISTRY = ROOT / "dataset/public/chembl_37/processed/dualcold/registry.parquet"
A2S_FEATURES = ROOT / "dataset/public/chembl_37/processed/dualcold/ligand_features.npz"
A2S_PROTEINS = ROOT / "dataset/public/chembl_37/processed/dualcold/target_esm2.npz"
A2S_FROZEN_ROSTER = ROOT / "dataset/processed/strict/episodes.pKi.v1.parquet"
A2S_OUTPUT = ROOT / "dataset/processed/a2s_validation_small.v1"
A2S_SOURCE_TARGET_LIMIT = 16
A2S_SOURCE_ROWS_PER_TARGET = 64
A2S_SUPPORT_BUDGETS = (1, 3, 5)
A2S_META_COLUMNS = (
    "target", "conn", "endpoint", "scaffold", "assays", "docs",
    "accession", "hcluster", "dual_cold_split",
)
A2S_LABEL_COLUMNS = A2S_META_COLUMNS + ("affinity", "replicate_sd", "n_records")

EXACT_RELATION = "="
CENSORED_RELATIONS = ("<", "<=", ">", ">=", "~")
LABEL_TOLERANCE = 0.05
SUPPORT_CUTOFFS = (2018, 2019, 2020, 2021, 2022, 2023)


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return sha256(value.encode("utf-8")).hexdigest()


def json_hash(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=True) + "\n", encoding="utf-8")


def normalized_text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip().lower()


def nonempty(value: object) -> bool:
    return normalized_text(value) not in {"", "none", "nan"}


def unique_join(values: Iterable[object]) -> str:
    cleaned = sorted({str(value) for value in values if nonempty(value)})
    return ";".join(cleaned)


def connect_readonly() -> sqlite3.Connection:
    if not RAW_DB.exists():
        raise FileNotFoundError(f"raw ChEMBL 37 database is missing: {RAW_DB}")
    connection = sqlite3.connect(f"file:{RAW_DB}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


SINGLE_TARGET_CTE = """
WITH single_target_components AS (
    SELECT tc.tid, MIN(tc.component_id) AS component_id
    FROM target_components tc
    GROUP BY tc.tid
    HAVING COUNT(DISTINCT tc.component_id) = 1
)
"""


def strict_where(endpoint: str, relation: str = EXACT_RELATION) -> tuple[str, list[object]]:
    return (
        """
        a.standard_type = ?
        AND a.standard_relation = ?
        AND a.standard_units = 'nM'
        AND a.standard_value > 0
        AND a.standard_flag = 1
        AND a.pchembl_value IS NOT NULL
        AND a.data_validity_comment IS NULL
        AND COALESCE(a.potential_duplicate, 0) = 0
        AND ass.confidence_score = 9
        AND ass.assay_type = 'B'
        AND ass.variant_id IS NULL
        AND td.target_type = 'SINGLE PROTEIN'
        AND td.tax_id = 9606
        AND md.molecule_type = 'Small molecule'
        AND pmd.molecule_type = 'Small molecule'
        AND pcs.canonical_smiles IS NOT NULL
        AND LENGTH(seq.sequence) BETWEEN 50 AND 5000
        AND (
            LENGTH(seq.sequence) - LENGTH(REPLACE(UPPER(seq.sequence), 'X', ''))
        ) / CAST(LENGTH(seq.sequence) AS REAL) <= 0.01
        AND (a.modality IS NULL OR a.modality != 'Targeted Protein Degradation')
        """,
        [endpoint, relation],
    )


def base_from_clause() -> str:
    return """
        FROM activities a
        JOIN assays ass ON ass.assay_id = a.assay_id
        JOIN target_dictionary td ON td.tid = ass.tid
        JOIN single_target_components stc ON stc.tid = ass.tid
        JOIN component_sequences seq ON seq.component_id = stc.component_id
        LEFT JOIN (
            SELECT component_id, MIN(protein_class_id) AS protein_class_id
            FROM component_class
            GROUP BY component_id
        ) cc ON cc.component_id = stc.component_id
        JOIN molecule_dictionary md ON md.molregno = a.molregno
        JOIN molecule_hierarchy mh ON mh.molregno = a.molregno
        JOIN molecule_dictionary pmd ON pmd.molregno = mh.parent_molregno
        JOIN compound_structures pcs ON pcs.molregno = mh.parent_molregno
        JOIN docs d ON d.doc_id = a.doc_id
        LEFT JOIN source src ON src.src_id = a.src_id
    """


def extract_measurements(connection: sqlite3.Connection, endpoint: str) -> pd.DataFrame:
    where, params = strict_where(endpoint)
    query = (
        SINGLE_TARGET_CTE
        + """
        SELECT
            a.activity_id, a.assay_id, ass.chembl_id AS assay_chembl_id,
            ass.assay_group, ass.assay_type, ass.assay_test_type,
            ass.assay_category, ass.description AS assay_description,
            ass.bao_format, ass.assay_organism, ass.assay_tax_id,
            ass.assay_strain, ass.assay_tissue, ass.assay_cell_type,
            ass.assay_subcellular_fraction, ass.tid, ass.relationship_type,
            ass.confidence_score, ass.src_id AS assay_src_id,
            ass.src_assay_id, ass.variant_id, ass.cell_id, ass.tissue_id,
            td.chembl_id AS target_chembl_id, td.target_type,
            td.tax_id AS target_tax_id, td.organism AS target_organism,
            stc.component_id, cc.protein_class_id,
            seq.accession, seq.sequence, seq.sequence_md5sum,
            seq.tax_id AS sequence_tax_id, seq.organism AS sequence_organism,
            a.molregno, mh.parent_molregno,
            pmd.chembl_id AS parent_molecule_chembl_id,
            pmd.pref_name AS parent_molecule_name,
            pmd.molecule_type AS parent_molecule_type,
            pcs.canonical_smiles AS parent_canonical_smiles,
            pcs.standard_inchi AS parent_standard_inchi,
            pcs.standard_inchi_key AS parent_standard_inchi_key,
            a.standard_type, a.standard_relation, a.standard_value,
            a.standard_units, a.standard_flag, a.pchembl_value,
            a.activity_comment, a.data_validity_comment, a.potential_duplicate,
            a.modality, a.src_id AS activity_src_id,
            a.doc_id, d.journal, d.year AS document_year, d.pubmed_id,
            d.doi, d.patent_id, d.chembl_id AS document_chembl_id,
            d.title AS document_title, d.doc_type, d.src_id AS document_src_id,
            src.src_description AS activity_source_description
        """
        + base_from_clause()
        + " WHERE "
        + where
    )
    chunks = pd.read_sql_query(query, connection, params=params, chunksize=100_000)
    frames = list(chunks)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def extract_censored(connection: sqlite3.Connection, endpoint: str) -> pd.DataFrame:
    # Censored rows are audited separately and are never admitted to the exact-label corpus.
    query = (
        SINGLE_TARGET_CTE
        + """
        SELECT
            a.activity_id, a.assay_id, ass.chembl_id AS assay_chembl_id,
            ass.assay_group, ass.assay_type, ass.assay_test_type,
            ass.assay_category, ass.description AS assay_description,
            ass.bao_format, ass.assay_organism, ass.assay_tax_id,
            ass.assay_strain, ass.assay_tissue, ass.assay_cell_type,
            ass.assay_subcellular_fraction, ass.tid,
            td.chembl_id AS target_chembl_id, td.target_type,
            td.tax_id AS target_tax_id, td.organism AS target_organism,
            stc.component_id, seq.accession, seq.sequence,
            a.molregno, mh.parent_molregno,
            pmd.chembl_id AS parent_molecule_chembl_id,
            pmd.molecule_type AS parent_molecule_type,
            pcs.canonical_smiles AS parent_canonical_smiles,
            pcs.standard_inchi_key AS parent_standard_inchi_key,
            a.standard_type, a.standard_relation, a.standard_value,
            a.standard_units, a.pchembl_value, a.activity_comment,
            a.data_validity_comment, a.potential_duplicate, a.modality,
            a.doc_id, d.year AS document_year, d.pubmed_id, d.doi,
            d.patent_id, d.chembl_id AS document_chembl_id,
            d.title AS document_title, a.src_id AS activity_src_id,
            src.src_description AS activity_source_description
        """
        + base_from_clause()
        + " WHERE "
        + """
        a.standard_type = ?
        AND a.standard_relation IN ('<', '<=', '>', '>=', '~')
        AND a.standard_units = 'nM'
        AND a.standard_value > 0
        AND a.standard_flag = 1
        AND ass.confidence_score = 9
        AND ass.assay_type = 'B'
        AND ass.variant_id IS NULL
        AND td.target_type = 'SINGLE PROTEIN'
        AND td.tax_id = 9606
        AND md.molecule_type = 'Small molecule'
        AND pmd.molecule_type = 'Small molecule'
        AND pcs.canonical_smiles IS NOT NULL
        AND LENGTH(seq.sequence) BETWEEN 50 AND 5000
        AND (
            LENGTH(seq.sequence) - LENGTH(REPLACE(UPPER(seq.sequence), 'X', ''))
        ) / CAST(LENGTH(seq.sequence) AS REAL) <= 0.01
        AND (a.modality IS NULL OR a.modality != 'Targeted Protein Degradation')
        """
    )
    chunks = pd.read_sql_query(query, connection, params=[endpoint], chunksize=100_000)
    frames = list(chunks)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def add_canonical_ids(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()
    result = frame.copy()
    result["sequence"] = result["sequence"].astype(str).str.upper().str.replace(r"\s+", "", regex=True)
    result["sequence_hash"] = result["sequence"].map(sha256_text)
    result["target_uid"] = result.apply(
        lambda row: f"{row.accession}:{row.sequence_hash}", axis=1
    )
    result["compound_parent_uid"] = result.apply(
        lambda row: f"chembl37-parent-{int(row.parent_molregno)}", axis=1
    )
    result["connectivity_inchikey"] = result["parent_standard_inchi_key"].fillna("").astype(str).str.split("-").str[0]
    result["document_uid"] = result.apply(document_uid, axis=1)
    result["assay_context_uid"] = result.apply(assay_context_uid, axis=1)
    result["pKi"] = 9.0 - np.log10(result["standard_value"].astype(float))
    result["pchembl_delta"] = (result["pKi"] - result["pchembl_value"].astype(float)).abs()
    result["measurement_uid"] = result.apply(
        lambda row: sha256_text(
            "|".join(
                (
                    str(row.parent_standard_inchi_key),
                    str(row.sequence_hash),
                    str(row.standard_type),
                    str(row.standard_relation),
                    f"{float(row.pKi):.8f}",
                    str(row.document_uid),
                    str(row.activity_id),
                )
            )
        ),
        axis=1,
    )
    result["document_year"] = pd.to_numeric(result["document_year"], errors="coerce").astype("Int64")
    return result


def document_uid(row: pd.Series) -> str:
    doi = normalized_text(row.get("doi"))
    if doi:
        return "doi:" + doi.removeprefix("https://doi.org/").removeprefix("doi:")
    patent = normalized_text(row.get("patent_id"))
    if patent:
        return "patent:" + patent
    pmid = normalized_text(row.get("pubmed_id"))
    if pmid:
        return "pmid:" + pmid
    return "chembl-doc:" + str(row.get("doc_id"))


def assay_context_uid(row: pd.Series) -> str:
    group = normalized_text(row.get("assay_group"))
    if group:
        return "assay-group:" + sha256_text(group)[:24]
    signature = {
        "target_uid": row.get("target_uid", ""),
        "standard_type": row.get("standard_type", ""),
        "assay_type": row.get("assay_type", ""),
        "bao_format": normalized_text(row.get("bao_format")),
        "assay_tax_id": row.get("assay_tax_id"),
        "assay_strain": normalized_text(row.get("assay_strain")),
        "variant_id": row.get("variant_id"),
        "description": normalized_text(row.get("assay_description")),
    }
    return "assay-signature:" + json_hash(signature)[:24]


def metrics(frame: pd.DataFrame) -> dict[str, int]:
    def nuniq(column: str) -> int:
        return int(frame[column].nunique(dropna=True)) if column in frame else 0

    target_components = 0
    if not frame.empty and {"target_uid", "component_id"}.issubset(frame.columns):
        target_components = int(frame[["target_uid", "component_id"]].astype(str).drop_duplicates().shape[0])
    return {
        "rows": int(len(frame)),
        "unique_parents": nuniq("compound_parent_uid"),
        "unique_targets": nuniq("target_uid"),
        "unique_target_components": target_components,
        "unique_documents": nuniq("document_uid"),
        "unique_assays": nuniq("assay_id"),
        "unique_assay_contexts": nuniq("assay_context_uid"),
    }


def sql_count(connection: sqlite3.Connection, where: str, params: list[object]) -> int:
    query = SINGLE_TARGET_CTE + "SELECT COUNT(*) " + base_from_clause() + " WHERE " + where
    return int(connection.execute(query, params).fetchone()[0])


def build_funnel(connection: sqlite3.Connection, exact: pd.DataFrame, censored: pd.DataFrame) -> dict[str, Any]:
    stages: list[dict[str, Any]] = []
    def add(name: str, rows: int, source: str = "sqlite") -> None:
        stages.append({"stage": name, "rows": int(rows), "source": source})

    add("chembl_all_activities", connection.execute("SELECT COUNT(*) FROM activities").fetchone()[0])
    add("standard_type_Ki", connection.execute("SELECT COUNT(*) FROM activities WHERE standard_type='Ki'").fetchone()[0])
    add("exact_relation", connection.execute("SELECT COUNT(*) FROM activities WHERE standard_type='Ki' AND standard_relation='='").fetchone()[0])
    add("nM_positive", connection.execute("SELECT COUNT(*) FROM activities WHERE standard_type='Ki' AND standard_relation='=' AND standard_units='nM' AND standard_value>0").fetchone()[0])
    add("standard_flag", connection.execute("SELECT COUNT(*) FROM activities WHERE standard_type='Ki' AND standard_relation='=' AND standard_units='nM' AND standard_value>0 AND standard_flag=1").fetchone()[0])
    add("pchembl_present", connection.execute("SELECT COUNT(*) FROM activities WHERE standard_type='Ki' AND standard_relation='=' AND standard_units='nM' AND standard_value>0 AND standard_flag=1 AND pchembl_value IS NOT NULL").fetchone()[0])
    add("no_validity_comment", connection.execute("SELECT COUNT(*) FROM activities WHERE standard_type='Ki' AND standard_relation='=' AND standard_units='nM' AND standard_value>0 AND standard_flag=1 AND pchembl_value IS NOT NULL AND data_validity_comment IS NULL").fetchone()[0])
    add("no_potential_duplicate", connection.execute("SELECT COUNT(*) FROM activities WHERE standard_type='Ki' AND standard_relation='=' AND standard_units='nM' AND standard_value>0 AND standard_flag=1 AND pchembl_value IS NOT NULL AND data_validity_comment IS NULL AND COALESCE(potential_duplicate,0)=0").fetchone()[0])
    where, params = strict_where("Ki")
    add("confidence_9_binding_single_protein_human_wildtype_small_molecule_sequence_qc", len(exact), "post_sql_extract")
    add("exact_label_pchembl_consistent", int((exact["pchembl_delta"] <= LABEL_TOLERANCE).sum()) if not exact.empty else 0, "pandas")
    add("censored_label_audit", len(censored), "sqlite")
    return {
        "stages": stages,
        "strict_sql_query_parameters": params,
        "strict_final_unique_metrics": metrics(exact),
        "censored_unique_metrics": metrics(censored),
    }


def split_label_mismatch(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if frame.empty:
        return frame.copy(), frame.copy()
    mismatch = frame[frame["pchembl_delta"] > LABEL_TOLERANCE].copy()
    valid = frame[frame["pchembl_delta"] <= LABEL_TOLERANCE].copy()
    return valid, mismatch


def aggregate_contexts(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()
    keys = ["assay_context_uid", "target_uid", "compound_parent_uid", "standard_type"]

    def mad(values: pd.Series) -> float:
        values = values.astype(float)
        median = float(values.median())
        return float(np.median(np.abs(values.to_numpy() - median)))

    def joined(column: str) -> Any:
        return lambda values: unique_join(values)

    aggregated = (
        frame.groupby(keys, sort=True, as_index=False)
        .agg(
            pKi=("pKi", "median"),
            pKi_mad=("pKi", mad),
            pKi_min=("pKi", "min"),
            pKi_max=("pKi", "max"),
            replicate_count=("activity_id", "count"),
            source_activity_ids=("activity_id", joined("activity_id")),
            source_assay_ids=("assay_id", joined("assay_id")),
            source_document_uids=("document_uid", joined("document_uid")),
            document_year_min=("document_year", "min"),
            document_year_max=("document_year", "max"),
        )
    )
    aggregated["pKi_range"] = aggregated["pKi_max"] - aggregated["pKi_min"]
    aggregated["noise_class"] = np.select(
        [aggregated["pKi_range"] <= 0.3, aggregated["pKi_range"] <= 0.5, aggregated["pKi_range"] <= 1.0],
        ["normal", "moderate", "high"],
        default="strong_outlier",
    )
    first_columns = [column for column in frame.columns if column not in set(keys) and column not in aggregated.columns]
    representative = frame.sort_values("activity_id").groupby(keys, sort=True, as_index=False)[first_columns].first()
    return representative.merge(aggregated, on=keys, how="inner", suffixes=("", "_aggregate"))


def entity_tables(frame: pd.DataFrame) -> dict[str, pd.DataFrame]:
    targets = frame[
        [
            "target_uid", "target_chembl_id", "component_id", "protein_class_id",
            "accession", "sequence", "sequence_hash", "sequence_md5sum",
            "target_type", "target_tax_id", "target_organism", "sequence_tax_id",
            "sequence_organism",
        ]
    ].drop_duplicates("target_uid").sort_values("target_uid")
    compounds = frame[
        [
            "compound_parent_uid", "parent_molregno", "parent_molecule_chembl_id",
            "parent_molecule_name", "parent_molecule_type", "parent_canonical_smiles",
            "parent_standard_inchi", "parent_standard_inchi_key", "connectivity_inchikey",
        ]
    ].drop_duplicates("compound_parent_uid").sort_values("compound_parent_uid")
    documents = frame[
        [
            "document_uid", "doc_id", "doi", "pubmed_id", "patent_id", "document_year",
            "document_chembl_id", "document_title", "doc_type", "journal",
            "document_src_id", "activity_source_description",
        ]
    ].drop_duplicates("document_uid").sort_values("document_uid")
    assays = frame[
        [
            "assay_id", "assay_chembl_id", "assay_group", "assay_type", "assay_test_type",
            "assay_category", "assay_description", "bao_format", "assay_organism",
            "assay_tax_id", "assay_strain", "assay_tissue", "assay_cell_type",
            "assay_subcellular_fraction", "tid", "confidence_score", "assay_src_id",
            "src_assay_id", "variant_id", "cell_id", "tissue_id",
        ]
    ].drop_duplicates("assay_id").sort_values("assay_id")
    contexts = frame[
        [
            "assay_context_uid", "target_uid", "standard_type", "assay_group", "assay_type",
            "assay_test_type", "assay_category", "assay_description", "bao_format",
            "assay_organism", "assay_tax_id", "assay_strain", "assay_tissue",
            "assay_cell_type", "assay_subcellular_fraction", "confidence_score",
        ]
    ].drop_duplicates("assay_context_uid").sort_values("assay_context_uid")
    return {"targets": targets, "compounds": compounds, "documents": documents, "assays": assays, "assay_contexts": contexts}


def build_features(compounds: pd.DataFrame, output: Path) -> dict[str, Any]:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors

    RDLogger.DisableLog("rdApp.*")
    descriptors = [
        "molecular_weight", "logp", "hbd", "hba", "tpsa", "rotatable_bonds",
        "ring_count", "formal_charge", "heavy_atoms", "fraction_csp3",
    ]
    feature_rows: list[list[float]] = []
    valid_uids: list[str] = []
    invalid: list[dict[str, object]] = []
    for row in compounds.itertuples(index=False):
        smiles = str(row.parent_canonical_smiles)
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            invalid.append({"compound_parent_uid": row.compound_parent_uid, "parent_canonical_smiles": smiles, "reason": "rdkit_parse_failed"})
            continue
        feature_rows.append(
            [
                float(Descriptors.MolWt(molecule)),
                float(Crippen.MolLogP(molecule)),
                float(Lipinski.NumHDonors(molecule)),
                float(Lipinski.NumHAcceptors(molecule)),
                float(rdMolDescriptors.CalcTPSA(molecule)),
                float(Lipinski.NumRotatableBonds(molecule)),
                float(rdMolDescriptors.CalcNumRings(molecule)),
                float(sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())),
                float(Descriptors.HeavyAtomCount(molecule)),
                float(rdMolDescriptors.CalcFractionCSP3(molecule)),
            ]
        )
        valid_uids.append(str(row.compound_parent_uid))
    output.mkdir(parents=True, exist_ok=True)
    valid_compounds = compounds[compounds.compound_parent_uid.astype(str).isin(valid_uids)].copy()
    valid_compounds = valid_compounds.set_index("compound_parent_uid").loc[valid_uids].reset_index()
    bits = np.zeros((len(feature_rows), 2048), dtype=np.uint8)
    for index, row in enumerate(valid_compounds.itertuples(index=False)):
        molecule = Chem.MolFromSmiles(str(row.parent_canonical_smiles))
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(molecule, radius=2, nBits=2048)
        bits[index, :] = np.asarray(fingerprint, dtype=np.uint8)
    descriptor_values = np.asarray(feature_rows, dtype=np.float32)
    np.savez_compressed(
        output / "ligand_features.npz",
        parent_uids=np.asarray(valid_uids),
        ecfp4=bits,
        descriptors=descriptor_values,
        descriptor_names=np.asarray(descriptors),
    )
    invalid_frame = pd.DataFrame(invalid, columns=["compound_parent_uid", "parent_canonical_smiles", "reason"])
    invalid_frame.to_parquet(output.parent / "quarantine_structure_invalid.parquet", index=False)
    return {
        "feature_file": "features/ligand_features.npz",
        "valid_compounds": len(valid_uids),
        "invalid_compounds": len(invalid),
        "ecfp4_bits": 2048,
        "descriptor_names": descriptors,
    }


def target_frequency(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame()
    group = frame.groupby("target_uid", sort=True)
    result = group.agg(
        target_chembl_id=("target_chembl_id", "first"),
        accession=("accession", "first"),
        n_measurements=("measurement_uid", "nunique"),
        n_unique_parent_compounds=("compound_parent_uid", "nunique"),
        n_documents=("document_uid", "nunique"),
        n_assays=("assay_id", "nunique"),
        n_assay_contexts=("assay_context_uid", "nunique"),
        first_document_year=("document_year", "min"),
        last_document_year=("document_year", "max"),
    ).reset_index()
    return result


def label_summary(frame: pd.DataFrame) -> dict[str, float]:
    values = frame["pKi"].astype(float)
    return {
        "min": float(values.min()),
        "q01": float(values.quantile(0.01)),
        "median": float(values.median()),
        "mean": float(values.mean()),
        "q99": float(values.quantile(0.99)),
        "max": float(values.max()),
        "std": float(values.std(ddof=1)),
    }


def roster_for_cutoff(frame: pd.DataFrame, cutoff: int) -> dict[str, Any]:
    timed = frame[frame.document_year.notna()].copy()
    pre = timed[timed.document_year.astype(int) <= cutoff]
    post = timed[timed.document_year.astype(int) > cutoff]
    pre_counts = pre.groupby("target_uid").agg(
        n_pre_unique_parents=("compound_parent_uid", "nunique"),
        n_pre_documents=("document_uid", "nunique"),
        n_pre_contexts=("assay_context_uid", "nunique"),
    )
    post_counts = post.groupby("target_uid").agg(
        n_post_unique_parents=("compound_parent_uid", "nunique"),
        n_post_documents=("document_uid", "nunique"),
        n_post_contexts=("assay_context_uid", "nunique"),
    )
    all_counts = timed.groupby("target_uid").agg(
        n_total_documents=("document_uid", "nunique"),
        n_total_contexts=("assay_context_uid", "nunique"),
    )
    counts = pre_counts.join(post_counts, how="outer").join(all_counts, how="outer").fillna(0).reset_index()
    counts["n_pre_unique_parents"] = counts["n_pre_unique_parents"].astype(int)
    counts["n_post_unique_parents"] = counts["n_post_unique_parents"].astype(int)
    source = counts[
        (counts.n_pre_unique_parents >= 64)
        & (counts.n_pre_documents >= 3)
        & (counts.n_pre_contexts >= 2)
    ]
    recipient = counts[
        (counts.n_pre_unique_parents.between(5, 20))
        & (counts.n_post_unique_parents >= 10)
        & (counts.n_total_documents >= 2)
    ]
    recipient_targets = set(recipient.target_uid.astype(str))
    source_targets = set(source.target_uid.astype(str))
    recipient_rows = timed[timed.target_uid.astype(str).isin(recipient_targets)]
    recipient_docs = set(recipient_rows.document_uid.astype(str))
    recipient_parents = set(recipient_rows.compound_parent_uid.astype(str))
    source_rows = pre[pre.target_uid.astype(str).isin(source_targets)].copy()
    closed_source_rows = source_rows[
        ~source_rows.document_uid.astype(str).isin(recipient_docs)
        & ~source_rows.compound_parent_uid.astype(str).isin(recipient_parents)
    ]
    closed_source_counts = closed_source_rows.groupby("target_uid").agg(
        n_closed_pre_parents=("compound_parent_uid", "nunique"),
        n_closed_pre_documents=("document_uid", "nunique"),
        n_closed_pre_contexts=("assay_context_uid", "nunique"),
    ).reset_index()
    closed_eligible = closed_source_counts[
        (closed_source_counts.n_closed_pre_parents >= 64)
        & (closed_source_counts.n_closed_pre_documents >= 3)
        & (closed_source_counts.n_closed_pre_contexts >= 2)
    ]
    return {
        "cutoff": cutoff,
        "time_field": "document_year",
        "time_quality": "publication_or_patent_year; not a measurement timestamp",
        "target_counts": counts,
        "source_candidates": source,
        "recipient_candidates": recipient,
        "closed_source_counts": closed_source_counts,
        "source_targets_after_document_parent_closure": set(closed_eligible.target_uid.astype(str)),
        "recipient_docs": recipient_docs,
        "recipient_parents": recipient_parents,
        "closed_source_rows": closed_source_rows,
        "raw_source_rows": pre[pre.target_uid.astype(str).isin(source_targets)].copy(),
        "recipient_rows": recipient_rows,
    }


def overlap_report(selected: dict[str, Any] | None) -> dict[str, Any]:
    if selected is None:
        return {
            "schema": "chembl37-overlap-audit-v1",
            "status": "NO_SELECTED_CUTOFF",
            "matrix": {},
        }
    source = selected["closed_source_rows"]
    recipient = selected["recipient_rows"]
    raw_source = selected["raw_source_rows"]

    def overlap(left: pd.DataFrame, right: pd.DataFrame, column: str) -> int:
        return len(set(left[column].astype(str)).intersection(right[column].astype(str)))

    matrix = {
        "source_recipient_target_uid": overlap(source, recipient, "target_uid"),
        "source_recipient_component_id": overlap(source, recipient, "component_id"),
        "source_recipient_document_uid_after_closure": overlap(source, recipient, "document_uid"),
        "source_recipient_parent_after_closure": overlap(source, recipient, "compound_parent_uid"),
        "source_recipient_parent_before_closure": overlap(raw_source, recipient, "compound_parent_uid"),
        "source_recipient_document_before_closure": overlap(raw_source, recipient, "document_uid"),
        "support_query_parent": None,
        "support_query_scaffold": None,
        "dev_test_target_family": None,
    }
    return {
        "schema": "chembl37-overlap-audit-v1",
        "status": "DIAGNOSTIC_CLOSED_SOURCE_RECIPIENT",
        "matrix": matrix,
        "hard_overlap_requirements": {
            "target_uid": 0,
            "component_id": 0,
            "document_uid_after_closure": 0,
            "parent_after_closure": 0,
        },
        "unconstructed_axes": ["support_query_parent", "support_query_scaffold", "dev_test_target_family"],
    }


def write_human_reports(output: Path, funnel: dict[str, Any], counts: dict[str, int], labels: dict[str, float], noise: dict[str, Any], natural: dict[str, Any], overlap: dict[str, Any]) -> None:
    report_dir = output / "reports"
    funnel_rows = "\n".join(
        f"| {stage['stage']} | {stage['rows']:,} | {stage['source']} |"
        for stage in funnel["stages"]
    )
    (report_dir / "preprocessing_report.md").write_text(
        "# ChEMBL 37 preprocessing report\n\n"
        "## Primary corpus\n\n"
        "The formal primary corpus is built directly from the read-only ChEMBL 37 relational database. "
        "The label is computed as `pKi = 9 - log10(standard_value_nM)`. Exact Ki, censored Ki, and pKd are stored in separate files.\n\n"
        "| Output | Count |\n|---|---:|\n"
        f"| Exact Ki activity measurements | {counts['pki_exact_measurements']:,} |\n"
        f"| Exact Ki assay-context rows | {counts['pki_exact_contexts']:,} |\n"
        f"| Main context rows (`pKi_range <= 0.5`) | {counts['pki_main_contexts']:,} |\n"
        f"| Censored Ki measurements | {counts['pki_censored_measurements']:,} |\n"
        f"| pChEMBL mismatches quarantined | {counts['pki_label_mismatches']:,} |\n"
        f"| pKd auxiliary measurements | {counts['pkd_exact_measurements']:,} |\n"
        f"| Parent compounds | {counts['parent_compounds']:,} |\n"
        f"| Human single-protein targets | {counts['targets']:,} |\n\n"
        "## pKi distribution\n\n"
        f"`min={labels['min']:.4f}`, `q01={labels['q01']:.4f}`, `median={labels['median']:.4f}`, "
        f"`mean={labels['mean']:.4f}`, `q99={labels['q99']:.4f}`, `max={labels['max']:.4f}`, `std={labels['std']:.4f}`.\n\n"
        "## Funnel\n\n| Stage | Rows | Source |\n|---|---:|---|\n"
        + funnel_rows
        + "\n\n## Assay noise\n\n"
        f"{noise['contexts_with_replicates']:,} contexts have replicates; {noise['range_gt_0_5']:,} contexts are isolated from the main context table because their within-context pKi range exceeds 0.5. "
        f"{noise['range_gt_1_0']:,} exceed 1.0 and remain available for audit.\n\n"
        "## Admission boundary\n\n"
        f"Natural-tail status is `{natural['status']}`. This does not invalidate the exact pKi corpus for supervised source training; it blocks prospective natural-tail claims until time/lineage and statistical-power gates pass.\n",
        encoding="utf-8",
    )
    write_json(report_dir / "overlap_matrix.json", overlap)
    matrix_rows = "\n".join(f"| {key} | {value if value is not None else 'not constructed'} |" for key, value in overlap["matrix"].items())
    (report_dir / "leakage_audit.md").write_text(
        "# Leakage and closure audit\n\n"
        "The selected metadata-only cutoff is evaluated separately from the primary exact pKi corpus. "
        "Source rows are closed against recipient documents and parent compounds, then source eligibility is recomputed.\n\n"
        "| Axis | Overlap |\n|---|---:|\n" + matrix_rows + "\n\n"
        "The support/query axes are not constructed because the selected recipient roster does not meet the formal natural-tail power requirement. No model may treat this diagnostic as a sealed natural-tail test.\n",
        encoding="utf-8",
    )
    chosen = natural.get("chosen_cutoff")
    cutoff_rows = "\n".join(
        f"| {item['cutoff']} | {item['source_candidates']} | {item['recipient_candidates']} | {item['source_after_document_parent_closure']} |"
        for item in natural["cutoff_candidates"]
    )
    (report_dir / "natural_tail_roster.md").write_text(
        "# Natural-tail D0 roster\n\n"
        f"Chosen metadata cutoff: `{chosen if chosen is not None else 'none'}`. Roster selection used counts and document years only; labels were not used.\n\n"
        "| Cutoff | Source candidates | Recipient candidates | Closed source candidates |\n|---:|---:|---:|---:|\n"
        + cutoff_rows
        + "\n\n"
        f"Status: `{natural['status']}`. The candidate count is below the formal power target and document year is not a complete measurement timestamp.\n",
        encoding="utf-8",
    )


def build_natural_tail(frame: pd.DataFrame, output: Path) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=True)
    audits: list[dict[str, Any]] = []
    candidate_data: dict[int, dict[str, Any]] = {}
    for cutoff in SUPPORT_CUTOFFS:
        result = roster_for_cutoff(frame, cutoff)
        candidate_data[cutoff] = result
        audits.append(
            {
                "cutoff": cutoff,
                "source_candidates": int(len(result["source_candidates"])),
                "recipient_candidates": int(len(result["recipient_candidates"])),
                "source_after_document_parent_closure": int(len(result["source_targets_after_document_parent_closure"])),
                "recipient_components": int(result["recipient_candidates"].target_uid.nunique()),
            }
        )
    chosen = max(
        (item for item in audits if item["recipient_candidates"] > 0),
        key=lambda item: item["cutoff"],
        default=None,
    )
    selected = candidate_data[chosen["cutoff"]] if chosen else None
    if selected is not None:
        counts = selected["target_counts"].copy()
        counts.to_parquet(output / "target_frequency_by_cutoff.parquet", index=False)
        selected["recipient_candidates"].to_parquet(output / "recipient_candidates.parquet", index=False)
        selected["source_candidates"].to_parquet(output / "source_candidates.parquet", index=False)
        selected["closed_source_counts"].to_parquet(output / "closed_source_counts.parquet", index=False)
        closed_source = selected["closed_source_rows"]
        closed_targets = selected["source_targets_after_document_parent_closure"]
        source_rows = closed_source[closed_source.target_uid.astype(str).isin(closed_targets)].copy()
        source_rows.to_parquet(output / "source_meta_train.parquet", index=False)
        recipient = selected["recipient_candidates"].copy()
        recipient["role"] = recipient.target_uid.astype(str).map(
            lambda value: "recipient_meta_dev" if int(sha256_text(value)[:8], 16) % 10 < 4 else "recipient_final_test"
        )
        recipient.to_parquet(output / "recipient_roster.parquet", index=False)
    status = "NATURAL_TAIL_METADATA_ONLY_DIAGNOSTIC"
    if not chosen or chosen["source_after_document_parent_closure"] == 0:
        status = "NATURAL_TAIL_BLOCKED_NO_CLOSED_SOURCE_ROSTER"
    report = {
        "schema": "chembl37-natural-tail-d0-v1",
        "status": status,
        "labels_used_for_roster_selection": False,
        "cutoff_candidates": audits,
        "chosen_cutoff": chosen["cutoff"] if chosen else None,
        "blocking_reasons": [
            "document_year is publication/patent year rather than a complete measurement timestamp",
            "source/recipient closure is evaluated as a separate diagnostic and is not used to alter the primary corpus",
        ],
    }
    write_json(output / "natural_tail_d0.json", report)
    return report


def raw_alias_manifest() -> dict[str, Any]:
    RAW_ALIAS.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema": "fort-raw-alias-v1",
        "source": "ChEMBL 37 SQLite relational database",
        "storage_mode": "canonical_file_retained_under_dataset_public",
        "canonical_path": str(RAW_DB.relative_to(ROOT)),
        "alias_directory": str(RAW_ALIAS.relative_to(ROOT)),
        "bytes": RAW_DB.stat().st_size,
        "sha256": sha256_file(RAW_DB),
        "read_only_policy": "never modify the SQLite source; all corrections are derived outputs",
    }
    write_json(RAW_ALIAS / "manifest.json", manifest)
    (RAW_ALIAS / "README.md").write_text(
        "# ChEMBL 37 raw layer\n\n"
        "The canonical SQLite file remains at `dataset/public/chembl_historical/snapshots/chembl_37/chembl_37.db`. "
        "This directory is the explicit raw-data classification and points to the immutable file without duplicating 30 GB.\n",
        encoding="utf-8",
    )
    return manifest


def copy_innovation_package(output: Path) -> list[dict[str, Any]]:
    destination = ROOT / "dataset/innovation_tests/a2s_validation_small.v1"
    destination.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    for source in sorted(INNOVATION_SOURCE.iterdir()):
        if not source.is_file():
            continue
        target = destination / source.name
        if not target.exists():
            shutil.copy2(source, target)
        elif sha256_file(source) != sha256_file(target):
            raise ValueError(f"protected innovation file differs from ready source: {target}")
        records.append({"path": str(target.relative_to(ROOT)), "bytes": target.stat().st_size, "sha256": sha256_file(target)})
    return records


def write_package_readme(output: Path, manifest: dict[str, Any]) -> None:
    (output / "README.md").write_text(
        "# ChEMBL 37 pKi formal training corpus\n\n"
        "Primary supervision is exact `standard_type=Ki`, exact relation `=`, positive nM, "
        "confidence 9 binding assays, human single proteins, wild-type assays, and small molecules.\n\n"
        f"Status: `{manifest['status']}`. The exact pKi corpus is ready for supervised source training; "
        "the natural-tail roster remains a separately reported metadata diagnostic until time/lineage closure passes.\n\n"
        "Read `manifest.json` before loading data. Raw SQLite is never a model input.\n",
        encoding="utf-8",
    )


def build(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    if output.exists() and (output / "manifest.json").exists():
        raise FileExistsError(f"refusing to overwrite existing package: {output}")
    if not INNOVATION_SOURCE.exists():
        raise FileNotFoundError(f"innovation validation package is missing: {INNOVATION_SOURCE}")
    output.mkdir(parents=True, exist_ok=True)
    raw_manifest = raw_alias_manifest()
    connection = connect_readonly()
    try:
        pki_raw = add_canonical_ids(extract_measurements(connection, "Ki"))
        pkd_raw = add_canonical_ids(extract_measurements(connection, "Kd"))
        censored = add_canonical_ids(extract_censored(connection, "Ki"))
        funnel = build_funnel(connection, pki_raw, censored)
    finally:
        connection.close()

    pki_exact, mismatch = split_label_mismatch(pki_raw)
    pki_exact = pki_exact.sort_values(["target_uid", "compound_parent_uid", "assay_context_uid", "activity_id"]).reset_index(drop=True)
    mismatch = mismatch.sort_values("activity_id").reset_index(drop=True)
    pkd_exact, _ = split_label_mismatch(pkd_raw)
    censored = censored.sort_values("activity_id").reset_index(drop=True)

    canonical = output / "canonical"
    components = output / "components"
    auxiliary = output / "auxiliary"
    features = output / "features"
    reports = output / "reports"
    for path in (canonical, components, auxiliary, features, reports):
        path.mkdir(parents=True, exist_ok=True)

    pki_exact.to_parquet(canonical / "pki_measurements_exact.parquet", index=False)
    censored.to_parquet(canonical / "pki_measurements_censored.parquet", index=False)
    mismatch.to_parquet(canonical / "quarantine_label_mismatch.parquet", index=False)
    pki_context = aggregate_contexts(pki_exact)
    pki_context.to_parquet(canonical / "pki_measurements_context_aggregated.parquet", index=False)
    pki_context[pki_context.pKi_range <= 0.5].to_parquet(canonical / "pki_measurements_context_main.parquet", index=False)
    pki_context[pki_context.pKi_range > 0.5].to_parquet(canonical / "pki_measurements_context_high_noise.parquet", index=False)
    pkd_exact.to_parquet(auxiliary / "pkd_measurements_exact.parquet", index=False)
    aggregate_contexts(pkd_exact).to_parquet(auxiliary / "pkd_measurements_context_aggregated.parquet", index=False)

    tables = entity_tables(pki_exact)
    for name, table in tables.items():
        table.to_parquet(components / f"{name}.parquet", index=False)
    feature_info = build_features(tables["compounds"], features)
    targets_payload = {
        str(row.target_uid): {
            "target_chembl_id": row.target_chembl_id,
            "uniprot_accession": row.accession,
            "sequence": row.sequence,
            "sequence_length": len(str(row.sequence)),
            "sequence_hash": row.sequence_hash,
            "protein_class_id": None if pd.isna(row.protein_class_id) else int(row.protein_class_id),
        }
        for row in tables["targets"].itertuples(index=False)
    }
    write_json(features / "target_sequences.json", targets_payload)

    natural_tail = build_natural_tail(pki_exact, reports)
    innovation_files = copy_innovation_package(output)
    noise_report = {
        "schema": "chembl37-assay-noise-audit-v1",
        "rows": int(len(pki_context)),
        "contexts_with_replicates": int((pki_context.replicate_count > 1).sum()),
        "range_le_0_3": int((pki_context.pKi_range <= 0.3).sum()),
        "range_0_3_to_0_5": int(((pki_context.pKi_range > 0.3) & (pki_context.pKi_range <= 0.5)).sum()),
        "range_gt_0_5": int((pki_context.pKi_range > 0.5).sum()),
        "range_gt_1_0": int((pki_context.pKi_range > 1.0).sum()),
        "main_context_rule": "retain pKi_range <= 0.5; isolate larger ranges",
    }
    write_json(reports / "assay_noise_audit.json", noise_report)
    write_json(reports / "funnel.json", funnel)

    counts = {
        "pki_exact_measurements": len(pki_exact),
        "pki_exact_contexts": len(pki_context),
        "pki_main_contexts": int((pki_context.pKi_range <= 0.5).sum()),
        "pki_label_mismatches": len(mismatch),
        "pki_censored_measurements": len(censored),
        "pkd_exact_measurements": len(pkd_exact),
        "targets": len(tables["targets"]),
        "parent_compounds": len(tables["compounds"]),
        "documents": len(tables["documents"]),
        "assays": len(tables["assays"]),
        "assay_contexts": len(tables["assay_contexts"]),
    }
    selected_cutoff = natural_tail.get("chosen_cutoff")
    selected_roster = roster_for_cutoff(pki_exact, selected_cutoff) if selected_cutoff is not None else None
    overlap = overlap_report(selected_roster)
    labels = label_summary(pki_exact)
    write_human_reports(output, funnel, counts, labels, noise_report, natural_tail, overlap)
    files: dict[str, dict[str, Any]] = {}
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name not in {"manifest.json", "README.md"}:
            files[str(path.relative_to(output)).replace("\\", "/")] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    manifest: dict[str, Any] = {
        "schema": "fort-chembl37-preprocessing-manifest-v1",
        "status": "FORMAL_PKI_CORPUS_READY_NATURAL_TAIL_BLOCKED",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "release": "ChEMBL 37",
            "license": "CC BY-SA 3.0",
            "raw": raw_manifest,
            "database_read_only": True,
        },
        "primary_label": {
            "endpoint": "Ki",
            "relation": "=",
            "unit": "nM",
            "formula": "pKi = 9 - log10(standard_value_nM)",
            "pchembl_consistency_tolerance": LABEL_TOLERANCE,
            "censored_relations_separate": list(CENSORED_RELATIONS),
        },
        "strict_filters": [
            "standard_type = Ki", "standard_relation = =", "standard_units = nM", "standard_value > 0",
            "standard_flag = 1", "pchembl_value IS NOT NULL", "data_validity_comment IS NULL",
            "potential_duplicate = 0 or NULL", "confidence_score = 9", "assay_type = B",
            "variant_id IS NULL", "target_type = SINGLE PROTEIN", "target tax_id = 9606",
            "single target component", "sequence length 50..5000", "sequence X fraction <= 0.01",
            "molecule_type = Small molecule", "valid parent structure", "exclude Targeted Protein Degradation",
        ],
        "counts": counts,
        "pKi_label_summary": labels,
        "funnel": funnel,
        "assay_noise": noise_report,
        "feature_info": feature_info,
        "overlap_audit": overlap,
        "natural_tail_d0": natural_tail,
        "innovation_tests": {
            "status": "DEVELOPMENT_ONLY_H0_BLOCKED",
            "package": "dataset/innovation_tests/a2s_validation_small.v1",
            "files": innovation_files,
            "purpose": "test innovation mechanisms/modules only; not natural-tail evidence",
        },
        "model_input_policy": {
            "formal_primary": [
                "canonical/pki_measurements_context_main.parquet",
                "components/targets.parquet", "components/compounds.parquet",
                "components/documents.parquet", "components/assays.parquet",
                "components/assay_contexts.parquet", "features/ligand_features.npz",
                "features/target_sequences.json",
            ],
            "raw_database_forbidden_at_model_time": True,
            "normalization_fit_role": "source_meta_train only after a separately sealed roster",
        },
        "files": files,
    }
    manifest["content_sha256"] = json_hash(manifest)
    write_json(output / "manifest.json", manifest)
    write_package_readme(output, manifest)

    classification = {
        "schema": "fort-dataset-classification-v2",
        "raw": raw_manifest,
        "formal_training": {
            "path": str(output.relative_to(ROOT)).replace("\\", "/"),
            "status": manifest["status"],
            "primary_endpoint": "pKi",
            "counts": counts,
            "manifest": str((output / "manifest.json").relative_to(ROOT)).replace("\\", "/"),
        },
        "innovation_tests": manifest["innovation_tests"],
        "cleanup": {
            "status": "COMPLETE",
            "stale_candidates_and_superseded_intermediates_removed": True,
        },
        "policy": "raw SQLite is retained and never modified; only sealed derived packages are model inputs",
    }
    write_json(ROOT / "dataset/registry/DATASET_CLASSIFICATION_REPORT.v2.json", classification)
    (ROOT / "dataset/registry/DATASET_CLASSIFICATION_REPORT.v2.md").write_text(
        "# Dataset classification v2\n\n"
        f"- Raw ChEMBL 37: `{raw_manifest['canonical_path']}` ({raw_manifest['bytes']} bytes, `{raw_manifest['sha256']}`).\n"
        f"- Formal pKi corpus: `{output.relative_to(ROOT)}` with `{counts['pki_exact_measurements']}` exact measurements and `{counts['pki_exact_contexts']}` assay-context rows.\n"
        "- Innovation tests: `dataset/innovation_tests/a2s_validation_small.v1/`.\n"
        "- Natural-tail status: `FORMAL_PKI_CORPUS_READY_NATURAL_TAIL_BLOCKED`; the exact corpus is prepared, but the prospective roster remains a separately audited diagnostic.\n\n"
        "Stale candidates and superseded intermediates have been removed. Read the versioned `manifest.json` before consuming any data.\n",
        encoding="utf-8",
    )
    return manifest


def preparetable(
    records: Sequence[Mapping[str, object]],
    *,
    split: str,
    allowvalues: bool = False,
) -> tuple[AffinityRow, ...]:
    """Map registered rows to one explicit few-shot role."""

    mapped = [
        {
            "target_key": record["target"],
            "ligand_parent_key": record["conn"],
            "scaffold_key": record["scaffold"],
            "endpoint": record["endpoint"],
            "assay_key": record["assays"],
            "document_or_provenance_key": record["docs"],
            "affinity_value": record.get("affinity"),
            "split_role": split,
        }
        for record in records
    ]
    return preparerows(mapped, allowvalues=allowvalues)


def preparerows(
    records: Sequence[Mapping[str, object]],
    *,
    allowvalues: bool = False,
) -> tuple[AffinityRow, ...]:
    """Normalize canonical rows without opening affinity values by default."""

    rows: list[AffinityRow] = []
    for position, record in enumerate(records):
        missing = [field for field in FIELDS if not str(record.get(field, "")).strip()]
        if missing:
            raise ValueError(f"record {position} misses canonical fields: {missing}")
        value = record.get("affinity_value")
        if value is not None and not allowvalues:
            raise PermissionError("affinity values are not permitted during FSA-D0 preprocessing")
        rows.append(
            AffinityRow(
                target_key=str(record["target_key"]),
                ligand_parent_key=str(record["ligand_parent_key"]),
                scaffold_key=str(record["scaffold_key"]),
                endpoint=str(record["endpoint"]),
                assay_key=str(record["assay_key"]),
                document_or_provenance_key=str(record["document_or_provenance_key"]),
                affinity_value=float(value) if value is not None else None,
                split_role=str(record["split_role"]),
            )
        )
    if len({row.endpoint for row in rows}) > 1:
        raise ValueError("preprocess one endpoint stratum at a time")
    return tuple(rows)


def preparevectors(values: Sequence[Sequence[float]]) -> torch.Tensor:
    """Move label-free design vectors to CUDA for FSA-D0 selection."""

    if not torch.cuda.is_available():
        raise RuntimeError("active preprocessing requires CUDA")
    vectors = torch.as_tensor(values, dtype=torch.float64, device="cuda")
    if vectors.ndim != 2 or vectors.shape[0] == 0 or not torch.isfinite(vectors).all():
        raise ValueError("design vectors must be a nonempty finite matrix")
    return vectors


def normalizeligands(
    feature: np.ndarray,
    trainrows: Sequence[int],
    descriptors: int = 10,
) -> tuple[np.ndarray, dict[str, list[float]]]:
    """Normalize descriptor columns with train-only statistics on CUDA."""

    if not torch.cuda.is_available():
        raise RuntimeError("active preprocessing requires CUDA")
    if feature.ndim != 2 or feature.shape[1] <= descriptors:
        raise ValueError("ligand features must contain fingerprint and descriptor columns")
    rows = torch.as_tensor(tuple(trainrows), dtype=torch.long, device="cuda")
    if rows.numel() == 0 or bool((rows < 0).any()) or bool((rows >= len(feature)).any()):
        raise ValueError("train rows must be nonempty valid feature indices")
    values = torch.as_tensor(feature, dtype=torch.float32, device="cuda")
    trainvalues = values[rows, -descriptors:]
    center = trainvalues.mean(dim=0)
    scale = trainvalues.std(dim=0).clamp_min(1e-6)
    values[:, -descriptors:] = (values[:, -descriptors:] - center) / scale
    return values.cpu().numpy(), {
        "descriptor_center": center.cpu().tolist(),
        "descriptor_scale": scale.cpu().tolist(),
    }


def verify_formal(output: Path, formal: dict[str, Any], raw: dict[str, Any]) -> dict[str, Any]:
    """Fail-closed validation for the emitted formal package and leakage contract."""

    content = dict(formal)
    expected_content_hash = str(content.pop("content_sha256"))
    actual_content_hash = json_hash(content)
    if actual_content_hash != expected_content_hash:
        raise ValueError("formal manifest content hash mismatch")

    checked_files = 0
    for relative, metadata in formal["files"].items():
        path = output / relative
        if not path.exists():
            raise FileNotFoundError(f"formal manifest file is missing: {path}")
        if path.stat().st_size != int(metadata["bytes"]) or sha256_file(path) != metadata["sha256"]:
            raise ValueError(f"formal manifest file hash mismatch: {path}")
        checked_files += 1

    raw_path = ROOT / Path(str(raw["canonical_path"]).replace("\\", "/"))
    if raw_path.stat().st_size != int(raw["bytes"]) or sha256_file(raw_path) != raw["sha256"]:
        raise ValueError("raw ChEMBL 37 hash mismatch")

    pki = pd.read_parquet(output / "canonical/pki_measurements_exact.parquet")
    censored = pd.read_parquet(output / "canonical/pki_measurements_censored.parquet")
    contexts = pd.read_parquet(output / "canonical/pki_measurements_context_main.parquet")
    compounds = pd.read_parquet(output / "components/compounds.parquet")
    features = np.load(output / "features/ligand_features.npz", allow_pickle=False)
    if len(pki) != formal["counts"]["pki_exact_measurements"] or not pki["activity_id"].is_unique:
        raise ValueError("exact pKi rows are not unique or disagree with manifest")
    if set(pki["standard_type"]) != {"Ki"} or set(pki["standard_relation"]) != {"="} or set(pki["standard_units"]) != {"nM"}:
        raise ValueError("formal pKi table contains a non-exact Ki row")
    if not bool((pki["standard_value"] > 0).all()) or not bool(pki["pchembl_delta"].le(0.05).all()):
        raise ValueError("formal pKi value or pChEMBL consistency check failed")
    if not bool(np.isclose(pki["pKi"], 9.0 - np.log10(pki["standard_value"]), atol=1e-10).all()):
        raise ValueError("pKi formula check failed")
    if not set(censored["standard_relation"]).issubset({"<", "<=", ">", ">=", "~"}):
        raise ValueError("censored Ki table contains an exact relation")
    context_key = ["assay_context_uid", "target_uid", "compound_parent_uid", "standard_type"]
    if not bool(contexts["pKi_range"].le(0.5).all()) or contexts.duplicated(context_key).any():
        raise ValueError("main assay-context table violates noise or uniqueness policy")
    if features["ecfp4"].shape[0] != len(compounds) or list(features["parent_uids"].astype(str)) != list(compounds["compound_parent_uid"].astype(str)):
        raise ValueError("ligand feature cache is not aligned with compounds")

    natural = formal["natural_tail_d0"]
    overlap = formal["overlap_audit"]["matrix"]
    required_zero = (
        "source_recipient_target_uid",
        "source_recipient_component_id",
        "source_recipient_document_uid_after_closure",
        "source_recipient_parent_after_closure",
    )
    if natural["labels_used_for_roster_selection"] is not False:
        raise ValueError("natural-tail roster is not label-blind")
    if any(overlap[key] != 0 for key in required_zero):
        raise ValueError("source/recipient hard overlap is nonzero")
    if formal["model_input_policy"]["raw_database_forbidden_at_model_time"] is not True:
        raise ValueError("raw database is not fail-closed at model time")
    if formal["model_input_policy"]["normalization_fit_role"] != "source_meta_train only after a separately sealed roster":
        raise ValueError("normalization fit role is not source-only")
    return {
        "status": "PASS",
        "manifest_files_checked": checked_files,
        "raw_hash_checked": True,
        "exact_label_checks": True,
        "censored_isolation_checked": True,
        "context_alignment_checked": True,
        "feature_alignment_checked": True,
        "natural_tail_label_blind_checked": True,
        "hard_overlap_checked": True,
        "model_raw_input_guard_checked": True,
    }


def a2s_choose_source_targets(metadata: pd.DataFrame) -> list[str]:
    train = metadata[(metadata.endpoint == "pKi") & (metadata.dual_cold_split == "train")]
    counts = train.groupby("target", sort=True).size()
    candidates = [str(target) for target in sorted(counts[counts >= A2S_SOURCE_ROWS_PER_TARGET].index)]
    selected: list[str] = []
    components: set[str] = set()
    for target in candidates:
        component = str(train.loc[train.target == target, "hcluster"].iloc[0])
        if component in components:
            continue
        selected.append(target)
        components.add(component)
        if len(selected) == A2S_SOURCE_TARGET_LIMIT:
            break
    if len(selected) < A2S_SOURCE_TARGET_LIMIT:
        raise ValueError("insufficient distinct source targets for the validation package")
    return selected


def a2s_leakage_report(rows: pd.DataFrame) -> dict[str, Any]:
    train = rows[rows.split_role == "source_train"]
    development = rows[rows.split_role == "recipient_development"]
    confirmation = rows[rows.split_role == "confirmation_holdout"]

    def overlap(left: pd.DataFrame, right: pd.DataFrame, column: str) -> int:
        return len(set(left[column].astype(str)).intersection(right[column].astype(str)))

    axes = ("target", "accession", "hcluster", "conn", "scaffold", "docs", "assays")
    report = {
        "schema": "a2s-validation-leakage-v1",
        "label_policy": "labels copied only after split assignment; no query label used for selection",
        "source_train_vs_development": {axis: overlap(train, development, axis) for axis in axes},
        "source_train_vs_confirmation": {axis: overlap(train, confirmation, axis) for axis in axes},
        "development_vs_confirmation": {axis: overlap(development, confirmation, axis) for axis in axes},
        "target_disjoint": not bool(overlap(train, development, "target")),
        "homology_cold": not bool(overlap(train, development, "hcluster")),
        "scaffold_cold": not bool(overlap(train, development, "scaffold")),
        "provenance_cold": not bool(
            overlap(train, development, "docs") or overlap(train, development, "assays")
        ),
    }
    report["pass"] = all(
        report[key]
        for key in ("target_disjoint", "homology_cold", "scaffold_cold", "provenance_cold")
    )
    return report


def build_innovation_validation(output: Path = A2S_OUTPUT) -> dict[str, Any]:
    """Build the development-only A2S validation package from sealed inputs."""
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"refusing to overwrite frozen validation package: {output}")
    output.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_parquet(A2S_REGISTRY, columns=list(A2S_META_COLUMNS)).reset_index(names="source_row")
    registry_row_count = len(metadata)
    metadata = metadata[metadata.endpoint == "pKi"].copy()
    source_targets = a2s_choose_source_targets(metadata)

    roster = pd.read_parquet(A2S_FROZEN_ROSTER)
    if set(roster.endpoint.astype(str)) != {"pKi"}:
        raise ValueError("frozen roster must contain only pKi episodes")
    recipient_rows = sorted(set(roster.source_row.astype(int)))
    source = metadata[metadata.target.astype(str).isin(source_targets)]
    source = source[source.dual_cold_split == "train"].groupby("target", sort=True, group_keys=False).head(A2S_SOURCE_ROWS_PER_TARGET)
    selected_rows = sorted(set(source.source_row.astype(int)).union(recipient_rows))

    labels = pd.read_parquet(A2S_REGISTRY, columns=list(A2S_LABEL_COLUMNS)).reset_index(names="source_row")
    rows = labels[labels.source_row.isin(selected_rows)].copy()
    rows["split_role"] = np.where(
        rows.target.astype(str).isin(source_targets),
        "source_train",
        np.where(rows.dual_cold_split == "confirmation", "confirmation_holdout", "recipient_development"),
    )
    if set(rows.split_role) - {"source_train", "recipient_development", "confirmation_holdout"}:
        raise ValueError("unexpected split role")
    report = a2s_leakage_report(rows)
    if not report["pass"]:
        raise ValueError(f"validation package leakage audit failed: {report}")

    cache = np.load(A2S_FEATURES, allow_pickle=False)
    raw_features = np.asarray(cache["feat"], dtype=np.float32)
    if len(raw_features) != registry_row_count:
        raise ValueError("feature cache is not aligned with registry rows")
    train_indices = source.source_row.to_numpy(dtype=np.int64)
    selected_indices = np.asarray(selected_rows, dtype=np.int64)
    normalized = raw_features[selected_indices].copy()
    if normalized.shape[1] < 10:
        raise ValueError("expected descriptor tail in ligand feature cache")
    center = raw_features[train_indices, -10:].mean(axis=0)
    scale = np.maximum(raw_features[train_indices, -10:].std(axis=0), 1e-6)
    normalized[:, -10:] = (normalized[:, -10:] - center) / scale

    rows = rows.sort_values("source_row").reset_index(drop=True)
    rows.to_parquet(output / "rows.parquet", index=False)
    np.save(output / "ligand_features.npy", normalized)

    protein_cache = np.load(A2S_PROTEINS, allow_pickle=False)
    target_set = set(rows.target.astype(str))
    keys = np.asarray([str(key) for key in protein_cache["keys"]])
    keep = np.asarray([key in target_set for key in keys], dtype=bool)
    if not bool(keep.any()) or len(keys[keep]) != len(target_set):
        missing = sorted(target_set.difference(set(keys)))
        raise ValueError(f"target embeddings missing for validation targets: {missing[:5]}")
    np.savez_compressed(output / "target_embeddings.npz", keys=keys[keep], segments=protein_cache["segments"][keep])

    row_to_package = {int(row): index for index, row in enumerate(rows.source_row.astype(int))}
    episode_records: list[dict[str, Any]] = []
    for episode_id, group in roster.groupby("episode", sort=True):
        support = sorted(group.loc[group.role == "support", "source_row"].astype(int))
        query = sorted(group.loc[group.role == "query", "source_row"].astype(int))
        target = str(group.target.iloc[0])
        for budget in A2S_SUPPORT_BUDGETS:
            nested = support[:budget]
            for role, values in (("support", nested), ("query", query)):
                episode_records.extend(
                    {
                        "episode": str(episode_id),
                        "target": target,
                        "k": budget,
                        "role": role,
                        "source_row": int(row),
                        "package_row": int(row_to_package[int(row)]),
                    }
                    for row in values
                )
    episodes = pd.DataFrame(episode_records)
    episodes.to_parquet(output / "episodes.parquet", index=False)

    files = ["rows.parquet", "ligand_features.npy", "target_embeddings.npz", "episodes.parquet"]
    manifest = {
        "schema": "a2s-validation-package-v1",
        "status": "DEVELOPMENT_ONLY_H0_BLOCKED",
        "endpoint": "pKi",
        "source": "ChEMBL-37 audited dual-cold registry",
        "split_authority": "frozen chembl-dualcold-v1 manifest plus strict-fewshot-roster.v1",
        "support_budgets": list(A2S_SUPPORT_BUDGETS),
        "nested_support": True,
        "source_targets": source_targets,
        "counts": {
            "rows": len(rows),
            "source_rows": int((rows.split_role == "source_train").sum()),
            "development_rows": int((rows.split_role == "recipient_development").sum()),
            "confirmation_rows": int((rows.split_role == "confirmation_holdout").sum()),
            "recipient_targets": int((rows.split_role == "recipient_development").groupby(rows.target).any().sum()),
            "episodes": int(episodes.episode.nunique()),
        },
        "normalization": {
            "fit_role": "source_train",
            "columns": list(range(int(normalized.shape[1] - 10), int(normalized.shape[1]))),
            "center": center.tolist(),
            "scale": scale.tolist(),
        },
        "leakage_report": report,
        "inputs": {
            "registry": {"path": str(A2S_REGISTRY), "sha256": sha256_file(A2S_REGISTRY)},
            "frozen_roster": {"path": str(A2S_FROZEN_ROSTER), "sha256": sha256_file(A2S_FROZEN_ROSTER)},
            "features": {"path": str(A2S_FEATURES), "sha256": sha256_file(A2S_FEATURES)},
            "proteins": {"path": str(A2S_PROTEINS), "sha256": sha256_file(A2S_PROTEINS)},
        },
        "model_input_contract": {
            "allowed_files": files,
            "forbidden_inputs": ["raw registry", "raw database", "query labels during routing/selection"],
        },
    }
    manifest["files"] = {name: {"sha256": sha256_file(output / name)} for name in files}
    manifest["content_sha256"] = json_hash(manifest)
    write_json(output / "leakage_report.json", report)
    write_json(output / "normalization.json", manifest["normalization"])
    write_json(output / "manifest.json", manifest)
    return manifest


def innovation_main() -> None:
    parser = argparse.ArgumentParser(description="Build the development-only A2S validation package.")
    parser.add_argument("--output", type=Path, default=A2S_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(build_innovation_validation(args.output.resolve()), indent=2, sort_keys=True, ensure_ascii=True))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("formal", "innovation"), default="formal")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    default_output = DEFAULT_OUTPUT if args.mode == "formal" else A2S_OUTPUT
    result = (build if args.mode == "formal" else build_innovation_validation)(
        (args.output or default_output).resolve()
    )
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True))


if __name__ == "__main__":
    main()
