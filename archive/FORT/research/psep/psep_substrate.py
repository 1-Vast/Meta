"""PSEP: a provenance-separated adaptation substrate built from Papyrus 05.7++.

The A2S programme closed with `NO_CHEMICAL_ADAPTATION_OBJECT_SURVIVES_SEPARATION`
at 92 homology components, and registered a single reopening condition: not a new
architecture, but *more independent components* -- approximately 445 -- with
Papyrus 05.7 named as the realistic route.  This module builds that corpus.

Design decisions that matter, and why:

``unit``        The adaptation problem is one (accession, endpoint) pair.  Ki and
                IC50 for the same protein are different estimands (IC50 carries a
                substrate-concentration term), so they are never mixed inside a
                unit.  They may share a component; the component is the
                independent statistical unit and collapses them correctly.

``component``   Homology cluster, single-linkage over 4-mer Jaccard at tau=0.20.
                The threshold is *calibrated against the existing ChEMBL
                `hcluster` partition*, not chosen freely: on 863 shared
                accessions the within-cluster minimum Jaccard is 0.224 and the
                between-cluster maximum is 0.248, so tau=0.20 never splits an
                existing cluster and merges slightly more.  Merging more means
                *fewer* components and *less* power, which is the conservative
                direction.

``firewall``    Strictly stronger than the original.  Every accession carrying a
                `probe` or `locked` role in the A2S source lock, and every
                accession on the recipient roster, is sealed -- and so is every
                accession in the same homology component as a sealed one.  The
                ChEMBL `locked` set and the recipient roster are never read for
                labels here; only their accessions and sequences are used, to
                exclude them.

``role``        The surviving components are split discover/validate/confirm at
                60/20/20 by component hash.  Mechanism search may read
                `discover`.  `confirm` stays sealed for a mechanism that passes
                the admission gate.

No labels from any sealed role are read.  The base model is target-agnostic and
cross-fitted by component, so it cannot memorise a target.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, save_npz
from scipy.sparse.csgraph import connected_components
from sklearn.feature_extraction.text import CountVectorizer

ROOT = Path(__file__).resolve().parents[2]
PAPYRUS = ROOT / "dataset" / "public" / "papyrus_05_7" / "raw"
ACTIVITIES = PAPYRUS / "05.7++_combined_set_without_stereochemistry.tsv.xz"
PROTEINS = PAPYRUS / "05.7_combined_set_protein_targets.tsv.xz"
CHEMBL_TARGETS = (
    ROOT / "dataset" / "formal_training" / "chembl37_pki_formal.v4" / "components" / "targets.parquet"
)
SOURCE_LOCK = ROOT / "reports" / "active" / "a2s_source_information_gate_lock_v2_2026-08-01.json"
RECIPIENTS = ROOT / "dataset" / "formal_training" / "a2s_d0r_roster.v3" / "recipients.parquet"
DEFAULT_OUTPUT = ROOT / "dataset" / "processed" / "psep.v1"

SEED = 20260802
TAU = 0.20                 # 4-mer Jaccard, calibrated against ChEMBL hcluster
MORGAN_BITS = 1024
MORGAN_RADIUS = 2
MIN_EVAL_ROWS = 8          # a document must hold this many rows to be splittable
MIN_TRAIN_ROWS = 20
ROLE_FRACTIONS = (("discover", 0.60), ("validate", 0.20), ("confirm", 0.20))
ENDPOINTS = {"pKi": "type_Ki", "pIC50": "type_IC50", "pKd": "type_KD"}

DESCRIPTORS = (
    "MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
    "NumRotatableBonds", "RingCount", "FractionCSP3", "HeavyAtomCount",
    "NumAromaticRings",
)


# --------------------------------------------------------------------------- #
# Activities
# --------------------------------------------------------------------------- #


def load_activities() -> pd.DataFrame:
    columns = [
        "Quality", "source", "connectivity", "SMILES", "target_id", "accession",
        "Protein_Type", "AID", "doc_id", "all_doc_ids", "Year", "relation",
        "pchembl_value_Mean", "pchembl_value_N", *ENDPOINTS.values(),
    ]
    frame = pd.read_csv(ACTIVITIES, sep="\t", usecols=columns, low_memory=False)
    for column in ENDPOINTS.values():
        frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0)
    frame["pchembl_value_Mean"] = pd.to_numeric(frame.pchembl_value_Mean, errors="coerce")
    frame = frame[
        (frame.Protein_Type == "WT")
        & (frame.relation == "=")
        & frame.pchembl_value_Mean.notna()
        & frame.SMILES.notna()
    ]
    parts = []
    for endpoint, flag in ENDPOINTS.items():
        # Papyrus never mixes activity types inside a row, so a positive type
        # flag is an exact endpoint assignment rather than a majority vote.
        part = frame.loc[frame[flag] > 0].copy()
        part["endpoint"] = endpoint
        parts.append(part)
    activities = pd.concat(parts, ignore_index=True)
    activities["unit"] = activities.accession.astype(str) + "|" + activities.endpoint
    activities["docs"] = activities.doc_id.astype(str)
    activities["assays"] = activities.AID.astype(str)
    activities["affinity"] = activities.pchembl_value_Mean.astype(np.float64)
    return activities


def eligible_units(activities: pd.DataFrame) -> set[str]:
    """Units deep enough for a document-disjoint within-unit split."""

    keep: set[str] = set()
    for unit, group in activities.groupby("unit", sort=True):
        if len(group) < MIN_TRAIN_ROWS + MIN_EVAL_ROWS:
            continue
        counts = group.docs.value_counts()
        if int((counts >= MIN_EVAL_ROWS).sum()) >= 2:
            keep.add(str(unit))
    return keep


# --------------------------------------------------------------------------- #
# Homology components and the firewall
# --------------------------------------------------------------------------- #


def sealed_accessions() -> tuple[set[str], dict[str, object]]:
    lock = json.loads(SOURCE_LOCK.read_text(encoding="utf-8"))
    assignment = lock["components"]["assignment"]
    role_by_component = lock["components"]["role_by_component"]
    targets = pd.read_parquet(CHEMBL_TARGETS)[["target_chembl_id", "accession"]]
    to_accession = targets.groupby("target_chembl_id").accession.first()

    sealed_roles = {"probe", "locked"}
    chembl_sealed = [
        target for target, component in assignment.items()
        if role_by_component.get(component) in sealed_roles
    ]
    accessions = set(to_accession.reindex(chembl_sealed).dropna().tolist())
    recipients = set(pd.read_parquet(RECIPIENTS).accession.astype(str).tolist())
    provenance = {
        "chembl_sealed_targets": len(chembl_sealed),
        "chembl_sealed_accessions": len(accessions),
        "recipient_accessions": len(recipients),
        "sealed_roles": sorted(sealed_roles),
    }
    return accessions | recipients, provenance


def sequence_table() -> dict[str, str]:
    proteins = pd.read_csv(PROTEINS, sep="\t", low_memory=False)
    mapping: dict[str, str] = {}
    for target_id, sequence in zip(proteins.target_id, proteins.Sequence):
        if isinstance(sequence, str):
            mapping.setdefault(str(target_id).split("_")[0], sequence)
    chembl = pd.read_parquet(CHEMBL_TARGETS)[["accession", "sequence"]].dropna()
    for accession, sequence in zip(chembl.accession, chembl.sequence):
        mapping.setdefault(str(accession), str(sequence))
    return mapping


def homology_components(
    accessions: list[str], sealed: list[str], sequences: dict[str, str]
) -> tuple[np.ndarray, set[int], dict[str, object]]:
    """Single-linkage over 4-mer Jaccard.  Sealed accessions join the graph so
    that homology to a sealed protein is itself disqualifying."""

    nodes = accessions + sealed
    vectoriser = CountVectorizer(analyzer="char", ngram_range=(4, 4), binary=True, lowercase=False)
    matrix = vectoriser.fit_transform([sequences[node] for node in nodes]).astype(np.float32)
    matrix.data[:] = 1.0
    intersection = (matrix @ matrix.T).toarray()
    sizes = np.asarray(matrix.sum(axis=1)).ravel()
    jaccard = intersection / np.maximum(sizes[:, None] + sizes[None, :] - intersection, 1.0)
    count, labels = connected_components(csr_matrix(jaccard >= TAU), directed=False)
    sealed_labels = set(labels[len(accessions):].tolist())
    stats = {
        "tau": TAU,
        "kmer": 4,
        "linkage": "single",
        "nodes": len(nodes),
        "components_total": int(count),
        "sealed_components": len(sealed_labels),
        "calibration": "tau chosen so no existing ChEMBL hcluster is split (within-min 0.224 > tau, between-max 0.248)",
    }
    return labels[: len(accessions)], sealed_labels, stats


def assign_role(component: int) -> str:
    digest = int(sha256(f"{SEED}:role:{component}".encode()).hexdigest()[:8], 16)
    draw = (digest % 10_000) / 10_000.0
    cumulative = 0.0
    for name, fraction in ROLE_FRACTIONS:
        cumulative += fraction
        if draw < cumulative:
            return name
    return ROLE_FRACTIONS[-1][0]


# --------------------------------------------------------------------------- #
# Chemistry
# --------------------------------------------------------------------------- #


def featurise(smiles: list[str]) -> tuple[csr_matrix, np.ndarray, list[str], np.ndarray]:
    """Morgan bits, descriptors and Murcko scaffolds for unique structures."""

    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")
    generator = GetMorganGenerator(radius=MORGAN_RADIUS, fpSize=MORGAN_BITS)
    functions = [getattr(Descriptors, name) for name in DESCRIPTORS]

    indptr, indices = [0], []
    descriptors = np.zeros((len(smiles), len(DESCRIPTORS)), dtype=np.float64)
    scaffolds: list[str] = []
    valid = np.zeros(len(smiles), dtype=bool)
    for position, text in enumerate(smiles):
        molecule = Chem.MolFromSmiles(text)
        if molecule is None:
            scaffolds.append("")
            indptr.append(len(indices))
            continue
        valid[position] = True
        indices.extend(generator.GetFingerprint(molecule).GetOnBits())
        indptr.append(len(indices))
        descriptors[position] = [function(molecule) for function in functions]
        try:
            core = MurckoScaffold.GetScaffoldForMol(molecule)
            scaffolds.append(Chem.MolToSmiles(core) if core is not None else "")
        except Exception:
            scaffolds.append("")
        if position % 25_000 == 0:
            print(f"  featurised {position}/{len(smiles)}", flush=True)
    bits = csr_matrix(
        (np.ones(len(indices), dtype=np.float32), np.asarray(indices), np.asarray(indptr)),
        shape=(len(smiles), MORGAN_BITS),
    )
    return bits, descriptors, scaffolds, valid


# --------------------------------------------------------------------------- #
# Target-agnostic base, cross-fitted by component
# --------------------------------------------------------------------------- #


def component_fold(component: int, folds: int = 5) -> int:
    return int(sha256(f"{SEED}:oof:{component}".encode()).hexdigest()[:8], 16) % folds


def fit_base(
    bits: csr_matrix, descriptors: np.ndarray, affinity: np.ndarray, folds: np.ndarray
) -> tuple[np.ndarray, dict[str, object]]:
    """Ridge on the complete ligand design vector, cross-fitted by component.

    The base is deliberately target-agnostic: it sees chemistry only.  Every
    target-specific and context-specific effect therefore lands in the residual,
    which is exactly the object D0 interrogates.
    """

    from scipy.sparse import hstack
    from sklearn.linear_model import Ridge

    scale = descriptors.std(axis=0)
    scale[scale < 1e-9] = 1.0
    standardised = (descriptors - descriptors.mean(axis=0)) / scale
    design = hstack([bits, csr_matrix(standardised.astype(np.float32))]).tocsr()

    prediction = np.full(len(affinity), np.nan, dtype=np.float64)
    per_fold = []
    for fold in sorted(set(folds.tolist())):
        held = folds == fold
        train = ~held
        model = Ridge(alpha=10.0, solver="sparse_cg", max_iter=3000, tol=1e-4)
        model.fit(design[train], affinity[train])
        prediction[held] = model.predict(design[held])
        per_fold.append({"fold": int(fold), "train_rows": int(train.sum()), "held_rows": int(held.sum())})
    if not np.isfinite(prediction).all():
        raise RuntimeError("component cross-fit left missing base predictions")
    residual = affinity - prediction
    stats = {
        "model": "ridge(alpha=10) on 1024 Morgan bits + 10 standardised descriptors",
        "cross_fit": "by homology component, 5 folds",
        "target_agnostic": True,
        "folds": per_fold,
        "rmse": float(np.sqrt(np.mean(residual ** 2))),
        "affinity_sd": float(np.std(affinity)),
    }
    return prediction, stats


# --------------------------------------------------------------------------- #
# Build
# --------------------------------------------------------------------------- #


def build(output: Path) -> dict[str, object]:
    started = time.time()
    print("loading activities ...", flush=True)
    activities = load_activities()
    print(f"  {len(activities)} rows, {activities.unit.nunique()} units", flush=True)

    keep = eligible_units(activities)
    activities = activities[activities.unit.isin(keep)].copy()
    print(f"eligible: {len(activities)} rows, {len(keep)} units", flush=True)

    sealed, firewall_stats = sealed_accessions()
    sequences = sequence_table()
    accessions = sorted(set(activities.accession.astype(str)) & set(sequences))
    sealed_present = sorted(sealed & set(sequences))
    labels, sealed_labels, homology_stats = homology_components(accessions, sealed_present, sequences)
    component_of = dict(zip(accessions, labels.tolist()))

    activities = activities[activities.accession.astype(str).isin(component_of)].copy()
    activities["component"] = activities.accession.astype(str).map(component_of)
    before = activities.unit.nunique()
    activities = activities[~activities.component.isin(sealed_labels)].copy()
    firewall_stats.update({
        "units_before_firewall": int(before),
        "units_after_firewall": int(activities.unit.nunique()),
        "components_after_firewall": int(activities.component.nunique()),
    })
    print(f"firewall: {before} -> {activities.unit.nunique()} units, "
          f"{activities.component.nunique()} components", flush=True)

    # Eligibility is re-checked after the firewall so every retained unit is
    # still splittable.
    keep = eligible_units(activities)
    activities = activities[activities.unit.isin(keep)].copy()
    activities["role"] = activities.component.map(assign_role)

    print("featurising unique structures ...", flush=True)
    structures = activities[["connectivity", "SMILES"]].drop_duplicates("connectivity")
    bits, descriptors, scaffolds, valid = featurise(structures.SMILES.tolist())
    position = {key: index for index, key in enumerate(structures.connectivity)}
    activities["structure_row"] = activities.connectivity.map(position)
    activities = activities[valid[activities.structure_row.to_numpy()]].copy()

    # Re-index structures to those actually retained, then re-check eligibility.
    keep = eligible_units(activities)
    activities = activities[activities.unit.isin(keep)].reset_index(drop=True)
    rows = activities.structure_row.to_numpy()
    activities["scaffold"] = [scaffolds[row] for row in rows]

    print("fitting target-agnostic base ...", flush=True)
    folds = activities.component.map(component_fold).to_numpy()
    affinity = activities.affinity.to_numpy(dtype=np.float64)
    base, base_stats = fit_base(bits[rows], descriptors[rows], affinity, folds)
    activities["base"] = base
    activities["oof_fold"] = folds

    output.mkdir(parents=True, exist_ok=True)
    columns = [
        "unit", "accession", "endpoint", "component", "role", "oof_fold",
        "connectivity", "structure_row", "scaffold", "docs", "assays", "Year",
        "affinity", "base",
    ]
    table = activities[columns].rename(columns={"Year": "year"})
    table.to_parquet(output / "rows.parquet", index=False)
    save_npz(output / "morgan.npz", bits.tocsr())
    np.save(output / "descriptors.npy", descriptors)
    structures.reset_index(drop=True).to_parquet(output / "structures.parquet", index=False)

    manifest = {
        "schema": "psep.v1",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "seed": SEED,
        "source": {
            "activities": str(ACTIVITIES),
            "proteins": str(PROTEINS),
            "release": "Papyrus 05.7++ (high quality, no stereochemistry)",
        },
        "counts": {
            "rows": int(len(table)),
            "units": int(table.unit.nunique()),
            "accessions": int(table.accession.nunique()),
            "components": int(table.component.nunique()),
            "structures": int(len(structures)),
            "documents": int(table.docs.nunique()),
            "assays": int(table.assays.nunique()),
            "by_endpoint": table.groupby("endpoint").size().to_dict(),
            "by_role": table.groupby("role").component.nunique().to_dict(),
        },
        "eligibility": {
            "min_train_rows": MIN_TRAIN_ROWS,
            "min_eval_rows": MIN_EVAL_ROWS,
            "rule": "unit needs >=28 rows and >=2 documents holding >=8 rows",
        },
        "homology": homology_stats,
        "firewall": firewall_stats,
        "base": base_stats,
        "roles": {name: fraction for name, fraction in ROLE_FRACTIONS},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the PSEP provenance-separated substrate")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    manifest = build(arguments.output)
    print(json.dumps(manifest["counts"], indent=2))


if __name__ == "__main__":
    main()
