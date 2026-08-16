"""Gate D0-R: natural-tail roster under the per-recipient document-ordered estimand.

Metadata-only. Reads target/compound/document/assay identity, document years and
counts. Never reads ``pKi``, ``standard_value`` or ``pchembl_value``.

Design (registered in task.md, section "2026-08-01 A2S-SDO Preregistration And
D0 Protocol Revision"):

  For a scarce recipient r, order its documents by (document_year, document_uid).
  Support S_r is drawn label-blind from documents with year <= tau_r; query Q_r is
  the parent compounds appearing only in documents with year > tau_r and in no
  support document. tau_r is the earliest split year admitting k=5 support parents
  and >= QUERY_FLOOR closed query parents.

Outputs an immutable versioned package. Does not overwrite the sealed v4 corpus.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset" / "formal_training" / "chembl37_pki_formal.v4"
DEFAULT_OUT = ROOT / "dataset" / "formal_training" / "a2s_d0r_roster.v3"

# ---- frozen protocol constants (preregistered; not tuned on any score) -------
SCARCE_MIN, SCARCE_MAX = 5, 50      # recipient total unique parent compounds
SOURCE_MIN = 100                     # source unique parent compounds after closure
SUPPORT_K = (1, 3, 5)                # nested support budgets
QUERY_FLOOR = 10                     # closed query parents per recipient
N_DRAWS = 5                          # label-blind nested support draws
HOMOLOGY_IDENTITY = 0.40             # sequence-identity clustering threshold
RECIPIENT_FLOOR = 50                 # D0-R PASS floor on recipients
COMPONENT_FLOOR = 25                 # D0-R PASS floor on independent components
MAX_COMPONENT_SHARE = 0.15           # D0-R PASS ceiling on largest component
SIGMA_DELTA_GRID = (0.2, 0.3, 0.4, 0.5)
SEED = 1729

META_COLUMNS = [
    "target_uid", "accession", "protein_class_id", "compound_parent_uid",
    "connectivity_inchikey", "document_uid", "document_year",
    "assay_context_uid", "assay_id", "document_src_id", "measurement_uid",
]
FORBIDDEN = {"pKi", "pKd", "standard_value", "pchembl_value", "pchembl_delta"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


# --------------------------------------------------------------------- scaffolds
def murcko_scaffolds(smiles_by_parent: dict[str, str]) -> dict[str, str]:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")
    out: dict[str, str] = {}
    for uid, smi in smiles_by_parent.items():
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                out[uid] = "INVALID"
                continue
            core = MurckoScaffold.GetScaffoldForMol(mol)
            out[uid] = Chem.MolToSmiles(core) if core is not None else "ACYCLIC"
        except Exception:
            out[uid] = "INVALID"
    return out


# --------------------------------------------------------------------- homology
def homology_components(sequences: dict[str, str], identity: float) -> dict[str, int]:
    """Union-find clustering at a sequence-identity threshold.

    Stage 1 prefilters candidate pairs by 4-mer Jaccard; stage 2 confirms with a
    parasail Smith-Waterman identity over the shorter sequence.
    """
    import parasail

    uids = sorted(sequences)
    kmers = {u: {sequences[u][i:i + 4] for i in range(len(sequences[u]) - 3)} for u in uids}

    parent = {u: u for u in uids}

    def find(a: str) -> str:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    matrix = parasail.blosum62
    for i, a in enumerate(uids):
        ka = kmers[a]
        for b in uids[i + 1:]:
            kb = kmers[b]
            inter = len(ka & kb)
            if inter == 0:
                continue
            jac = inter / len(ka | kb)
            if jac < 0.03:                       # loose prefilter, recall-oriented
                continue
            sa, sb = sequences[a], sequences[b]
            res = parasail.sw_stats_striped_16(sa, sb, 11, 1, matrix)
            denom = min(len(sa), len(sb))
            if denom and res.matches / denom >= identity:
                union(a, b)
    roots = {}
    comp: dict[str, int] = {}
    for u in uids:
        r = find(u)
        comp[u] = roots.setdefault(r, len(roots))
    return comp


# --------------------------------------------------------------------- roster
def pick_split_year(frame: pd.DataFrame, k_max: int, query_floor: int) -> dict[str, Any] | None:
    """Earliest split year admitting k_max support parents and query_floor query parents."""
    years = sorted(frame.document_year.unique())
    for split in years[:-1]:
        sup = frame[frame.document_year <= split]
        qry = frame[frame.document_year > split]
        sup_par = set(sup.compound_parent_uid)
        if len(sup_par) < k_max:
            continue
        sup_doc = set(sup.document_uid)
        qry = qry[~qry.document_uid.isin(sup_doc)]           # document-disjoint
        qry_par = set(qry.compound_parent_uid) - sup_par      # parent-disjoint
        if len(qry_par) < query_floor:
            continue
        return {
            "tau": int(split),
            "support_pool": sorted(sup_par),
            "query": sorted(qry_par),
            "support_docs": sorted(sup_doc),
            "query_docs": sorted(qry[qry.compound_parent_uid.isin(qry_par)].document_uid.unique()),
        }
    return None


def nested_draws(pool: list[str], scaffolds: dict[str, str], rng: np.random.Generator,
                 n_draws: int, budgets: tuple[int, ...]) -> list[dict[int, list[str]]]:
    """Label-blind nested support draws, scaffold-diverse, deterministic per seed."""
    draws: list[dict[int, list[str]]] = []
    seen: set[tuple[str, ...]] = set()
    kmax = max(budgets)
    for _ in range(n_draws * 20):
        if len(draws) >= n_draws:
            break
        order = rng.permutation(len(pool))
        chosen: list[str] = []
        used_scaffolds: set[str] = set()
        for idx in order:                       # pass 1: distinct scaffolds
            uid = pool[idx]
            scaf = scaffolds.get(uid, "INVALID")
            if scaf in used_scaffolds:
                continue
            used_scaffolds.add(scaf)
            chosen.append(uid)
            if len(chosen) == kmax:
                break
        for idx in order:                       # pass 2: top up if scaffold-poor
            if len(chosen) == kmax:
                break
            uid = pool[idx]
            if uid not in chosen:
                chosen.append(uid)
        if len(chosen) < kmax:
            continue
        # uniqueness is on the SET, not the order: a support pool of exactly kmax
        # admits only one k=kmax support set no matter how it is permuted
        key = frozenset(chosen)
        if key in seen:
            continue
        seen.add(key)
        draws.append({k: chosen[:k] for k in budgets})       # nested by construction
    return draws


def build(output: Path) -> dict[str, Any]:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite sealed package: {output}")
    main_context = CORPUS / "canonical" / "pki_measurements_context_main.parquet"
    rows = pd.read_parquet(main_context, columns=META_COLUMNS)
    assert not (set(rows.columns) & FORBIDDEN), "outcome firewall violated"
    rows = rows[rows.document_year.notna()].copy()
    rows["document_year"] = rows.document_year.astype(int)

    per_target = rows.groupby("target_uid").compound_parent_uid.nunique()
    cutflow: list[dict[str, Any]] = [
        {"stage": "corpus_targets_with_document_year", "targets": int(per_target.size)}
    ]

    # ---- scaffolds ---------------------------------------------------------
    compounds = pd.read_parquet(
        CORPUS / "components" / "compounds.parquet",
        columns=["compound_parent_uid", "parent_canonical_smiles"])
    scaffolds = murcko_scaffolds(dict(zip(compounds.compound_parent_uid,
                                          compounds.parent_canonical_smiles)))

    # ---- recipient candidates ---------------------------------------------
    scarce = set(per_target[per_target.between(SCARCE_MIN, SCARCE_MAX)].index)
    cutflow.append({"stage": f"scarce_band_{SCARCE_MIN}_{SCARCE_MAX}_parents",
                    "targets": len(scarce)})

    recipients: list[dict[str, Any]] = []
    episodes: list[dict[str, Any]] = []
    query_rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(SEED)

    for target, frame in rows[rows.target_uid.isin(scarce)].groupby("target_uid", sort=True):
        split = pick_split_year(frame, max(SUPPORT_K), QUERY_FLOOR)
        if split is None:
            continue
        draws = nested_draws(split["support_pool"], scaffolds, rng, N_DRAWS, SUPPORT_K)
        if len(draws) < N_DRAWS:
            continue
        support_observations = (
            frame[
                (frame.document_year <= split["tau"])
                & frame.compound_parent_uid.isin(split["support_pool"])
            ]
            .sort_values(
                ["document_year", "document_uid", "assay_context_uid", "measurement_uid"],
                kind="stable",
            )
            .drop_duplicates("compound_parent_uid", keep="first")
            .set_index("compound_parent_uid")
        )
        query_observations = (
            frame[
                (frame.document_year > split["tau"])
                & ~frame.document_uid.isin(split["support_docs"])
                & frame.compound_parent_uid.isin(split["query"])
            ]
            .sort_values(
                ["document_year", "document_uid", "assay_context_uid", "measurement_uid"],
                kind="stable",
            )
            .drop_duplicates("compound_parent_uid", keep="first")
            .set_index("compound_parent_uid")
        )
        if set(split["support_pool"]) - set(support_observations.index):
            raise RuntimeError("support observation identity is incomplete")
        if set(split["query"]) - set(query_observations.index):
            raise RuntimeError("query observation identity is incomplete")
        recipients.append({
            "target_uid": target,
            "accession": frame.accession.iloc[0],
            "protein_class_id": int(frame.protein_class_id.iloc[0]),
            "tau": split["tau"],
            "n_support_pool": len(split["support_pool"]),
            "n_query": len(split["query"]),
            "n_support_docs": len(split["support_docs"]),
            "n_query_docs": len(split["query_docs"]),
            "n_total_parents": int(per_target[target]),
        })
        for draw_id, draw in enumerate(draws):
            for k, members in draw.items():
                for uid in members:
                    observation = support_observations.loc[uid]
                    episodes.append({"target_uid": target, "draw_id": draw_id,
                                     "k": k, "compound_parent_uid": uid,
                                     "scaffold": scaffolds.get(uid, "INVALID"),
                                     "measurement_uid": str(observation.measurement_uid),
                                     "document_uid": str(observation.document_uid),
                                     "document_year": int(observation.document_year),
                                     "assay_context_uid": str(observation.assay_context_uid)})
        for uid in split["query"]:
            observation = query_observations.loc[uid]
            query_rows.append({"target_uid": target, "compound_parent_uid": uid,
                               "scaffold": scaffolds.get(uid, "INVALID"),
                               "measurement_uid": str(observation.measurement_uid),
                               "document_uid": str(observation.document_uid),
                               "document_year": int(observation.document_year),
                               "assay_context_uid": str(observation.assay_context_uid)})

    cutflow.append({"stage": "document_ordered_split_and_5_distinct_draws",
                    "targets": len(recipients)})
    recipient_frame = pd.DataFrame(recipients)
    if recipient_frame.empty:
        raise SystemExit("D0-R: zero recipients; DATA STOP")

    recipient_ids = set(recipient_frame.target_uid)
    recipient_rows = rows[rows.target_uid.isin(recipient_ids)]
    recipient_docs = set(recipient_rows.document_uid)
    recipient_parents = set(recipient_rows.compound_parent_uid)
    recipient_accessions = set(recipient_frame.accession)

    # ---- source pool with document/parent closure --------------------------
    source_candidates = set(per_target[per_target >= SOURCE_MIN].index) - recipient_ids
    src = rows[rows.target_uid.isin(source_candidates)]
    src = src[~src.accession.isin(recipient_accessions)]
    cutflow.append({"stage": f"source_candidates_{SOURCE_MIN}_parents_accession_disjoint",
                    "targets": int(src.target_uid.nunique())})
    closed = src[~src.document_uid.isin(recipient_docs)
                 & ~src.compound_parent_uid.isin(recipient_parents)]
    closed_counts = closed.groupby("target_uid").compound_parent_uid.nunique()
    source_ids = set(closed_counts[closed_counts >= SOURCE_MIN].index)
    closed = closed[closed.target_uid.isin(source_ids)]
    cutflow.append({"stage": "source_after_document_parent_closure",
                    "targets": len(source_ids)})

    # ---- homology components over recipients + sources ---------------------
    targets_tbl = pd.read_parquet(CORPUS / "components" / "targets.parquet",
                                  columns=["target_uid", "sequence", "accession"])
    keep = recipient_ids | source_ids
    seqs = {r.target_uid: r.sequence for r in
            targets_tbl[targets_tbl.target_uid.isin(keep)].itertuples()}
    comp = homology_components(seqs, HOMOLOGY_IDENTITY)
    recipient_frame["component_id"] = recipient_frame.target_uid.map(comp)

    # a recipient sharing a homology component with any source is homology-warm
    source_components = {comp[t] for t in source_ids if t in comp}
    recipient_frame["homology_warm"] = recipient_frame.component_id.isin(source_components)

    comp_sizes = recipient_frame.component_id.value_counts()
    n_components = int(comp_sizes.size)
    largest_share = float(comp_sizes.iloc[0] / len(recipient_frame))

    # ---- query strata ------------------------------------------------------
    query_frame = pd.DataFrame(query_rows)
    source_conn = set(closed.connectivity_inchikey)
    conn_by_parent = dict(zip(rows.compound_parent_uid, rows.connectivity_inchikey))
    query_frame["connectivity_inchikey"] = query_frame.compound_parent_uid.map(conn_by_parent)
    query_frame["source_seen"] = query_frame.connectivity_inchikey.isin(source_conn)
    source_scaffolds = set(pd.Series(sorted(set(closed.compound_parent_uid)))
                           .map(scaffolds).dropna())
    query_frame["scaffold_cold"] = ~query_frame.scaffold.isin(source_scaffolds)

    episode_frame = pd.DataFrame(episodes)
    support_scaffold_by_target = episode_frame[episode_frame.k == max(SUPPORT_K)] \
        .groupby("target_uid").scaffold.apply(set).to_dict()
    query_frame["support_scaffold_overlap"] = [
        row.scaffold in support_scaffold_by_target.get(row.target_uid, set())
        for row in query_frame.itertuples()]

    # ---- power -------------------------------------------------------------
    n_units = n_components
    mde80 = {str(s): round(2.802 * s / np.sqrt(n_units), 4) for s in SIGMA_DELTA_GRID}

    # ---- hard overlap audit ------------------------------------------------
    overlap = {
        "source_recipient_target_uid": len(source_ids & recipient_ids),
        "source_recipient_accession": len(set(closed.accession) & recipient_accessions),
        "source_recipient_document_uid": len(set(closed.document_uid) & recipient_docs),
        "source_recipient_parent": len(set(closed.compound_parent_uid) & recipient_parents),
        "source_recipient_assay_id": len(set(closed.assay_id) & set(recipient_rows.assay_id)),
        "source_recipient_homology_component": len(
            source_components & set(recipient_frame.component_id)),
        "support_query_parent": 0,      # disjoint by construction, verified below
        "support_query_document": 0,
        "support_query_measurement": 0,
    }
    sup_pairs = set(zip(episode_frame.target_uid, episode_frame.compound_parent_uid))
    qry_pairs = set(zip(query_frame.target_uid, query_frame.compound_parent_uid))
    overlap["support_query_parent"] = len(sup_pairs & qry_pairs)
    sup_docs = set(zip(episode_frame.target_uid, episode_frame.document_uid))
    qry_docs = set(zip(query_frame.target_uid, query_frame.document_uid))
    overlap["support_query_document"] = len(sup_docs & qry_docs)
    sup_measurements = set(episode_frame.measurement_uid)
    qry_measurements = set(query_frame.measurement_uid)
    overlap["support_query_measurement"] = len(sup_measurements & qry_measurements)

    passed = (len(recipient_frame) >= RECIPIENT_FLOOR
              and n_components >= COMPONENT_FLOOR
              and largest_share <= MAX_COMPONENT_SHARE
              and overlap["source_recipient_target_uid"] == 0
              and overlap["source_recipient_accession"] == 0
              and overlap["source_recipient_document_uid"] == 0
              and overlap["source_recipient_parent"] == 0
              and overlap["support_query_parent"] == 0
              and overlap["support_query_document"] == 0
              and overlap["support_query_measurement"] == 0
              and int(query_frame.groupby("target_uid").size().min()) >= QUERY_FLOOR)

    output.mkdir(parents=True, exist_ok=True)
    recipient_frame.to_parquet(output / "recipients.parquet", index=False)
    episode_frame.to_parquet(output / "support_draws.parquet", index=False)
    query_frame.to_parquet(output / "query.parquet", index=False)
    pd.DataFrame({"target_uid": sorted(source_ids)}).assign(
        component_id=lambda d: d.target_uid.map(comp)).to_parquet(
        output / "sources.parquet", index=False)
    closed[["target_uid", "compound_parent_uid", "measurement_uid", "document_uid",
            "document_year", "assay_context_uid", "connectivity_inchikey"]].to_parquet(
                output / "source_rows.parquet", index=False)

    report = {
        "schema": "a2s-d0r-roster-v2",
        "status": "PASS" if passed else "DATA_STOP",
        "estimand": "per-recipient document-ordered support->query; "
                    "support from documents <= tau_r, query from strictly later documents",
        "labels_used_for_roster_selection": False,
        "label_table": "canonical/pki_measurements_context_main.parquet",
        "frozen_observation_key": "measurement_uid",
        "protocol": {
            "scarce_band_parents": [SCARCE_MIN, SCARCE_MAX],
            "source_min_parents": SOURCE_MIN,
            "support_budgets": list(SUPPORT_K),
            "query_floor": QUERY_FLOOR,
            "support_draws": N_DRAWS,
            "homology_identity": HOMOLOGY_IDENTITY,
            "seed": SEED,
        },
        "cutflow": cutflow,
        "counts": {
            "recipients": int(len(recipient_frame)),
            "recipient_accessions": int(recipient_frame.accession.nunique()),
            "independent_components": n_components,
            "largest_component_share": round(largest_share, 4),
            "homology_warm_recipients": int(recipient_frame.homology_warm.sum()),
            "sources": len(source_ids),
            "source_rows": int(len(closed)),
            "source_parents": int(closed.compound_parent_uid.nunique()),
            "query_rows": int(len(query_frame)),
            "query_source_seen_fraction": round(float(query_frame.source_seen.mean()), 4),
            "query_scaffold_cold_fraction": round(float(query_frame.scaffold_cold.mean()), 4),
            "query_support_scaffold_overlap_fraction":
                round(float(query_frame.support_scaffold_overlap.mean()), 4),
        },
        "query_depth": {
            "min": int(query_frame.groupby("target_uid").size().min()),
            "q25": float(query_frame.groupby("target_uid").size().quantile(0.25)),
            "median": float(query_frame.groupby("target_uid").size().median()),
            "max": int(query_frame.groupby("target_uid").size().max()),
        },
        "tau": {
            "min": int(recipient_frame.tau.min()),
            "median": float(recipient_frame.tau.median()),
            "max": int(recipient_frame.tau.max()),
        },
        "overlap": overlap,
        "power": {"statistical_unit": "independent homology component",
                  "n_units": n_units, "mde80_by_sigma_delta": mde80},
        "pass_rule": {
            "recipient_floor": RECIPIENT_FLOOR,
            "component_floor": COMPONENT_FLOOR,
            "max_component_share": MAX_COMPONENT_SHARE,
        },
    }
    (output / "d0r_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    files = {p.name: {"sha256": sha256_file(p), "bytes": p.stat().st_size}
             for p in sorted(output.iterdir()) if p.is_file()}
    manifest = {"schema": "a2s-d0r-manifest-v2", "status": report["status"],
                "source_package": str(CORPUS.relative_to(ROOT)).replace("\\", "/"),
                "label_table": "canonical/pki_measurements_context_main.parquet",
                "labels_read_during_build": False,
                "files": files,
                "content_sha256": sha256_text(json.dumps(files, sort_keys=True))}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Gate D0-R natural-tail roster (metadata only)")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    if args.out.exists() and (args.out / "manifest.json").exists():
        raise SystemExit(f"refusing to overwrite sealed package: {args.out}")
    report = build(args.out)
    print(json.dumps({k: report[k] for k in
                      ("status", "counts", "query_depth", "tau", "overlap", "power")}, indent=2))
    print("\ncutflow:")
    for stage in report["cutflow"]:
        print(f"  {stage['stage']:<55} {stage['targets']:>6}")


if __name__ == "__main__":
    main()
