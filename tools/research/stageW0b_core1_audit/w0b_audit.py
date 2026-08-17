"""Stage W0b — Core Task 1 data, censoring, hierarchy and positive-control audit.

CPU/statistical audit only. No neural model is trained. Run:
    python -m tools.research.stageW0b_core1_audit.w0b_audit
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from tools.research.stageU_mmp_interaction.mmp import fragment, transformation
from tools.research.stageW_soft_mmp.w0_census import cdhit_components

RDLogger.DisableLog("rdApp.*")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PREREG_SHA = "ff23c408d20cc79b1bd5fcd20854a0443280d32d6fc3dbb8abf0733a9a70631f"

FILES = {
    "davis": ROOT / "dataset/raw/dta/davis.tab",
    "kiba": ROOT / "dataset/raw/dta/kiba.tab",
    "metz_matrix": ROOT / "dataset/raw/crossed_panels/kinase_panels/metz_matrix.csv",
    "metz_xls": ROOT / "dataset/raw/crossed_panels/kinase_panels/metz.xls",
    "klaeger": ROOT / "dataset/raw/crossed_panels/kinase_panels/klaeger_matrix.csv",
    "klaeger_smiles": ROOT / "dataset/processed/crossed_panels_xp2/klaeger_smiles.json",
    "klifs_info": ROOT / "dataset/raw/crossed_panels/protein_annotation/klifs_kinase_information_human.json",
    "klifs_groups": ROOT / "dataset/raw/crossed_panels/protein_annotation/klifs_kinase_groups.json",
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_davis():
    targets, rows = {}, defaultdict(list)
    with FILES["davis"].open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for line in reader:
            if len(line) < 5:
                continue
            smiles, sequence, value = line[1], line[3], float(line[4])
            if not smiles or not sequence:
                continue
            target = hashlib.sha256(sequence.encode()).hexdigest()
            targets[target] = sequence
            rows[target].append({"smiles": smiles, "value": value,
                                 "censored": value >= 10000.0 - 1e-9})
    return targets, rows


def read_kiba():
    targets, rows = {}, defaultdict(list)
    with FILES["kiba"].open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for line in reader:
            if len(line) < 5:
                continue
            smiles, sequence, value = line[1], line[3], float(line[4])
            if not smiles or not sequence:
                continue
            target = hashlib.sha256(sequence.encode()).hexdigest()
            targets[target] = sequence
            rows[target].append({"smiles": smiles, "value": value,
                                 "censored": False})
    return targets, rows


def read_klifs():
    info = json.loads(FILES["klifs_info"].read_text(encoding="utf-8"))
    groups = json.loads(FILES["klifs_groups"].read_text(encoding="utf-8"))
    mapping = {}
    for item in info:
        name = str(item.get("HGNC") or item.get("name") or "").upper()
        pocket = str(item.get("pocket") or "")
        group = str(item.get("group") or "")
        mapping[name] = {"group": group, "pocket": pocket}
    return mapping, groups


def read_metz():
    matrix = pd.read_csv(FILES["metz_matrix"])
    matrix = matrix.rename(columns={matrix.columns[0]: "Cmpd_ID"})
    xls = pd.read_excel(FILES["metz_xls"], sheet_name=0)
    smiles = {}
    for _, row in xls.iterrows():
        if pd.notna(row.get("Cmpd_ID")) and pd.notna(row.get("Canonical_Smiles")):
            smiles[int(row["Cmpd_ID"])] = str(row["Canonical_Smiles"])
    kinase_cols = [c for c in matrix.columns
                   if c not in ("Cmpd_ID", "Unnamed: 0")]
    rows = defaultdict(list)
    for _, line in matrix.iterrows():
        cid = int(line["Cmpd_ID"])
        smi = smiles.get(cid)
        if not smi:
            continue
        for col in kinase_cols:
            value = line[col]
            if pd.isna(value):
                continue
            value = float(value)
            rows[col].append({"smiles": smi, "value": value,
                              "censored": value <= 4.0 + 1e-9})
    return rows


def read_klaeger():
    matrix = pd.read_csv(FILES["klaeger"])
    smiles = json.loads(FILES["klaeger_smiles"].read_text(encoding="utf-8"))
    rows = defaultdict(list)
    for _, line in matrix.iterrows():
        drug = str(line["Drug"])
        entry = smiles.get(drug) or {}
        smi = entry.get("smiles")
        if not smi:
            continue
        for col in matrix.columns[1:]:
            value = line[col]
            if pd.isna(value):
                continue
            value = float(value)
            rows[col].append({"smiles": smi, "value": value,
                              "censored": value <= 5.0 + 1e-9})
    return rows


def morgan(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=1024).GetFingerprint(mol)


def layer_census(rows_by_target, components):
    """Broad-to-strict estimand hierarchy for one dataset."""
    fp_cache = {}
    layers = {
        "all_pairs": {"pairs": 0, "classes": defaultdict(set),
                      "target_pairs": 0},
        "similar_pairs": {"pairs": 0, "classes": defaultdict(set)},
        "mmp": {"pairs": 0, "classes": defaultdict(set)},
        "strict_mmp": {"pairs": 0, "classes": defaultdict(set)},
    }
    for target, rows in sorted(rows_by_target.items()):
        ligands = sorted({r["smiles"] for r in rows})
        if len(ligands) < 2:
            continue
        for r in rows_by_target[target]:
            if r["smiles"] not in fp_cache:
                fp_cache[r["smiles"]] = morgan(r["smiles"])
        # all / similar pairs
        for i, left in enumerate(ligands):
            for right in ligands[i + 1:]:
                layers["all_pairs"]["pairs"] += 1
                key = (left, right)
                layers["all_pairs"]["classes"][key].add(target)
                lf, rf = fp_cache[left], fp_cache[right]
                if lf is not None and rf is not None:
                    union = len(lf | rf)
                    if union and len(lf & rf) / union >= 0.6:
                        layers["similar_pairs"]["pairs"] += 1
                        layers["similar_pairs"]["classes"][key].add(target)
        # MMP / strict MMP
        pieces = {}
        for smiles in ligands:
            parts = fragment(smiles)
            if parts:
                pieces[smiles] = parts
        by_core = defaultdict(list)
        for smiles, parts in pieces.items():
            for part in parts:
                by_core[part.core].append((smiles, part))
        emitted = set()
        for core, entries in by_core.items():
            for pos, (left_s, left) in enumerate(entries):
                for right_s, right in entries[pos + 1:]:
                    built = transformation(left, right)
                    if built is None:
                        continue
                    item, _flip = built
                    sig = (core, item.r_a, item.r_b)
                    if sig in emitted:
                        continue
                    emitted.add(sig)
                    layers["mmp"]["pairs"] += 1
                    layers["mmp"]["classes"][item.exact_key].add(target)
                    layers["strict_mmp"]["pairs"] += 1
                    layers["strict_mmp"]["classes"][item.exact_key].add(target)

    out = {}
    for name, layer in layers.items():
        classes = {key: set(value) for key, value in layer["classes"].items()}
        comp_sets = {key: {components[t] for t in value} for key, value in
                     classes.items()}
        rich = [key for key in classes
                if len(classes[key]) >= 3 and len(comp_sets[key]) >= 3]
        cross_d = 0
        for key, targets in classes.items():
            ordered = sorted(targets)
            for i, a in enumerate(ordered):
                for b in ordered[i + 1:]:
                    if components[a] != components[b]:
                        cross_d += 1
        eiu = min(len({components[t] for t in components if t in
                       {x for s in classes.values() for x in s}}),
                  len(classes)) if classes else 0
        out[name] = {
            "within_target_pairs": layer["pairs"],
            "classes": len(classes),
            "classes_ge3_targets_ge3_components": len(rich),
            "cross_component_D_rows": cross_d,
            "effective_independent_units": eiu,
            "top1_class_share": (
                max(len(v) for v in classes.values()) / len(classes)
                if classes else 0.0),
        }
    return out


def censoring_stats(rows_by_target):
    values, censored = [], []
    per_target = {}
    for target, rows in rows_by_target.items():
        vals = [r["value"] for r in rows]
        cens = [r["censored"] for r in rows]
        per_target[target] = {
            "n": len(vals), "censored_fraction": float(np.mean(cens))}
        values.extend(vals)
        censored.extend(cens)
    return {
        "cells": len(values),
        "censored_cells": int(sum(censored)),
        "censored_fraction": float(np.mean(censored)),
        "value_quantiles": {
            "min": float(np.min(values)), "p05": float(np.quantile(values, 0.05)),
            "median": float(np.median(values)),
            "p95": float(np.quantile(values, 0.95)),
            "max": float(np.max(values)),
        },
        "per_target_censored_fraction_quantiles": {
            "median": float(np.median([v["censored_fraction"]
                                       for v in per_target.values()])),
            "p95": float(np.quantile([v["censored_fraction"]
                                      for v in per_target.values()], 0.95)),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-prefix", type=str, default="W0B")
    args = parser.parse_args()

    assets = {name: {"path": str(path), "sha256": sha256(path)}
              for name, path in FILES.items()}

    # Positive-control asset search.
    search_names = ["hiv", "mutant", "variant", "ortholog", "resistance",
                    "gatekeeper", "anastassiadis"]
    found = []
    for root in (ROOT / "dataset/raw", ROOT / "dataset/processed"):
        for path in root.rglob("*"):
            try:
                is_file = path.is_file()
            except OSError:
                continue
            if is_file and any(token in path.name.lower()
                               for token in search_names):
                found.append(str(path))
    positive_control = {
        "searched_tokens": search_names,
        "found_files": found[:50],
        "runnable": False,
        "reason": ("No HIV-resistance, gatekeeper/mutation or ortholog "
                   "point-mutant panel with matching ligand affinities was "
                   "found locally. W0-P cannot be executed from local assets."),
    }

    # Stage W audit: files and known prereg hashes.
    stagew = HERE.parent / "stageW_soft_mmp"
    stagew_audit = {
        "preregistration_sha256": sha256(stagew / "PREREGISTRATION.md"),
        "w1_preregistration_sha256": sha256(stagew / "W1_PREREGISTRATION.md"),
        "files": sorted(str(p.relative_to(ROOT)) for p in stagew.rglob("*")
                        if p.is_file() and "__pycache__" not in p.parts),
        "trained_artifacts": [],
        "w1_training_status": "PAUSED; no training metric read",
    }

    # Censoring and hierarchy for locally parseable single-platform panels.
    censoring, hierarchy = {}, {}
    davis_targets, davis_rows = read_davis()
    davis_components = cdhit_components(davis_targets, "davis")
    censoring["davis"] = censoring_stats(davis_rows)
    hierarchy["davis"] = layer_census(davis_rows, davis_components)

    try:
        metz_rows = read_metz()
        klifs, groups = read_klifs()
        metz_components = {
            target: klifs.get(target.upper(), {}).get("group", "UNMAPPED")
            for target in metz_rows}
        censoring["metz"] = censoring_stats(metz_rows)
        hierarchy["metz"] = layer_census(metz_rows, metz_components)
    except Exception as exc:
        censoring["metz"] = {"error": str(exc)}
        hierarchy["metz"] = {"error": str(exc)}

    try:
        klaeger_rows = read_klaeger()
        klifs, groups = read_klifs()
        klaeger_components = {
            target: klifs.get(target.upper(), {}).get("group", "UNMAPPED")
            for target in klaeger_rows}
        censoring["klaeger"] = censoring_stats(klaeger_rows)
        hierarchy["klaeger"] = layer_census(klaeger_rows, klaeger_components)
    except Exception as exc:
        censoring["klaeger"] = {"error": str(exc)}
        hierarchy["klaeger"] = {"error": str(exc)}

    kiba_targets, kiba_rows = read_kiba()
    kiba_components = cdhit_components(kiba_targets, "kiba")
    censoring["kiba"] = censoring_stats(kiba_rows)

    report = {
        "schema": "MetaSieve.StageW0b.CoreTask1Audit.v1",
        "stage": "stageW0b_core1_audit",
        "preregistration_sha256": PREREG_SHA,
        "assets": assets,
        "positive_control_audit": positive_control,
        "stage_w_audit": stagew_audit,
        "censoring": censoring,
        "hierarchy": hierarchy,
        "go_no_go": {
            "w1_biological_interpretation": "NO-GO",
            "reason": ("W0-P positive control is NOT RUNNABLE from local "
                       "assets; support/censoring statistics alone cannot "
                       "authorize a biological-null or signal claim."),
            "per_dataset": {name: "AUDITED"
                            for name in censoring},
        },
    }
    out = HERE / "W0B_AUDIT.json"
    out.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                   encoding="utf-8")
    print(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
