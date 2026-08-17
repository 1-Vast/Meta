"""Stage W0-P panel construction from local BindingDB near-identical sequences."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.qpsmp_data import QPSMPData

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
SPLIT_VIEW = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views"
CDHIT = ROOT / "tools/runtime/cdhit/4.8.1/cd-hit.exe"
PREREG_SHA = "ba0b51ec419b0275a129e69e4cb45db1bccbdd138000893ee7daf881e7bacbf1"


def cluster98(sequences):
    ordered = sorted(sequences)
    short = {target: f"s{index}" for index, target in enumerate(ordered)}
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        fasta = tmp / "p.fasta"
        with fasta.open("w", encoding="utf-8") as handle:
            for target in ordered:
                handle.write(f">{short[target]}\n{sequences[target]}\n")
        subprocess.run([str(CDHIT), "-i", str(fasta), "-o", str(tmp / "o"),
                        "-c", "0.98", "-n", "5", "-G", "1", "-g", "1",
                        "-M", "8000", "-T", "1"], check=True,
                       capture_output=True)
        clusters, members = [], []
        for line in (tmp / "o.clstr").read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                if members:
                    clusters.append(members)
                members = []
            elif ">" in line:
                name = line.split(">", 1)[1].split()[0].replace("...", "").strip()
                members.append(ordered[int(name[1:])])
        if members:
            clusters.append(members)
    return [m for m in clusters if len(m) > 1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "W0P_PANEL.json")
    args = parser.parse_args()

    data = QPSMPData(
        CORPUS,
        ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank",
        ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank",
        ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank_compact",
        split_directory=ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1",
        split_view=SPLIT_VIEW)
    proteins = [json.loads(line) for line in
                (CORPUS / "proteins.jsonl").open(encoding="utf-8")]
    sequences = {row["sequence_sha256"]: row["sequence"] for row in proteins}

    cells = defaultdict(dict)
    for cell in data.cells:
        if cell["split"] in ("meta_train", "meta_val"):
            current = cells[(cell["target_id"], cell["ligand_id"])]
            current["pK"] = max(current.get("pK", float("-inf")), float(cell["pK"]))
            current["panel_count"] = current.get("panel_count", 0) + int(cell["panel_count"])

    clusters = cluster98(sequences)
    panel = []
    for members in clusters:
        for i, left in enumerate(members):
            for right in members[i + 1:]:
                sa, sb = sequences[left], sequences[right]
                if len(sa) != len(sb):
                    continue
                mismatches = [j for j, (x, y) in enumerate(zip(sa, sb))
                              if x != y]
                if not (1 <= len(mismatches) <= 5):
                    continue
                shared = set()
                for ligand_id in {k[1] for k in cells if k[0] == left}:
                    if (right, ligand_id) in cells:
                        shared.add(ligand_id)
                if len(shared) < 3:
                    continue
                rows = []
                for ligand_id in sorted(shared):
                    ya = cells[(left, ligand_id)]["pK"]
                    yb = cells[(right, ligand_id)]["pK"]
                    rows.append({"ligand_id": ligand_id, "y_a": ya,
                                 "y_b": yb, "delta_y": yb - ya})
                deltas = np.asarray([r["delta_y"] for r in rows])
                panel.append({
                    "target_a": left, "target_b": right,
                    "len": len(sa),
                    "mutation_positions": mismatches,
                    "n_shared_ligands": len(rows),
                    "rows": rows,
                    "delta_mean": float(deltas.mean()),
                    "delta_sd": float(deltas.std(ddof=1)) if len(deltas) > 1 else 0.0,
                    "sign_consistency": float(
                        max((deltas > 0).mean(), (deltas < 0).mean())),
                })

    total_rows = sum(p["n_shared_ligands"] for p in panel)
    report = {
        "schema": "MetaSieve.StageW0P.Panel.v1",
        "stage": "stageW0P_positive_control",
        "preregistration_sha256": PREREG_SHA,
        "meta_test": data.seal_record(),
        "pairs": panel,
        "summary": {
            "pairs": len(panel),
            "total_rows": total_rows,
            "sufficient": len(panel) >= 3 and total_rows >= 20,
        },
    }
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(report["summary"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
