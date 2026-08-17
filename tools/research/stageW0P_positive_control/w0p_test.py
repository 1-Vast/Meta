"""Stage W0-P positive-control test: low-capacity bilinear model, LOOCV."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np
import torch
import torch.nn as nn

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from scripts.qpsmp_data import QPSMPData

RDLogger.DisableLog("rdApp.*")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
SPLIT_VIEW = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views"
PREREG_SHA = "ba0b51ec419b0275a129e69e4cb45db1bccbdd138000893ee7daf881e7bacbf1"

AAS = "ARNDCQEGHILKMFPSTWYV"
BLOSUM62 = [
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
]


def blosum(seq, old, new):
    i, j = AAS.index(old), AAS.index(new)
    return BLOSUM62[i][j]


def load_panel_and_smiles():
    panel = json.loads((HERE / "W0P_PANEL.json").read_text(encoding="utf-8"))
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
    fingerprints = {}
    for pair in panel["pairs"]:
        for row in pair["rows"]:
            if row["ligand_id"] not in fingerprints:
                smiles = data._ligand_smiles.get(row["ligand_id"])
                mol = Chem.MolFromSmiles(smiles) if smiles else None
                if mol is None:
                    fingerprints[row["ligand_id"]] = np.zeros(1024, dtype=np.float32)
                else:
                    fp = rdFingerprintGenerator.GetMorganGenerator(
                        radius=2, fpSize=1024).GetFingerprint(mol)
                    vec = np.zeros(1024, dtype=np.float32)
                    for bit in fp.GetOnBits():
                        vec[bit] = 1.0
                    fingerprints[row["ligand_id"]] = vec
    return panel, sequences, fingerprints, data


def build_features(panel, sequences, fingerprints, data, control, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for pair in panel["pairs"]:
        seq = sequences[pair["target_a"]]
        mutations = list(zip(pair["mutation_positions"],
                             [seq[p] for p in pair["mutation_positions"]],
                             [sequences[pair["target_b"]][p]
                              for p in pair["mutation_positions"]]))
        for row in pair["rows"]:
            ligand = fingerprints[row["ligand_id"]]
            if control == "correct":
                positions = mutations
            elif control == "random_positions":
                positions = []
                for pos, old, new in mutations:
                    while True:
                        q = int(rng.integers(0, len(seq)))
                        if q not in pair["mutation_positions"]:
                            positions.append((q, seq[q], new))
                            break
            elif control == "blosum_positions":
                positions = []
                for pos, old, new in mutations:
                    scores = []
                    for q in range(len(seq)):
                        if q in pair["mutation_positions"]:
                            continue
                        scores.append((abs(blosum(seq, seq[q], old)
                                          - blosum(seq, old, old)), q))
                    scores.sort()
                    positions.append((scores[0][1], seq[scores[0][1]], new))
            elif control == "global_pooled":
                positions = []
            elif control == "random_protein":
                positions = []
            elif control == "ligand_only":
                positions = []
            else:
                raise ValueError(control)
            rows.append({
                "ligand": ligand,
                "positions": positions,
                "len": len(seq),
                "delta_y": row["delta_y"],
                "pooled_a": None if control != "global_pooled" else
                data.protein_for_target(pair["target_a"])[0].numpy(),
                "pooled_b": None if control != "global_pooled" else
                data.protein_for_target(pair["target_b"])[0].numpy(),
                "control": control,
                "pair": (pair["target_a"], pair["target_b"]),
            })
    return rows


class Bilinear(nn.Module):
    def __init__(self, ligand_dim=1024, hidden=32):
        super().__init__()
        self.ligand = nn.Sequential(nn.Linear(ligand_dim, hidden), nn.ReLU())
        self.aa = nn.Linear(40, hidden)
        self.pos = nn.Linear(1, hidden)
        self.combine = nn.Sequential(nn.Linear(hidden, hidden), nn.ReLU())

    def protein_feature(self, row):
        if row["control"] in ("global_pooled",):
            diff = torch.tensor(row["pooled_b"] - row["pooled_a"],
                                dtype=torch.float32)
            # project 640 -> 32 via a small shared head
            return self.global_proj(diff) if hasattr(self, "global_proj") else diff
        if row["control"] == "random_protein":
            return self.random_feature
        if not row["positions"]:
            return torch.zeros(1, 32)
        vectors = []
        for pos, old, new in row["positions"]:
            aa = torch.zeros(40)
            aa[AAS.index(old)] = 1.0
            aa[20 + AAS.index(new)] = 1.0
            p = torch.tensor([pos / max(row["len"], 1)], dtype=torch.float32)
            vectors.append(self.aa(aa) + self.pos(p))
        return self.combine(torch.stack(vectors).mean(0)).unsqueeze(0)


class GlobalBilinear(Bilinear):
    def __init__(self):
        super().__init__()
        self.global_proj = nn.Linear(640, 32)

    def protein_feature(self, row):
        diff = torch.tensor(row["pooled_b"] - row["pooled_a"],
                            dtype=torch.float32)
        return self.global_proj(diff).unsqueeze(0)


class RandomBilinear(Bilinear):
    def __init__(self, seed):
        super().__init__()
        gen = torch.Generator().manual_seed(seed)
        self.random_feature = torch.randn(32, generator=gen) * 0.1

    def protein_feature(self, row):
        return self.random_feature.unsqueeze(0)


def make_model(control, seed):
    if control == "global_pooled":
        return GlobalBilinear()
    if control == "random_protein":
        return RandomBilinear(seed)
    return Bilinear()


def run_control(rows, control, seed):
    torch.manual_seed(seed)
    ligand_t = torch.stack([torch.tensor(r["ligand"]) for r in rows])
    truth_t = torch.tensor([r["delta_y"] for r in rows], dtype=torch.float32)
    pairs = sorted({r["pair"] for r in rows})

    def protein_batch(model):
        if control == "global_pooled":
            diffs = torch.stack([torch.tensor(r["pooled_b"] - r["pooled_a"])
                                 for r in rows])
            return model.global_proj(diffs.float())
        if control == "random_protein":
            return model.random_feature.expand(len(rows), -1)
        if control == "ligand_only":
            return torch.zeros(len(rows), model.ligand[0].out_features)
        out = []
        for r in rows:
            if not r["positions"]:
                out.append(torch.zeros(1, model.ligand[0].out_features))
                continue
            vectors = []
            for pos, old, new in r["positions"]:
                aa = torch.zeros(40)
                aa[AAS.index(old)] = 1.0
                aa[20 + AAS.index(new)] = 1.0
                p = torch.tensor([pos / max(r["len"], 1)], dtype=torch.float32)
                vectors.append(model.aa(aa) + model.pos(p))
            out.append(model.combine(torch.stack(vectors).mean(0)))
        return torch.stack(out)

    predictions = np.zeros(len(rows), dtype=np.float32)
    for test_pair in pairs:
        test_idx = [i for i, r in enumerate(rows) if r["pair"] == test_pair]
        train = [i for i in range(len(rows)) if i not in test_idx]
        model = make_model(control, seed)
        optimizer = torch.optim.AdamW(model.parameters(), lr=0.01,
                                      weight_decay=1e-4)
        for _ in range(300):
            optimizer.zero_grad()
            lig = model.ligand(ligand_t[train])
            prot = protein_batch(model)[train]
            pred = (lig * prot).sum(-1)
            loss = (pred - truth_t[train]).pow(2).mean()
            loss.backward()
            optimizer.step()
        model.eval()
        lig = model.ligand(ligand_t[test_idx])
        prot = protein_batch(model)[test_idx]
        predictions[test_idx] = (lig * prot).sum(-1).detach().numpy()
    truth = truth_t.numpy()
    pred = predictions
    return {
        "pearson": float(np.corrcoef(truth, pred)[0, 1])
        if pred.std() > 1e-12 else 0.0,
        "spearman": float(np.corrcoef(np.argsort(np.argsort(truth)),
                                      np.argsort(np.argsort(pred)))[0, 1])
        if pred.std() > 1e-12 else 0.0,
        "sign_accuracy": float(np.mean(np.sign(truth) == np.sign(pred))),
        "mse": float(np.mean((truth - pred) ** 2)),
        "scheme": "leave-one-pair-out",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    args = parser.parse_args()
    panel, sequences, fingerprints, data = load_panel_and_smiles()
    controls = ["correct", "random_positions", "blosum_positions",
                "global_pooled", "random_protein", "ligand_only"]
    out = {"schema": "MetaSieve.StageW0P.Result.v1",
           "stage": "stageW0P_positive_control",
           "preregistration_sha256": PREREG_SHA,
           "rows": panel["summary"]["total_rows"],
           "seeds": args.seeds, "controls": {}}
    for control in controls:
        metrics = []
        for seed in args.seeds:
            rows = build_features(panel, sequences, fingerprints, data,
                                  control, seed)
            metrics.append(run_control(rows, control, seed))
        out["controls"][control] = {
            "pearson": float(np.mean([m["pearson"] for m in metrics])),
            "spearman": float(np.mean([m["spearman"] for m in metrics])),
            "sign_accuracy": float(np.mean([m["sign_accuracy"]
                                            for m in metrics])),
            "mse": float(np.mean([m["mse"] for m in metrics])),
            "per_seed": metrics,
        }
    (HERE / "W0P_RESULT.json").write_text(
        json.dumps(out, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(out["controls"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
