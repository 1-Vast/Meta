"""Compute frozen ChemBERTa-77M pooled embeddings for all governed ligands.

External representation lane (ligand-side language model, local snapshot
DeepChem/ChemBERTa-77M-MLM). SMILES -> mean-pooled last hidden state.
Output: tools/runtime/chemberta_ligand_pooled/embeddings.npz + manifest.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

MODEL_DIR = "C:/Users/59964/.cache/huggingface/hub/models--DeepChem--ChemBERTa-77M-MLM/snapshots/ed8a5374f2024ec8da53760af91a33fb8f6a15ff"
OUT = ROOT / "tools/runtime/chemberta_ligand_pooled"


def main():
    if OUT.exists():
        print("already exists")
        return 0
    OUT.mkdir(parents=True)
    ligands = []
    with (ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0/ligands.jsonl").open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                ligands.append(json.loads(line))
    print("ligands:", len(ligands))

    from transformers import AutoModel, AutoTokenizer
    tokenizer = AutoTokenizer.from_pretrained(MODEL_DIR, local_files_only=True)
    model = AutoModel.from_pretrained(MODEL_DIR, local_files_only=True)
    model.to("cuda").eval()
    hidden = int(model.config.hidden_size)
    print("hidden:", hidden)
    keys, vectors = [], []
    batch = 64
    with torch.inference_mode():
        for start in range(0, len(ligands), batch):
            chunk = ligands[start:start + batch]
            smiles = [row["smiles"] for row in chunk]
            tokens = tokenizer(smiles, return_tensors="pt", padding=True,
                               truncation=True, max_length=512)
            tokens = {k: v.to("cuda") for k, v in tokens.items()}
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                out = model(**tokens).last_hidden_state
            mask = tokens["attention_mask"].unsqueeze(-1).float()
            pooled = (out * mask).sum(1) / mask.sum(1).clamp_min(1.0)
            for row, vector in zip(chunk, pooled):
                keys.append(row["drug_key"])
                vectors.append(vector.float().cpu().numpy())
            if (start // batch) % 25 == 0:
                print(start, "/", len(ligands), flush=True)
    np.savez_compressed(OUT / "embeddings.npz", keys=np.asarray(keys),
                        pooled=np.stack(vectors))
    (OUT / "manifest.json").write_text(json.dumps({
        "schema": "MetaSieve.StageM.ChemBERTaLigand.v1",
        "model_dir": MODEL_DIR,
        "hidden_dim": hidden,
        "pooling": "attention_masked_mean",
        "ligands": len(keys),
    }, indent=1), encoding="utf-8")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
