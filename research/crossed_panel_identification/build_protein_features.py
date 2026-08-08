"""Build protein representations for the XP1-B zero-shot arms.

Includes the representation the MetaSieve pipeline itself uses (ESM-2 t30 150M,
cached locally as `facebook/esm2_t30_150M_UR50D`), so that the zero-shot verdict
applies to the production encoder and not only to hand-made pocket features.
"""
from __future__ import annotations

import json
import os
import ssl
import sys
import urllib.request

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panels import load_klifs, load_metz, map_kinases  # noqa: E402

CACHE = r"D:\MetaSieve\dataset\processed\crossed_panels"
os.makedirs(CACHE, exist_ok=True)
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE


def fetch_sequences(uniprot_ids):
    path = os.path.join(CACHE, "uniprot_sequences.json")
    seqs = json.load(open(path)) if os.path.exists(path) else {}
    todo = [u for u in uniprot_ids if u and u not in seqs]
    for i, u in enumerate(todo):
        url = f"https://rest.uniprot.org/uniprotkb/{u}.fasta"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            txt = urllib.request.urlopen(req, timeout=60, context=CTX).read().decode()
            seqs[u] = "".join(l.strip() for l in txt.splitlines() if not l.startswith(">"))
        except Exception as e:
            print(f"  !! {u}: {type(e).__name__}")
        if (i + 1) % 20 == 0:
            print(f"  fetched {i+1}/{len(todo)}")
    json.dump(seqs, open(path, "w"))
    return seqs


def esm_embed(records, tag, max_len=1024, batch=4):
    """Mean-pooled ESM-2 t30 150M embedding for a list of (key, sequence)."""
    import torch
    from transformers import AutoTokenizer, EsmModel

    path = os.path.join(CACHE, f"esm2_t30_{tag}.npz")
    if os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: v for k, v in zip(d["keys"], d["emb"])}
    name = "facebook/esm2_t30_150M_UR50D"
    tok = AutoTokenizer.from_pretrained(name)
    mdl = EsmModel.from_pretrained(name).eval()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    mdl = mdl.to(dev)
    keys, embs = [], []
    with torch.no_grad():
        for i in range(0, len(records), batch):
            chunk = records[i:i + batch]
            enc = tok([s[:max_len] for _, s in chunk], return_tensors="pt",
                      padding=True, truncation=True, max_length=max_len).to(dev)
            out = mdl(**enc).last_hidden_state
            mask = enc["attention_mask"].unsqueeze(-1).float()
            pooled = (out * mask).sum(1) / mask.sum(1)
            for (k, _), e in zip(chunk, pooled.cpu().numpy()):
                keys.append(k)
                embs.append(e)
            print(f"  esm[{tag}] {min(i+batch, len(records))}/{len(records)}", end="\r")
    print()
    embs = np.stack(embs)
    np.savez(path, keys=np.array(keys), emb=embs)
    return {k: v for k, v in zip(keys, embs)}


def build(density=0.60):
    Y, M, cid, kin = load_metz(density)
    hit, miss = map_kinases(kin, load_klifs())
    assert not miss, miss
    uni = [hit[k]["uniprot"] for k in kin]
    pockets = [hit[k]["pocket"] for k in kin]
    print(f"{len(kin)} kinases, {len(set(uni))} uniprot ids")
    seqs = fetch_sequences(uni)
    have = [u for u in uni if u in seqs]
    print(f"sequences: {len(have)}/{len(uni)}")

    full = esm_embed([(u, seqs[u]) for u in dict.fromkeys(have)], "kinase_full")
    pock = esm_embed([(k, p) for k, p in zip(kin, pockets)], "kinase_pocket85")

    Efull = np.stack([full[u] if u in full else np.zeros(640) for u in uni])
    Epock = np.stack([pock[k] for k in kin])
    out = os.path.join(CACHE, f"metz{int(density*100)}_protein_features.npz")
    np.savez(out, kinases=np.array(kin), uniprot=np.array(uni),
             esm_full=Efull, esm_pocket=Epock)
    print("wrote", out, Efull.shape, Epock.shape)


if __name__ == "__main__":
    build(0.60)
    build(0.70)
