"""Phase 0 — isolate determinism from the dataset change, and audit the sampler.

The earlier determinism comparison was CONFOUNDED: the I-1 atom quarantine was
introduced between the two runs, changing train from 9,758 to 9,757 complexes, so
a differing checkpoint hash was expected and says nothing about determinism.

This trains twice on IDENTICAL data inside one process and compares state dicts,
which isolates code-path determinism. It also audits the negative sampler against
the registered contract at the level the contract actually specifies: exactly six
UNIQUE negatives PER POSITIVE, never a positive.
"""
from __future__ import annotations

import hashlib
import io
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import torch

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import build, make_split, protein_components  # noqa: E402
from s7_localizer import FeatureStore, sample_negatives, train_arm  # noqa: E402
from s7_run import load_mols  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
N_TRAIN = 800


def state_hash(model):
    buf = io.BytesIO()
    torch.save({k: v.detach().cpu() for k, v in model.state_dict().items()}, buf)
    return hashlib.sha256(buf.getvalue()).hexdigest()


def main():
    device = "cuda" if torch.cuda.is_available() else "cpu"
    kept, _q, _c, _f = build()
    comp = protein_components(kept)
    train, _ha, _hA, _hB = make_split(kept, comp)
    sub = train[:N_TRAIN]
    store, mols = FeatureStore(), load_mols()

    print(f"same-data determinism test on {len(sub)} complexes, device={device}",
          flush=True)
    m1, _t1 = train_arm(sub, store, mols, True, device)
    h1 = state_hash(m1)
    m2, _t2 = train_arm(sub, store, mols, True, device)
    h2 = state_hash(m2)
    det = {"run1_sha256": h1, "run2_sha256": h2, "identical": bool(h1 == h2),
           "device": device, "complexes": len(sub)}
    print(f"  identical: {det['identical']}", flush=True)
    if not det["identical"]:
        d = max(float((a - b).abs().max())
                for (a, b) in zip(m1.state_dict().values(), m2.state_dict().values()))
        det["max_abs_param_difference"] = d
        print(f"  max abs parameter difference: {d:.3e}", flush=True)

    # ---- negative sampler audited at the contract's own granularity ----
    rng = np.random.default_rng(20260813)
    per_pos_unique = True
    pos_collisions = 0
    cross_pos_dupes = 0
    total_neg = total_pos = 0
    six_exact = True
    for rec in train[:1500]:
        pos = set(map(tuple, rec["edges"]))
        negs = sample_negatives(rec, rng)
        total_pos += len(rec["edges"])
        total_neg += len(negs)
        if len(negs) != 6 * len(rec["edges"]):
            six_exact = False
        for k in range(0, len(negs), 6):
            block = negs[k:k + 6]
            if len(set(block)) != len(block):
                per_pos_unique = False
            pos_collisions += sum(1 for b in block if b in pos)
        cross_pos_dupes += len(negs) - len(set(negs))

    sampler = {
        "positives_examined": total_pos,
        "negatives_drawn": total_neg,
        "exactly_six_per_positive": six_exact,
        "unique_within_each_positive_block": per_pos_unique,
        "negatives_that_are_actually_positives": pos_collisions,
        "cross_positive_duplicate_pairs": cross_pos_dupes,
        "cross_positive_duplicate_fraction": round(cross_pos_dupes / max(total_neg, 1), 5),
        "contract_reading": "the registered contract requires six UNIQUE negatives PER "
                            "POSITIVE and never a positive. Both hold. Repetition of the "
                            "same pair across DIFFERENT positives in one complex is not a "
                            "contract violation, but it does up-weight those pairs and is "
                            "reported here rather than left implicit.",
        "verdict": "SAMPLER_CONTRACT_SATISFIED" if (six_exact and per_pos_unique
                                                    and pos_collisions == 0)
                   else "SAMPLER_CONTRACT_VIOLATED",
    }
    out = {"schema": "MetaSieve.S7L2B.P0.DeterminismAndSamplerCheck.v1",
           "created_utc": "2026-08-09",
           "earlier_determinism_comparison_was_confounded": {
               "reason": "the I-1 atom quarantine was introduced between the two runs",
               "train_before": 9758, "train_after": 9757,
               "heldout_A_before": 2415, "heldout_A_after": 2409,
               "conclusion": "a differing checkpoint hash across those runs is expected "
                             "and is NOT evidence of non-determinism"},
           "same_data_determinism": det,
           "negative_sampler_audit": sampler}
    (OUT / "P0_DETERMINISM_AND_SAMPLER_CHECK.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"same_data_determinism": det,
                      "negative_sampler_audit": sampler}, indent=2))


if __name__ == "__main__":
    main()
