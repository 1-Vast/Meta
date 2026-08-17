"""Build the KIBA ESM-2 150M protein token bank for Stage W W1."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.cache_structure_proteins import cache_structure_proteins
from tools.research.stageW_soft_mmp.w0_census import KIBA, read_dataset

HERE = Path(__file__).resolve().parent
MODEL_ID = "facebook/esm2_t30_150M_UR50D"
REVISION = "a695f6045e2e32885fa60af20c13cb35398ce30c"
SNAPSHOT = Path.home() / ".cache/huggingface/hub/models--facebook--esm2_t30_150M_UR50D/snapshots/a695f6045e2e32885fa60af20c13cb35398ce30c"
OUTPUT = HERE / "kiba_esm2_t30_slots128"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    _rows, targets, _ligands = read_dataset(KIBA)
    records = HERE / "kiba_targets_records.jsonl"
    with records.open("w", encoding="utf-8") as handle:
        for target_id in sorted(targets):
            handle.write(json.dumps({
                "sequence_sha256": target_id,
                "sequence": targets[target_id],
            }, sort_keys=True) + "\n")
    manifest = cache_structure_proteins(
        records, args.output, model_id=MODEL_ID, revision=REVISION,
        snapshot_path=SNAPSHOT, shard_size=128, token_budget=2048,
        max_batch=8, device="cuda")
    print(json.dumps({
        "output": str(args.output),
        "records": manifest["records"],
        "residue_slots": manifest["residue_slots"],
        "hidden_dim": manifest["hidden_dim"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
