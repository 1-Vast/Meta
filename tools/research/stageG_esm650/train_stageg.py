"""Stage G: ESM-650M residue-input trunk screen.

The incumbent similarity_only recipe trained unchanged except that the protein
bank inputs are the local ESM-2 650M embeddings (1280-dim pooled +
128-slot residues) instead of the governed 150M bank (640-dim). Everything
else — seed, budget, losses, partition, leak-free checkpoint selection, GPU
verification — is byte-identical to the Stage D T2 arm, so the only variable
is the external protein representation.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, TrainConfig,
)
from tools.research.stageD_level_panel.train_staged import (
    StageEConfig, train as train_staged_train,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
BANK_650 = ROOT / "tools/runtime/esm2_t33_650M_protein_bank"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--arm", default="G", choices=("G",))
    parser.add_argument("--seed", type=int, default=20260815)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--device",
                        default="cuda" if __import__("torch").cuda.is_available()
                        else "cpu")
    args = parser.parse_args()
    progress_path = args.output / "progress.jsonl"
    if progress_path.exists() and not args.force:
        raise SystemExit(f"{progress_path} exists; pass --force to overwrite")
    data = QPSMPData(CORPUS, BANK_650, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    config = StageEConfig(
        base=TrainConfig(arch="similarity_only", steps=args.steps,
                         seed=args.seed, split_directory=str(SPLIT),
                         device=args.device, amp=False),
        arm="T2")
    args.output.mkdir(parents=True, exist_ok=True)
    train_staged_train(data, config, args.output, progress_path=progress_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
