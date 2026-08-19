"""Q2d-1e pre-flight: span-init + L2 keeps null energy low (600-step
short runs, 4 restarts, correct arm, M1 seed 0, CPU)."""
import sys, json
from pathlib import Path
import numpy as np
import torch

HERE = Path(__file__).resolve().parent
X0C = HERE.parent / "stageX0c_measurement_qualification_20260818"
STAGE_D = HERE.parent / "stageQ2d1d_spanrestricted_interaction_20260818"
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(STAGE_D))
import truth_d as truth
import runner_e as R
import q2

fz = np.load(STAGE_D / 'q2d1d_features.npz', allow_pickle=False)
P_t = fz['P_t'].astype(np.float32)
L_t = fz['L_t'].astype(np.float32)
splits = json.loads(open(STAGE_D / 'Q2D1D_SPLITS.json', encoding='utf-8').read())
for k in ('train_cells', 'pc', 'lc', 'dc'):
    splits[k] = np.asarray(splits[k], dtype=np.int64)
for k in ('cold_row', 'cold_lig', 'train_row', 'train_lig'):
    splits[k] = np.asarray(splits[k], dtype=bool)
device = 'cpu'
t = truth.generate_truth('M1', 0, P_t, L_t, splits)
Vsp, r = truth._span_projection(P_t, splits)
Lt_dev = torch.from_numpy(L_t).float()
for restart in (0, 1, 2, 3):
    model, val, _ = R.train_level(P_t, 'correct', t, 'A', 0, splits, device, restart,
                                  Lt_dev, Vsp=Vsp, max_steps=600)
    A = model.A.weight.detach().numpy().T
    null_frac = 1 - float(np.linalg.norm(Vsp @ (Vsp.T @ A)) ** 2 / np.linalg.norm(A) ** 2)
    print('restart', restart, 'null_frac', round(null_frac, 4),
          'mon', round(val, 3), flush=True)
