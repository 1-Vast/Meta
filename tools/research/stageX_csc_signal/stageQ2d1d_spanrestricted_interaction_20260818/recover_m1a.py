"""Q2d-1d M1-A-only recovery after the original ladder crashed at M2
truth generation (NameError: PCA_VT not defined in truth_d.py).

Runs the EXACT frozen code path of runner_d.main restricted to M1 level A
(seeds 0,1,2, all 8 arms, correct 8 restarts, negatives 1) and writes
Q2D1D_LADDER.json (schema LADDER.v1) carrying only M1:A plus provenance
flags. The recovered dz values must match the original runner_d.log
printouts (cross-checked by recover_m1a_check.py). Verdict evidence only;
no truth/split/gate change.
"""
import json
import sys
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
X0C = HERE.parent / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(HERE))
import runner_d as R  # noqa: E402
import truth_d as truth  # noqa: E402
import q2  # noqa: E402
import x0_i1  # noqa: E402
from x0_common import stable_rng  # noqa: E402


def main() -> int:
    device = "cuda" if torch.cuda.is_available() else "cpu"
    rows, compounds, prot_feats, lig_feats, scaffolds, _meta = x0_i1.load_features()
    fz = np.load(HERE / "q2d1d_features.npz", allow_pickle=False)
    P_t = fz["P_t"].astype(np.float32)
    L_t = fz["L_t"].astype(np.float32)
    splits = json.loads(open(HERE / "Q2D1D_SPLITS.json", encoding="utf-8").read())
    for k in ("train_cells", "pc", "lc", "dc"):
        splits[k] = np.asarray(splits[k], dtype=np.int64)
    for k in ("cold_row", "cold_lig", "train_row", "train_lig"):
        splits[k] = np.asarray(splits[k], dtype=bool)
    n_rows = P_t.shape[0]
    Lt_dev = torch.from_numpy(L_t).float().to(device)
    rng_arm = stable_rng("stageQ2d1d", "arms")
    shuf = rng_arm.permutation(n_rows)
    fams = np.asarray([q2.family_of_parent(x0_i1._parent_of(r)) for r in rows])
    fam_perm = np.arange(n_rows)
    for f in set(fams.tolist()):
        idx = np.where(fams == f)[0]
        fam_perm[idx] = idx[rng_arm.permutation(len(idx))]
    rand_p = rng_arm.normal(0, 1, size=(n_rows, P_t.shape[1])).astype(np.float32)
    results = {"M1": {"A": {}}}
    cens = {"M1": {"A": 0}}
    for seed in R.SEEDS:
        t0 = truth.generate_truth("M1", seed, P_t, L_t, splits)
        arm_inputs = R.build_arm_inputs(P_t, t0, shuf, fam_perm, rand_p)
        res, n_cens_total = R.run_seed_arms(P_t, L_t, splits, device, Lt_dev,
                                            arm_inputs, "M1", "A", seed)
        results["M1"]["A"][str(seed)] = res
        cens["M1"]["A"] += n_cens_total
    out = {
        "schema": "MetaSieve.StageQ2d1d.LADDER.v1",
        "preregistration_sha256": R.PREREG_SHA,
        "results": results,
        "censored_counts": cens,
        "repro_A_value_level": {},
        "recovered": True,
        "recovery_note": ("M1-A-only recovery run with the exact frozen "
                          "runner_d code path after the original ladder "
                          "crashed at M2 truth generation (NameError: PCA_VT "
                          "not defined in truth_d.py). M2/M3/NC1/NC2 levels "
                          "and the value-level reproduction checks were never "
                          "produced by the original ladder."),
    }
    json.dump(out, open(HERE / "Q2D1D_LADDER.json", "w"), indent=1)
    print("wrote Q2D1D_LADDER.json (recovered M1:A only)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
