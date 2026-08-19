# One-shot generator: truth_e.py = frozen truth_d.py + documented repairs.
from pathlib import Path

SIG = Path(r"D:\MetaSieve\tools\research\stageX_csc_signal")
SRC = SIG / "stageQ2d1d_spanrestricted_interaction_20260818" / "truth_d.py"
DST = SIG / "stageQ2d1e_spaninit_interaction_20260818" / "truth_e.py"

t = SRC.read_text(encoding="utf-8")


def rep(s, old, new, label):
    n = s.count(old)
    if n != 1:
        raise SystemExit(f"{label}: expected 1 occurrence, got {n}")
    return s.replace(old, new)


t = rep(t,
'''"""Stage Q2d-1c truth generator + feature-space oracle precheck.
Prereg SHA: 25b8b9129120d0a770ba353cd56a8a388dd847778e4e5cb2b488ee0cbfee7106.
Fixes vs Q2d-1b: PCA-32 protein features, resolved-only ligand pool (157),
feature-smoothed double centring. ALS is DIAGNOSTIC ONLY.
"""''',
'''"""Stage Q2d-1e truth generator: byte-copy of the frozen Q2d-1d truth_d.py
(source prereg 25b8b912..., stage SHA baf4bb72...) with documented repairs
(addendum AD1_M2_PCA_VT_FIX.md). M2 gets PCA_VT loaded from the frozen
label-free feature artifact; NC1 realizes its frozen description (I = 0,
main effects only) incl. the undefined 0/0 normalization; NC2 returns
A=None/B=None (no feature map). The 1d copy raised UnboundLocalError on
these never-before-executed branches. M1/M2/M3 streams are bit-identical
to truth_d. Everything else untouched.
"""''', "docstring")

t = rep(t,
'''PREREG_SHA = "baf4bb72df02e9411d6b8d4815302ec91c7526cc15447b6e80cd06383d546991"
TAU = 1.0
RANK = 4''',
'''STAGE_D = HERE.parent / "stageQ2d1d_spanrestricted_interaction_20260818"
PREREG_SHA = "baf4bb72df02e9411d6b8d4815302ec91c7526cc15447b6e80cd06383d546991"
TAU = 1.0
RANK = 4''', "STAGE_D const")

t = rep(t,
"LIG_PROJ_DIM = 48\n",
"LIG_PROJ_DIM = 48\n\n# AD1 repair 1: the frozen M2 definition needs the pre-compression PCA\n# basis; it lives in the frozen label-free feature artifact.\nPCA_VT = np.load(STAGE_D / \"q2d1d_features.npz\",\n                 allow_pickle=False)[\"PCA_VT\"].astype(np.float64)\n",
"PCA_VT load")

t = rep(t,
'''    if mechanism == "NC2":
        I_raw_all = F_r @ F_l.T
    else:
        I_raw_all = (P_t @ A) @ (L_t @ B).T
        if mechanism == "M3":
            I_raw_all = np.tanh(I_raw_all / np.sqrt(RANK))
''',
'''    if mechanism == "NC2":
        I_raw_all = F_r @ F_l.T
    elif mechanism == "NC1":
        # AD1 repair 2: NC1 = no interaction (frozen description: I = 0,
        # main effects only); no feature-conditioned map exists.
        I_raw_all = np.zeros((n_rows, n_lig))
        A = None
        B = None
    else:
        I_raw_all = (P_t @ A) @ (L_t @ B).T
        if mechanism == "M3":
            I_raw_all = np.tanh(I_raw_all / np.sqrt(RANK))
''', "NC1 dispatch")

t = rep(t,
'''    elif mechanism == "NC2":
        F_r = rng.normal(0, 1, size=(n_rows, RANK))
        F_l = rng.normal(0, 1, size=(n_lig, RANK))
''',
'''    elif mechanism == "NC2":
        F_r = rng.normal(0, 1, size=(n_rows, RANK))
        F_l = rng.normal(0, 1, size=(n_lig, RANK))
        # AD1 repair 2: NC2 has no feature-conditioned map.
        A = None
        B = None
''', "NC2 repair")

t = rep(t,
'''    sd_tr = I_c[tr_i, tr_j].std()
    I = I_c / sd_tr * TAU
''',
'''    sd_tr = I_c[tr_i, tr_j].std()
    if sd_tr > 0:
        I = I_c / sd_tr * TAU
    else:
        # AD1 repair 3: NC1 has zero interaction; the frozen normalization
        # is 0/0 undefined there and the description says I = 0.
        I = np.zeros_like(I_c)
''', "NC1 sd guard")

DST.write_text(t, encoding="utf-8")
print("truth_e.py regenerated,", len(t), "chars")
