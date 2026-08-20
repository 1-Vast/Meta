# stageCIIP2_olr_potential_20260820

OLR-Potential (Orthogonalized Ligand-Routed interaction potential) — CIIP-2
successor stage. Base prereg SHA a7b17e8a...; ADD-1 (SPB construction)
SHA aa8d06af.... This stage never touches production model/ or scripts/.

## Artifacts

| file | role |
|---|---|
| PREREGISTRATION.md | frozen rules (base) |
| ADDENDUM_ADD1_SPB.md | frozen SPB split correction (pre-training) |
| ADDENDUM_ADD2_EVAL.md | frozen evaluation-definition correction (pre-real-data) |
| olr.py | model, data, metrics, controls (deployed path: SGD only) |
| runner.py | arms A0-A5 + C-* controls, qualification, smoke, phase5 |
| gen_erased.py | X-erased ESM state generator (me_/we_ keys per pair) |
| ERASED_ESM.npz | erased residue states, 49/49 pairs, sha a8f36905... |
| tests/test_structure.py | structural tests T1-T10 + split checks (12 pass) |
| QUALIFICATION.json | instrument qualification (INSTRUMENT_UNDERPOWERED) |
| PHASE4_SMOKE.json | real-data single-seed smoke + gates (gate b FAILED) |
| PHASE_REPORT.md | phase 3-5 consolidated record + verdict ladder |
| diag_fit*.py, diag_sweep.py, diag_lastepoch.py | synthetic instrument diagnostics |

## Implementation amendments (documented, none result-driven)

AM-1: the prereg's "C-erased: A5 evaluated on X-erased sequences" is
mathematically degenerate for potential models because erased WT and erased
MT states are identical, so the contrast is identically zero. Implemented
instead as C-erased-site: the VARIANT construct is replaced by its
X-erased states (site residue masked, context preserved), WT kept real.
Full-erasure equality (we_ == me_) is asserted by structural test T9.
Additionally the informative train-side analogue is implicitly covered by
C-randprot and C-wrongmut.

AM-2: A3-oid is provably an identity on panel-centered contrasts (any
separable main-effect component added to s cancels; asserted in T3). A3 is
therefore retained as an A2-equivalence/determinism check plus raw-endpoint
main-effect accounting, not as a distinct modeling arm.

AM-3: ligand_prior/family_prior build a padded panel-width profile before
averaging (per-pair ligand subsets have width 179-183, not 183).

## Reproduction

    /d/anaconda/envs/drug/python.exe tests/test_structure.py
    /d/anaconda/envs/drug/python.exe gen_erased.py
    /d/anaconda/envs/drug/python.exe runner.py qualify
    /d/anaconda/envs/drug/python.exe runner.py smoke 11
    /d/anaconda/envs/drug/python.exe runner.py phase5
