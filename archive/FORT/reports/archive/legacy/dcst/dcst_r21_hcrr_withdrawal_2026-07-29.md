# DCST-R21 HCRR withdrawal

Date: 2026-07-29  
Decision: `WITHDRAW_HCRR_BEFORE_LABEL_LOAD`

R21-T1 was preregistered but never implemented or executed. No ChEMBL TRAIN
affinity was loaded for R21, and no CUDA model was trained.

The new source-identifiability audit supplied by the user changes the
admissible question. A high-confidence affinity curriculum still lacks the
required Stage-1 information object:

- exact same-target paired complexes;
- real ligand-conditioned structural deltas;
- repeated directed edits across targets/families;
- independent provenance and target-domain support.

R18-R20 already show that PDB-link counts and combinatorial affinity
rectangles do not supply those conditions. Continuing R21 would therefore
test another outcome curriculum on the same information-deficient substrate,
not the proposed RBSDD mechanism.

The frozen R21 preregistration is retained as an audit artifact. It is
withdrawn, not failed, and may not be cited as an empirical negative result.

