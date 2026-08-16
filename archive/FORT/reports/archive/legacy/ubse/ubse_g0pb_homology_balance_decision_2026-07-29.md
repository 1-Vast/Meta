# UBSE-G0PB homology-balance decision

Date: 2026-07-29  
Decision:
`REQUEST_UBSE_G1_CENTERED_CONTACT_STUDENT_PREREGISTRATION`

## Outcome

The one permitted removal-only correction passes every unchanged G0P gate.
Two homology components exceeded the pre-removal 5% share cap; removing their
200 panels left:

| Metric | Result |
| --- | ---: |
| Balanced panels / pair contrasts | 1,412 / 3,383 |
| Exact targets / homology components | 957 / 492 |
| Scaffolds / PubMed IDs | 897 / 871 |
| Conflict components | 385 |
| Largest conflict component | 17.6346% |
| Resource ceiling | 492 |
| Deterministic conflict-free packing | 450 |
| Largest homology / scaffold / PubMed share | 3.7535% / 0.9207% / 2.4788% |
| Frozen audit manifest | 88 panels, all three resources unique |
| Residual closed training substrate | 1,324 panels / 3,243 contrasts |
| Exact ChEMBL-TRAIN accession support | 239 / 957 = 24.9739% |

No second balancing iteration or result-selected pruning was performed.
Binding residues, affinity fields/values, coordinates, protected
features/labels, and sealed outcomes remained unread.

## Authorized continuation

G0PB establishes that a strict same-scaffold, within-publication residual
experiment is topologically executable. It does not establish that sequence
and 2D chemistry can predict the residual.

The only authorized next step is a separately frozen G1 pilot:

- use the 88-panel manifest as the untouched homology-, scaffold-, and
  PubMed-cold audit set;
- derive a disjoint validation set from the remaining train candidates before
  reading contacts;
- train a frozen-protein-encoder residue-contact student on CUDA;
- compare a multiplicative residue-ligand interaction arm against the exact
  same-parameter additive two-tower null;
- score only within-panel centered contact contrasts;
- require ligand-feature derangement and protein-free destruction;
- retain affinity, coordinates, development/confirmation labels, and sealed
  outcomes behind the firewall.

Passing G1 may request a ChEMBL-TRAIN metadata/support transfer gate. It does
not directly authorize Stage-2 affinity fitting or a causal
ligand-intervention claim.

## Artifacts

- Preregistration:
  `reports/active/ubse_g0pb_homology_balance_preregistration_2026-07-29.md`
- Result: `reports/active/ubse_g0pb_seed1729.json`
- Manifest:
  `dataset/public/biolip2/processed/ubse_g0pb_panels.parquet`
- Implementation: `research/ubse_g0pb.py`
- Tests: `tests/test_ubse_g0pb.py`
