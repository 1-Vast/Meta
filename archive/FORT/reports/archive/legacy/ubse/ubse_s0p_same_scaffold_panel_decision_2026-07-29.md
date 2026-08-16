# UBSE-S0P whole-component split decision

Date: 2026-07-29  
Decision:
`STOP_UBSE_S0P_SAME_SCAFFOLD_INDEPENDENT_TOPOLOGY_INADEQUATE`

## Result

The earlier whole-component panel program retained 2,170 same-scaffold
panels after removing five generic scaffold tokens, but failed:

- largest homology-scaffold-PubMed conflict component:
  1,266 / 2,170 = `58.3410%`;
- balanced whole-component fold counts:
  `1,266 / 226 / 226 / 226 / 226`;
- fit side:
  904 panels / 718 exact targets, below the frozen 1,000 / 600 joint floor;
- largest homology share:
  `5.3456%`, above the frozen 5% cap.

Resource ceiling 686 and deterministic conflict-free packing 482 show that
the failure was not marginal entity count. It was the inability to assign the
entire observed panel population to balanced independent folds.

The audit decoded no affinity and used no coordinates, protected features or
labels, or sealed outcomes. It did inspect whether residue lists were
nonempty/malformed; later G0P improved the topology firewall by avoiding the
binding-residue field entirely.

## Relationship to G0PB

S0P remains a binding STOP for population-wide five-fold claims. The later
G0P/G0PB program asks a narrower question: whether one can preselect an
88-panel conflict-free mechanism audit and remove all direct homology,
scaffold, and PubMed conflicts from its training substrate. G0PB may admit
that pilot, but cannot be cited as reversing S0P's full-population result.

Artifacts:

- `reports/active/ubse_s0p_same_scaffold_panel_preregistration_2026-07-29.md`
- `reports/active/ubse_s0p_seed1729.json`
- `research/ubse_s0p.py`
- `tests/test_ubse_s0p.py`
