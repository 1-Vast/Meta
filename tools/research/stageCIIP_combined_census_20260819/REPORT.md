# Stage CIIP-0c combined local-panel census — SUFFICIENT (2026-08-19)

Prereg SHA 6952ed1a...; artifact CENSUS.json; read-only audit.

## Verdict

**SUFFICIENT: the Duong-Ly panel ALONE supports both CIIP-1A and
CIIP-1B** under the frozen stop rules:

- usable single-mutant WT-variant pairs: **70**
- identical ligands per pair: median **183** (min 179, max 183) — a
  near-complete within-panel screen
- parents: 21 with pairs; **12 with >= 2 single-mutant pairs**
  (ABL1, ALK, CKIT, CMET, EGFR, FGFR3, FGFR4, FLT3, LRRK2, PDGFRA, RET,
  TIE2) — exceeds the frozen CIIP-1B threshold of 10 held-out parents
- multi-mutant rows excluded from the single-mutant estimand and
  counted separately: 6 (CKIT V559D/T670I, V559D/V654A; EGFR
  d747-749/A750P, d747-752/P753S, d746-750/T790M, L858R/T790M)
- endpoint: % inhibition (larger = stronger inhibition); observed
  range -12.5 .. 191.3 (out-of-bounds values recorded; windsorization/
  interval treatment is a future-prereg decision, never silent)
- NA fraction 0.23%; 0 unparsed rows
- ligand-dependent mutation-effect variance (pairs with >= 2 ligands):
  median 138.3 (%²) — large, real mutation-effect heterogeneity (H1
  admission relevant)

## Admissibility classification (frozen)

| Surface | Classification |
|---|---|
| Davis Kd panel | same-platform: CIIP-1A yes, CIIP-1B no (7 held-out parents) |
| Duong-Ly % inhibition panel | same-platform: CIIP-1A yes, CIIP-1B **yes** (12 parents) |
| Anastassiadis 2011 | 0 within-panel WT-variant pairs; cross-endpoint replication reference only |
| Cross-panel pooling (Davis+Duong-Ly) | **FORBIDDEN** (Kd vs % inhibition are different endpoints; never relabeled) |

## Consequences

- CIIP-1A/1B can now be preregistered on the Duong-Ly panel (single
  platform, single endpoint, functional % inhibition — never called
  pK/Ki/Kd). Q2d terminal archival must be complete first (diagnostic
  running; not a training blocker but a governance ordering).
- KiRHub remains DATA BLOCKER; any future KiRHub use requires a
  first-hand, legal, locally SHA-pinned raw table.
- No return to BindingDB exact-MMP retraining.
