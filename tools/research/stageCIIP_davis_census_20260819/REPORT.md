# Stage CIIP-0b Davis panel census — INSUFFICIENT ALONE (2026-08-19)

Prereg SHA 2dd8b708...; artifact CENSUS.json; read-only audit.

## Verdict

**INSUFFICIENT for the full CIIP program on Davis alone** (frozen stop
rule): held-out parents = 10 with pairs / 7 with >=2 mutants < 10. The
panel is otherwise strong: 67 usable WT-variant pairs, median 33 common
ligands per pair (min 4, max 43), all 67 pairs have >=3 common ligands,
and the log10(mutant/WT) variance median is 0.254 (p10-p90
0.054-1.31) — real mutation-effect heterogeneity exists.

| Item | Value | State |
|---|---|---|
| 1 usable WT-variant pairs | 67 | OK |
| 2 identical ligands/pair | median 33 (min 4, max 43) | OK |
| 3 parents/mutants/fusions | 379 genes; 388 WT / 54 mutant rows; 0 fusion-flagged | OK |
| 4 condition completeness | single assay condition class; measured-cell fraction 29.6% | PARTIAL |
| 5 duplicates/saturation | duplicate (gene,mutant) keys listed (ABL1 13 YES-rows etc.); NA 70.4%; cap-9900 0.04%; NA semantics (untested vs >10 uM) UNKNOWN | PARTIAL |
| 6 endpoint | Kd [nM] competition binding; larger = weaker | OK (never relabeled) |
| 7 connectivity | 10 parents with pairs, median 4.5 pairs/parent | OK |
| 8 held-out parents | 7 parents with >=2 mutants (<10) | **FAIL** |
| 9 coverage/variance | 67 pairs with >=3 ligands; logratio var median 0.254 | OK |
| 10 availability | local, SHA-pinned, nbt.1990 SI | OK |

## Consequences (frozen)

- CIIP-1A/1B training on Davis ALONE is NOT authorized.
- The same 10-item census must run on the other local panels
  (Anastassiadis 2011, Duong-Ly) and the COMBINED surface re-checked
  against the frozen thresholds before any CIIP preregistration.
- This is a data-admission negative, not a biological statement.
