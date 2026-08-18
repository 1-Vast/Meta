# Stage X0 correction plan (round 2, 2026-08-18)

Authority: user directive + `report/STAGE_X_ROUND1_REVIEW_20260818.md`. The frozen
preregistration (`STAGE_X0_PREREGISTRATION.md`, SHA-256
`03cdc907df3e778f5fe79fb1a238d35ebb6ece5e9e743db181728ba6b25e9683`) is never
modified; thresholds are never changed. Round-1 artifacts that encode the
defective implementation (`X0_INSTRUMENTS.json`, `X0_PLANTED.json`,
`x0_planted.py`, `x0_instruments.py`) are preserved as negative evidence and
are not overwritten. Corrected implementations are new files with the `x0_iN_`
prefix and new artifact names.

## Step 0 — audit baseline (DONE in this round)
- [x] SHA-256 of frozen prereg matches X0_PREREGISTRATION_SHA256.txt.
- [x] No Python training processes running.
- [x] File ownership established (round-1 agent committed files vs reviewer
      uncommitted files) — see history.md.
- [x] Q9DE49 = Danio rerio pdgfra (wrong species for human PDGFRalpha);
      P15056 residue 599 is T (so V599E is historical numbering for canonical
      V600E); KLIFS API resolves PDGFRa -> P16234.

## Step 1 — reproduce the reported defects as failing tests
(see plan details in history.md; tests: tests/test_x0_corrections.py)

## Step 2 — data foundations: P16234 fetch, BRAF renumbering evidence,
PubChem SMILES/scaffolds for 183 compounds, KLIFS pocket sequences.

## Step 3 — corrected WT-mutant pair table (x0_pair_table.py -> X0_PAIR_TABLE.json)

## Step 4 — corrected I2 (x0_i2.py -> X0_I2.json, X0_I2_REPORT.md)

## Step 5 — corrected I1 planted-signal control (x0_i1.py -> X0_I1.json, X0_I1_REPORT.md)

## Step 6 — I6 production-dataflow suite (tests/test_x0_production_dataflow.py)

## Step 7 — I3 (x0_i3.py), I4 (x0_i4.py), I5 (x0_i5.py)

## Step 8 — X0_RESULT.json + X0_REPORT.md + audits + ledger + task/history updates

## Step 9 — final pytest runs

Full detail (verification criteria per step) is in the round-2 history entry.
