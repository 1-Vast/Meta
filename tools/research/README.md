# Research workspace

All unadmitted model and training experiments start here. Each family must have
a preregistration, isolated implementation, structural tests and explicit stop
gates. `main.py` must not import this package.

After a family passes its synthetic and real-data admission gates, move reusable
model code to `model/`, workflow code to `scripts/`, and regression contracts to
`tools/tests/`. Then remove the research copy so there is one implementation.
Failed families leave only a compact decision in `report/` and are deleted from
this directory.

## Stage index (2026-08-17/18 cycle, all with PREREGISTRATION.md + REPORT.md)

- `stageD_level_panel/` - D0 audits (decomposition, level identifiability,
  level anatomy, occupancy strata, doc transfer), Stage E (panel level head +
  orthogonal routing, REJECTED), shared evaluation/contrast tooling.
- `stageF_pairwise/` - pairwise learned transport (REJECTED).
- `stageG_esm650/` - ESM-650M residue-input lane (single-seed promising,
  NOT CONFIRMED across three seeds); multiseed comparison tooling.
- `stageI_lm/` - live ESM-150M LoRA lane (REJECTED; engineering note on
  chunked-LoRA memory bounding).
- `stageJ_assay/` - D0c journal identifiability + assay-aware level head
  (REJECTED by resolved k=2/3 ranking degradation).
- `stageK_contrastive/` - contrastive coembedding; K-REG = first all-k
  resolved MSE improvement across three seeds, NOT CONFIRMED on centered.
- `stageL_gated/` - support-gated level head (REJECTED; best k=0
  calibration on record).
- `stageM_chemberta/` - ChemBERTa-77M ligand LM probes (REJECTED at
  identifiability).
- `stageN_audit/` - final boundary audit (bitwise verification of all
  load-bearing numbers; seal audit).
- `stageP_go/` - ProteinKG25 GO annotation probes (REJECTED at
  identifiability).
- `stageQ_frozenhead/` - decoupled frozen-feature level head (REJECTED;
  closes the level-head composition axis).
- `stageR_daviskiba/` - FROZEN Davis/KIBA boundary-check plan
  (PREREGISTRATION.md; NOT AUTHORIZED, NOT RUN; inventories only).
- `GOAL_ACTIVE.md` - durable goal record (goal tools unavailable in the
  originating session).
