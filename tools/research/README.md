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
- `stageS_sar_field/` - cross-target protein-conditioned SAR field, a
  conservative potential-difference formulation trained on signed within-target
  dpK (REJECTED: 1 of 6 gates for the protein arm, 0 of 6 for the
  counterfactual arm; a shuffled-protein control reproduces the entire measured
  protein gain). First stage in the record mounted on the physical split view
  with no checkpoint selection at all. Retains the ligand-only positive result
  (Pearson +0.160 out of component, +0.400 on the most novel tercile).
- `stageT_mmp/` - true-MMP transformation-space test, **SCOPE CORRECTED**. T0
  measurement reliability (100% provenance recovery, hash-verified; same-panel
  difference variance ~0.858 pK^2 vs cross-panel ~1.221; several quantities
  recorded as not identifiable) stands. T2 **FAILED** its frozen gate 3/10, and
  that rejection of the **coarsened-key pooled-protein discriminator** stands.
  But its `exact_key` omitted the shared core, so 40.4% of fit `D` rows compared
  targets with disjoint cores (uncancelled residual median 0.269 pK vs a truth
  sd of 0.804): the "cancels mu_tau exactly" claim and the **global closure of
  the family are WITHDRAWN**. See `CORRECTION_20260817_CORE_KEY.md`; the
  1,112-key figure must not be reused (core-inclusive count 1,001).
- `stageU_mmp_interaction/` - **STOPPED and superseded**; no number used as
  evidence. Correct chemistry (core-inclusive key, core-consuming descriptor,
  interaction-variance gate, local region operator) but frozen 4 minutes after
  Stage T's metrics were read and missing four load-bearing controls. Its
  preregistration is retained unedited. Audit in
  `stageV_core_mmp/STAGE_U_GOVERNANCE_AUDIT.md`.
- `stageV_core_mmp/` - corrected core-inclusive Phase-1 test, **STOPPED BEFORE
  TRAINING; no model built**. Inherits every Stage U threshold verbatim, adds
  the four missing controls and two repairs, loosens nothing. V0 admission
  passes on all five size thresholds (1,001 core-inclusive rich keys) but
  **fails two concentration caps** (one target = 29.63% of fit observations);
  V0b leaves the primary internal surface at **32 rows / 4 components** with
  **0 internal rich keys**; V1 interaction variance `MS_effect` 0.4517 vs
  supervision noise 0.8576, `theta` -0.4059 [-0.6889, -0.0577]. Verdict:
  the requested estimand is **NOT ESTIMABLE** on this corpus - a statement about
  support, **not** biological absence. 31 structural/leakage tests pass.
- `GOAL_ACTIVE.md` - durable goal record (goal tools unavailable in the
  originating session).
