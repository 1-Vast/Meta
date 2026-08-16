# Gate PA preregistration - dense-panel interaction identifiability and feature alignment

Registered 2026-07-25, before any PA statistic was computed. Train-only. No panel development cell,
no panel confirmation cell and no ChEMBL confirmation label is read.

## Why a new substrate rather than a new head

Every gate on the sparse observational ChEMBL graph failed on the same object, not on the same
model. `BM2_PIRR_GATE_P0_FAIL_STOP` could pack only 160 entity-disjoint 2x2 blocks;
`P0_CYCLE_A_BIOLOGICAL_FAIL_STOP` measured a nuisance-orthogonal projected label SD of 0.35596 pK
against the frozen 0.5 minimum with 51.451% of residual energy on 1% of ligands; and
`P0_HOAB_X0/X1` found neither bilinear nor RBF-kernel alignment in that residual
(`p=0.4386`, `p=0.5292`). task.md names the only admissible new research input: an
endpoint-consistent dense cross-target panel that directly identifies protein-conditioned ligand
reordering.

That panel already exists inside the frozen local ChEMBL-37 extract. Document `CHEMBL1201862`
(Metz et al. 2011 kinase Ki profiling) contributes 36,241 pKi records in one publication with one
endpoint, giving a (20 ligands, 5 targets)-dense bipartite core. Endpoint and document are therefore
*constant*, not nuisances to firewall. The registered substrate is
`tools/kinase_panel_registry.py` -> `dataset/public/chembl_37/processed/panel_metz/`, registry
sha256 `94da6bb5a59c2911672fde982530c8dd6a673c194b2b2d7b4638df7768c8173e`: 12,574 train cells over
112 targets in 101 homology components and 619 anchor ligands in 344 scaffolds; 5,212 development
cells over 271 scaffold-disjoint query ligands (max anchor Tanimoto 0.9091); 1,427 sealed
confirmation cells over 30 targets and 195 ligands. 21 panel targets whose homology component
appears in the sealed confirmation split of the main ChEMBL registry were dropped, so this substrate
cannot contaminate that sealed set.

## What PA asks

PA separates two questions that every previous gate confounded.

1. **Is the target-ligand interaction a real, reproducible, systematic object on this substrate?**
   This does not involve protein features at all.
2. **Does the frozen protein representation align with it?** This is the identical HOAB-X0 test on a
   substrate where question 1 has an answer.

## Criteria (frozen before running)

| id | criterion | threshold | provenance |
|---|---|---|---|
| PA1 | numerical validity of the exact fixed-effect projection | relative KKT < `1e-8`, idempotence < `1e-7`, LSMR/LSQR disagreement < `1e-6` | P0-Cycle-A, unchanged |
| PA2 | nuisance-orthogonal projected label SD | `>= 0.5` pK | P0-Cycle-A, unchanged |
| PA3 | top 1% of any nuisance axis carries at most half the residual energy | `<= 0.50` | P0-Cycle-A, unchanged |
| PA4 | interaction reproducibility: held-out-cell rank-8 completion of the residual beats the value-permuted control, target-grouped bootstrap | `LCB95 > 0` | zero threshold; no tuned constant |
| PA5 | frozen-feature alignment: adaptive score against 4,096 graph-exposure-matched target-feature permutations | `p_adapt <= 0.01` | HOAB-X0, unchanged |

The nuisance space is the target and ligand incidence design with unit edge weights. Assay is
exactly 1:1 with target on this panel and document is constant, so neither adds a column; this is
stated in advance rather than discovered afterwards.

PA4 fits a rank-8 alternating-least-squares factorization of the projected residual on 4 of 5
deterministic cell folds and scores the held-out fold, against the identical procedure applied to a
seed-fixed permutation of the residual values. The paired statistic is the per-target reduction in
held-out squared error; inference is a bootstrap over targets. Rank 8 is the same function-space
dimension the BM0/BM1 posterior uses; it is not selected from PA outcomes.

PA5 reuses the frozen bases without modification: 32 centered PCA coordinates of pooled ESM-2 650M
train-target embeddings, and 64 centered PCA coordinates of the 64-bin count-Morgan plus ten
physicochemical descriptor basis over unique panel ligands. Permutations are restricted to blocks
matched on label-free graph exposure. No affinity label selects either basis.

## Explicitly not applied

P0's `entity_disjoint_blocks >= 100` and `mde80_sign_accuracy_delta <= 0.10` packing criteria are
**not** PA criteria and P0's failure is not revisited. Those thresholds exist to defend against
uncontrolled shared provenance across documents. A single-document panel controls provenance by
construction, and no single panel of 112 targets can ever supply 196 mutually target-disjoint
blocks, so importing that criterion here would be a structurally unpassable test rather than a
scientific one. Entity-disjoint packing counts and the cross-document sign agreement are reported as
diagnostics with P0's `0.65` reference value, and are non-gating.

## What a pass authorizes

A PA pass authorizes exactly one next step: the separately preregistered panel power audit, which
freezes a Gate PB threshold from B0 retraining noise before any target-conditioned arm is scored. It
does not authorize multi-seed runs, Hierarchical MoT, long training, confirmation access or any
change to the frozen thresholds above. A single failed criterion stops the panel route and is
recorded in `history.md` with its cause.
