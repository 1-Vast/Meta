# Pan-Family Single-Cold-start (PFSC) — preregistration; gate PFSC-0

Date: 2026-07-26. Route code **PFSC**. Authorized by the scope decision taken after the TR reframe
closed (TR0_PREMISE_FAIL_STOP) and the hidden-dense-panel hunt showed strict dual cold-start remains
data-blocked on all audited open data. PFSC is a **deliberate, declared relaxation of the ligand-cold
axis**: it tests **single cold-start** (novel target, *known* ligands), not the strict dual-cold task.
No dual-cold claim is made from PFSC. It is a mechanism probe, not confirmation.

## 1. Substrate (frozen, train-only, no dev/confirmation label read)

The pan-family promiscuous-ligand block discovered inside the ChEMBL-37 dual-cold **train** partition
(registry sha256 `0e754f73...6b784`), endpoint **pKi** only:

* ligands measured across `>= 10` distinct homology clusters ("promiscuous"): 771 ligands;
* their targets: 187, spanning **172 homology clusters** (genuinely pan-family, vs within-kinome
  panels Metz/KirHub);
* frozen ESM-2 pooled 1280-d features cover all 187 targets
  (`dataset/public/chembl_37/processed/dualcold/target_esm2.npz`).

The block registry and its sha256 are frozen by the PFSC-0 runner before scoring. `docs` and `assays`
are pipe-delimited ID sets per edge; both are parsed to sets for isolation.

## 2. Scientific hypothesis

Within the kinome, KLIFS group was not a resolvable transfer resolution (TR-0 G0 failed: groups too
similar). **Across families**, near-family protein neighbours are genuinely different from distant
ones, so protein-similarity conditioning should reorder a novel target's *known* ligands better than
ligand potency alone, and that gain should **require near-family proximity** and **survive
cross-document isolation**. If so, protein/family transfer is load-bearing for pan-family single
cold-start — the positive that within-kinome could not produce — and the dual-cold failure is
localised specifically to the novel-ligand axis.

## 3. Estimand and predictor

Single cold-start, leave-homology-cluster-out. For a held target `t` (homology cluster held out) and
one of its block-ligands `L`, predict `t`'s pKi for `L` from OTHER targets that measured the **same**
`L`, restricted to:

* **cross-family**: different homology cluster from `t` (automatic under leave-cluster-out);
* **cross-document (mandatory)**: the neighbour's `docs` set for `L` is **disjoint** from `t`'s `docs`
  set for `L` (removes the same-paper/same-protocol shortcut).

Arms (all use the identical cross-doc cross-family neighbour set per `(t,L)`):

* `ligand_only` (B0): uniform mean of neighbour pKi for `L` (protein-agnostic potency).
* `esm_kernel`: neighbour mean weighted by `clip(cos(ESM_t, ESM_neighbor),0,inf)^2` (near-favouring).
* `protein_shuffle`: ESM vectors permuted across targets (train-fold).
* `random_protein`: random Gaussian target vectors.
* `far_kernel` (proximity-cold): neighbour mean weighted by `clip(max_sim - sim,0,inf)^2`
  (distant-favouring); same neighbour set as `esm_kernel`.

ESM features are PCA-whitened on the fold's training targets (64 components, as in A1). Nonparametric;
no learned parameters, no training, no tuning. Seed 1729.

## 4. Scoring, unit, power

Within-target Spearman between predicted and true pKi over `t`'s **scorable** ligands. A ligand is
scorable if it has `>= 3` cross-doc cross-family neighbour measurements; a target is scored if it has
`>= 5` scorable ligands (feasibility audit: 90 targets / **83 homology components**). Statistical unit
= homology component (component-macro Spearman); grouped bootstrap 10,000 draws. MDE80 at paired
SD 0.10 over 83 components ~ 0.016; threshold `max(0.03, MDE80) = 0.03`.

## 5. Frozen gates

* **G1 (protein value):** paired `esm_kernel - ligand_only` `>= 0.03` with grouped LCB95 `> 0`.
* **G2 (specificity):** paired `esm_kernel - protein_shuffle` and `esm_kernel - random_protein` both
  LCB95 `> 0`.
* **G0 (proximity is load-bearing — the identifying control):** paired `esm_kernel - far_kernel`
  LCB95 `> 0`. This is the pan-family analogue of the TR-0 own-group-cold control; passing here is the
  result within-kinome could not produce.

A G1, G2 or G0 failure returns `PFSC0_FAIL_STOP`. A full pass returns
`PFSC0_PASS_AUTHORIZE_PFSC1` and authorises PFSC-1: a proper single-cold-start predictor with
calibration/selective prediction and a matched-capacity control, plus a robustness rerun under
**cross-assay** (not only cross-document) isolation. No dual-cold claim, no confirmation/sealed
access, no multi-seed, no long training is authorised by PFSC-0.

## 6. Leakage, confounding, risks

* **Target axis:** homology-cluster held out (accession/homology-disjoint). ✓
* **Ligand axis:** WARM by construction (single cold-start). Declared; this is not dual-cold. ✓
* **Document axis:** strict cross-doc disjointness enforced (primary control against protocol
  shortcut); cross-assay reported as robustness in PFSC-1. ✓
* **Promiscuity bias:** block ligands bind many targets, so they are *less* selective, biasing the
  target-specific reordering signal DOWN — a conservative bias; a positive despite it is real.
* **Pseudoreplication:** unit = homology component; grouped bootstrap; neighbours never treated as
  independent pairs.
* **Not a dual-cold result:** PFSC cannot and does not claim strict dual cold-start; it localises where
  the dual-cold signal lives and delivers a single-cold-start capability if the gates pass.

`sealed_test_consumed=false`; `confirmation_labels_read=true` (pre-existing; PFSC reads train only).
