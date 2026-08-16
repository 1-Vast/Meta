# AMOB-DTA evaluation — Assay-Monotone Ordinal Bridge for Dual-Cold DTA

Date: 2026-07-26. Evaluated as a scientific hypothesis against the binding evidence in `history.md`
/`task.md`. AMOB is **admissible** (new information source: dense low-fidelity ordinal panels; new
observation structure: within-target within-assay monotone-ordinal constraints; new identification:
non-matched-molecular-pair ordinal contrasts, unlike the closed MISO four-state double-difference).
This report delivers the nine required outputs and one final verdict.

---

## 1. Falsification-first literature and novelty audit

**Internal prior work is binding and highly relevant.** An extensive prior program (MISO/AXIS/DICE/
PRISM, `history.md` L554–L1231) already built exactly on dense low-fidelity kinase panels (PKIS1/PKIS2
% inhibition, DAVIS pKd) to identify protein-conditioned interaction. Its verdicts:

* `DATA_LIMITED` / `R1B_EXTERNAL_DATA_REJECTED`. Low-fidelity continuous bridge reliability
  `rho≈0.336` (E1), max Pearson ≈0.58. **The present O0 residual agreement `+0.434` is the same
  phenomenon** — the low-fidelity signal is real — measured here on non-MMP within-target pairs.
* The decisive structural obstacle (L1157, L1173): *dense matched-molecular-pair coverage and
  continuous high-SNR endpoints are anti-correlated across public kinase data* (dense-MMP PKIS = %
  inhibition; continuous DAVIS/Metz = MMP-sparse), and the four-state protein-conditioned
  double-difference was single-source-dominated (PKIS2) with external validation unattainable.

AMOB partly *sidesteps* the MMP obstacle: within-target **ordinal** constraints do not require matched
molecular pairs, so they use far more of the dense panel than the double-difference did. This is the
one genuine mechanistic difference from MISO, and it is why O0 is measurable and positive here.

**External prior work substantially pre-empts the data hypothesis.** The staged Harmonic-Discovery
paper (Theisen, Wang, Ravikumar, Rahman, Cichońska, *bioRxiv* 2024, doi:10.1101/2024.03.07.583951,
`activity-integration`) **already integrates single-dose (low-fidelity) with IC50 (exact) to improve
compound-kinase bioactivity**, with pairwise kernel ridge regression (single vs multi-stage) evaluated
on compound-kinase, **scaffold** and homology-**cluster** splits. Related: pQSAR 2.0 (Martin et al.,
*JCIM* 2017) profile-QSAR; large-scale sparse-kinase modelling (Bosc et al./CDDLeiden, *JCIM* 2023,
Papyrus); ActFound few-shot bioactivity (GPL-3). Mathematical guidance: 1-bit/ordinal matrix
completion (Davenport et al. 2014), matrix completion under unknown monotone transforms (Ganti et al.
2015, single-index MC), inductive MC with side information (Jain & Dhillon 2013), multi-fidelity
co-kriging (Kennedy & O'Hagan 2000).

**Novelty verdict (falsification-first): AMOB's central data hypothesis is NOT new.** "Dense
low-fidelity single-dose identifies compound-kinase interaction that sparse exact data misses" is
published externally (Harmonic 2024) and tested internally to `DATA_LIMITED` (MISO). AMOB's admissible
residual novelty is narrow and specific: (i) an explicit per-(target,assay) **unknown-monotone** ordinal
likelihood (vs Harmonic's regression), (ii) a **target-centered interaction component that cannot
represent a target intercept** (isolating reordering from potency), (iii) **strict simultaneous
dual-cold** (target *and* ligand) under an **assay/document-provenance firewall** that neither MISO nor
Harmonic enforced. This qualifies as *one* innovation mechanism (the assay-monotone ordinal bridge),
but the honest expectation must be tempered: the underlying signal is already known to be real, and the
unresolved question is transfer under strict leakage control — exactly where MISO stopped.

**Why AMOB ≠ ordinary contrastive/pairwise/matrix-completion/multi-task learning.** Ordinary pairwise
ranking or contrastive pretraining pools all comparisons and would relearn generic potency (B0);
AMOB restricts supervision to *within-target within-assay* order and *target-centers* the interaction,
so the gradient carries only reordering, not potency or a target intercept. It is not standard matrix
completion because the observation is an unknown monotone transform of the latent affinity with
per-assay saturation/censoring, not the latent value itself; and it is not generic multi-task learning
because the tasks share one identifiability constraint (monotone-invariant order) rather than a shared
readout head.

## 2. Mechanism definition and identifiability

Latent interaction `z_td = f_θ(p_t, m_d)`; low-fidelity observation
`s_td^(a) = h_{t,a}(z_td) + ε`, with `h_{t,a}` an unknown per-(target,assay) monotone map. Supervision
is the within-target within-assay order only:
`sign(s_ti^(a) − s_tj^(a)) ≈ sign(f_θ(p_t,m_i) − f_θ(p_t,m_j))`, via a tie/censoring-aware
Bradley–Terry/ordinal likelihood; comparisons within a frozen saturation band are dropped or
reliability-weighted (here 5% < %inh < 95%). High-fidelity bridge `ŷ(t,d)=b(d)+γ·f̃_θ(t,d)`, with
`b(d)` the cross-fitted ligand-only B0 and `f̃_θ` target-centered (`Σ_d f̃_θ(t,d)=0` per target) so it
cannot encode a target intercept. pKi and pKd endpoint-separated; exact fine-tuning limited to a small
calibration readout.

**Identifiability assumptions.** (A1) monotone invariance: `h_{t,a}` monotone ⇒ within-assay order of
`s` equals order of `z` outside saturation. (A2) protein-conditioning is carried by `p_t` to unseen
targets — the load-bearing and historically-failed assumption. (A3) assay/document independence: the
low-fidelity order is not a same-campaign protocol artifact of the exact reference. **O0 tests A1
(passes) and the target-specificity of the signal (passes); it cannot test A2 (transfer) or A3 (assay
firewall) on the available data — see §3, §9.**

## 3. Open-data and licensing audit (Stage O0 step 1)

* **`activity-integration`** (Harmonic Discovery): single-dose `prepped_labeled_single_dose_data.csv`
  (423 targets, 8,359 compounds, 16,514 rows) + exact `prepped_activities_all_chembl_pubchem.csv`
  (465 targets, 210,862 rows; `measured` 141,193 / `inferred` 69,669). **LICENSE file empty** → *not
  admissible* as a training/redistribution source under OPEN_DATA_ONLY. Underlying data is ChEMBL +
  PubChem (open), so it is **reconstructible at source**, but the staged CSVs may be used only as a
  read-only feasibility diagnostic.
* **No assay or document identifiers** in either CSV → the mandatory assay/document-provenance firewall
  (O0 step 2) **cannot be evaluated**; the single-dose→exact agreement cannot be shown independent of
  same-campaign protocol correlation with this data.
* **Leakage vs the FORT substrate:** activity-integration accessions overlap my ChEMBL dual-cold splits
  — 136/559 train, **52/217 development, 26/163 sealed confirmation** — so the staged data is not
  firewalled against my dev/confirmation; admissible use requires accession+homology+scaffold+document
  +assay firewalling against those splits.
* **KCGS2.0.xlsx** = metadata only (compound/kinase lists), not the % inhibition matrix.
* **kinase-modelling** (CDDLeiden, Papyrus): empty LICENSE; Papyrus is aggregated ChEMBL → high
  leakage risk; not used.
* Prior OPEN-S already recorded PKIS2/KCGS as open (Zenodo, CC BY 4.0) but percent-inhibition — barred
  from *continuous-affinity confirmation*; AMOB's ordinal use is not barred by that clause, but the
  above license/provenance/leakage gaps are binding for a certified run.

## 4. Cheapest discriminating feasibility experiment (Stage O0 step 4) — DONE

Read-only diagnostic on staged data (`scratchpad/amob_o0_diagnostic.py`, seed 1729): does the
low-fidelity within-target order recover **protein-conditioned** reordering (residual after removing
generic ligand potency) that matches independent **measured** exact affinity, target-specifically and
non-circularly? Single-dose %inh at 1 µM (densest, n=8072), reliable band 5–95%, leave-target-out
generic-potency removal on both sides, statistical unit = target, grouped bootstrap.

| statistic | value | interpretation |
|---|---|---|
| raw within-target ρ (%inh vs exact) | +0.602 [+0.555,+0.645] | potency + protein, dominated by potency |
| generic potency recovery ρ(b_lo,b_ex) | +0.476 | low-fi recovers B0 potency |
| **residual (protein-conditioned) ρ, mixed** | **+0.347 [+0.291,+0.401]** (n=88) | signal beyond potency |
| **residual ρ, MEASURED-only exact** | **+0.434 [+0.375,+0.492]** (n=64) | non-circular; stronger |
| within-target shuffle null | −0.005 [−0.012,+0.002] | significance floor |
| **wrong-target control** | **+0.064 [−0.024,+0.154]** (n.s.) | **target-specific** |

Circularity ruled out: at fixed %inh∈[45,55) the exact value spans 2.9–9.7 (sd 0.95); single-dose
pairs are 92% `measured` (not model-`inferred`), and single-dose `activity_value` matches the
*measured* exact value (0% match to inferred). **O0 ordinal-identifiability sub-gate: PASS** — well
above the program's +0.03 materiality floor, target-specific, non-circular. This is the first
identifiable, protein-conditioned within-target reordering signal obtained in the program.

## 5. Leakage-safe model and data-flow design (for a certified O1)

Encoders (mature, non-innovative): ligand = pose-free Morgan/pharmacophore/graph features (transfer
across Bemis–Murcko scaffolds); protein = aligned functional-site coordinate (KLIFS-85 for kinases;
family-specific catalytic-site alignment otherwise) — **not** pooled whole-protein embeddings, per the
program's evidence that pooled ESM blurs reordering. Innovation module = the assay-monotone ordinal
bridge (§2). Data flow: reconstruct single-dose + dose-response from ChEMBL/PubChem **with assay and
document IDs**; firewall every evaluation target (+homologues) and compound (+scaffolds/neighbours) out
of *all* low- and high-fidelity training; train `f_θ` on within-target within-assay ordinal
constraints (saturation-filtered, frozen rule); freeze `b(d)`=B0; fit only `γ` + a small calibration
readout on exact; keep pKi/pKd separate; k≤1 falls back exactly to B0. No docking/poses/geometric
fields/expert-banks/free-form Bayesian adaptation.

## 6. Matched baselines and destructive controls (O1)

ligand-only B0; exact-affinity-only interaction; **AMOB**; raw low-fidelity regression (no ordinal);
ligand-only ordinal pretraining; pooled-embedding vs aligned-site; shuffled target reps; random
matched target reps; shuffled ligand assignment; within-target low-fidelity label permutation;
**assay/campaign permutation** (the key new firewall). Primary metric: paired target-component macro
Spearman under strict dual-cold; secondary: target-centered RMSE, negative-transfer rate, worst
component, calibration. Freeze required effect at empirical MDE before scoring; require positive
grouped LCB, no material RMSE loss, and collapse under both target and ligand destruction.

## 7. Preregistered O1 protocol (one-seed; conditional three-seed)

O1 is authorized only after O0 step-2 firewall is satisfiable, which requires the reconstructed,
assay-annotated, leakage-firewalled dataset (§3, §5). Then one short seed, matched budgets, identical
eval rows, all §6 arms. One-seed pass requires every registered condition (effect ≥ MDE with LCB>0;
protein-shuffle, random-protein, wrong-support, label- and **assay**-permutation all collapse the gain;
aligned-site beats pooled; RMSE safe). Any single failure stops. Three seeds only after a full
one-seed pass; confirmation/sealed access never in O0/O1.

## 8. Effect sizes, grouped CIs, empirical power

O0 (above): residual protein-conditioned ρ **+0.434 [+0.375,+0.492]** measured-only (n=64 targets),
wrong-target +0.064 [−0.024,+0.154], shuffle null ≈0. Power: 88/110/169 targets have ≥8/≥8/... usable
compounds; ample for the O0 unit test. **O1 power is unestablished** and is the historically binding
limitation (MISO: single-source dominance, external validation unattainable). No O1 effect size,
dual-cold transfer number, or multi-seed result exists — none is claimed.

## 9. Final verdict

**`OPEN_DATA_INSUFFICIENT_FOR_AMOB`** — with an explicit, load-bearing qualification.

The O0 ordinal-identifiability sub-gate **passed strongly and non-circularly** (residual ρ +0.434,
target-specific), refuting `LOW_FIDELITY_ORDINAL_BRIDGE_NOT_IDENTIFIABLE` and
`NO_ADMISSIBLE_NEW_MECHANISM_FOUND`: the dense low-fidelity ordinal source *does* supply identifiable,
protein-conditioned within-target reordering — the first such signal in the program. But the mechanism
**cannot be certified**, for three converging reasons:

1. **The mandatory assay/document-provenance firewall (O0 step 2) is unrunnable on the available
   admissible data** — the staged CSVs carry no assay/document IDs, so the +0.434 cannot be shown
   independent of same-campaign single-dose→dose-response protocol correlation (the exact confound that
   inflated CROSSDOC/PFSC signals earlier in this program).
2. **Licensing/leakage:** the only staged source is unlicensed and overlaps 26 of my sealed
   confirmation accessions; a certified run needs reconstruction from ChEMBL/PubChem with assay/document
   provenance and full firewalling against dev/confirmation.
3. **Transfer (A2/O1) is untested and faces the program's binding wall:** dual-cold deployment predicts
   an *unseen* target's reordering from protein features `p_t` alone, and four independent mechanisms
   plus TR-0/PFSC-0 have shown protein features do not carry transferable reordering. O0's per-target,
   feature-free identifiability does not establish that `f_θ(p_t,·)` generalizes across homology
   clusters.

I decline `AMOB_IDENTIFIABLE_AND_TRANSFERS` (transfer untested; no multi-seed; no firewall) and
`LOW_FIDELITY_SIGNAL_REAL_BUT_NO_HIGH_FIDELITY_TRANSFER` (transfer not tested-and-failed — asserting it
would overclaim a negative). The honest state is: **admissible new mechanism, strong O0 identifiability,
but the open data in hand is insufficient to complete the mandatory firewall and the transfer/power
gates.**

**Concrete path to resolve (highest-value next step).** Reconstruct single-dose % inhibition + exact
pKi/pKd from open ChEMBL/PubChem *with assay and document identifiers*, firewalled against the FORT
dev/confirmation splits; re-run O0 step 4 under **assay/document isolation** (the decisive test — if
+0.434 survives cross-document, AMOB is genuinely new identifying information; if it collapses like
CROSSDOC/PFSC, it was protocol correlation); only then run the §6–§7 O1 dual-cold transfer gate. Do not
claim dual-cold generalization, SOTA, or mechanism recovery before those gates pass.

`sealed_test_consumed=false`; `confirmation_labels_read=true` (pre-existing; AMOB read no FORT
confirmation labels — the staged diagnostic used external Harmonic-Discovery CSVs only).
