# UBSE-A1 typed-event OT successor design review

**Binding correction (2026-07-30):** this v1 design must not be used directly
as an A1 preregistration. Unbalanced OT can improve by changing
pair-conditioned row/column marginals without residue-FG coupling; the v1
rectangle lacks a functional-group axis; its dustbin is underdefined; and the
previously read G1 audit is not a fresh confirmation set. See
`ubse_a1_coupling_identifiability_correction_2026-07-30.md`. Current status:
`REVISE_UBSE_A1_BEFORE_PREREGISTRATION`.

Date: 2026-07-29  
Status: design review only; not an A1 preregistration or execution
authorization  
Prerequisites: corrected A0C coordinate gate and P0A proposal gate must both
pass

## Decision

The smallest successor that adds information not present in stopped UBSE-G1
has two non-skippable stages:

\[
\mathrm{A1\!-\!R:\ independent\ typed\ event\ reliability}
\quad\longrightarrow\quad
\mathrm{A1\!-\!S:\ deployable\ monomer+ligand\ OT\ student}.
\]

Neither stage reads affinity, confirmation, or sealed outcomes. P0A, if it
passes, is only a frozen target-marginal proposal model.

## Novelty boundary

[LINKER](https://pubs.acs.org/doi/10.1021/acs.jcim.6c00527) already derives
residue-by-functional-group-by-seven-type supervision from experimental
complexes and predicts that tensor from protein sequence and ligand SMILES.
Therefore typed maps, functional-group abstraction, PLIP supervision, and
structure-free interaction-map inference are not novel.

The defensible combined contribution is narrower:

> independently reliability-certified and source-closed 3D event
> supervision, distilled into a deployment-matched predicted-monomer plus
> 2D-ligand unbalanced-OT student, with a typed pair-burden exact null,
> fixed-training-measure purification, and placement-identifiable rectangle
> certificates under strict dual-cold closure.

The teacher supplies residue-functional-group correspondence, event type and
geometry, and within-target ligand-dependent spatial redistribution. It does
not supply affinity, binding free energy, or a causal ligand-induced
transition.

## A1-R: independent event reliability

The 3,467-row A0 main manifest contains only four cross-PubMed/PDB repeats of
an exact `(target_key, conn)` unit. It cannot certify its own event extractor.

Build A1-R from the P0A-legal 62,849-row fit complement, then remove every A0
target, PDB, PubMed, and scaffold. Label-blind accounting leaves:

- 719 cross-PubMed/PDB repeat units over 611 targets;
- 458 units over 350 targets with a legal same-target wrong ligand.

Before coordinates are read, select 160 distinct targets by seed-1729
SHA-256 order. For each target:

1. select one exact `(target_key, conn)` unit;
2. use two different PubMed/PDB structures as the correct repeat;
3. choose a same-target, different-connectivity, different-PubMed/PDB wrong
   ligand;
4. prefer same scaffold, then minimum heavy-atom difference, then frozen
   lexical order;
5. do not use binding-residue or event labels in selection.

A1-R must remain disjoint from A0 in exact target, PDB, PubMed, and scaffold,
and its events never enter A1-S fit.

### Locator and extractor contract

Use the corrected A0C fields:

- BioLiP column 7 only as `biolip_filename_serial`;
- BioLiP column 20 as `mmcif_auth_seq_id`.

Never infer one from the other. Bind sequence, ligand heavy-atom graph,
chain, CCD, model, altloc, occupancy, water, metal, covalent-ligand, missing
residue, and symmetry rules before event extraction.

The current official PLIP release is
[3.0.1](https://github.com/pharmai/plip/releases); a future preregistration
must freeze the exact package/commit/dependency hashes. Extract seven
LINKER-comparable channels:

- hydrogen bond;
- hydrophobic;
- pi stacking;
- pi cation;
- salt bridge;
- water bridge;
- halogen bond.

Metal and covalent events are recorded but excluded from v1. Ligand
functional-group instances that are graph-automorphism equivalent must be
collapsed into a symmetry class.

### A1-R hard gates

All must pass:

- at least 128/160 complete target units with two correct repeats and one
  legal wrong ligand;
- locator, graph, and sequence unique-resolution rate at least 90%;
- median correct residue-by-type Jaccard at least 0.50;
- median correct residue-by-FG-by-type Jaccard at least 0.40;
- correct-minus-wrong residue-by-type Jaccard at least 0.10;
- 2,000-target-bootstrap lower 95% bound for that margin above 0.05;
- at least three event types each covering 40 targets;
- for every adequately supported type, repeat Jaccard at least 0.30 and
  wrong-ligand margin lower bound above zero;
- zero A0 target/PDB/PubMed/scaffold overlap and zero forbidden outcome
  access.

Failure returns
`STOP_UBSE_A1R_TYPED_EVENT_TEACHER_UNRELIABLE`. Types or thresholds may not
be removed after seeing the result.

## A1-S input and source closure

Only after A1-R passes, extract A0 fit and validation events. Freeze the
extractor, model, controls, epochs, and gates before validation extraction;
extract audit events once only after validation passes.

The student may receive:

- exact target sequence;
- frozen three-seed P0A ensemble logits;
- sequence-matched predicted-monomer residue coordinates and confidence;
- ligand 2D graph and deterministic functional-group instances.

It may not receive holo coordinates, ligand pose, binding-residue labels,
source IDs, observed event counts, other test ligands, or affinity.
Predicted monomers require a frozen source/version/hash, sequence identity at
least 0.95, and target coverage at least 0.90. A stripped holo receptor is
not a deployment substitute.

Use the P0A ensemble only to propose
\(K_t=\min(256,L_t)\) residues. Teacher-positive proposal recall must be at
least 0.85 on fit/validation and 0.80 on the one-time audit.

## Exact null and residual model

Use frozen protein features plus two invariant graph layers on the predicted
monomer, and a two-layer ligand MPNN pooled to functional groups. Keep the
learned model below three million parameters.

For event type \(k\), learn deployment-side marginals:

\[
m_{tik}=\operatorname{softmax}_i f_k(h_{ti}),\qquad
q_{lgk}=\operatorname{softmax}_g g_k(z_{lg}),
\]

and a position-free typed pair burden:

\[
c_{tlk}=\operatorname{softplus}b_k(\operatorname{pool}h_t,
\operatorname{pool}z_l).
\]

The exact null is:

\[
\lambda^{(0)}_{tligk}=c_{tlk}m_{tik}q_{lgk}.
\]

No held holo event count may replace \(c_{tlk}\).

The residual uses rank-32 compatibility and dustbin-enabled unbalanced
entropic OT:

\[
s_{igk}=\langle U_kh_i,V_kz_g\rangle/\sqrt{32},
\]

\[
\pi_k=\arg\max_{\pi\ge0}
\langle\pi,s_k\rangle
-\epsilon KL(\pi\Vert m_k\otimes q_k)
-\tau_r KL(\pi\mathbf1\Vert m_k)
-\tau_c KL(\pi^\top\mathbf1\Vert q_k).
\]

After normalizing real-real mass:

\[
\lambda^{(1)}_{tligk}
=c_{tlk}\pi_{tligk}/\sum_{ig}\pi_{tligk}.
\]

Zero compatibility recovers the exact null. The residual can only move typed
mass across residue-FG locations; it cannot win by changing event count.

Train the null for 20 fixed epochs, freeze it, then train only OT
compatibility and dustbin parameters for 30 fixed epochs. P0A never receives
A1 gradients.

## Fixed-measure purification and identifiable topology

Purification is a fit-only regularizer/evaluation operator, never a test-batch
forward operation. Freeze weights \(W\) before labels and penalize the
deviation between the log-rate residual and its weighted two-way purified
projection only on fit panels.

A placement-identifiable event rectangle requires two ligands and two
residues with an observed checkerboard:

\[
Y_{l_1ik}=1,\;Y_{l_1jk}=0,\quad
Y_{l_2ik}=0,\;Y_{l_2jk}=1.
\]

Score:

\[
\Delta=
[\log\Lambda_{l_1i}-\log\Lambda_{l_1j}]
-[\log\Lambda_{l_2i}-\log\Lambda_{l_2j}].
\]

The typed burden and target marginal cancel under the exact null. Require at
least 500 fit, 32 validation, and conditionally 40 audit
placement-identifiable panels spanning at least three event types before GPU
training.

## Required controls

- exact \(m\otimes q\otimes c\) null;
- same-panel wrong ligand;
- matched wrong protein;
- residue-feature/monomer-position permutation;
- FG identity erase/derangement with FG count preserved;
- cyclic event-channel permutation;
- P0A/target-marginal only;
- equal-budget margin-preserving 2-by-2 event switches;
- parameter-matched LINKER-form sequence-plus-SMILES baseline trained on the
  same closed fit source.

## Validation and audit gates

Use the same gates first on validation and then once on audit:

- group-resolved macro AP at least 0.15 and at least twice prevalence;
- full-minus-exact-null AP at least 0.03, bootstrap lower bound above zero,
  all three seeds positive;
- rectangle directional accuracy at least 0.60, lower bound above 0.55;
- purified log-rate cosine at least 0.10;
- count-normalized assignment Recall@1 at least 0.65 and at least 0.10 above
  every control;
- wrong-ligand, wrong-protein, and structure-position decrements each at
  least 0.05 with positive lower bounds;
- FG and event shuffle AP decrements each at least 0.03 with positive lower
  bounds;
- beat the margin-preserving shuffled-label arm by AP 0.03 and directional
  0.05;
- three-seed AP range at most 0.05 and all values finite;
- beat the parameter-matched LINKER-form baseline by at least 0.03 in AP and
  directional accuracy, with positive bootstrap lower bounds;
- CUDA execution, no CPU model offload, and no forbidden access.

If typed signal passes but not the LINKER-form comparison, report
`A1_SIGNAL_REAL_BUT_MODEL_NOVELTY_NOT_ESTABLISHED`.

Full success returns only
`FREEZE_UBSE_A1_TYPED_EVENT_STUDENT_FOR_DOMAIN_SUPPORT_ONLY`. It requests a
later affinity-blind target-domain support gate; it does not unlock affinity,
Stage-2, confirmation, or sealed access.

## RTX 4060 8-GB execution budget

- top-256 proposal;
- FP16 GNN/heads, float32 Sinkhorn;
- at most 16 complexes per microbatch;
- \(\sum K_tG_l\le65,536\);
- at most 1,024 ligand atoms per microbatch;
- gradient accumulation 4;
- cached ESM/P0A/monomer/ligand features;
- median GPU utilization target at least 75%;
- peak framebuffer below 7.2 GB;
- telemetry for utilization, memory, power, and temperature; no overclock.

Once external coordinates, monomers, and fixed extractor dependencies are
available, the end-to-end estimate is 8--24 hours. A0C or P0A not passing
returns `WAIT_UBSE_A1_PREREQUISITES`.
