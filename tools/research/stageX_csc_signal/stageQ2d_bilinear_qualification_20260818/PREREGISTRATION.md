# Stage Q2d preregistration — isomorphic bilinear qualification (2026-08-18)

Frozen BEFORE any Q2d computation. SHA-256 recorded in PREREGISTRATION_SHA256.txt.
This stage is synthetic-only; it reads no external biological labels.
Baseline state carried from Q2c (unchanged): pipeline qualification FAILED at
Q2 (X0c); biological conclusion UNRESOLVED; B1/B2/C/D NOT AUTHORIZED.

## Mechanism claim under test (bounded)

A transferable statistical interaction field: pocket-residue
physicochemical/positional properties x ligand substructure, learned as a
gradient-trained low-rank bilinear map. No claim of atomic contacts,
conformational change or binding free energy is made from these data.

## Logical anchor

The Q2 planted interaction is generated as a low-rank bilinear of the
KLIFS-aligned pocket one-hot and ECFP4: I(p,l) = (P@U)(L@V)^T, double-centred,
standardized to sd = tau*. Therefore the one-hot pocket is part of the truth
generation, not an information-poor stand-in. If a gradient-trained model of
the SAME functional form cannot recover I, the suspects are the learner,
loss routing, endpoint transform, split power or optimization - not the
one-hot representation. Q2d-1 isolates these steps one at a time.

## Q2d-1 isomorphic bilinear positive control

Model (gradient-trained; no ridge, no closed form, no test-time gradients):
  p = aligned pocket one-hot (1700) [or z-scales in Q2d-2]
  l = ECFP4 (2048)
  I_hat = inter_scale * ((p A) . (l B) + inter_bias),  rank 4
  y_hat = mu + p_b(row) + l_b(lig) + I_hat

Arms: exact_bilinear; additive_only (inter parameters absent); ligand_only
(p=0); shuffled_protein (row-permuted p); random_protein (Gaussian p);
oracle_latent (p = P@U, 4-dim); no_interaction_head (full model, inter_scale
frozen at 0 with requires_grad=False - the Q2c-0 drift bug is not repeated).

Phases (each on truth seeds 0/1/2, median reported):
  A  z-scale identity link, no sigmoid, no censoring, full panel;
  B  + sigmoid endpoint (determinate logit), no censoring;
  C  + panel missingness (70% observed cells, frozen MCAR seed);
  D  + interval censoring at 0/100 (floor semantics);
  E  + main-effect competition (pm/lm learned through linear encoders shared
     with the interaction projections A/B; no nonlinearity).

Frozen gate (identical thresholds, never moved retroactively): exact_bilinear
median over seeds must reach Spearman >= 0.30 AND dead-zone sign accuracy
>= 0.70 AND gap vs ligand_only >= 0.05. Every negative arm must FAIL the
same gate. Protocol: exact_bilinear 8 restarts with validation-loss selection;
negative arms 1 restart x 6000 steps; same splits/graph/noise as Q2.

Decision rule: if exact_bilinear passes phase A, the ladder attributes the
Q2 failure to the later step(s) where the pass disappears; if it fails phase
A, the harness/optimization is unqualified and Q2d-2 is not started.

## Q2d-2 representation replacement (only after Q2d-1 passes)

Candidates: aligned pocket one-hot; pocket amino-acid z-scales;
pair-centered local ESM; KLIFS-aligned per-position ESM (kept as
[position, residue embedding], NOT pooled into one vector, interacting with
ligand fragments through the shared low-rank head); ESM + z-scales; global
pooled ESM; shuffled/local random. The Q1 pair_centered_local_esm pass only
certifies pocket-membership readability; ESM must re-pass this gate.

## Q2d-3 matched-delta training (only after Q2d-2)

Add training contrasts: endpoint loss + within-protein ligand contrast
d(p,l1,l2) = y(p,l1) - y(p,l2) + WT-mutant same-ligand contrast
d(p_wt,p_mut,l) = y(p_mut,l) - y(p_wt,l), plus interaction orthogonality and
modality gradient coverage. Contrasts remove protein/ligand main effects in
the loss itself rather than by post-hoc projection.

## B1 execution rules (authorized only after Q2d passes)

Primary panel: Saifudeen 2026 SAME-STUDY (409 WT, 349 variants/fusions,
92 inhibitors, 1 uM, Km ATP, duplicates). Duong-Ly is B1-R cross-study
replication with explicit study/batch covariates (its WT values come from
Anastassiadis and carry cross-study batch). B1 rules: single mutants only;
exact matched WT required; same or compatible construct/substrate;
same-ligand WT->mutant differences; parallel estimands (responsive-window,
sign-only, censoring-aware); held-out parent and pocket component;
protein/component cluster bootstrap; known resistance pairs only as a
preregistered positive-control subset; primary result over ALL legal pairs.

## Literature corrections recorded (from independent review)

- Merget et al. 2017 (J. Med. Chem.) trains ligand-side random-forest models
  per kinase; it is NOT a 234x380 Kronecker PCM and does not demonstrate
  unseen-kinase transfer. Withdrawn as PCM evidence.
- PCM evidence base: 317 WT/mutant kinases x 38 inhibitors x 12,046 Kd with
  protein-ligand cross descriptors (PMC2910025); HIV variant PCM with
  mutation-site physicochemical x drug descriptors (PMC3002298, external
  antivirogram validation PMC3578754); active-site sequence DTA on BindingDB
  (PMC9516689).
- The claim that the literature uniformly uses delta-delta targets is
  WITHDRAWN: many PCM works predict absolute values; delta-delta is chosen
  here because it structurally zeroes ligand/protein main effects.
- MdrDB is an external stress test only, not a clean cold-protein benchmark
  (GDSC/DepMap mixing, structure-derived entries); audit before any use.

## Governance

SHA-256 seeds; no Python hash(); artifacts carry schema / prereg SHA / input
SHA; commands.jsonl appended; restricted data never committed; stageX0c and
stageQ2c artifacts are read-only inputs; separate commits for implementation
and documentation.
