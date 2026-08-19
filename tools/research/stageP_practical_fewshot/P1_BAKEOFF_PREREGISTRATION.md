# Stage P1 bake-off preregistration — practical few-shot baselines (2026-08-19)

Frozen BEFORE any bake-off computation. Scope: BindingDB-Ki only (the
admitted governed corpus). Davis/KIBA/PKIS/Saifudeen follow the same
protocol in their own stages with dataset-specific heads. This document
does not authorize changes to model/ or scripts/; all code lives under
tools/research/stageP_practical_fewshot/.

## 1. Corpus and split

- Corpus: dataset/processed/meta_fewshot/bindingdb_ki_main_v0 (read-only,
  manifest content SHA frozen in the first bake-off artifact). pKi =
  9 - log10(Ki[nM]); exact positive uncensored Ki only; panel-aggregated
  (see corpus manifest).
- New P-line split (built by scripts in this stage, frozen before first
  label use): cdhit40 clusters (from the corpus's own cdhit40 artifact)
  assigned to p_train / p_val / p_test with 60/20/20 split of TARGETS,
  cluster-balanced by target count, stable SHA-256 seed; whole clusters
  move together (protein-component cold). Ligands are NEVER split by
  scaffold or series in P1/P2 (same-series support/query allowed by
  design). Assert: target sets disjoint; split recorded per cell id.
- p_val = model/hyperparameter/checkpoint selection surface (cold targets);
  p_test = final reporting surface; p_test labels are never used for
  selection. k-eligibility: a target is eligible at k iff it has >= k+query
  unique ligands; reported per k with the eligible-target census in the
  artifact.

## 2. Frozen episode bank

- One bank, built once, SHA-pinned, shared by ALL arms:
  - For each eval target (p_val, p_test) and draw d in 0..D-1: ligand-unique
    rng ordering (rng keyed by stage|split|target|draw, never by arm);
    support(k) = first k ligands (nested: support(k1) subset of support(k2)
    for k1<k2, k in {0,1,2,3,5,10,20,40}); query = next Q ligands after
    max_k (Q frozen per layer; P1 Q=8 when available, else all remaining).
  - k=0 episodes: same rows, support empty, same query; the model code path
    for k=0 must not change because a few-shot module exists — the few-shot
    module is bypassed with zero support (frozen assertion test).
  - D (draws), Q, and the max_k list frozen in the artifact manifest.
- Each record: split, component (cdhit40 cluster), target, support cell ids
  (per k), query cell ids, draw, donor target (for foreign-protein
  controls).

## 3. Arms (bake-off list, all measured under the SAME bank/split/budget)

1. ligand_only: per-target support-mean of query-independent labels;
   k=0 -> training-set global mean (and per-family variants reported).
2. fixed Tanimoto: Morgan ECFP4 2048-bit (RDKit, frozen radius/config)
   Tanimoto-weighted support label average (kNN with k_sup neighbours);
   k=0 -> global mean. Closed-form; admissible as BASELINE only, never as
   the final method.
3. ordinary fine-tuning: shared pretrained-by-episodic-training backbone,
   support-only SGD adaptation at test time; checkpoint = best support
   loss (query labels never enter adaptation or selection).
4. first-order MAML: inner loop on support, outer meta-objective on query
   during TRAINING split only; at eval, same support adaptation as (3).
5. CNP: support encoder -> latent -> query decoder (per-target adaptation
   through the latent only).
6. FS-CAP-style: ligand-only support encoder (no protein features),
   matching the published comparison against Tanimoto/kNN.
7. ActFound-style: within-task pairwise (support pairs, difference
   supervision) during training; eval = adapted level prediction.
8. AdaMBind-style: task sampling + loss per the published paper; the exact
   spec (sampling schedule, gradient-consistency term, hypernetwork
   topology) is appended as a SHA-frozen addendum after FULL-text
   inspection and BEFORE its first run. Not run before the addendum.
9. current admitted baseline: QPSMP/BPSF with its existing governed
   protocol and checkpoints (read-only reference); if its split differs it
   is reported as parameter/split delta, never silently matched.

Backbone for arms 3-8: one frozen interaction trunk (protein feature =
corpus protein-bank feature; ligand feature = ECFP4-2048; bilinear
low-rank interaction + additive heads; details frozen in the
implementation artifact). Any arm-specific extra module must be recorded
as a parameter-count delta and ablated. Ridge / pseudoinverse / closed-form
regression is FORBIDDEN as any arm's final method (Tanimoto arm is an
explicitly allowed baseline exception).

## 4. Budgets, seeds, checkpoint rule

- Single-seed screening on p_val for every arm; promoted arms rerun with
  >= 3 fixed SHA-256 seeds on p_val+p_test (seeds 20260861/2/3 frozen).
- Training budget per arm: identical wall-clock-equivalent step x batch
  budget (frozen constant per arm family, recorded in the artifact);
  identical optimizer (AdamW, lr/wd frozen) unless the arm's published
  method mandates a different optimizer (recorded as a delta).
- Checkpoint rule for training: best p_val aggregated metric (frozen
  metric = mean over k of MSE); for test-time adaptation: best SUPPORT
  loss. Query labels never enter adaptation or selection (asserted by
  tests).
- RNG: SHA-256 keyed streams, shared minibatch order across arms
  (keyed by seed|phase|step, never by arm). Python hash() banned.

## 5. Metrics and stratification (P1)

- Primary: MSE, RMSE, CI (concordance), Spearman, Pearson, centered MSE
  per k in {5,10,20,40}; stress k in {1,2,3}; k=0 reported for continuity.
- Stratification (computed once from the bank, label-free where possible):
  support-query mean Tanimoto bands (low/mid/high), Bemis-Murcko scaffold
  novelty (query scaffold unseen in training ligands / seen), activity-cliff
  cells (|delta pKi| >= 2.0 to nearest support neighbour, frozen threshold).
- Uncertainty: paired target-level (component) bootstrap, cluster-level
  intervals; comparisons are paired per target per draw.

## 6. P1 promotion gates (never moved)

A candidate arm promotes iff on p_test:
- paired improvement over the matched ligand_only AND Tanimoto baselines
  (MSE/RMSE improved with component-bootstrap interval not covering 0);
- CI and Spearman not significantly degraded (paired interval excludes
  degradation beyond frozen tolerance);
- improvement is not explainable by exact-ligand recall (ablation: remove
  query ligands whose exact InChIKey appears in support/training from the
  improvement computation; the residual must remain);
- direction-consistent across >= 3 fixed seeds;
- component-level intervals support the main claim.
AdaMBind/MAML adoption is decided by these measurements only.

## 7. P2 screening addendum (same bank, k in {0,5})

- Regression metrics always; screening metrics (EF 1/5/10%, BEDROC a=20,
  PR-AUC) reported only when a frozen active/inactive threshold label
  exists for the corpus (pKi >= 6.0 = active, frozen; sensitivity
  reported).
- Ligand-novelty stratification: query scaffolds unseen in p_train.

## 8. M1 mechanism layer (unchanged stress protocol)

- Reuses the frozen double-cold artifact
  (dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1) read-only
  and the existing counterfactual suite (correct/ligand-only/shuffled/
  family-preserving/matched-wrong/residue-permuted/random/no-interaction,
  component bootstrap), k in {0,1,2,3,5}. M1 decides ONLY mechanism
  claims; its failure never blocks P1/P2 performance conclusions.
