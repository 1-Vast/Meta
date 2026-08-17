# Stage E preregistration — panel-set level head + orthogonal level/shape training

Frozen before any Stage E arm trained. Date: 2026-08-17. Development evidence
only; meta_test is sealed and never constructed in this stage.

## The two preregistered innovations (maximum two)

- **I1 (framework): panel-set level readout.** The additive zero-shot protein
  level (protein_head(P), the incumbent's only level channel) is replaced by
  a level head that consumes the protein summary plus order-invariant
  mean/max pooling over the *query ligand encodings*. Query ligands are legal
  model inputs at every k, including k=0. Rationale (measured, not assumed):
  D0_LEVEL_IDENTIFIABILITY shows panel-composition features carry level signal
  (MLP level MSE 1.887 vs 2.155 for the meta_train constant; shuffled-panel
  control 5.075), and D0_LEVEL_ANATOMY shows panel composition transfers 23.9%
  of meta_train between-target level variance across components versus 11.9%
  for the protein sequence embedding and -1.1% for component identity.
- **I2 (training): orthogonal level/shape routing.** The level head is trained
  only by the per-episode level term (panel_level - mean(query_y))^2; the
  interaction path (interaction head, contact dictionary, ligand baseline) is
  trained only by the centered/ranking terms. The centered term's gradient
  into a per-episode constant is identically zero, so the routing is enforced
  by construction, not by detached tricks. The incumbent's smooth_l1 on the
  full post-adaptation prediction — a 68% level-dominated objective at k=0 —
  is dropped for the candidate.

No closed-form solver, no ridge/pseudoinverse, no inner loop, no query-label
adaptation, no cross-dataset information. Single-stage end-to-end training.

## Arms (identical seed/budget/partition/optimizer; one code path)

| arm | model | loss |
|---|---|---|
| T2 | incumbent similarity_only trunk | Stage B recipe: smooth_l1(post) + 1.0·smooth_l1(pre) + 0.5·ranknet + 0.5·centered + 0.05·dictionary |
| T2-LEVEL | T2 model | recipe + 1.0·level term routed to the episode mean of the zero-shot endpoint (loss-only ablation) |
| LSP | PanelLevelShapeModel | 1.0·smooth_l1(pre) + 0.5·ranknet + 0.5·centered + 0.05·dictionary + 1.0·level term on the panel level head (routed) |
| LSP-NOROUTE | PanelLevelShapeModel | Stage B recipe (framework-only ablation, no routing) |

Budget: 1,200 optimization steps, 3 episodes/step, seed 20260815, AdamW
(lr 3e-4 backbone 0.25×, transport 3e-4), grad clip 1.0, amp off, float32.
Leak-free protocol (Stage B): meta_train components partitioned once
(seed 20260818) into 227 fit / 31 internal-validation; training episodes from
fit components only; checkpoint selection on internal-validation MSE only;
meta_val read exactly once after freezing. Uniform component→target sampling.

## Preregistered gates (single-seed screen, on the frozen meta_val banks)

G1. LSP beats T2 on k=0 MSE (restored pK²) with a resolved paired component
    bootstrap interval (lower bound > 0), and no k degrades with a resolved
    interval.
G2. Ranking never trades: for every k, LSP's Spearman, Pearson and CI are not
    lower than T2's by a resolved interval; activity-cliff sign does not
    degrade by more than 0.03.
G3. Correct-support dependence: at k>=1, permuted-support and matched-wrong
    MSE for LSP is higher than correct by a resolved interval (label-bound
    mechanism), and the incremental dependence is not less than T2's by more
    than the preregistered slack.
G4. Wrong-protein perturbation degrades LSP by at least as much as T2 (no
    protein-blindness regression), within preregistered slack.
G5. Level attribution: LSP's k=0 level² is lower than T2's by a resolved
    interval AND the panel level head's out-of-sample level MSE on the frozen
    meta_val bank is below the Stage C ESM-MLP record (1.6357) — otherwise the
    mechanism is declared to carry no new level information.
G6. Cost: trainable parameters ≤ 2.0× T2; peak VRAM and wall time within 1.5×
    T2; gradient coverage recorded for every parameter group.

Stop rules (any one fires → no multi-seed, no meta_test):
S1. G1 fails (no resolved k=0 gain, or any k degrades resolved).
S2. G2 fails (resolved ranking regression at any k).
S3. G5 fails AND the panel level head's level MSE exceeds the meta_train-only
    constant — the framework innovation is then declared inert, not rescued.
S4. Any control inverts (wrong support improves MSE with resolved interval).

Promotion path: all gates pass on the single seed → ≥3 fixed seeds (nested k,
component bootstrap) → freeze architecture/hyperparameters/checkpoints →
meta_test opened exactly once with written authorization.

## Required controls and ablations (this stage's scope)

- leak-free uniform-training baseline: T2.
- level-only / shape-only / zero-shot / few-shot decomposition: reported per k
  per arm (level², centered MSE, zero_shot condition, no_adaptation).
- support-label permutation + matched-wrong support: k>=1 conditions.
- wrong-protein: full-system perturbation per episode.
- exact-ligand / scaffold novelty / low recall strata: the double-cold split
  makes all meta_val ligands, scaffolds and documents novel by construction;
  novelty is quantified per episode (max/mean Tanimoto to meta_train).
- assay/target/panel covariate shuffle: D0 shuffled-panel control is the
  no-training evidence; the trained candidate's level head is additionally
  audited on a panel-covariate-shuffled forward (features permuted across
  episodes) in the single-seed report.
- gradient-routing ablation: LSP-NOROUTE; loss-only ablation: T2-LEVEL.
- parameter count / peak VRAM / training time / gradient coverage: recorded in
  each arm's RESULT.json and GPU_PROBE.json.

## Authorities and hygiene

- D0 evidence: D0_AUDIT_DECOMPOSITION.json, D0_LEVEL_IDENTIFIABILITY.json,
  D0_LEVEL_ANATOMY.json, D0_OCCUPANCY_STRATA.json (this directory).
- The five re-examination questions from the governing task are answered in
  D0_REPORT.md.
- meta_val figures remain development evidence (single seed); they must not be
  quoted as independent confirmation.