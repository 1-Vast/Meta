# Stage 0 audit of the retained no-warmup baseline

Date: 2026-08-15. Scope: `model/qpsmp_meta.py`, `model/bpsf.py`,
`model/encoders.py`, `scripts/qpsmp_data.py`, `scripts/train_qpsmp.py`,
`scripts/evaluate_qpsmp.py`, `scripts/evaluate_checkpoint_nested.py`.

Numerical authority for this audit:

```text
RECHECK_seed{20260812,20260813,20260814}.json   checkpoint re-evaluation
DIAG_seed20260812.json                           gradient/activation/parameter
REFERENCE_frozen_bank.json                       label-only reference predictors
THROUGHPUT.json                                  per-architecture cost
```

Nothing in this directory modifies a retained result. No historical
`RESULT.json` was edited.

## 1. Baseline reproduction

`scripts/evaluate_checkpoint_nested.py` reloaded each retained checkpoint and
reproduced the retained frozen-bank metrics exactly (`RECHECK_*.json` versus
`stageE_no_warmup_seed*_100step/RESULT.json`). Checkpoint save/load and the
fixed evaluation bank are sound.

Retained three-seed mean MSE (pK^2): k0 2.12, k1 1.70, k2 1.34, k3 1.27,
k5 1.25 — confirmed.

## 2. What the retained model actually is

Measured on `stageE_no_warmup_seed20260812`:

| quantity | value |
|---|---:|
| zero-shot spread across queries inside one episode | **0.065 pK** |
| zero-shot change when the protein is swapped for a cross-component donor | **0.0093 pK** |
| wrong-protein zero-shot MSE gap | 0.0119 |
| wrong-protein full-prediction gap, k=1/2/3/5 | 0.0055 / 0.0029 / 0.0023 / 0.0014 |
| primitive response, abs mean (cross-episode std) | 0.4913 (0.00029) |

The zero-shot endpoint is, to three decimal places, a constant. It is
insensitive to the ligand and insensitive to the protein. Every reported
few-shot gain is therefore the shrunken support mean around that constant.

### Label-only reference predictors on the same frozen bank

| predictor | k0 | k1 | k2 | k3 | k5 |
|---|---:|---:|---:|---:|---:|
| meta-train global mean (constant) | 2.063 | 2.063 | 2.063 | 2.063 | 2.063 |
| **retained model** | **2.047** | **1.657** | **1.306** | **1.243** | **1.217** |
| raw support mean | — | 2.361 | 1.873 | 1.545 | 1.288 |
| meta-train ligand-average prior | 1.591 | 1.591 | 1.591 | 1.591 | 1.591 |
| ligand prior + support offset | — | 2.684 | 1.937 | 1.742 | 1.275 |
| oracle target mean (level ceiling) | 1.096 | 1.096 | 1.096 | 1.096 | 1.096 |

The trained 3.79M-parameter interaction model beats a single global constant by
0.016 MSE at k=0, and loses to a trivial ligand-average lookup by 0.456 MSE.
Only 32% of frozen-bank query ligands even appear in `meta_train`, so the
ligand prior is a weak reference, not a strong one.

## 3. Module-by-module findings

Labels: **BUG** = confirmed bug, **ARCH** = architectural limitation,
**OPT** = optimization limitation, **DIAG** = missing diagnostic.

### 3.1 Protein encoder and residue refinement

* **BUG — dead trainable branch.** `ProteinEncoder.proj` and
  `ProteinEncoder.norm` (160,512 parameters) feed only `pooled`, which is
  consumed only by `QPSMPMetaLearner.protein_level`, whose scalar reaches only
  the `additive` diagnostic. `prediction` never sees it. Gradient is exactly
  zero at k=0,1,2,3,5.
* **ARCH.** The only live protein path is `bank_proj`, a single `Linear`
  applied to 128 pooled ESM slots. There is no residue refinement on the live
  path.
* **OPT.** Protein-side gradient norm is 0.16–0.90 against 2.0–10.8 for the
  ligand tower at every k.

### 3.2 Ligand encoder and graph masks

* No correctness defect found. Adjacency is derived from non-zero bond
  features, the zero-atom guard is present, masked max-pooling is guarded with
  `nan_to_num`, and padding is consistent through `compact_episode`.
* **OPT.** `ligand_only_mse_pk` = 2.082 is *worse* than the global constant
  (2.063). After 100 steps the ligand tower carries no affinity information.

### 3.3 Localization and padding

* **OPT + ARCH.** `residue_query`/`residue_key` gradient norms are
  0.0013–0.0077, three orders of magnitude below the ligand tower. The
  top-`k` gather passes gradient only through the 32 selected slots, so the
  pocket selector is effectively frozen at initialisation.
* Padding/trimming in `compact_episode` and `forward` is consistent; no bug.

### 3.4 Atom–residue pair field

* **ARCH.** `BipartitePairBlock` is the most expensive module in the model
  (peak 5.95–6.75 GB, 5.3–5.6 s per 4-episode step) and its endpoint varies by
  0.065 pK across queries. Cost is decoupled from output variation.
* `_LatentReadout` softmax-pools every A×R pair position into 24 latents and
  then averages the latents (`slots.mean(1)`), followed by a `LayerNorm` inside
  `cross_head`. Three successive averaging/normalising steps remove per-ligand
  contrast before the scalar readout.

### 3.5 Endpoint and slot pooling

* **BUG — dead computation.** `qpsmp_meta.py:413` binds `query_section` and
  never uses it. `SectionLatentEncoder.section` (4,608 parameters) is computed
  on every forward with zero gradient.
* **BUG — task-dead dictionary.** `_LatentReadout.response_weight`
  (2,304 parameters) receives gradient only from `support_match_loss`
  (weight 0.05), which is a pure self-consistency regulariser containing no
  task term. The resulting responses are a constant 0.4913 with cross-episode
  std 0.00029.

### 3.6 Support-zero and residual construction

* **ARCH.** In `QPSMPMetaLearner.infer`,
  `prediction = zero_shot + level_adjustment + (coefficients - level_adjustment)`.
  The two `level_adjustment` terms cancel identically, and `coefficients` is
  built from `transport_residual.detach()`. Consequently the few-shot loss
  delivers **no gradient at all** to the endpoint evaluated on support ligands.
  The `level_adjustment` output field looks gradient-carrying but is inert.
* **DIAG.** Nothing in the repository measured this cancellation.

### 3.7 Level / local / final output decomposition

* **ARCH (confirmed at gradient level).** At k=1 the kernel softmax is over one
  element, so `local == level`. `sar_cut_mse_pk == full_mse_pk ==
  level_only_mse_pk` in all three retained seeds, and
  `local_scale_logit`, `log_temperature` and `interaction_key.*` have exactly
  zero gradient at k=1.
* Learned scalars after training: temperature 5.003, `local_scale` 0.5008,
  shrinkage prior 1.9969 — all within 0.4% of their initial values. The kernel
  never moved.

### 3.8 Gradient coverage for every k

Dead (exactly zero-gradient) trainable tensors under the query loss:

| k | 0 | 1 | 2 | 3 | 5 |
|---|---:|---:|---:|---:|---:|
| dead tensors | 23 | 22 | 17 | 17 | 17 |

Dead at **every** k: `meta.ligand_baseline.*`, `meta.protein_level.*`,
`meta.term.state.weight`, `pair_section.latent.interaction.response_weight`,
`pair_section.latent.section.weight`, `protein_encoder.proj.*`,
`protein_encoder.norm.*` — **303,362 of 3,788,937 trainable parameters (8.0%)
are unreachable from the query loss.** The project's own admission gate
("no dead trainable branch") is violated by the retained baseline.

### 3.9 Checkpoint save/load and resume semantics

* Correct. `{model_state, config}` round-trips, `evaluate_checkpoint_nested.py`
  rebuilds strictly, and re-evaluation reproduced the retained metrics.
* **DIAG.** No optimizer/scheduler state is saved, so a long run cannot be
  resumed; only `--pretrained-checkpoint` warm-starting exists.

### 3.10 Controls

* Permuted-label gap is large only at k=1 (1.267), where it is produced by the
  magnitude-matched flip, not by ligand-specific structure. It collapses to
  0.021–0.078 at k=2/3/5.
* Matched wrong support (`foreign_code_state`) gaps of 1.5–3.2 are produced by
  a donor target's different absolute level, not by binding specificity.
* Wrong protein changes the prediction by 0.0014–0.0055 MSE. **The model is
  protein-blind.**

### 3.11 Parameters, activations, memory, utilization

| quantity | value |
|---|---:|
| trainable parameters | 3,788,937 |
| pair trunk share | 2,485,440 (66%) |
| peak CUDA memory (100-step runs) | 5.95 – 6.75 GB of 8.0 GB |
| seconds per 4-episode step | 5.3 – 5.6 |
| total training exposure | 400 episodes, 3,002 of 12,633 meta-train query cells |

**OPT.** The retained baseline performs 100 optimizer updates and sees roughly
one quarter of the training corpus once. That is not a converged model.

### 3.12 Reproducibility and evaluation power

* **DIAG.** Re-running the identical seed and configuration for 400 steps
  produced a different validation trajectory from step 40 onward (1.708 versus
  1.624). CUDA training here is not bitwise deterministic.
* **DIAG.** `eval_targets_per_component=1` yields **6 episodes per k** on the
  frozen bank. The same checkpoint scores zero-shot MSE 4.19 on randomly drawn
  meta-test episodes versus 2.05 on the frozen bank. All retained stage-to-stage
  deltas of order 0.03–0.05 MSE are inside this bank's noise.

## 4. Budget probe and wide-bank re-measurement

### 4.1 Budget probe (one changed variable: 100 -> 400 steps)

`stage0_budget_probe_seed20260812_400step/` repeats the retained configuration
with four times the steps (1,600 episodes, 6,888 of 12,633 meta-train query
cells). Frozen-bank test MSE:

| k | 0 | 1 | 2 | 3 | 5 |
|---|---:|---:|---:|---:|---:|
| retained 100 steps | 2.047 | 1.657 | 1.306 | 1.243 | 1.217 |
| 400 steps | 2.308 | 1.774 | 1.367 | 1.263 | 1.240 |

Validation admission score improved (1.595 against 1.624) while every test
metric worsened, and peak CUDA memory reached **8,399 MB on an 8,188 MB
device**. **More budget alone does not fix the retained recipe**, and
validation selection on a 6-episode bank does not transfer.

### 4.2 Trunk capacity probe (synthetic, controlled)

Both trunks trained for 1,200 steps on a noiseless protein-by-ligand
contact-type bilinear task (`TRUNK_CAPACITY*.json`). Held-out MSE relative to
label variance:

| learning rate | 3e-4 | 1e-3 | 3e-3 | 1e-2 | s / 1200 steps | parameters |
|---|---:|---:|---:|---:|---:|---:|
| retained pair trunk | 0.008 | **0.003** | 0.991 | 1.001 | 171 | 447,496 |
| interaction grammar | 0.014 | 0.017 | **0.009** | 0.400 | 39 | 183,005 |

This **falsifies a capacity explanation**. The retained trunk can represent a
protein-conditioned interaction; it collapses to a constant only above a
learning-rate threshold. The failure on the real corpus is an optimization and
budget failure, not an expressivity failure. The grammar trunk reaches the same
fit **4.4x faster per step with 2.4x fewer parameters** and stays stable over a
wider learning-rate band.

### 4.3 Wide-bank re-measurement (42 episodes, 7 components, all eligible targets)

The frozen 6-episode bank is a favourable subsample. On the wide bank:

| predictor | k0 | k1 | k2 | k3 | k5 |
|---|---:|---:|---:|---:|---:|
| meta-train global mean | 3.441 | 3.441 | 3.441 | 3.441 | 3.441 |
| meta-train ligand-average prior | 3.119 | 3.119 | 3.119 | 3.119 | 3.119 |
| raw support mean | — | 2.346 | 2.180 | 1.918 | 1.523 |
| **retained baseline, 100 steps** | **3.589** | **2.386** | **1.993** | **1.831** | **1.581** |
| retained baseline, 400 steps | 3.954 | 2.555 | 2.121 | 1.913 | 1.646 |
| oracle target mean (level ceiling) | 1.100 | 1.100 | 1.100 | 1.100 | 1.100 |

On the full cold-target meta-test set the retained model is **worse than a
constant at k=0** and **worse than the plain support mean at k=5**. It beats the
support mean only at k=2 and k=3, by 0.09-0.19 MSE. The retained frozen-bank
numbers (2.12 / 1.70 / 1.34 / 1.27 / 1.25) are real but are not representative
of cold-target performance.

## 5. Consequences for the next hypothesis

1. The zero-shot endpoint is a constant in practice, but **not** because the
   trunk cannot represent an interaction. The binding constraint is cost:
   5.3-5.6 s per 4-episode step and 6-8.4 GB peak on an 8 GB device cap the
   achievable number of updates, and the recipe degrades rather than improves
   when those updates are simply increased.
2. Adding another scalar output correction cannot help, because the quantity it
   would correct is a constant.
3. Any new mechanism must keep every trainable tensor on the path to
   `prediction`; the retained baseline does not (8.0% are unreachable).
4. k=1 must gain a query-specific channel; it is identically scalar today.
5. Decisions must not be taken from the 6-episode frozen bank alone. A wide
   bank over all eligible meta-test targets is reported alongside it from here
   on, with the frozen bank retained as the protocol comparator.
6. Selection must not be taken from a 6-episode validation bank either; the
   400-step probe improved validation while worsening every test metric.
