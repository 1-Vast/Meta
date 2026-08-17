# Stage V preregistration — core-inclusive MMP × local protein-region interaction

Frozen **before any Stage V census statistic, any interaction-variance statistic
and any trained-arm evaluation metric was read.**

## 0. Provenance of this document, stated plainly

This stage is an **adaptive correction informed by Stage T**, not an independent
confirmation. Stage T's full gate result was read at 17:08 before Stage U was
frozen at 17:12, and Stage V supersedes Stage U (see
`STAGE_U_GOVERNANCE_AUDIT.md`). No claim made here may be described as an
independent replication of anything.

**Every numeric threshold in sections 3, 4 and 6 is inherited verbatim from
Stage U's `PREREGISTRATION.md` (SHA-256 `fdc0a830…a645c11`), which was frozen
before any core-inclusive census number existed.** That is what keeps them
non-retrofitted. Stage V adds controls and stricter evaluability conditions; it
**loosens nothing**.

### Knowledge state at freezing time

Stage 0 forensics (`STAGE0_FORENSICS.json`) were mandated before this document
and are therefore known. They are recorded here so no later reader has to guess
what was known:

* within-target across-core `delta_y` gap: median 0.269 pK, p95 1.268, max 3.401;
* 40.4% of Stage T fit `D` rows and 28.9% of internal rows had disjoint core sets;
* under the core-inclusive key: fit 12,740 `D` rows / 99 components,
  internal 546 `D` rows / 10 components, of which **32 rows over 4 components**
  have a key repeated in fit and **514 rows over 7 components** do not;
* fit keys with >=3 targets and >=3 components: 1,001.

Because these bank sizes were known, **no threshold in this document was chosen
with reference to them.** All are inherited from Stage U.

## 1. Estimand

    tau = (shared core, R_a -> R_b, attachment environment, stereochemistry,
           formal-charge change)
    delta_y(t,tau) = mu_tau + delta(t,tau) + noise
    D(tau,t1,t2)   = delta_y(t1,tau) - delta_y(t2,tau)

`D` cancels the target affinity level **and** `mu_tau` exactly — and, unlike
Stage T, that statement is now true, because `tau` carries the complete chemical
context and two targets are compared **only** when they realise the identical
`tau`.

The question: can protein sequence-derived **region** tokens predict `D` on
unseen protein components, better than zero, better than a global protein
summary, better than a shuffled protein, better than a similarity-matched wrong
protein, better than permuted residue tokens, and better than a capacity-matched
random protein representation?

## 2. Governance

Identical to Stage U §1 and inherited: governed BindingDB-Ki `main_v0`, double-
cold `v1`, **physically isolated split view**; the development-validation split
used for nothing and enforced by a parsed-AST test; frozen
`scripts/internal_validation.py` partition (227 fit / 31 internal); SHA-256
seeds only, Python `hash()` forbidden and AST-tested; labels never touch feature
construction, key definition, split construction or hyperparameter selection.

Supervision reliability is inherited from Stage T's T0
(`tools/research/stageT_mmp/T0_RELIABILITY.json`), which is unaffected by the
key defect: `sigma2_same = 2 x within-assay variance = 0.858 pK^2`,
`sigma2_cross = 1.221 pK^2`, with L2 recorded non-identifiable. It is a
supervision-reliability estimate on a small selected subset and **never** an MSE
floor.

## 3. V0 — core-inclusive census

MMP construction is Stage U's (`rdMMPA.FragmentMol`, Hussain–Rea, single cut,
core = larger fragment, isomeric SMILES throughout, canonical direction by
R-group SMILES sort, deduplication by lower cell index).

    exact key  = sha256( core_isomeric | repr(attachment environment) |
                         R_a_isomeric >> R_b_isomeric )
    coarse key = sha256( core_stereo_stripped | element | aromatic |
                         R_a_stereo_stripped >> R_b_stereo_stripped )

Attachment environment = (element, aromatic, in_ring, degree, formal_charge,
hybridization) of the core atom bearing the cut.

Primary bank: same target, Ki endpoint, and the two cells share an identical
governed `panel_id`. Strata S1/S2 (same panel) and S3 (cross panel, weak, never
pooled) as in Stage T.

### Frozen V0 admission gate — inherited verbatim from Stage U §2.5

1. same-panel fit observations >= **2,000**;
2. fit targets >= **50**;
3. at least **30** exact keys each spanning >= **3 targets** and >= **3 components**;
4. internal same-panel observations >= **300**;
5. internal components >= **10**;
6. no domination: top-1 exact-key observation share <= **0.05**; top-10 <= **0.20**;
   top-1 target share <= **0.25**; top-5 <= **0.75**; top-1 component share
   <= **0.25**; top-5 <= **0.75**.

If any fails: stop, negative report, no neural model.

### Frozen V0b evaluability gate — inherited from Stage U §4.3/§4.6 gate 10

Stage U froze the rule that a `D` evaluation surface with **< 100 rows** is
recorded `not_evaluable` and the route cannot pass on it. That rule is applied
here to **every** evaluation surface, not only the disjoint one:

* `internal_repeated` — internal `D` rows whose exact key also occurs in fit;
* `internal_disjoint` — internal `D` rows whose **coarse** key is absent from fit;
* `internal_all` — all internal same-panel `D` rows.

A surface with < 100 rows is `not_evaluable` and no gate may be declared passed
on it.

## 4. V1 — interaction variance, inherited verbatim from Stage U §3

Fit components only; internal reported as descriptive consistency and never
gating.

1. aggregate same-panel `delta_y` per `(exact key, target)` by median;
2. restrict to exact keys with >= 2 target effects;
3. per key `SS_tau = sum_t (dy_t - mean_tau)^2`, `df = k_tau - 1`; pool
   `MS_effect = sum SS_tau / sum df`;
4. `sigma2_noise = sigma2_same` from T0 for one observation;
5. `theta = MS_effect - sigma2_noise`.

**Hierarchical bootstrap:** 2,000 draws, seed 20260820. (a) resample exact keys
and protein components with replacement, multiplicity = product of the two draw
counts, recompute the weighted pooled between-key MS; (b) resample T0 L1
repeated-measure groups with replacement, recompute pooled residual variance,
double it; (c) pair by position, `theta_b = MS_b - sigma2_noise_b`.

**Frozen V1 gate: the 2.5th percentile of `theta_b` must be > 0.**
If V1 fails: **stop before any neural training.** No larger network is a rescue.

## 5. V2 — the corrected operator (runs only if V0 and V1 pass)

### 5.1 Architecture contract

* **Edit branch**: label-blind transformation token from core / `R_a` / `R_b`
  counts, attachment environment, charge and stereo flags, plus 256-bit folded
  Morgan fingerprints of core, `R_a`, `R_b` (Stage U `edit_features`).
* **Protein branch**: the governed ESM-2 150M bank's 128 **ordered** residue
  region tokens (640-d), projected to model width, fixed sinusoidal slot
  encoding, mask handled.
* The edit token **queries** region tokens through two multi-head cross-attention
  layers (4 heads, width 128, FFN 256, dropout 0.1). **No pooled protein
  summary, no target embedding, no target index, no component ID, no
  assay/document ID** reaches the operator.
* `D_hat(tau,p1,p2) = R(tau,p1) - R(tau,p2)`: identity, protein-pair
  antisymmetry and cycle consistency hold for every parameter setting.
* **NEW — nested zero predictor.** `R` is multiplied by a scalar gate
  `g = softplus(a)` with `a` initialised so `g ≈ 0`. The zero predictor is
  therefore *inside* the hypothesis class and the candidate can abstain; it
  cannot lose to `A_zero` on error merely because its output amplitude is
  uncontrolled. This repairs Stage U defect 6.

### 5.2 Matched arms

| arm | description |
|---|---|
| A `A_zero` | constant response; `D_hat` identically 0 |
| B `B_global` | global ESM pooled summary + edit token (Stage S/T reference) |
| C `C_local` | **candidate**: edit token cross-attending ordered region tokens |
| D `D_local_shuffled` | arm C trained on stable cross-component shuffled proteins |
| E `E_label_shuffled` | arm C, correct protein, trained on within-key permuted `D` |

### 5.3 Primary causal controls — substitutions inside the trained candidate C

All on identical rows, only the protein input replaced:

1. **correct protein**;
2. **stable shuffled protein** (cross-component permutation within the
   recipient's own population);
3. **similarity-matched wrong protein** (different CD-HIT40 component, most
   similar admissible by cosine on frozen pooled ESM, drawn from the recipient's
   own population);
4. **NEW — residue-token permutation** (the 128 ordered region tokens permuted
   by a stable seed; the protein's content is preserved and only its order is
   destroyed). Repairs Stage U defect 1;
5. **NEW — capacity-matched random protein representation** (a fixed random
   token matrix of identical shape and per-feature moments, one stable draw per
   target). Repairs Stage U defect 2;
6. **protein-blind reference** (`E_mean_tokens`-style: region tokens replaced by
   the target-independent masked mean over fit components, which makes `D_hat`
   structurally 0).

Arms D and E are **secondary** controls. The primary correct-vs-shuffled gate is
the paired substitution (control 2), never the separately trained arm.

### 5.4 Matched-arm discipline — NEW, repairs Stage U defect 4

* **identical parameter initialization** for every arm whose architecture
  matches (C, D, E share an architecture and are initialised from the same
  seeded state dict);
* **identical minibatch order** across all arms (the batch seed excludes the arm
  name);
* identical optimizer, schedule and budget;
* **checkpoint-free fixed-budget** evaluation of the final parameters.

### 5.5 Shortcut diagnostic — NEW, repairs Stage U defect 3

Stage U's `fit_unsampled` (a random 10% of fit `D` rows) retains the same
targets and keys as training and cannot detect target-key memorisation. It is
**replaced** by a **fit-target-held-out** bank: a stable, label-blind subset of
**fit components** is withheld from training entirely, and `D` rows whose *both*
targets lie in the withheld components form `fit_heldout`. Memorising a target
key cannot solve it.

Both `fit_heldout` and the internal surfaces are reported.

### 5.6 Hyperparameters — inherited verbatim from Stage U §4.5

Seed `20260821` (screen); confirmation seeds `20260822`, `20260823`. Steps
`4000`; batch `256`; AdamW lr `3e-4`, weight decay `1e-4`, cosine to 0; Huber
delta `1.0`; gradient clip `5.0`; row sampling weight `1/sqrt(deg_exact_key)`.

## 6. Frozen V2 success conditions — inherited from Stage U §4.6, with defect 5 repaired

Evaluated on the primary evaluable internal surface; paired contrasts use
identical rows. Two-way cluster bootstrap over components and keys, 2,000 draws,
seed `20260820`; effective independent units `min(#components, #keys)` reported
with every interval.

1. `C_local` minus `B_global` Pearson >= **+0.05**, **and** `C_local` beats
   `A_zero` on **MSE and MAE** (a Pearson contrast against a constant predictor
   is undefined — this is the repair of Stage U defect 5, and it is strictly
   stricter, since C must now beat the zero predictor on error rather than being
   compared to an undefined quantity);
2. both differences have clustered 95% lower bounds **> 0**;
3. C correct-input minus C **shuffled**-input Pearson >= **+0.05**, lower bound > 0;
4. C correct-input minus C **matched-wrong**-input Pearson >= **+0.05**, lower bound > 0;
5. C vs `D_local_shuffled`: MSE delta < 0 **and** Spearman delta > 0 **and**
   CI delta > 0 **and** sign-accuracy delta > 0;
6. C minus `E_label_shuffled` Pearson >= **+0.05**, lower bound > 0;
7. protein-induced shift (C correct minus C shuffled) correlates with truth
   >= **+0.10** — large movement without alignment is failure;
8. on `fit_heldout`, C correct-input minus C shuffled-input Pearson >= **+0.05**
   with lower bound > 0;
9. leave-one-out influence: no single transformation key and no single protein
   component accounts for more than **50%** of the C-vs-shuffled effect;
10. on a second evaluable surface, C correct minus C shuffled Pearson >= **0**
    (transformation/scaffold-cold does not reverse); a surface with < 100 rows is
    `not_evaluable` and the route cannot pass on it;
11. **NEW** — C correct-input beats **residue-permuted**-input and
    **capacity-matched random**-protein input on Pearson by >= **+0.05** each,
    lower bounds > 0.

**A single seed may reject; a single seed may not confirm.** Confirmation
requires all conditions in each of three fixed seeds (20260821/22/23).

**Disclosure required in the report:** bootstrap intervals are computed from one
checkpoint per arm and therefore capture sampling variability of the evaluation
rows, **not** training stochasticity. Only the three-seed protocol speaks to the
latter.

## 7. Stop and claim rules

* **V0 census fails** -> close only the *estimability of the exact-MMP route on
  this dataset*. This is a statement about support, not about biology.
* **V0b leaves no evaluable surface** -> same conclusion: the requested estimand
  is not identifiable here. Insufficient support is **not** biological absence.
* **V1 fails** -> stop before neural training; the interaction variance does not
  exceed supervision noise on this corpus.
* **V2 single-seed fails** -> reject **this operator and this protocol**. Do not
  claim that all protein-conditioned interaction representations are impossible.
* **V2 passes** -> update Phase 1 only; no few-shot adaptation until the
  three-seed confirmation and all controls pass; nothing promoted to `model/` or
  `scripts/`.

## 8. Verification required before any training

Deterministic core-inclusive MMP decomposition; inverse/sign consistency;
attachment, stereochemistry and charge preservation; core presence in the key
and in the descriptor; no cross-target, cross-panel or cross-split
contamination; physical meta-test seal; bank stability across `PYTHONHASHSEED`;
identity, protein-pair antisymmetry and cycle consistency; nested-zero-gate
behaviour; no dead trainable parameters; no target-ID or pooled-protein bypass
in the local operator; no label path into inputs, keys, splits or selection.
