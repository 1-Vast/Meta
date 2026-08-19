# Stage CIIP-1A collapse audit — read-only diagnosis (2026-08-19)

Schema: MetaSieve.StageCIIP1A.CollapseAudit.v1. Every frozen input is
pinned by SHA-256 in STAGE1_COLLAPSE_AUDIT.json; this audit performs no
training, no retraining, and writes nothing to frozen artifacts.

## 1. Output collapse table (13 test mutation pairs)

| metric | unified_local | free_pairwise | ligand_only |
|---|---|---|---|
| N_total | 13 | 13 | 13 |
| N_target_informative (median, |c|>=10) | 69 | 69 | 69 |
| N_prediction_nonconstant | **3** | **3** | 0 |
| N_rank_evaluable | **3** | **3** | 0 |
| finite_pairs / total_pairs | 3/13 | 3/13 | 0/13 |
| pair-mean Spearman | -0.0409 | +0.2605 | undefined (zero floor) |
| dead-zone sign accuracy | 0.4869 | 0.5242 | 0.4909 |
| centered MSE (pair-mean) | 232.1 | 232.3 | 236.7 |

Per-pair details (unified_local), full tables in the JSON:

| pair | parent | var_true | mse | R2 | slope | sign_acc | sp | nonconstant |
|---|---|---|---|---|---|---|---|---|
| ABL1 Q252H | 174.8 | 184.4 | -0.055 | -9.98 | 0.176 | -0.492 | yes |
| ABL1 T315I | 286.7 | 216.3 | +0.245 | +0.85 | 0.708 | +0.326 | yes |
| FGFR4 V550L | 118.4 | 119.9 | -0.013 | +0.44 | 0.594 | +0.042 | yes |
| ALK R1275Q, KIT A829P, KIT V560G, MET D1228H, MET D1228N, EGFR L861Q, LRRK2 R1441C, RET G691S, RET M918T, TEK Y897S (10 pairs) | var_true 72..571 | mse == var_true exactly | 0.0 | — | 0.40..0.57 | undefined | **constant (var_pred = 0 exactly)** |

NaN Spearman on the 10 constant pairs is mathematically undefined, NOT
a zero correlation; those pairs carry no ranking information.

## 2. The collapse is structurally forced — the decisive finding

The set of constant-prediction test pairs is EXACTLY the set of test
pairs whose KLIFS pocket one-hot inputs are IDENTICAL for WT and
variant (ΔP = 0):

- 38 of 65 pairs (20 train / 8 val / 10 test) have **ΔP = 0**: the
  frozen KLIFS pocket one-hot does not encode the mutation at all for
  these rows (mutation site outside the aligned 85-position pocket, or
  the mutation application produced no feature change).
- For ΔP = 0, the potential contrast is identically zero by
  antisymmetry: s(P_v,L) − s(P_wt,L) = s(P,L) − s(P,L) = 0. No
  training can change this.
- Every nonzero ΔP has norm exactly √2 (a single one-hot flip; the old
  residue is not even cleared), effective rank 18 of 1700 dimensions.

Therefore the "output collapse" on 10/13 test pairs is a
**representation bottleneck of the frozen input features**, not an
optimization failure and not evidence about the potential architecture.
The same applies to free_pairwise (its 3/13 nonconstant pairs are the
same three pairs).

## 3. Potential variance sources (read-only)

| quantity | value |
|---|---|
| KLIFS mutation-diff matrix | 65×1700, rank 20, effective rank 18.0, 38 zero rows, nonzero norms ≡ √2 |
| local-ESM mutation-diff matrix (Q1 cache, 49/65 pairs, radius-6 window at verified site) | rank 49/49, effective rank 31.9, norms 0.256..0.807 (all nonzero) |
| Q1 frozen evidence: klifs_pocket | selectivity −0.086 [−0.229, +0.019] — not significant |
| Q1 frozen evidence: pair_centered_local_esm | selectivity +0.189 [+0.033, +0.363] — the only genuine protein representation that passed |
| delta_alpha at init (ReLU encoder on one-hot) | median 0.0; 38/65 pairs map to IDENTICAL alpha (forced by ΔP=0) |
| Cov(psi) over 183 ligands at init | trace 0.204, effective rank 6.1/8 |
| potential output var at init | 0.0066 (vs target c var ≈ 72..571) |
| f output var at init | 0.044 |

At initialization the encoder already erases 38/65 mutation pairs
(delta_alpha = 0) because the input difference is zero; for the
remaining pairs the potential's output variance is orders of magnitude
below the target variance.

## 4. Gradient competition audit (initialization; no trajectory was persisted)

Reproducing the trainer's exact epoch-0 first batch (same SHA-256 keyed
rng streams, same frozen row masks):

| quantity | value |
|---|---|
| L_abs (raw % inhibition) | 8413.5 |
| L_contrast (centered mutation) | 228.7 |
| g_abs on s-params (alpha,psi) | 163.8 |
| g_ctr on s-params | 0.152 |
| **R_g = ||g_abs|| / ||g_ctr||** | **1081** |
| **C_g = cosine(g_abs, g_ctr)** | **−0.016** |
| g_ctr on b_P / b_L | 0 / 0 (centering cancels main effects — correct) |
| g_abs on b_P / b_L / enc | 316 / 242 / 2203 |

Per the frozen interpretation rules: R_g >> 1 with C_g < 0 — the
absolute objective's gradient on the interaction parameters is ~1080×
the contrast objective's, so the joint step is dominated by the
absolute task; the contrast signal is a ~0.1% perturbation inside the
update. Consistent end-state evidence: best val contrast MSE = 239.8 ≈
the zero-prediction floor (ligand_only = 251.8), i.e. the contrast
objective was never effectively optimized, while the absolute fit
succeeded (f variance at checkpoint = 488). Even the three informative
test pairs show no learned structure (sp −0.49 / +0.33 / +0.04).

## 5. What the free-pairwise result really means

- Nonconstant on only 3/13 pairs (same three pairs as unified, all
  ΔP ≠ 0); parent coverage = 2 parents (ABL1, FGFR4).
- Paired gap vs unified on jointly evaluable pairs: +0.55 (Q252H),
  −0.02 (T315I), +0.37 (V550L); bootstrap 2.5% lower bound of the
  free−unified gap = +0.265 — driven by two pairs in two parents.
- The five conditions for claiming "the integrable potential's
  expression is insufficient" are NOT met: (1) not stable nonconstant
  across parents (2/8 parents); (2) R2/slope advantage confined to 2
  pairs; (3) driven by few pairs (3/13); (4) same inputs/split/seed
  hold, but (5) unified's gradients are healthy only in flow, while its
  inputs are zero for most pairs.
- The free-pairwise advantage is real but small-sample and cannot be
  attributed to non-integrability. The scalar field y(P,L) is still a
  legitimate model class; the observed gap reflects input
  degeneracy + objective competition, not an integrability limit.

## 6. Diagnosis

Evidence-weighted primary causes:

1. **Representation bottleneck — PROVEN.** 38/65 pairs (10/13 test)
   have zero input difference; Q1 independently shows klifs_pocket not
   significant (−0.086) while local ESM passed (+0.189). The frozen
   KLIFS one-hot input is insufficient for this estimand.
2. **Objective competition — PROVEN.** R_g ≈ 1081, C_g ≈ −0.02 at
   initialization; the contrast objective never improved from the zero
   floor even for informative pairs.
3. Potential capacity — NOT implicated (collapse is explained without
   capacity arguments).
4. Evaluation — sound (undefined Spearman handled as undefined; the
   zero floor is the honest baseline).

Per the frozen decision tree, representation AND objective are both
implicated → the authorized successor is the **preregistered 2×2**
{KLIFS one-hot, local ESM} × {joint loss, centered-only loss}, keeping
EVERYTHING else frozen (ECFP4, rank 8, hidden 64, potential formula,
AdamW 1e-3/wd 1e-4, 200 epochs, batch 512, split, seeds, sampling,
gates, metrics, checkpoint rule). A new preregistration must be frozen
before any successor training; no single-factor attribution will be
claimed from any 2×2 cell.

## 7. Final status

```text
tested one-hot potential: FAIL
biological protein-conditioned signal: UNRESOLVED
primary cause: representation + objective (representation dominant; objective competition proven co-cause)
authorized successor: 2x2 (KLIFS/ESM x joint/centered-only), new preregistration required before training
CIIP-1B: NOT AUTHORIZED
BindingDB Potential Bridge: NOT AUTHORIZED
Production integration: NOT AUTHORIZED
```
