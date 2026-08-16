# Stage 6: chemistry-grounded support weighting — audited result

Numerical authority: `SIGNAL_meta_val.json`, `SIGNAL_meta_val_production1024.json`,
`SIGNAL_meta_train.json`, `MULTISEED_meta_val.json` (+`.rows.jsonl`),
`PAIRED_ANALYSIS.json`, `REPORTING_meta_test.json` (+`.rows.jsonl`),
`PAIRED_ANALYSIS_meta_test.json`, `RESIDUAL_AUDIT_{A,F}_seed20260812.json`,
`arm_*/RESULT.json`, `seed*/RESULT.json`. Gates fixed in
`PREREGISTRATION.md`; audit findings and corrections in `AUDIT_ADDENDUM.md`.

## Decision

**H1 is ACCEPTED as a validated few-shot mechanism at k>=2, and REJECTED as a
demonstrated improvement over the incumbent `grammar` transport.**

The within-checkpoint contrast — the same trained model with and without the
support weighting — passes on **both** splits, three seeds, with component-level
bootstrap lower bounds above zero for MSE, concordance and Spearman. The
cross-arm contrast between independently trained trunks **contradicts itself
across splits** and no lower bound excludes zero, so no superiority claim is
made.

## 1. Corrected bottleneck diagnosis

Stage 5 concluded "the objective is the blocker". That was over-generalised from
two failed parameterisations. A label-and-chemistry-only audit falsifies it: with
a **fixed** Morgan/Tanimoto kernel, similarity-weighted support labels beat the
support mean on both splits.

Production 1024-bit contract (`SIGNAL_meta_val_production1024.json`), which
matches `scripts/qpsmp_data.py` exactly:

| k | support mean | Tanimoto-softmax | oracle best support | level ceiling |
|---|---:|---:|---:|---:|
| 2 | 1.163 | **0.976** | 0.570 | 0.615 |
| 3 | 1.096 | **0.886** | 0.357 | 0.615 |
| 5 | 1.002 | **0.747** | 0.186 | 0.615 |

The earlier 2048-bit audit gave 0.976 / 0.887 / 0.747 and Pearson -0.347; the
production 1024-bit contract gives 0.976 / 0.886 / 0.747 and Pearson -0.349.
**The bit-width mismatch was immaterial**, and all figures quoted here are the
production-compatible ones.

### The gain is in the weighting, not the residual decomposition

The audit above weights raw labels; the model transports `r_k = y_k - f0(L_k)`.
Applying Tanimoto weighting to residuals of a **frozen `grammar` checkpoint that
never saw the mechanism** (`RESIDUAL_AUDIT_A_seed20260812.json`):

| k | mean residual | Tanimoto residual | nearest residual | oracle residual |
|---|---:|---:|---:|---:|
| 2 | 1.303 | **1.175** | 1.033 | 0.821 |
| 3 | 1.269 | **1.115** | 0.907 | 0.581 |
| 5 | 1.201 | **0.965** | 0.700 | 0.276 |

So the weighting helps residuals even without co-adaptation of `f0`: the
bottleneck was the similarity representation, not the residual decomposition.
Two correctness checks fall out of this audit — `mean_residual` reproduces the
`grammar` arm's measured `level_only` exactly (1.3026 / 1.2685 / 1.2015), and
`tanimoto_residual` on the F checkpoint reproduces F's measured `full` exactly
(1.0648 / 1.0194 / 0.9018).

**Unexploited headroom:** `nearest_residual` beats `tanimoto_residual` at every
k, so the learned scale (gamma ~ 7.99) is too soft. `oracle_residual` at k=5 is
0.246 against the mechanism's 0.902. Not pursued in this stage.

## 2. What the mechanism actually is

A **fixed Morgan/Tanimoto support kernel with learned scalar calibration** — not
a non-learned estimator and not "similarity only" in the sense of nothing being
trained:

```text
f    = f0(P,Lq) + s(n) * sum_k w_qk * r_k
w_qk = softmax_k( gamma * Tanimoto1024(Lq, L_k) )
```

`gamma` and the shrinkage `s(n)` are trained; the kernel is not. Measured
`gamma` after training: 7.982 / 7.985 / 7.990 against an initialisation of 8.0,
so almost all behaviour is carried by the fixed kernel.

By construction the weighting is **inactive at k=0** and **degenerate at k=1**
(softmax over one support). Confirmed empirically: `full - level_only` is exactly
`-0.000000` at k=0 and k=1 in all three seeds.

## 3. Primary evidence: within-checkpoint, three seeds, full banks

Same checkpoint, same episodes, differing only by whether the weighting is
applied. `meta_val` = 44 episodes, `meta_test` = 42 episodes, all eligible
targets.

### Per-seed (positive = the mechanism helps)

| split | seed | k | dMSE | dCI | dSpearman | permutation gap |
|---|---|---:|---:|---:|---:|---:|
| val | 20260812 | 2 | +0.142 | +0.056 | +0.145 | +0.431 |
| val | 20260812 | 3 | +0.201 | +0.085 | +0.235 | +0.425 |
| val | 20260812 | 5 | +0.298 | +0.129 | +0.343 | +0.513 |
| val | 20260813 | 2 | +0.137 | +0.058 | +0.134 | +0.428 |
| val | 20260813 | 3 | +0.200 | +0.071 | +0.204 | +0.399 |
| val | 20260813 | 5 | +0.308 | +0.126 | +0.336 | +0.485 |
| val | 20260814 | 2 | +0.147 | +0.037 | +0.135 | +0.423 |
| val | 20260814 | 3 | +0.202 | +0.076 | +0.246 | +0.412 |
| val | 20260814 | 5 | +0.292 | +0.133 | +0.373 | +0.489 |

**9/9 positive for MSE, CI and Spearman; permutation gaps +0.40 to +0.51.**
The incumbent `grammar` arm on the same episodes gives dCI of only +0.009 to
+0.033 and permutation gaps of +0.062 to +0.199.

### Aggregate paired bootstrap, 9,999 draws

| split | k | metric | component-level mean | 95% CI | LB>0 |
|---|---:|---|---:|---|---|
| val | 2 | MSE | +0.119 | [+0.052, +0.192] | yes |
| val | 2 | CI | +0.038 | [+0.016, +0.061] | yes |
| val | 2 | Spearman | +0.106 | [+0.049, +0.164] | yes |
| val | 3 | MSE | +0.154 | [+0.076, +0.241] | yes |
| val | 3 | CI | +0.053 | [-0.015, +0.106] | **no** |
| val | 3 | Spearman | +0.156 | [-0.005, +0.289] | **no** |
| val | 5 | MSE | +0.235 | [+0.126, +0.347] | yes |
| val | 5 | CI | +0.104 | [+0.030, +0.174] | yes |
| val | 5 | Spearman | +0.275 | [+0.084, +0.455] | yes |
| test | 2 | MSE | +0.097 | [+0.002, +0.202] | yes |
| test | 2 | CI | +0.047 | [+0.018, +0.070] | yes |
| test | 2 | Spearman | +0.086 | [+0.023, +0.140] | yes |
| test | 3 | MSE | +0.141 | [+0.036, +0.255] | yes |
| test | 3 | CI | +0.106 | [+0.039, +0.198] | yes |
| test | 3 | Spearman | +0.250 | [+0.105, +0.445] | yes |
| test | 5 | MSE | +0.118 | [+0.038, +0.196] | yes |
| test | 5 | CI | +0.143 | [+0.096, +0.190] | yes |
| test | 5 | Spearman | +0.326 | [+0.209, +0.436] | yes |

**Component-level lower bounds are NOT all positive.** 16 of 18 exclude zero; the
two that do not are `meta_val` k=3 CI ([-0.015, +0.106]) and `meta_val` k=3
Spearman ([-0.005, +0.289]). Their point estimates remain clearly positive and
their target-level lower bounds are positive, but the component-level interval
does not exclude zero and the ranking claim at `meta_val` k=3 is therefore not
established. Target-level bootstrap is positive for every metric on `meta_val`
and for CI/Spearman on `meta_test`.

**These intervals are conditional on the three trained seeds.** The analysis
averages the three seeds per (component, target) *before* resampling components,
so seed-to-seed variance is not resampled and does not enter the interval width.
The intervals describe uncertainty over homology components given these three
checkpoints, not uncertainty over retraining. A seed-resampling or
seed-as-random-effect analysis would give wider intervals; with three seeds it
would also be very low powered, which is why it was not attempted.

**This is the first mechanism in this project to improve MSE and ranking
together, reproducibly, across seeds and across both splits.**

## 4. What is NOT established

### Cross-arm superiority contradicts across splits

F `full` minus A `full`, component-level bootstrap:

| k | meta_val | meta_test |
|---|---|---|
| 0 | +0.121 [-0.041, +0.301] | +0.077 [-0.328, +0.544] |
| 1 | +0.058 [-0.075, +0.279] | **-0.224** [-0.545, +0.044] |
| 2 | +0.115 [-0.004, +0.269] | **-0.124** [-0.598, +0.186] |
| 3 | +0.149 [+0.028, +0.340] | **-0.091** [-0.503, +0.233] |
| 5 | +0.253 [+0.071, +0.488] | **-0.079** [-0.398, +0.155] |

On `meta_val` F beats A; on `meta_test` A beats F at every k>=1. No lower bound
excludes zero on `meta_test`. The incumbent's own within-checkpoint gains are
much larger on `meta_test` (+0.20 to +0.27) than on `meta_val` (+0.08 to +0.13),
i.e. the `rho` gate is strongly split-dependent while the Tanimoto kernel is
consistent. **Neither arm is established as the better model.**

### k=0 and k=1 are not attributable

The mechanism is inactive at k=0 and degenerate at k=1. Any k=0/k=1 difference
between arms is training-dynamics variation, and both cross-arm lower bounds
cross zero. No k=0 or k=1 claim is made.

### Wrong-protein is a full-system perturbation

Tanimoto weighting consumes no protein information. The wrong-protein
intervention changes `f0` and therefore recomputes every residual, so a positive
wrong-protein gap does **not** demonstrate protein-specific transport. It is
reported as a full-system perturbation only. Protein specificity of the
zero-shot trunk is measured separately by the zero-shot-only wrong-protein
control.

### The small automatic meta_test evaluation was misleading

The 6-episode evaluation written by each training run showed no F advantage and
even F `full` worse than F `level_only` at k=2/3/5 for seed 20260813. The
full-power 42-episode re-evaluation reverses that: F `full` beats F `level_only`
in 9/9 cells. The 6-episode bank is underpowered and must not be used for any
decision.

## 5. Verified facts

| check | result |
|---|---|
| ligand rows / unique drug keys | 9,880 / 9,880 |
| RDKit parse failures | 0 |
| unique canonical molecules | 9,880 |
| canonical duplicate groups | 0 |
| cells | 17,717 |
| fingerprints derived only from corpus SMILES | yes |
| support/query ligand separation | enforced in `materialize`, tested |
| query-label leakage | none; query labels are loss/metric targets only |
| `tests/test_similarity_grammar.py` | 15 passed (13 + 2 added by this audit) |
| all six checkpoints load `strict=True` after the freeze fix | yes |

## 6. Corrections applied by this audit

1. **Wording.** "Similarity only" / "zero-learning" replaced by *fixed
   Morgan/Tanimoto kernel with learned scalar calibration*; measured `gamma`
   reported.
2. **Dead parameters.** With `use_learned_key=False`, `transport.key.weight` and
   `transport.log_temperature` had `grad=None`. AdamW skips `grad=None`
   parameters, so **no completed run is invalidated**. They are now frozen
   rather than removed, preserving state-dict keys and initialisation order;
   both facts are covered by new tests. Trainable count for `similarity_only`
   drops from 7,244,891 to 7,232,602.
3. **Zero-fingerprint docstring** corrected: a zero fingerprint receives the
   lowest similarity logit but still a finite softmax weight; it is not
   excluded. No ligand in the active corpus exercises this path.
4. **Fingerprint contract.** Signal audit rerun at the production 1024-bit
   width; old 2048-bit numbers are labelled as such.
5. **Attribution.** Primary evidence moved to the within-checkpoint contrast;
   cross-arm comparisons demoted and reported with their split disagreement.
6. **Evaluation power.** All decisions now use the complete 44-episode
   `meta_val` bank, never the 6-episode automatic bank, and never `meta_test`.

## 7. Unresolved risks

* Cross-arm superiority is unresolved and split-contradictory.
* `meta_val` was used for checkpoint selection (17 episodes) and for the
  decision (44 episodes); those overlap, so the within-checkpoint effect size on
  `meta_val` is optimistic. The `meta_test` replication is the mitigation, and it
  is consumed for tuning history.
* Only 11 `meta_val` and 6-7 `meta_test` homology components bound every
  bootstrap.
* `gamma` barely trains, so the mechanism is close to a fixed prior; the
  `nearest_residual` result implies a sharper kernel would do better.
* k=0 remains the dominant error term (3.5-3.7 on `meta_test`) and is untouched.

## 8. Next decision

Stage 6 is frozen as specified. The next experiment must address **one** thing
and must not be combined with another change:

* resolve the cross-arm contradiction by training A and F under identical
  initialisation and comparing on a bank neither selected on; or
* sharpen the kernel (the `nearest_residual` headroom), which is a
  single-parameter change to `gamma`'s parameterisation.

k=0/k=1 improvements stay out of this stage.
