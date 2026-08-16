# Hierarchical protein-conditioned Bayesian meta-learner — decision (2026-07-26)

Preregistration: `reports/active/hier_bayes_meta_preregistration.md` (thresholds frozen before any
result was read; identical to the global-prior diagnostic; no relaxation). Single seed 1729, CUDA
(`drug` env, RTX 4060). Development diagnostic on **spent** Metz panel rows; confirmation, Davis and
sealed labels untouched. This is not a confirmatory gate.

## What was tested

One structural change to `model/bayesian_meta.py::TransformerBayesianMetaLearner`: the prior over the
target function became **hierarchical and protein-conditioned** — a shared global covariance factor
`L0` plus a low-parameter frozen-ESM-2-predicted correction `Delta_t` to the covariance *factor*, so
the protein orients (rotates/reshapes) the full prior covariance. Prior mean stayed structurally zero,
so exact `k<=1` fallback to B0 and the "no signed zero-shot correction" contract were preserved. The
`Delta_t` map was initialised at ~0 (starts identical to the global prior) and had to earn any protein
dependence from the unchanged episodic query objective. `--protein-prior` was the only change; folds,
episodes, seeds, optimiser, capacity elsewhere and evaluation rows are byte-identical to the control.

The control (global full-covariance prior) was re-run in the same session and **reproduced the recorded
baseline exactly** (`reports/active/transformer_bayes_meta_repro.json` == `..._short.json`), confirming
the pipeline and that the model edit did not touch the global path.

## Result (paired, same folds/episodes/GPU, 104 targets / 96 paired components)

| arm | control ρ | control RMSE | treatment ρ | treatment RMSE |
|---|---:|---:|---:|---:|
| B0 | 0.2737 | 0.7984 | 0.2737 | 0.7984 |
| model | 0.2983 | 0.8076 | 0.2917 | 0.8153 |
| protein_shuffle | 0.2929 | 0.8099 | 0.2954 | 0.8114 |
| protein_random | 0.2819 | 0.8122 | 0.2950 | 0.8109 |
| label_permuted | 0.2452 | 0.8510 | 0.2578 | 0.8246 |
| wrong_support | 0.2675 | 0.8350 | 0.2742 | 0.8211 |

Paired contrasts (grouped component bootstrap):

| contrast | control | treatment |
|---|---|---|
| model − B0 | +0.0253 [−0.0078, +0.0597] | +0.0191 [−0.0080, +0.0493] |
| model − protein_shuffle | +0.0048 [−0.0031, +0.0135] | **−0.0039 [−0.0148, +0.0068]** |
| model − protein_random | +0.0158 [−0.0065, +0.0452] | **−0.0037 [−0.0208, +0.0134]** |
| model − label_permuted | +0.0551 [+0.0075, +0.1053] | +0.0336 [−0.0069, +0.0758] |
| model − wrong_support | +0.0321 [−0.0131, +0.0776] | +0.0184 [−0.0151, +0.0546] |

Frozen checks (treatment): effect≥0.03 **fail** (0.0191); LCB95>0 **fail**; RMSE≤1.02·B0 **fail**
(0.8153 > 0.8144); protein specificity **fail** (both contrasts negative); wrong-support **fail**;
label specificity **fail** (LCB crosses 0); exact k=0 **pass**; support-permutation invariance
**pass** (4.8e−7); support-offset invariance **pass**; finite positive variance **pass**.
Negative-transfer rate rose from 0.356 to 0.423.

## Verdict

```text
HIER_BAYES_META_SHORT_FAIL_REVIEW
PROTEIN_CONDITIONED_PRIOR_NOT_LOAD_BEARING
```

The structural hypothesis is **falsified on this substrate**: giving the protein the ability to orient
the Bayesian prior did not make protein identity load-bearing. It did the opposite — shuffled and
random proteins scored at or above the true protein, the transferable point effect shrank, and
negative transfer rose. The extra protein→covariance capacity fit training-component idiosyncrasy that
does not generalise to held-out targets, so on held-out targets an arbitrary (shuffled/random)
orientation is no worse.

## What this localises

1. The bottleneck is **not** the prior's expressiveness. The global-prior model could already rotate
   via its full covariance; adding protein-conditioned orientation did not help and slightly hurt. The
   earlier note (task.md item 2, "cannot rotate or mix the SAR basis") is not the operative constraint
   for `TransformerBayesianMetaLearner`.
2. The modest ranking gain of the meta-adapter (~+0.02–0.025) is a **support-conditioned ligand
   kernel**, not protein-specific transfer: it survives protein shuffle in the control and is matched
   by shuffled/random protein in the treatment. This is the same redundancy diagnosed in
   `BM1_RR_FAIL_STOP` and `PANEL_GATE_PC_FAIL_STOP`, now shown to persist against a protein-conditioned
   prior.
3. Consistent with the standing analysis (task.md "Registered next decision point", item 1), the only
   arm with a positive point estimate above the frozen minimum remains the **unconstrained jointly
   trained zero-shot interaction I0** (+0.0376, LCB95 −0.0033 — under-powered, not null). The binding
   constraint is component-level **statistical power**, not adapter architecture. Another meta-adapter
   variant is not indicated; more independent endpoint-consistent panels (or a power-budgeted multi-seed
   protocol) are.

No confirmation, Davis or sealed label was read. Panel development rows were already spent before these
runs. No threshold was relaxed; no post-failure rescue (added capacity, contrast loss, or extra epochs)
was attempted, per the standing fail-fast rule.
