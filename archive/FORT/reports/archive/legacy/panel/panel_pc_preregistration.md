# Gate PC preregistration - Target-Conditional Bayesian Functional Meta-Adapter on the dense panel

Registered 2026-07-25 before any PC statistic was computed and before the k=4 power audit was run.
This is the blueprint's mandatory core module (BHR-MoT-DTA module 1) on the first substrate in this
program where the interaction it is supposed to infer has been shown to exist
(`PANEL_GATE_PA_PASS`, adaptive `p=0.000488`).

## What is unchanged

The posterior mathematics of `model/bayesian_meta.py` is reused without modification: a zero-mean
Gaussian prior over an 8-dimensional function space whose diagonal precision is predicted from the
frozen target representation, one exact FP32 Cholesky solve on support residuals, a task-level Bayes
factor mapped to a mixture weight by a fixed model prior, and an exact `k=0` fallback. Support
labels never enter an encoder, no target-id feature exists, the prior mean is structurally zero, and
`mu0(t,d) = b(d)` as the BM0 contract requires. Full-parameter MAML, support-label tokens, learned
free-form gates, Student-t likelihoods and test-time encoder updates remain rejected.

## What is new

Only the substrate and the evaluation protocol. Episodes come from the panel's development cells of
held-out homology components under the same leave-component-out folds as Gate PB, with `k=4`
scaffold-disjoint support ligands and a fixed query set of the remaining ligands at Tanimoto < 0.95
from the support. Meta-training uses the training components' own train cells and cross-fitted
out-of-fold base predictions; the adapter never sees a held-out component during its fit.

## Frozen threshold

`research/panel_power_k4.py` fits the identical ligand-only base at four seeds on the identical
fixed k=4 query rows before the adapter exists, and freezes
`max(0.03, empirical MDE80)` as the Gate PC threshold. The full-query panel audit already gave
MDE80 `0.0181`; the k=4 audit repeats it on the exact rows PC will score.

## Criteria

| id | criterion | threshold |
|---|---|---|
| PC1 | `BM0P - B0` paired component macro Spearman | `>= ` frozen k=4 threshold |
| PC2 | grouped component bootstrap of that difference | `LCB95 > 0` |
| PC3 | RMSE safety versus B0 | `<= 1.02 x B0` |
| PC4 | correct-target support beats wrong-target support | `LCB95 > 0` |
| PC5 | correct labels beat permuted support labels | `LCB95 > 0` |
| PC6 | real protein beats BOTH shuffled and random protein representations | both `LCB95 > 0` |
| PC7 | exact `k=0` fallback | maximum absolute deviation exactly `0` |
| PC8 | support-permutation invariance | maximum absolute deviation `< 1e-4` |
| PC9 | finite positive predictive variance | `min > 0` and finite |

PC1-PC5 and PC7-PC9 are the unchanged BM0 criteria. PC6 is the criterion that
`BM1_RR_FAIL_STOP` failed (`-0.0005` against protein shuffle, `+0.0035` against random protein);
it is included so that a support-conditioned ligand kernel masquerading as protein conditioning
cannot pass. A single failed criterion stops Gate PC.

## What a pass authorizes

Three seeds of Gate PC, and nothing else. The optional Hierarchical MoT, long training, pretrained
multimodal warm starts and confirmation access all remain unauthorized and require a separate
review, exactly as task.md states.
