# E0R1 Proposal: Identification by Deterministic Solver

Status: REVIEWED, NOT EXECUTED VERBATIM. The three-agent decision rejected
`StandardScaler + Ridge(alpha=10)` as the primary exact-identifiability witness.
The registered stage used raw-centered Moore-Penrose A/B/C plus conditional
corrected-objective D. See `E0R1_RESEARCH_SYNTHESIS.md` and the formal evidence
under `report/mechanism_refactor/p1r2b_e0r1_objective_design_solver_v1/`.
Written 2026-08-07 against `AGENT_HANDOFF.md`, `task.md`, `history.md` F-71/F-73/F-75.

## 1. Diagnosis: E0R0's learned heads did not fail, they did not finish

The `full_240` arm of `scripts/run_e0r0_typed_tensor.py` is an **exactly linear model
in frozen features**:

```python
def forward(self, tensor):            # FullTypedTensor
    return torch.einsum("nard,ard->n", tensor, self.energy)
```

`features` (`centered`) are precomputed and constant, `residuals` are the target, and
the analytic optimum reproduces the target to `4.51e-07`. The objective is therefore a
**convex least-squares problem with 240 parameters and 640 training rows**
(32 train tasks x 20 ligand states). It is solvable in closed form.

Every line of `full_240_trace.jsonl` says the run was travel-limited, not converged:

| epoch | ‖energy‖ | grad L2 (mean) | loss | train CI | holdout residual Pearson |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.731 | 0.0444 | 0.9794 | 0.4936 | 0.420 |
| 3 | 2.032 | 0.0439 | 0.9650 | 0.5446 | 0.459 |
| 55 | 24.284 | 0.0327 | 0.8087 | 0.7531 | 0.623 |
| 60 | 25.797 | 0.0328 | 0.8004 | 0.7535 | 0.630 |

- The parameter norm grows **linearly at ~0.43/epoch from a zero init** and shows no
  sign of levelling. The iterate is marching, not settling.
- The gradient norm is **flat across all 60 epochs** (0.044 -> 0.033). A converged
  convex optimum has a vanishing gradient.
- Loss, train CI and holdout Pearson are **monotone improving at the frozen endpoint**.
- Train CI is only `0.7535` where the analytic solution scores `1.00000`. The model is
  underfit **on the training folds**, which excludes generalization as the explanation.
- Learned-vs-teacher tensor cosine `0.21025` is what an unconverged iterate looks like
  after it has traversed only the high-variance feature directions.

With AdamW the per-step displacement is bounded by roughly `lr` per coordinate, so
`60 x 160 = 9,600` steps at `lr=1e-3` caps the reachable norm. The run consumed that
budget at a constant rate and stopped mid-flight. **`TYPED_TENSOR_LEARNED_HEADS_FAIL`
is a statement about the step budget, not about the typed tensor.**

The `cp_rank_6` arm has the same signature plus a second defect: factors initialized at
`randn * 0.1` make the reconstructed tensor entries ~`1e-3`, so the multiplicative
parameterization starts in a vanishing-gradient region far from any good basin.

### Why this matters upstream

The same optimizer, budget, init and loss produced the original E0 synthetic pre-gate
failure (F-71: correct CI `0.68553`, partner delta `+0.03618`). The learned full-240
head reaches `0.69013 / +0.05987` — the same place. E0S localized the loss to `T2 -> T3`
and could not separate realization from optimization because no trace existed
(F-73). **E0R0 persisted the trace, and the trace resolves it: optimization.**

## 2. Precedent: this is an engineering failure, not a scientific gate

The ledger already distinguishes these classes and repairs rather than stops:
F-30 (runtime non-viable), F-35 (serial CPU underuse), F-19 (backgrounding), F-38
(temp-dir ACL), F-44 (collection pollution). None consumed a scientific verdict.

E0R0 belongs in that class. Recommended reclassification, subject to user decision:

```text
P1R2B_E0R0_LEARNED_HEADS_NOT_CONVERGED_WITHIN_BUDGET
```

This does **not** promote E0 to PASS, does not unseal affinity labels, and does not
authorize E0-S, T, DAVIS or P2-P4. It records that no biological inference is available
from E0/E0R0 in either direction.

## 3. Proposed E0R1: identification by deterministic solver

The stop rule forbids "post-hoc epoch, learning-rate, loss or initializer retry", and it
should be honoured. Re-running the same estimator with better hyperparameters is exactly
the prohibited move and would repeat the F-29 threshold-substitution failure.

The clean alternative is to **replace the estimator with one that has no tunable knob**,
so the outcome cannot be selected. A closed-form solve is not a retry: it is the exact
minimizer of the already-registered objective class.

**E0R1-A (primary).** Solve the 240-parameter tensor head by closed-form ridge on the
frozen train folds, using the project's standing convention `StandardScaler + Ridge
alpha=10` (the same estimator used in P1R1, P1R2A, P1R2B0 and the E0 ligand prior).
Score once on the unchanged fold-4 holdout under the unchanged Gate
(`correct CI >= 0.80`, gains over ligand and deranged each `>= 0.10`, permutation error
`<= 1e-6`). Register alpha before running; do not scan it.

**E0R1-B (preflight, label-free, run first).** Report the **rank and condition number of
the 240-D train feature Gram matrix**, and the fraction of the analytic tensor's energy
lying in its null space. This is the genuine identifiability question and it is cheap:

- Full rank, small null-space energy -> the typed statistic is identifiable from this
  corpus and E0R1-A must recover it.
- Rank deficient with material null-space energy -> some `ligand-type x residue-type x
  distance-bin` cells never co-occur in EnergyPilot.v1. **No estimator can identify them**,
  which is a real scientific boundary and a corpus-coverage finding, not a model failure.

Run B before A. B alone may settle the stage.

### Decision table

| E0R1-B | E0R1-A | Verdict | Consequence |
|---|---|---|---|
| full rank | Gate PASS | `E0_OPTIMIZATION_DEFECT_CONFIRMED` | E0's synthetic pre-gate failure is an optimizer artifact. The E0 MAP training protocol must be repaired and the pre-gate rerun before E0-S is reconsidered. |
| full rank | Gate FAIL | `TYPED_STATISTIC_INSUFFICIENT` | Genuine: the frozen geometry plus typed statistic does not carry the teacher functional at the required precision. A real boundary. |
| rank deficient | either | `SYNTHETIC_CONTROL_NOT_IDENTIFIABLE_FROM_CORPUS` | The control is mis-specified for this corpus. Redesign the teacher over observed cells before any further E0 work. |

Note the asymmetry deliberately: a PASS here buys **nothing biological**. It only removes
the optimizer as an explanation and returns E0 to its fail-closed state pending a repaired
pre-gate. It must not be read as authorization for E0-S, real affinity, DAVIS or T.

## 4. Secondary defects to fix in the same registration

Fix these together so they cannot each cost a separate stage:

1. **No feature scaling.** `run_e0r0_typed_tensor.py` centers features
   (`correct - feature_center`) but never scales them. Contact-normalized typed
   co-occurrence counts span orders of magnitude across the 240 bins. Every other stage
   in this project standardizes before Ridge; E0/E0R0 dropped it.
2. **Convergence is not a gate precondition.** Register the general rule:
   *a gate may only be evaluated at a converged optimum; a run whose gradient norm has
   not decayed and whose loss is still falling is `NOT_RUN`, not `FAIL`.*
   This is the durable fix — it prevents the whole class of F-75-style verdicts.
   Add a regression test that fails a run whose final gradient norm exceeds a registered
   fraction of its initial value.
3. **CP init.** `_analytic_cp_tensor` already computes an SVD-based factorization.
   If CP training is kept at all, initialize its factors from the HOSVD of the ridge
   solution rather than `randn * 0.1`.
4. **Weight decay.** `WEIGHT_DECAY=1e-4` in AdamW shrinks toward zero while the target
   tensor's norm is set by `1/scale`. Register it as `0` for the tensor head or justify
   it; it is currently an unexamined bias against the known solution.
5. **Rank loss starvation.** `BATCH_SIZE=4` within one task yields `C(4,2)=6` pairs per
   step. Batching the full ~20 states per task yields 190. The rank term is currently
   near-inert relative to its weight of 1.
6. **Holdout too small to support the Gate.** Fold-4 holdout is **8 tasks / 8 closure
   components** (F-73 caveat). A `+-0.10` CI margin estimated on 8 macro units is very
   noisy, and one of the eight derangements has local identity `0.54545`, violating the
   `<40%` rule. Enlarge the synthetic holdout and repair that derangement, or a correct
   estimator can still fail by sampling noise.

Items 1-5 are defects in the *procedure*; none of them requires reading a label.

## 5. Strategic note, offered separately from the above

Across F-08 to F-75 the ledger contains roughly two dozen distinct attempts to show that
correct-protein information improves affinity ranking beyond a ligand-only baseline:
statistical, neural, cross-attention, GO-MF supervision at two scales, residue-site
supervision, fixed physicochemical pair features, a frozen published PSICHIC latent, and
a structural mechanistic bridge. P1B learns partner-specific geometry decisively
(AUPRC `0.43885` vs `0.05149` deranged). None of it transfers to affinity.

P1R2A gives the most economical statement of why: source candidate variance is
`77.58%` ligand main effect, `20.76%` protein main effect, `1.67%` interaction. The
quantity being chased is roughly one part in sixty of the signal, on benchmarks with
complete ligand overlap across target splits (noted as far back as the F-11 root-cause
assessment).

That is worth naming as a finding rather than a series of setbacks. Two honest options
exist alongside continuing: (a) power-analyze the target effect first — establish what
corpus size and measurement precision a `+0.03` CI contrast actually requires before
registering another gate; or (b) write up the negative result, which is well-governed,
extensively controlled, and more publishable than most positive DTA claims. Neither is a
recommendation to stop; both are cheaper than the next gate.
