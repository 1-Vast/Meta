# Stage 0/1 (R5): contract repairs and structural gates

Numerical authority for the gates: `pytest tests/test_reltransport_synthetic.py`
and `tests/test_stage0_contract_fixes.py` (this stage's tests). Gate outcomes
summarized in `RESULT.json`. No real-data performance claim is made here.
`meta_test` untouched.

## Part 1: the 2026-08-16 contract repairs

The audit (subagent report, filed as part of this stage) found six violations;
all six are fixed and covered by tests.

| # | violation | fix |
|---|---|---|
| 1 | evaluation wrong-protein donors drawn from `meta_train` (seen/unseen confound; docstring said the opposite) | `matched_donors(donor_pool, whitening_pool)`; evaluation uses `donor_pool="meta_val"`; training counterfactuals use `donor_pool="meta_train"` |
| 2 | whitening fitted on whatever pool the donors came from | `whitening_pool="meta_train"` always; verified by recomputation test |
| 3 | COMPARE payload described donors as meta_train/train-whitened while code used meta_val/meta_val-whitened; per-arm RESULT.json had no donor description at all | explicit `donors` block in every RESULT.json; COMPARE description corrected; rerun in R6 |
| 4 | gradient cosines saved from one episode of one step of one seed | per-step mean over all episodes + per-episode records; `gradient_summary` aggregates across steps (mean, conflict frequency); seed aggregation at stage level |
| 5 | `train_qpsmp.py`/`evaluate_qpsmp.py` auto-evaluated `meta_test` by default | physical seal: `QPSMPData(include_meta_test=False)`; explicit `--include-meta-test`/`--eval-meta-test` opt-in; `--split` default `meta_val`; historical A0 meta_test files stay sealed and unread |
| 6 | missing artifacts: checkpoint hash, component bootstrap, activation stats | checkpoint sha256 in RESULT.json; activation statistics per row; bootstrap in the R6 compare script; gradient-coverage census |

`test_stage0_contract_fixes.py`: 6 fast contract tests + 2 slow end-to-end
smokes (both pass). The old tests (test_level_shape, test_qpsmp_data, ...)
still pass (61 selected).

## Part 2: Core Innovation A design notes (post-hoc to the gates, but fixed before real training)

The gate suite forced three design decisions, each recorded here because they
changed the architecture:

1. **The relative potential is bilinear and the transport is an additive
   relative correction, after the multiplicative gate was falsified by the
   R6a screening.** The multiplicative form `(1+tanh(delta))*r_k` cannot
   express the optimal per-query correction (its functional form mismatches
   the optimal residual scaling) and measured inert (nogate gap 0.000,
   eliminated under S1/S3). The final transport implements the exact
   residual identity, label-free:
   `t(q) = shrink * sum_k a(q,k) * [r_k + delta_hat(q,k) - (f0(q) - f0(k))]`
   with the bilinear antisymmetric `delta_hat(P,i,j) =
   u(P)^T[psi(e_i)phi(e_j) - psi(e_j)phi(e_i)]`; when `delta_hat` agrees
   with the endpoint's implied relative the correction degenerates to the
   support residual (level shift), and at k=1 it is query-specific through
   both delta terms. Recorded honestly: the anchor set acts through two
   summary vectors rather than M independent reference ligands; every anchor
   still receives gradient, the exact anchor-mean-zero property is
   preserved, and no constant can enter the shape branch.
2. **The k=1 magnitude-matched label flip is evaluated, not trained.** The
   flip contrast destabilizes a query-specific gate (measured on the
   synthetic gates: it trains the gate to fit per-episode error ratios).
   The k>=2 permuted-label contrast and a new wrong-support-ligand contrast
   (same label, different ligand) remain in training; the flip is a
   Stage 3 evaluation control.
3. **The binding contrast target is the full squared error** with the
   endpoint frozen (deployment metric), while ranking remains the primary
   signal for the interaction trunk.
4. **An identifiability pin closes the routed shape branch's mean freedom**
   (`0.3 * mean_q(shape)^2`, label-free). The first R6 screening seed
   reproduced the Stage R3/R4 predrift defect (ligand_only MSE > 99 while
   the sum stayed near 3): under routing the shape mean is a null direction
   of every objective and the anchors drift while `target_level` compensates.
   The pin stops the drift; the same defect-fix pattern B3 needed.

## Part 3: gate outcomes

`tests/test_reltransport_synthetic.py`: **23/23 pass** (algebraic, gradient,
k=1 mechanism, and the three synthetic-training gates; synthetic gates
measured over multiple seeds/panels as preregistered).

The synthetic interaction gate deserves one sentence: three independently
trained tiny models order a held-out bilinear protein-by-ligand task at mean
CI 0.72 with the interaction branch and at 0.46 without it, positive in
23/24 (seed, protein) cells.

## Resources

- No GPU training beyond a 2-step pipeline smoke (`smoke/`, 238 MB peak,
  14.2 s, 1,866,105 parameters).
- Complete maintained test suite at this point: 294 historical + new tests
  all green in the touched files; full suite rerun at stage end.

## Verdict

Stage 0 contract repairs: **complete and tested**. Stage 1 gates: **all
pass**. The new model is authorized to enter Stage 2 (R6) four-arm screening.
No admission claim of any kind is made here.
