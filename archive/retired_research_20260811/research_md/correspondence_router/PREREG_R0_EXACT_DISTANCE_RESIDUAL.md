# R0 Development Preregistration: Exact Residue--Atom Distance Residual

Status: frozen before any R0 model fit or R0 model-result read.

## 1. Question and boundary

R0 asks one question only:

> Does an exact residue--ligand-atom interaction residual recover held-out ordered
> distance information beyond the frozen P1B slot prior and capacity-matched
> non-interaction alternatives?

The P1B checkpoint, its five distance bins `(0, 4, 6, 8, 10, 999)` Angstrom,
the 128-slot pooling rule, Meta-Section, support/query isolation, and the
uncentered positive dual-ridge section are frozen. R0 reads no affinity labels
and cannot authorize an affinity, `z`, Meta-Section, or production migration.

The development corpus is the previously consumed C0 structural corpus. It is
not an untouched confirmation set. A development PASS only authorizes one
future, separately preregistered confirmation on fresh held-out-B components.

## 2. Scientific axis

All learned arms use the same distance-bin residual rank, optimizer, number of
updates, sampling order, and early-stopping rule. The only candidate-axis change
is whether the residual receives exact residue and exact ligand-atom states and
may couple them bilinearly.

The frozen prior is

`q0(atom, residue, bin) = q_P1B(atom, slot(residue), bin)`.

The candidate posterior is

`q(atom, residue, bin) = softmax(log(q0 + eps) + delta(atom, residue, bin))`.

No pooled ligand or pooled protein path may enter `delta`. Exact indices must be
one-to-one on active rows. Only records with 100% sequence-to-coordinate mapping
are admissible.

## 3. Evaluation population and split

- A complex's complete mapped heavy-atom by residue distance tensor is atomic:
  it cannot be split.
- All mapped residue x ligand-heavy-atom cells are evaluated. There is no 6 A
  positive sampling and no far-cell negative sampling.
- The bootstrap and inference unit is the protein closure component: canonical
  receptor identity, registered sequence-homology edges, same PDB entry, repeat
  chains, and duplicate structures are closed together.
- Exact ligand graph and Bemis--Murcko scaffold are train-to-heldout exclusion
  filters. They are not unioned into protein components because that previously
  produced a non-informative giant component.
- No closure component, exact ligand graph, or scaffold may straddle train and
  heldout-A after exclusions. Heldout-B remains sealed.
- Three fixed seeds are fit and ensembled by arithmetic mean of their posterior
  probabilities. Seeds cannot be selected post hoc.

## 4. Primary score and guard

For `B=5` ordered bins, pair-level ranked probability score is

`RPS = mean_(k=1..B-1) (F_pred(k) - 1[y <= k])^2`.

Scores are averaged first within complex, then within protein closure component,
then equally across components. Lower is better. All confidence bounds use an
equal-weight component bootstrap and the fixed seed-ensemble prediction.

Full-bin negative log likelihood is the only guard. Expected-distance MAE,
threshold calibration, and contact AP at 6 A are diagnostics and cannot decide
PASS.

## 5. Frozen nulls and interventions

- `N0 PRIOR`: frozen P1B slot posterior copied to exact residues.
- `N1 SLOT_SHARED`: the candidate's parameter budget and training protocol,
  with residue states averaged within slot so residues in one slot are
  indistinguishable.
- `N2 ADDITIVE_EXACT`: exact atom and residue marginal residuals are allowed,
  but no atom--residue product or cross term is allowed.
- `N3 RES_DERANGE`: after fitting the full candidate, derange exact residue state
  and residue chemistry within each slot; keep atom inputs, prior, and labels
  fixed. Do not retrain.
- `N4 ATOM_DERANGE`: after fitting, derange atom state and atom chemistry within
  each ligand; keep residue inputs, prior, and labels fixed. Do not retrain.

Define `S_free = min(RPS(N0), RPS(N1), RPS(N2))` on heldout-A. This minimum,
not a preferred baseline chosen after inspection, is the comparison baseline.

Wrong-protein normalized-position substitution is reported only as a diagnostic
in R0 because its residue axis is counterfactual. It cannot rescue a failed
identity-derangement Gate.

## 6. Pre-fit ceiling and power audit

Before GPU fitting, record:

- `S_prior = RPS(N0)`;
- `S_exact_star = 0`, the one-hot exact-bin oracle;
- `S_add_star`, a label-only additive row + column + slot oracle evaluated by
  fixed checkerboard cell cross-fitting so the held-out cell is not read.

The minimum scientific effect is frozen as

`delta_star = 0.05 * (S_prior - S_exact_star) = 0.05 * S_prior`.

Training is forbidden and the verdict is `R0_DEV_NOT_RUN_FAIL_CLOSED` if any of
the following holds:

- `S_add_star - S_exact_star < delta_star`;
- the component-bootstrap 80% power minimum detectable effect at one-sided
  alpha 0.05 is greater than `delta_star`;
- the largest component share is at least 0.20 or component resampling is not
  meaningful;
- mapping is not exactly 100%, a split dependency straddles, active exact
  indices duplicate, or an ordered bin (including overflow) has zero support.

## 7. Gates

Let `S_full` be the seed-ensemble candidate RPS and `G = S_free - S_full`.

### G1: prior incremental

`S_prior - S_full >= delta_star`, its one-sided component-bootstrap 95% lower
confidence bound is positive, and all three individual seeds have positive
direction.

Failure verdict: `DISTANCE_RESIDUAL_NOT_INCREMENTAL`.

### G2: exact-pair incremental

`G >= delta_star`, its one-sided component-bootstrap 95% lower confidence bound
is positive, and all three individual seeds have positive direction.

Failure verdict after a G1 PASS: `MARGINAL_OR_SLOT_RECALIBRATION_ONLY`.

### G3: both exact identities are load-bearing

Both

`RPS(N3) - S_full >= 0.5 * G`

and

`RPS(N4) - S_full >= 0.5 * G`

must have one-sided component-bootstrap 95% lower confidence bounds above zero.

Failure verdict: `ONE_SIDED_IDENTITY_SHORTCUT`.

### NLL guard

`UCB95(NLL_full - NLL_free) <= 0.01 * NLL_free`.

Failure verdict: `ORDERED_SCORE_GAIN_MISCALIBRATED`. No post-hoc temperature or
threshold repair is allowed under this preregistration.

Passing every Gate yields only
`R0_DEV_EXACT_DISTANCE_RESIDUAL_IDENTIFIED`.

## 8. Literature-scoped borrowing

R0 borrows only the exact pair distance-distribution supervision evidenced by
Interformer (Nature Communications 2024) and the ordered radial continuity
motivation from CORDIAL (PNAS 2025). It does not transplant either architecture:

- https://www.nature.com/articles/s41467-024-54440-6
- https://github.com/tencent-ailab/Interformer
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12557521/
- https://github.com/bpBrownLab/CORDIAL

The affinity-oriented typed-energy ledger remains quarantined until a future
structural confirmation and a separate measured-affinity preregistration.
