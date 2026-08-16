# Stage R14 preregistration: a regression-compatible within-target ranking objective

Frozen before any R14 training run. Population, bank, seeds, budget,
selection rule and wrong-protein contract are identical to R7-R13.
`meta_test` remains sealed.

## One core innovation, and why only one

The mandate allows at most two mutually reinforcing core contributions and
prefers one model plus one training mechanism. **R14 claims one, a training
mechanism.** The model-side candidate (a protein-conditioned amplitude head)
was falsified by the Phase 2 diagnostics before implementation — see
`REPORT.md` Result 3 and `RESULT.json` finding F5. Substituting an unmeasured
alternative to preserve a two-innovation shape would be the auxiliary
decoration the mandate forbids, so the slot is left empty and the cycle is
narrower on purpose.

## Causal hypothesis

The within-target shape term splits exactly:

    shape = Var(y)(1 - r²) + (sd_p - r·sd_y)²

Phase 2 measured that **every ranking-primary arm in this project has a worse
ordering floor than the plain MSE-trained incumbent (8/8)**, and that the
cause is the training method rather than the architecture (G1, zero
architecture change: `r` 0.213 -> 0.134).

The proposed cause is the one the calibrated-ranking literature identifies
([RCR, arXiv:2211.01494](https://arxiv.org/abs/2211.01494)): a pointwise
regression loss is minimised at `p_i -> y_i`, while RankNet/hinge/softmax
ranking losses are minimised by any *scale-free* monotone arrangement. Those
are different global minima, so the ranking term spends capacity buying
coarse pair statistics — which is what CI and cliff-sign accuracy measure —
at the cost of the fine-grained correlation `r` that MSE needs. Under
gradient routing (R11) the conflict is moved, not removed; under loss
reweighting (R9, R10, R12) it is traded, not removed.

**Hypothesis H0.** If the within-target ranking term is replaced by one whose
global optimum *coincides* with the regression optimum, the conflict
disappears at its source, and `r` rises instead of falling.

## The mechanism

A within-target listwise cross-entropy whose normalizer uses the model's own
value link rather than `exp`:

    ListCE(s, y) = -(1/C) Σ_i w_i log[ T(s_i) / Σ_j T(s_j) ]

over the queries of one episode's target, with `T` the identity on a
positive-shifted pK score, `w_i = y_i - min_shift`. The alignment claim to be
derived and then verified numerically before any GPU run: the squared-error
optimum `s_i = y_i` also minimises `ListCE`, because then
`T(s_i)/Σ_j T(s_j) = y_i/Σ_j y_j` exactly.

The total objective is `L = MSE + α·ListCE`, with the incumbent architecture
unchanged, no routing, no gate, no cliff weighting, no variance term.

## Simpler alternatives, and why they are not the experiment

* **Just train MSE-only for longer.** That is A0. It holds the ordering
  record and the CI record already; the budget lever is separately recorded
  as untested and requires a matched A0 retrain.
* **Just delete the ranking term.** Also A0. The point of R14 is whether a
  ranking term can be added *without* the measured cost.
* **Tune α on the existing RankNet loss.** R9/R10/R12 are three dose/form
  sweeps of exactly that; all three failed. The claim here is that the loss
  *form's optimum* is wrong, not its weight.

## Information flow

Unchanged from the incumbent: pooled ESM protein slots and 2D ligand graphs
in, one scalar pK out. The objective sees only the query panel's predictions
and labels within a single episode's target. No query label enters the model,
no target memory, no retrieval, no test-time gradient, no closed-form solve.

## Expected failure modes

1. **The alignment identity does not hold for squared error.** The RCR proof
   is for a sigmoid link with binary relevance. If the derivation fails, the
   direction dies before implementation — this is checked first.
2. **Small panels.** 16 queries make the normalizer noisy; the gradient may
   be dominated by the largest-`y` ligand.
3. **The shift is a hidden hyperparameter.** `T` on a positive-shifted score
   re-weights low-affinity ligands. One value is fixed in advance; no sweep.
4. **`r` rises but calibration regresses**, reproducing the trade one level
   up. Gate O2 catches this.
5. **`r` does not move at all**, which would mean the loss form is not the
   cause and would close the objective axis entirely — the same verdict R12
   reached for the margin form, and an acceptable outcome.

## Cost

No new parameters (objective-only change). No change to peak memory or
wall time beyond the listwise term's `O(panel)` cost, which is negligible
against the forward pass. Incumbent: 7,294,171 parameters, ~6,500 MB peak,
~150 s per 1200-step run.

## The single variable

**Correction, before any run.** The incumbent A0 is *not* a ranking-free
control: its recorded config carries `ranking_loss_weight = 0.5`, a RankNet
term applied to the whole prediction alongside a dominant smooth-L1 term.
This makes the experiment cleaner than first drafted — the single variable is
a **loss-form swap in place**, and the incumbent is itself the matched
misaligned control.

The incumbent `similarity_only` configuration exactly, 1200 steps, seeds
20260815/16/17, changing only the form of the existing ranking term at its
existing weight 0.5: `RankNet -> regression-compatible ListCE`. No weight
change, no architecture change, no new parameter.

## Arms

* **A0** — the frozen incumbent checkpoints. RankNet @ 0.5. The misaligned
  control; already trained, not re-run.
* **R1** — identical, with the ranking term's *form* swapped to the
  regression-compatible ListCE @ 0.5. The candidate.
* **R3** — identical, with `ranking_loss_weight = 0.0`. **The necessity
  control.** Phase 2 showed the misaligned term costs `r`; R3 asks whether
  simply *deleting* it recovers the same ordering. If R3 matches R1 then the
  innovation is not a source — removal is — and the claim must be withdrawn.

R3 is the arm that makes O4 meaningful. Without it, any gain by R1 over A0 is
consistent with "the RankNet term was harmful", which is a deletion, not a
contribution.

## Gates (preregistered, primary metric first)

* **O1 (primary).** R1's k=0 within-target Pearson `r` exceeds A0's 0.213
  with a positive component-level paired bootstrap lower bound.
* **O2.** R1's k=0 calibration does not regress beyond A0's 1.236 by more
  than 0.02.
* **O3.** R1's k=0 MSE improves on A0's 2.149 (point estimate) and CI is at
  least A0's 0.580 - 0.01.
* **O4 (necessity of the innovation).** R1's `r` exceeds **R3's** with a
  positive component-level lower bound, so the *aligned ranking pressure* —
  not the mere removal of the misaligned one — is the source.
* **O5.** The direction of O1 holds in all three seeds.

## Failure conditions and what each closes

* **O1 fails** → the loss form is not the cause of the ordering deficit. That
  closes the objective axis for this family, and the next hypothesis must be
  the trunk's ligand representation (a separate, unstarted line).
* **O1 passes, O4 fails** → the gain is deletion of a harmful term, not the
  aligned construction. The core-innovation claim is **withdrawn** and the
  result is recorded as a negative-term finding. This is a real possible
  outcome and it is not a partial success.
* **O2 fails** → the conflict was relocated, not removed; same verdict as
  R11, recorded as such.
* **R3 beats both** → the honest conclusion is that this project's
  within-target ranking terms have been net-harmful throughout, which is a
  publishable negative and closes the axis.

## Advancement rule

Only if O1-O5 all pass do three fixed seeds with component-level paired
bootstrap become a formal run. `meta_test` opens exactly once, only for a
candidate that has passed every preregistered `meta_val` gate, and R14 does
not open it.

## Reporting requirements

Every run records config, split hash, seed, checkpoint sha256, per-target
predictions, donors, activation statistics, gradient coverage and resources,
per the R5 contract, and every arm is reported with the ordering/amplitude
decomposition of `scripts/r14_dispersion_audit.py` alongside MSE and CI —
because Phase 2 showed that MSE and CI alone hid the effect for seven stages.
