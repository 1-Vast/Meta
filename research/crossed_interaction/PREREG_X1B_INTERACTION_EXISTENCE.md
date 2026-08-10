# Preregistration — X1B

## Does a real crossed protein-by-ligand affinity interaction exist, and is the few-shot target section identifiable?

Stage identifier: `E-AFF-X1B_INTERACTION_EXISTENCE`

Written 2026-08-10, after `X1A_ICC_AUDIT.json`
(`X1_ICC_PRECONDITION_PASSED`, commit `f1cd61d`) and before any `DD` value has
been computed.

X1B is audit-only. It trains nothing, adds no module and introduces no
parameter. X2 requires X1B to pass **and** its own separate preregistration.

## 1. Why this stage is the meta-learning identifiability test

`README.md` fixes the intended few-shot predictor:

```text
m(P,L) = U^T phi(P,L),     d <= 5
y_hat  = f_L(L, endpoint) + w0^T m(P,L) + a_t^T m(P,L)
```

`a_t` is the target section estimated from the `k` support ligands. The frozen
theory makes the same object central: the task statistic is
`zeta = z(S, Q, gamma)`, so the learned map is support/query conditioned.

Now take the crossed contrast

```text
DD = y(P1,La) - y(P1,Lb) - y(P2,La) + y(P2,Lb)
```

Under **any** additive law `y(P,L) = mu + alpha(P) + beta(L) + noise`,
every term cancels and `DD = 0` identically. The ligand term `f_L(L,endpoint)`
cancels, and so does any target main effect, whatever their functional form.
What survives `DD` is exactly the part of the affinity surface that depends on
the target **and** the ligand jointly.

That is the same quantity `a_t^T m(P,L)` is meant to carry. The consequence is
sharp and is the reason this stage is run before any model:

> If the interaction variance is at or below measurement noise, then the
> support set can identify only a per-target scalar offset `alpha(P)`. The
> optimal `k`-shot rule degenerates to "shift by the support mean". There is
> then nothing for `a_t` to estimate, no task-adaptive direction to learn, and
> a meta-learner would be fitting noise. Conversely a real interaction variance
> is precisely what makes `a_t`, support rank and coverage well posed.

X1B therefore decides whether the few-shot core of MetaSieve is identifiable at
all. It is not a detour from the meta-learning goal; it is its precondition.

## 2. Frozen inputs

The X1A panel, unchanged: X0 `cells.jsonl` (sha256 `898df882…`),
`dependency_components.jsonl` (`8970d059…`), `panels.jsonl` (`f378cdd6…`), and
the pChEMBL values already read under the X1A contract. Ki and Kd stay
completely separate and are never merged, pooled or averaged.

### 2.1 Rectangle reconstruction, and its disclosure requirement

X0-B published cell-disjoint unit **counts** (Ki `11,168`, Kd `1,041`) and
per-cluster sizes, but not the packing itself. The packing is therefore
reconstructed under X0-B's own rule — a deterministic greedy cell-disjoint
packing, panels in ascending `panel_id`, targets and ligands in ascending key
order, a rectangle admitted only if all four cells are present and none of its
four cells has already been used.

X0-B recorded that its counts are "auditable greedy lower bounds, not maximum
packings", so an exact match is not guaranteed. The reconstructed count **must
be reported next to the X0-B count**, and any discrepancy disclosed rather than
reconciled by tuning the traversal order. The statistical unit is not replaced:
it remains the cell-disjoint rectangle inside a dependency cluster.

## 3. Estimand

For each unit, `DD` is formed from the four cell means on the pChEMBL scale.

The measurement noise carried into `DD` is the sum of the four cell-mean
variances,

```text
v_noise(unit) = sum over the four cells of  var_replicate / n_replicates(cell)
```

with `var_replicate` the endpoint's exact-assay replicate variance estimated in
X1A (`0.38190` for Ki, `3.45871` for Kd), and `n_replicates` the number of
measurements backing that cell.

The primary quantity is the noise-corrected interaction variance

```text
I_real^2 = max(0, E[DD^2] - E[v_noise])
I_real   = sqrt(I_real^2)                        noise-corrected interaction RMS
INR      = I_real / sqrt(E[v_noise])             interaction-to-noise ratio
```

Expectations are taken as the mean over units within a dependency cluster, then
the mean over clusters. **`mean(DD)` is not tested**: opposing selectivity
effects cancel, so a null mean is uninformative and a non-null mean is not
required.

## 4. Additive null

The null preserves panel structure, missingness, target degree, ligand degree,
document grouping and closure dependence **by construction**: the design is not
touched at all. Only the values are replaced.

For each endpoint, fit the additive model `mu + alpha(target) + beta(ligand)`
globally per endpoint by sparse least squares — the X1A amendment-01 estimator,
unchanged — then regenerate every cell value as

```text
y_null(cell) = fitted(cell) + Normal(0, var_replicate / n_replicates(cell))
```

and recompute `DD`, `E[DD^2]`, `I_real` and `INR` on the identical unit set.
Seed `20260904`, `200` null replicates. Under this null the true interaction is
exactly zero and only additive structure plus realistic heteroscedastic
measurement noise remains, so the null distribution of `I_real` is the
false-positive distribution of the whole pipeline.

## 5. Reported quantities

1. noise-corrected interaction RMS `I_real`, per endpoint;
2. interaction-to-noise ratio `INR`;
3. prevalence of `|DD| >= 1.0` log units, a label-blind biological floor fixed
   at one order of magnitude of crossed selectivity, alongside the same
   prevalence under the additive null;
4. one-sided 95% component-bootstrap bounds, resampling **dependency clusters**,
   `10,000` draws, seed `20260903`; rectangles and pair rows are never
   resampled as independent observations;
5. sensitivity to restriction to single-document and to single-assay units;
6. the rank-1 target-section diagnostic of section 6.

## 6. Rank-1 target-section diagnostic — reported, non-gating

The README's section starts at `d = 1`, so the natural question is whether `DD`
carries a *consistent per-target-pair direction*. For each target pair with at
least four units, compute the fraction of units whose `DD` shares the sign of
that pair's mean `DD`, and its noise-corrected excess over the same statistic
under the additive null.

A `d = 1` section is only well posed if a target pair has a stable preference
direction across ligand pairs. This is reported as a diagnostic and never as a
Gate, because a failed `d = 1` diagnostic with a passing interaction Gate would
mean the section needs a larger `d`, not that interaction is absent.

## 7. Gates — frozen

```text
X1B-A  I_real >= 0.30 log units, one-sided 95% component-bootstrap LCB > 0
X1B-B  INR   >= 0.50, one-sided 95% component-bootstrap LCB > 0
X1B-C  I_real exceeds the 95th percentile of the additive-null I_real
```

`0.30` is the interaction RMS X1A showed to be detectable for Ki at X0's frozen
interaction-to-noise ratio (`0.5 x 0.618 = 0.309`), rounded down to two
decimals; it is the smallest effect this source can resolve, not a chosen
effect size. `0.50` is X0's frozen ratio, reused unchanged. Kd's detectable RMS
is `0.930`, so Kd is expected to be the harder endpoint and its failure would
be a statement about Kd's measurement noise, not about biology.

Endpoints are judged separately. One may pass while the other fails.

## 8. Terminal verdicts

```text
REAL_CROSSED_AFFINITY_INTERACTION_NOT_IDENTIFIED   no endpoint passes A, B and C
REAL_CROSSED_AFFINITY_INTERACTION_IDENTIFIED       an endpoint passes A, B and C
                                                   with the effect stable under
                                                   the section 5.5 restrictions
INTERACTION_SIGNAL_FAMILY_OR_PANEL_RESTRICTED      an endpoint passes A, B and C
                                                   but the effect is carried by a
                                                   single document, assay family
                                                   or cluster and does not survive
                                                   the restrictions
```

## 9. Stopping rules and consequences

One run. No threshold, seed, floor, null design or cluster definition may
change after any `DD` statistic is read.

`REAL_CROSSED_AFFINITY_INTERACTION_NOT_IDENTIFIED` stops model development.
Under this outcome the honest statement is that the source data do not identify
a target-dependent ligand preference at the resolution available, so the
few-shot target section `a_t` is not identifiable from this source and X2 must
not be trained. That is a data and estimand result, **not** a neural-network
failure, and it may not be answered by increasing model capacity, adding
adaptation modules, or substituting same-target ligand differences, OOF
residuals, wrong-protein penalties or residue contact scores for `DD`.

Only `REAL_CROSSED_AFFINITY_INTERACTION_IDENTIFIED` authorizes writing the X2
preregistration for one minimal `q_theta` map with `d <= 4`.

## 10. Boundary

`model/`, production `scripts/`, `theory/`, CSMO, Band, the mesh, production
`z` and `A(F,z)=K(B(z)F(z))` are unmodified. All code stays under `research/`.
BindingDB, DAVIS, KIBA, PDBbind, recipient labels and every previously consumed
confirmation panel remain unread. X1B identifies at most the existence and
scale of a crossed interaction; it identifies no affinity energy, no causal
mechanism, no few-shot performance and no `z` admission.
