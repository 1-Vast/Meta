# Preregistration — P1R2B-PHASE2B-S5D

## Ligand-direction collapse and the symmetric-difference conditional estimand

Stage identifier: `P1R2B-PHASE2B-S5D_ESTIMAND_AND_COLLAPSE_DIAGNOSTICS`

Written 2026-08-10, after `PHASE2B_S4R_GATE.json`
(`REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED`, commit `f866e0f`) and before any
S5D code, statistic or Gate value exists.

## 1. What this stage is, and what it is not

S4R's registered stopping rule closes the pose-free **ligand representation**
repair route. This stage does not reopen it. It trains nothing, adds no
capacity, introduces no new representation and reuses the frozen S4R
checkpoints byte-for-byte.

It asks a different question, about the **estimand and its metric** rather than
about the inputs:

> S4R's gain survived a foreign ligand pair almost intact. Is that because the
> `AP_bidir` estimand rewards a sign-agnostic pocket-like field, and does an
> estimand that removes that confound exactly show ligand specificity that
> `AP_bidir` was diluting?

A pass here does not overturn the S4R verdict, which stands. A failure here
strengthens the S4R closure by explaining it.

## 2. The suspected defect

The registered score is

```text
s(P, La, Lb) = (I - Q_P Q_P^T) H_P W [g(La) - g(Lb)],   Q_P = onb{1, b^P}
```

and the primary metric is

```text
AP_bidir = mean( AP(s, gain), AP(-s, loss) ).
```

The gauge removes a two-dimensional nuisance space. The real nuisance is
larger: gains and losses both lie inside the binding pocket, so "which residues
are in the pocket at all" is a shared, ligand-independent structure that
`span{1, b^P}` only partially captures.

`AP` is convex near its floor and the per-pair chance level is small
(`0.0255`). A field that scores pocket residues high therefore raises
`AP(s, gain)` by more than it lowers `AP(-s, loss)`, so `AP_bidir` can exceed
chance **without the score carrying any information about the direction of
change**. That is precisely the observed S4R signature: foreign ligand pairs
cost `+0.000644` and the within-construct chemistry shuffle scores `0.051322`,
above the candidate's `0.046856`.

This is a hypothesis about the metric. It is falsifiable, and section 6 states
what would falsify it.

## 3. Frozen inputs

Byte-identical reuse, no retraining:

- the S4R execution views under `dataset/processed/s7_l2b_r0r/phase2b_s4r/`,
  verified against `PHASE2B_S4R_INPUT_AND_FIREWALL_MANIFEST.json`;
- the S4R `candidate`, `baseline41` and `permuted` checkpoints;
- the frozen control maps, sha256
  `e187a5f00f0b66328877bacd93b22471fe607e382e811f2674ecfc4a9dec9c33`;
- frozen ESM2 residue states, the `b_prior` and the `Q_P` gauge.

Heldout-B is not created and not read. No affinity value, ChEMBL, BindingDB,
DAVIS, KIBA, recipient or metaval value is opened. R6 stays closed.

**Heldout-A has been consumed twice already, by S3R and by S4R.** Every number
this stage produces is development evidence. No S5D result can confirm
anything, and a pass authorizes only a separately registered confirmation on a
panel that is not heldout-A.

## 4. D1 — ligand-steering collapse, label-free

For a fixed protein construct, the estimator maps each pair's ligand difference
to a residue field

```text
f_p = (I - Q_P Q_P^T) H_P W dg_p,     dg_p = g(L_a) - g(L_b).
```

If the estimator genuinely conditions on the ligand, the unit fields
`{f_p / ||f_p||}` across the pairs of one construct must vary. If instead one
direction dominates, the score is a per-protein field wearing a ligand-shaped
hat, and R3 follows automatically.

For every heldout-A construct with at least three eligible pairs, report the
**top principal energy fraction** `rho` of the mean-centred unit fields: the
share of total variance on the first principal direction. `rho -> 1` is total
collapse.

Report `rho` for three objects on the identical construct set:

```text
rho_dg     the unit ligand differences themselves, the data-side upper bound
rho_graph  the candidate's residue fields
rho_base   the baseline41 residue fields
```

`rho_dg` is decisive as a control. If the ligand differences within a construct
are already near-collinear, then no `W` could produce diverse residue fields and
the collapse is a property of the corpus, not of the estimator. Reporting
`rho_graph` without `rho_dg` would be an unfalsifiable claim.

Also report, per pair, the cosine between the candidate's true field and its
frozen foreign-pair field, and the median over the panel.

D1 has no Gate. It is interpreted by the pre-declared rule in section 6.

## 5. D2 — the symmetric-difference conditional estimand

Restrict every pair's comparison to the residues that actually changed,

```text
idx = gain union loss,
AP_cond = AP( s[idx], gain[idx] ).
```

Every residue in this comparison changed, so pocket membership is constant
across the two classes and cancels **exactly and non-parametrically**. The
statistic asks only: among the residues that changed, does the score know which
way they changed?

This estimator is not new. It is `ap_symdiff_conditional`, already implemented
and registered in `p2b_residue_residual.pair_metrics` and already aggregated by
the parent Phase 2B runner. S3R and S4R computed it per pair and did not
aggregate it. This stage aggregates it; it does not invent it.

Eligible pairs are those with `|idx| >= 2` and both classes non-empty, which is
the condition the existing implementation already applies. Eligibility depends
only on labels, so **every arm uses the identical eligible set**, verified by one
common-mask SHA-256. Aggregation is the frozen pair -> construct -> closure
component macro, and every contrast is a paired closure-component bootstrap,
`10,000` resamples, seed `20260903`, one-sided 95% lower bound.

Per-pair conditional chance is `AP` of a constant score on the same conditional
problem, evaluated in closed form by the same exact tied-AP estimator.

Arms, all reusing frozen `W`:

```text
candidate    S4R graph-aware W, true ligand pair
baseline41   S4R mean-pooled W, true ligand pair
foreign      S4R graph-aware W, frozen foreign ligand pair
permuted     S4R permuted-label learner, true ligand pair
chance       constant score on the conditional problem
```

## 6. Gates and thresholds

The conditional problem is roughly balanced, so its chance level is near `0.5`
rather than near `0.0255`. The margins below are the project's standing
practical-effect values reused unchanged; they are not re-tuned for this scale,
and that re-anchoring is declared here rather than discovered later.

| Gate | contrast | margin |
|---|---|---:|
| E1 | candidate - conditional chance | +0.05 |
| E2 | candidate - foreign ligand pair | +0.03 |
| E3 | candidate - trained permuted-label learner | +0.03 |

Non-gating, reported with full intervals:

```text
E4  candidate - baseline41      does the graph statistic help this estimand
E5  baseline41 - chance         the S3R representation under this estimand
```

E2 is the decisive Gate. E1 without E2 would only show that the conditional
estimand is easier, not that it is ligand-specific.

### D1 interpretation rule, pre-declared

```text
collapse confirmed        rho_graph median >= 0.80 and rho_graph >= rho_dg + 0.10
collapse not confirmed    otherwise
```

If collapse is not confirmed, the section 2 explanation of R3 is **wrong** and
the report must say so plainly rather than reaching for a second explanation.

## 7. Terminal verdicts

Exactly one, by earliest failed boundary:

```text
S5D_CONTRACT_FAIL_CLOSED
LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED
POSE_FREE_LIGAND_CONDITIONED_DIRECTION_ABSENT_UNDER_CONDITIONAL_ESTIMAND
CONDITIONAL_ESTIMAND_RECOVERS_LIGAND_SPECIFIC_DIRECTION_IN_DEVELOPMENT
```

`LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED` is reported when D1 falsifies the
mechanism; D2 is still run and reported, because the conditional estimand is
worth measuring regardless of why R3 failed.

`POSE_FREE_..._ABSENT_UNDER_CONDITIONAL_ESTIMAND` is the outcome when E1, E2 or
E3 fails. It **strengthens** the S4R closure: the residue direction is absent
even under an estimand built to remove the pocket confound exactly.

`CONDITIONAL_ESTIMAND_RECOVERS_..._IN_DEVELOPMENT` requires E1, E2 and E3 to
pass. It authorizes exactly one thing: writing a separate preregistration for a
conditional-estimand confirmation on a panel that is not heldout-A. It does not
reopen the representation route, does not open heldout-B or affinity, does not
admit anything to `z` and does not authorize any capacity increase.

## 8. Stopping rules

One run. No threshold, seed, margin, arm or eligibility rule may change after
any S5D statistic is read. No failed Gate may be rescued, and a failure may not
be followed by a third estimand variant on heldout-A.

## 9. Boundary

The frozen theory archive explicitly does not provide ranking, pairwise-ordering
or listwise guarantees, so both the S4R and the S5D estimands are measurement
frontends and neither inherits any theorem. The operator

```text
A(F, z) = K(B(z) F(z))
```

is unchanged, as are CSMO, Band, the positive ridge, the simplex and the fixed
mesh. This stage identifies at most a bounded conditional residue-direction
statistic in development. It does not identify residue-atom coupling,
interaction energy, affinity, selectivity, a few-shot section, a biological `z`
or a validated end-to-end DTA model.
