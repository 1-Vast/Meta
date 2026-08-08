# Theory-Biology Integration Design

Status: **design proposal, not validated evidence.** Nothing here is admitted to
`model/`, no Gate is changed, and no claim below is established. It is placed in
`research/` because it is an unvalidated hypothesis about how the frozen
operator and the validated bioinformatics should be joined.

## 1. What The Frozen Theory Actually Fixes

The frozen deployment is `D = (z_H^0, B(.), Delta_m, mu, h)` and the sole
operator is

```text
A(F,z) = K(B(z)F(z)),        B(z) = [ beta_0(z) | beta_1 | ... | beta_m ]
```

with `beta_0(z) = b_pop[kappa(z)]` a population band chosen by a **finite**
context map, `beta_1..beta_m` **fixed anchor bands**, and `p = F(z) in Delta_m`
the only continuously `z`-dependent object in the entire model.

Three consequences are usually skipped, and all three matter biologically.

**(a) The model is a mixture over a prior band and a fixed anchor ladder.** The
emitted object is `sum_k p_k beta_k`. Nothing else in the output chain varies
with the pair. Whatever biology is supposed to do, it must do it by moving `p`.

**(b) The theorem is about law-class distance, not order.** The controlled
quantity is `d_M = Hausdorff-W1` between emitted law classes.

> `W1(P, P + c) = |c|`, so `W1` is sensitive to common translations, while it may
> also respond to distribution-shape changes. The frozen theory controls
> law-class distances and declared Lipschitz losses but does not automatically
> derive pairwise or listwise affinity ranking.

`W1` is **not** only a location metric, and earlier wording in this document
that described it that way is superseded by the statement above. The scope
chapter separately declines "pairwise, listwise, or metric ranking" and
"derivation of ranking from affinity regression".

**(c) The output is a set of laws, so "I don't know" is expressible.** A wide
band is a valid output, not a failure mode.

## 2. Where The Disconnection Actually Is

The triage file records the interface as "compatible but not deeply integrated".
The three specific breaks:

**Break 1 — the readout measured the one functional the theory disclaims.**
Every affinity result (P1C, P1R1, P1R2A/B, E-AFF-P0, H0A, H0C) was scored by
within-task concordance macro-averaged over closure components. E-AFF-R0 checked
this against the repository's own `metrics.concordance` and found it **exactly**
invariant to any per-task shift or positive rescaling of predictions or labels,
with maximum deviation `0.0` on all four transforms. A simulated predictor that
knows a task's affinity level perfectly and nothing else scores exactly `0.5000`
at every variance share, including when the level holds `98.5%` of total
variance. The instrument assigns the entire task-level channel zero credit by
construction.

**Break 2 — the level channel was also removed upstream, twice.** In H0C the
geometry was shown `y - global_ligand_prior - task_local_ligand_nuisance`, where
the nuisance is fit on 20 labelled support examples **of the correct protein's
own task** ([run_eaff_h0c.py:241-245](run_eaff_h0c.py:241)). The resulting
`local_score` is then added to **both** the correct and the deranged arm
([run_eaff_h0c.py:264-270](run_eaff_h0c.py:264)). Both arms therefore receive
the correct protein's task level for free before the contrast is taken, which is
why "the same head retained nearly all of its value after replacing the
protein". The derangement control was structurally weakened, not merely
uninformative.

**Break 3 — `beta_0` and the anchors were never given biological jobs.** The
former production state used arbitrary bounded latent projections for `z`, and
the anchors carried no declared meaning. So the sign of every affinity increment
had to be **estimated** — a free direction `w in R^288` (E-AFF-P0), or a
per-task direction (H0A/H0C). P0 found no shared direction; H0A found the
per-task version non-transferable. That is the expected outcome when sign is a
parameter rather than a design property.

Note carefully what R0 does **not** say. It does not show that protein-specific
affinity lives in the location channel. It shows that if it does, nothing in the
evidence chain could have detected it.

## 3. Two Claims That Were Conflated

| | Claim A: level | Claim B: interaction |
|---|---|---|
| statement | the correct protein sets the affinity *level* of a chemical series better than ligand-only or sequence-family baselines | protein-by-ligand affinity is non-additive beyond assay noise |
| protein-specific | yes | yes |
| transferable | yes, if predicted from sequence/geometry | yes |
| what the frozen metric sees | this | partly |
| what within-task CI sees | **nothing** | this |
| data requirement | held-out proteins | crossed rectangles, `rho <= 0.0915` (X0-B) |
| status | **never tested** | conditionally testable |

The project has been trying to establish the harder claim with an instrument
blind to the easier one. Claim A is also the one the frozen operator is a
theorem about, and the one a DTA model actually needs.

## 4. The Proposed Integration

### 4.1 Sign becomes a design property, not a parameter

Order the anchors by stochastic dominance on the affinity range `V=[a,b]`. For a
band `beta = (L,U)` with `L <= F_P <= U`, write `beta <= beta'` when
`L' <= L` and `U' <= U` pointwise. Register the deployment with

```text
beta_1 <= beta_2 <= ... <= beta_m.
```

Because `E[X] = b - integral F`, the band's mean interval is

```text
[ b - integral U , b - integral L ],
```

both endpoints affine in the band. Band assembly is linear in `p`, so moving
weight from `beta_j` to `beta_k` with `j < k` **raises both endpoints of the
emitted affinity-mean interval, monotonically**.

The sign of the affinity increment is then fixed by the deployment, and biology
only has to say *how much* weight to move. This is a choice of anchors inside
the frozen `D`; it adds no theorem, strengthens no assumption, and changes no
proof. It removes exactly the object P0 and H0A failed to estimate.

State the limit of that claim explicitly:

```text
sign fixed by deployment convention
    !=
protein-specific affinity information identified from biological input
```

A model on ordered anchors that carries no protein information places weight
without regard to the protein and fails the Gate. The convention removes an
estimation problem; it manufactures no signal.

The repository's frozen `ANCHOR_SPEC` already realizes this ladder: six logistic
anchors at centres `0.15/0.30/0.45/0.60/0.75/0.90` plus one broad uniform
anchor. `audit_l0_operator_contract.py` verified dominance with maximum gap
`0.0`, strictly increasing band mean intervals, mixture monotonicity violation
`0.0`, and all `258` frozen theory files hash-matched. The uniform anchor sits
deliberately outside the order as a width and abstention channel.

### 4.2 The nuisance gets its own declared channel

Set `kappa(z) = (endpoint family, assay format stratum)`, a finite context, and
let `beta_0(z) = b_pop[kappa(z)]` be the population band for that stratum. Assay
and protocol nuisance then lives in `beta_0` **by declaration** instead of
leaking into the mechanism or being differenced away along with the signal.
This is what the context map is for and it has not been used.

### 4.3 `z_bio`: a small, bounded, pair-local mechanism statistic

`p = F(z_bio)` is driven by coordinates that are already validated as
correct-protein-dependent — the P1B contact/distance bridge (correct AUPRC
`0.43885` versus wrong-protein `0.05149`) and the T-BASIS-R0 radial basis
(partner gain `0.1561 [0.1070, 0.2007]`). Candidate coordinates, all normalized
to `[0,1]`, permutation-invariant over residues and ligand atoms, and computed
without any target or assay identifier:

1. interface extent — predicted contact mass, ligand-size normalized;
2. packing density — contact mass in short distance bins over total;
3. apolar complementarity — radial mass on apolar/apolar channels;
4. polar complementarity — radial mass on donor/acceptor channels at H-bond range;
5. charge complementarity — radial mass on charged channel pairs;
6. steric strain proxy — radial mass below contact distance;
7. support coverage — projection of the query onto the support row space.

Coordinate 7 is the identifiability coordinate, not a chemistry coordinate, and
must be declared as such. Ligand-size normalization is a gauge choice and must
be recorded as gauge-dependent rather than given a physical reading. None of
these is free energy and none may be described as such.

### 4.4 Few-shot adaptation becomes dimensionally legal

`p in Delta_m` has exactly `m` free parameters. With `k <= 5` supports the
frozen requirement `d_adapt <= k` becomes simply

```text
m <= k,
```

and `m = 4` leaves a residual degree of freedom. The previous attention-based
support mechanism had no such bound. Rank, conditioning and query coverage are
then reported on an `m`-dimensional problem, which is auditable.

### 4.5 Abstention is free

When the support does not identify `p` — rank deficiency, poor conditioning, or
a query outside the support row space — the correct output is mass on `beta_0`,
i.e. the population band for that assay stratum. That is a wide, valid, honest
law class meaning "prior only". Abstention is therefore a point of `Delta_m`
rather than a bolt-on mechanism, and it is legal because the operator emits a
**set** of laws. The user requirement "unidentified queries must permit
abstention" is satisfied by construction.

## 5. What This Does Not Fix

State plainly:

1. It creates no information. If the P1B geometry carries no protein-specific
   affinity level, the ordered anchors will faithfully report zero.
2. It fixes the sign-estimation problem, not the information problem. P0's
   negative result about a shared 288D direction stands as a statement about
   that parameterization.
3. It does not retroactively convert any past negative into a positive. H0A and
   H0C remain valid conclusions about within-task ranking.
4. It does not address Claim B. The crossed-interaction route remains gated on
   X0-B's `rho <= 0.0915` (Ki) and `rho <= 0.0164` (Kd).
5. A location-sensitive metric is easier to game than a rank metric. The gate in
   section 6 must therefore carry a sequence-only protein baseline, or a family
   shortcut will pass it.

## 6. The Gate That Would Test Claim A

Registered separately as `EAFF_L0_PREREGISTRATION.md`. Its essential features:

- **Held-out proteins.** Evaluation proteins appear in no training closure
  component, so the level must be predicted rather than memorized.
- **Location-sensitive primary metric** in the operator's own terms: the band
  loss `L(B(z)F(z), Y)` required by the frozen theory to be convex, bounded and
  Lipschitz in band sup norm, with `W1` to the observed task law as secondary.
- **Five-arm control ladder**, each arm emitting a band under the same operator:
  1. population band only (`p_0 = 1`);
  2. ligand-only;
  3. **protein-sequence-only** (ESM, no geometry) — defeats "kinases read high";
  4. correct protein with geometry;
  5. deranged protein with geometry.
- **No shared support nuisance across arms.** Each arm fits its own nuisance, or
  none does. This repairs Break 2; without it the derangement contrast is void.
- **Gate:** arm 4 must beat arms 2, 3 and 5, each by the frozen margin with a
  positive 95% lower bound. Beating 5 alone is insufficient, and beating 2 alone
  is insufficient — arm 3 is what makes it a mechanism claim rather than a
  family-level shortcut.
- **Assay stratification** through `kappa`, so a protocol offset cannot be
  mistaken for a protein effect.

Only if that Gate passes does `z_bio` become a candidate for admission to `z`,
and admission still additionally requires the sealed transfer Gate.

## 7. Reading Order For This Package

1. `EAFF_R0_REGISTRATION.md` and `EAFF_R0_RESULT.md` — why the instrument was
   blind.
2. this document — how the operator and the biology join.
3. `EAFF_L0_PREREGISTRATION.md` — the Gate that would decide Claim A.
4. `EAFF_X0_FEAS_RESULT.md` and `EAFF_X0B_RESULT.md` — the separate Claim B route.
