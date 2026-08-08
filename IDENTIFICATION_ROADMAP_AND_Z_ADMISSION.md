# Identification Roadmap and `z` Admission Protocol

Status: **PROPOSAL. Not a registered stage, not an artifact, not hashed.**
Written 2026-08-08 against `AGENT_HANDOFF.md`, `task.md`, `history.md` (F-01..F-92),
`research/e0_identifiability/`, and `model/`.

Nothing here authorizes execution. Every stage below must be separately
preregistered before selection or scoring. Each stage carries a freeze tag:

- **[E]** executable under current freezes, needs only a preregistration;
- **[R]** needs a new registration that changes a frozen quantity;
- **[A]** needs explicit user authorization to break an active freeze.

---

## Part I — One model that contains all four failures

### 1.1 The generative model

Every affinity observation in the governed corpus is a measurement of one cell
`(protein t, ligand l)` inside one assay stratum `s`:

```text
y[s,t,l]  =  mu[s]  +  alpha[t]  +  beta[l]  +  delta[t,l]  +  eps,
                                                 eps ~ N(0, sigma_s^2)
```

- `mu[s]` — assay/protocol offset. Target-independent by the L0 data contract
  (`kappa_raw` excludes accessions, variants and task IDs).
- `alpha[t]` — the protein's affinity **level**.
- `beta[l]` — the ligand's affinity level (promiscuity / general potency).
- `delta[t,l]` — the **non-additive protein-by-ligand interaction**.
- `sigma_assay = 0.47971 [0.47034, 0.48946]` log units, from `4,261` replicate
  cells and `4,840` df (F-91).

### 1.2 The four failures are three estimands, and one is a corollary

| Project failure code | Estimand in the model | Status |
|---|---|---|
| `WITHIN_TASK_RANKING_DIRECTION_NOT_IDENTIFIED` | is the ordering of `beta[l] + delta[t,l]` over `l` dependent on `t`? | **corollary of the next row** |
| `CROSSED_INTERACTION_EXISTENCE_NOT_YET_TESTED` | `Var(delta) > 0` ? | untested, conditionally testable |
| `PROTEIN_SPECIFIC_AFFINITY_LOCATION_NOT_YET_TESTED` | is `alpha[t]` predictable from protein features, given `mu[s]` and `beta[l]` ? | untested |
| `BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z` | interface + licensing question | blocked on the two above |

**The corollary, stated exactly.** Within-task concordance conditions on `t`,
which removes `mu[s]` and `alpha[t]`, and E-AFF-R0 proved it is *exactly*
invariant to per-task shift and positive rescaling (deviation `0.0` on all four
transforms). What survives is the ordering of `beta[l] + delta[t,l]` over `l`.
The `beta[l]` part is identical for every target. Therefore a
**correct-minus-deranged within-task CI contrast is a function of `delta` alone**:

```text
delta == 0   =>   E[ CI(correct) - CI(deranged) ]  =  0    exactly.
```

This is not an analogy. It means P1C, P1R1, P1R2A, P1R2B0, P1R2B1, E-AFF-P0,
H0A and H0C were all **estimating the same quantity**, `Var(delta)`, through a
low-power contrast — and H0A's result is the model's exact prediction:
correct-minus-ligand `+0.08821` (that is `beta`, large and real) alongside
correct-minus-deranged `+0.00864` (that is `delta`, small or absent).

One refinement, in the conservative direction. The deranged arm scores
out-of-distribution inputs through a direction fitted on correct-protein data,
so distribution shift alone can depress it — the F-07 lesson. Any small positive
correct-minus-deranged contrast (`+0.00864` in H0A, `+0.02359` in P1R1,
`+0.02373` in P1R2A) is therefore an **upper** bound on the `delta` signal, not
an estimate of it. The true `Var(delta)` evidence in the ledger is weaker than
the point estimates suggest.

Eight negative results therefore collapse into **one** unreplicated,
underpowered measurement of `Var(delta)`. E-AFF-X0 and X0-FEAS later established
that the governed corpus supports at most `36` Ki and `12` Kd independent
crossed clusters against a registered requirement of `245`. The ledger is
internally consistent: the chain was null because it was ~7x underpowered for
the only quantity it could see.

**Consequence for the roadmap: do not run another within-task ranking
experiment.** Failure 1 is not separately addressable. It is answered by, and
only by, the Claim B test. Running a ninth rank-metric stage would add a ninth
underpowered estimate of the same parameter.

### 1.3 Why the location channel was never measured

Two independent reasons, both now documented:

1. **The instrument was blind.** Within-task CI gives the entire `alpha[t]`
   channel zero credit by construction: a simulated predictor holding a task's
   level perfectly and nothing else scores exactly `0.5000` at every variance
   share up to `0.985` (E-AFF-R0).
2. **The channel was also removed upstream.** In H0C the geometry received
   `y - global_ligand_prior - task_local_ligand_nuisance`, with the nuisance
   fitted on 20 labelled supports of the correct protein's own task and the
   resulting score added to **both** arms. Both arms held the correct protein's
   level before the contrast was taken.

L0 and L0R were the first stages to target `alpha[t]` directly. Neither
produced a protein verdict.

---

## Part II — Two power results that are decidable from published numbers

These are derived from figures already in `EAFF_L0R_RESULT.md` and
`EAFF_X0_FEAS_RESULT.md`. They require no new data and should be re-derived
independently before any further Claim A spend.

### 2.1 The L0R positive control is unattainable at any panel size

The registered precondition is

```text
Delta  =  location_error(A0) - location_error(A1)  >=  0.1 * sigma_assay = 0.04797,
with a 95% closure-component bootstrap lower bound above zero.
```

Observed on `G = 78` components: `Delta_hat = 0.03421`, CI `[-0.03304, 0.10793]`,
so `SE_78 ≈ 0.070485 / 1.96 ≈ 0.03596`.

**Requirement (i), positive lower bound.** `SE(G) = SE_78 * sqrt(78/G)`. Solving
`Delta_hat > 1.96 * SE(G)` gives

```text
sqrt(78/G) < 0.48535   =>   G > 331.
```

**Requirement (ii), point estimate above threshold.** `Delta_hat = 0.03421 <
0.04797`. This is a threshold on the **point estimate**, so it does not improve
with `G`. If the true `Delta` equals the observed value, requirement (ii) fails
with probability tending to one as `G` grows.

**Supply.** X0-FEAS: the governed corpus has `245` closure components, only
`202` carry Ki rows, and `P0` already touched all `245`. At L0R registration
only `78` components still held eligible unconsumed Ki tasks.

```text
required G  >  331          available untouched G = 0
                            available G with task supply <= 78
                            Ki-bearing components in the corpus = 202
```

**Verdict available now:** the L0-family Gate is unattainable on this corpus
under its current estimand — short by a factor of ~4.2 on realised supply and
~1.6 even against the entire Ki-bearing component universe, and blocked outright
by requirement (ii). `L0R_RESULT`'s own first option, "enlarge the panel
substantially", is therefore **not** a viable repair. Its second and third
options — condition within assay strata, or change the estimand — are the only
live ones.

### 2.2 The Gate margin was imported from a different metric and is mis-scaled

`margin_L0 = 0.5 * sigma_assay = 0.23985` log units. Context:

- irreducible location-error floor at perfect level prediction
  `≈ sigma_assay * E|N(0,1)| = 0.47971 * 0.7979 ≈ 0.383`;
- observed arm location errors `1.09552 .. 1.18342`;
- therefore reducible range `≈ 1.13 - 0.38 ≈ 0.75` log units;
- the best available location signal in the corpus, ligand identity itself,
  delivers `0.034`, i.e. **4.5%** of the reducible range;
- the Gate demands `0.240`, i.e. **32%** of the reducible range — roughly **7x
  the strongest signal that exists in the data**.

`+0.03` was calibrated for a concordance-index contrast. Carrying `0.5 *
sigma_assay` across to an absolute-location error in log units re-uses the
*form* of the old margin without re-deriving its *scale*. A margin should be
derived from a decision-relevant requirement — what location accuracy makes the
emitted band useful — not inherited by analogy. **[R]**

---

## Part III — Roadmap

Ordering principle unchanged from the project's practice: stop at the earliest
unresolved boundary, and prefer a cheap label-free audit to an expensive panel.
New principle added: **every stage must be able to terminate**, so each Gate
below carries an equivalence arm as well as a superiority arm.

### Stage L1-FEAS — attainability audit for Claim A **[E]**

Purpose: do to L0 what X0-FEAS did to X0 — decide whether *any* location Gate is
attainable before consuming another panel. Cheap; mostly re-analysis.

1. Re-derive §2.1 and §2.2 independently.
2. **Variance decomposition of governed Ki `p_affinity`** into
   between-stratum `Var(mu[s])`, between-component, between-task-within-
   component, and within-task-between-ligand. Reads labels; fits no model;
   produces no arm contrast. This single number decides everything downstream:
   it is the fraction of location variance that `kappa` conditioning removes.
3. Label-blind component census over the **full** governed D0 closure
   (`37,783` tasks, `4,787` proteins, `459` components per X0-FEAS) for Ki
   components with `>= 20`-ligand tasks surviving D1 protected-DAVIS exclusion,
   partitioned into touched / untouched.
4. Recompute C1/C2/C3 of the L0 data contract **within** the stratified estimand.

Terminal verdicts: `A_GATE_ATTAINABLE_UNDER_STRATIFIED_ESTIMAND`,
`A_GATE_ATTAINABLE_ONLY_WITH_NEW_CORPUS`, `A_GATE_UNATTAINABLE`.

### Stage L1 — stratified location Gate **[R]** (only if L1-FEAS is green)

Four changes from L0R, each named and justified in advance.

**(a) Use `kappa`. This is the largest single defect in L0/L0R.**
`EAFF_L0_DATA_CONTRACT.md` defines the target-independent stratum

```text
kappa_raw = (assay_organism, bao_format, cell_id, tissue_id,
             subcellular_fraction, relationship_type, assay_parameters)
```

and verified it identifiable for Ki (`C1=79`, `C2=218`, `C3=0.630`). The runner
then ignored it: `run_eaff_l0r.py:219` computes a **single pooled** population
band, `pop = population_band(y_scaled[train], grid)`, for every stratum. The
frozen operator's own nuisance channel `beta_0(z) = b_pop[kappa(z)]` — the thing
`THEORY_BIOLOGY_INTEGRATION.md` §4.2 says "has not been used" — was not used.
All between-stratum location variance was left in the residual, competing with a
ligand signal worth `0.034` log units.

L1 must build `b_pop[kappa]` per stratum and evaluate **within-stratum**
location. This is using the frozen operator as designed, not modifying it.

**(b) Replace the positive control with an oracle injection.** L0R's control
asks the *ligand* to predict *absolute cross-target* location. Ligands do not do
that — a promiscuous inhibitor has a different absolute affinity against every
target — so the control conflates "the readout cannot detect location" with "the
ligand carries no location information". Register instead a synthetic control in
the project's own established style (the E0 synthetic trainability pre-gate):
inject a known per-task offset of size `c * sigma_assay` into a copy of the
labels and require the pipeline to recover it at the Gate margin, for
`c in {0.25, 0.5, 1.0}`. This cannot fail for biological reasons, and it
calibrates the readout's sensitivity in the same units as the Gate. Retain the
ligand contrast as a reported diagnostic, not a precondition.

**(c) Re-derive the margin** from a stated decision requirement, per §2.2.

**(d) Preserve the five-arm ladder** `A0..A4` unchanged, including the
sequence-only arm `A2` — without it a family-level shortcut ("kinases read
high") passes a location Gate trivially — and the rule that no nuisance is
shared across arms.

Terminal verdicts must include an **equivalence arm**:
`PROTEIN_LOCATION_IDENTIFIED_IN_SOURCE`,
`PROTEIN_LOCATION_ABSENT_AT_REGISTERED_SENSITIVITY` (upper bound below margin),
`INDETERMINATE`, plus the existing not-run codes.

### Stage X1-RHO — intra-cluster correlation, with abstention **[E]**

Mandated by F-90 and currently the binding precondition on all of Claim B.
Estimate `rho` for the double difference from exact-assay replicates; require

```text
UCB_95(rho)  <  rho*        rho*(Ki) = 0.0915,   rho*(Kd) = 0.0164.
```

Two cautions to register in advance. First, `rho` is estimated from `36` (Ki) or
`12` (Kd) clusters, so its own interval is wide; use a cluster-jackknife or REML
profile interval and **abstain on a wide interval, not just a high point
estimate**. Second, `rho*(Kd) = 0.0164` against `12` clusters is a very narrow
target; a Kd abstention should be treated as the expected outcome, not a
surprise.

### Stage X1 — crossed interaction variance test **[R]**

Only if X1-RHO passes. Unit: the cell-disjoint rectangle established by X0-B.

```text
DD  =  y[t1,l1] - y[t1,l2] - y[t2,l1] + y[t2,l2]
```

`DD` cancels `mu[s]`, `alpha[t]` and `beta[l]` exactly, leaving
`4 * delta`-contrast plus `4 * sigma^2` of noise. Three methodological upgrades
over the registered `chi^2` design:

1. **Studentize before aggregating.** `sigma_assay` almost certainly varies by
   stratum. A pooled `chi^2` against a single `sigma_assay` is anti-conservative
   under heteroscedasticity. Divide each `DD` by a cell-local replicate SD.
2. **Use a replicate-resampled null, not a parametric `chi^2`.** Build the null
   by resampling within-cell replicates. This carries the true noise law,
   including heteroscedasticity and any residual non-normality, and removes the
   test's sensitivity to the `+-2%` uncertainty in `sigma_assay` (which
   propagates to `+-4%` against a `25%` target variance ratio).
3. **Cluster bootstrap over the `36`/`12` clusters** for all intervals, never
   over rectangles.

Register the censoring limitation explicitly: D0-C admits exact-point Ki/Kd only,
so weak binders are excluded, which **truncates the affinity range and biases
`Var(delta)` downward**. A positive result is therefore conservative; a null is
not.

Terminal verdicts, with an equivalence arm so the line can close:
`INTERACTION_IDENTIFIED`,
`INTERACTION_ABSENT_AT_REGISTERED_SENSITIVITY` (`UCB(Var(delta)/sigma^2) < 0.25`),
`INDETERMINATE`, `NOT_RUN_RHO_ABSTAINED`.

### Stage X1-D — the crossed corpus the project already owns **[A]**

X0-FEAS established something structural: crossing and document-disjointness are
produced by opposite kinds of study, so no ChEMBL-shaped corpus will ever supply
clean crossed units. A genuinely crossed corpus does exist in this repository.

**DAVIS is a complete target-by-ligand matrix measured on one platform.** Its
source split alone is `220` targets crossed against the full compound panel, one
protocol, one document — zero assay-context confounding, and rectangles
available by construction rather than by lucky co-occurrence. It is the exact
instrument X0 could not find, and Claim B is an *existence* question that a
single well-crossed panel can answer definitively.

This is the highest-information experiment available for Claim B, and it is
blocked by governance rather than by data. Registering it requires an explicit,
one-way decision to move `DAVIS_LABEL_READS` off zero **for the source split
only**, with metaval and recipient sealing preserved and enforced by the existing
UID/permission mechanism (F-21). Two limitations must be registered with it:
DAVIS Kd is heavily censored at the assay ceiling, and a single platform bounds
external validity to "interaction exists", not "interaction is transferable".

I recommend this be put to an explicit decision rather than left implicit. It is
the only route on the board with a realistic chance of a decisive Claim B answer.

### Stage X1-B — independent replication **[R]**

`D0-B` (BindingDB) has been carried as separate and unexecuted since F-68. If
either L1 or X1 produces a positive, it must replicate on an independent source
before any admission. Registering the acquisition now, while the primary stages
run, avoids a repeat of the F-43 pattern where the replication corpus was
discovered to be insufficient only after it was needed.

---

## Part IV — The `z` admission protocol

### 4.0 What the frozen contract actually requires

From `model/config.py`: the frozen objects are `V`, mesh `h`, `Delta_m`, the
positive ridge `mu`, the band polytope, the assembly `B(z)` and the operator
`A(F,z) = K(B(z)F(z))`. `d_z` and the meaning of its blocks are **declared
engineering choices**. The constraint on `z` is that it be bounded,
finite-dimensional, deterministic, permutation-invariant in `S`, and free of
query labels. Its dimension is not fixed.

Admission therefore has two independent halves: a **mechanical conformance**
half that is testable today, and a **scientific licensing** half that is blocked
on Part III.

### 4.1 Gate Z0 — mechanical conformance **[E], runnable now**

All twelve are label-free and can be implemented and tested before any
scientific result exists.

| # | Check | Pass condition |
|---|---|---|
| Z0.1 | boundedness | `z in [0,1]^d_z` on a registered stress set; reuse `assert_bounded_observable` |
| Z0.2 | determinism | bitwise-identical `z` across 3 runs and 2 processes |
| Z0.3 | support permutation invariance | `max abs` deviation `<= 1e-12` over 100 random permutations of `S` |
| Z0.4 | query-label freedom | `z` unchanged under permutation and corruption of query labels |
| Z0.5 | atom/residue permutation invariance | `<= 1e-6`; this is the exact defect that produced `P1C_READOUT_NOT_PERMUTATION_INVARIANT` |
| Z0.6 | no identifier leakage | `z` unchanged under target-ID and assay-ID permutation; an ID-only probe must not reconstruct `z` |
| Z0.7 | `kappa` conformance | finite codomain, every context reachable in source, reads only declared context coordinates |
| Z0.8 | view-bank registration | `DEFAULT_VIEWS` re-registered against **named** biological coordinates |
| Z0.9 | deployment manifest | `validate_deployment_manifest` passes and **fails closed** under a mutation test on every bound hash |
| Z0.10 | operator integrity | re-run `audit_l0_operator_contract.py`: anchors in polytope, dominance gap `0.0`, strictly increasing band mean intervals, mixture monotonicity violation `0.0`, Hausdorff-`W1` stability violation `0.0`, `M=32`, all `258` theory files hash-matched |
| Z0.11 | few-shot dimensional legality | `d_adapt <= k` proven, not asserted |
| Z0.12 | abstention reachability | out-of-support queries place `>=` registered mass on `beta_0` |

Three of these currently **fail or cannot be evaluated**, and all three are
repairs that should happen regardless of how Part III resolves.

**Z0.7 is broken twice over.**

*First, the code path is dead.* `model/meta_operator.py:648`, inside
`build_band_operator`, executes `from .biological import kappa`.
`model/biological.py` was deleted in the F-79 consolidation. The function that
constructs `B(z)` for any admitted statistic will raise `ModuleNotFoundError` on
its first call. It is inside a function body, so import-time is clean and the
`70`-test suite does not reach it.

*Second, there are two incompatible definitions of `kappa` in the project.*

```text
production   context_index(z, cfg)  buckets z[:,12], z[:,16], z[:,26]
             = (mean support label bin) x (mass bin) x (unused continuous bin)
             n_context = 3 * 2 * 1 = 6

L0 contract  kappa_raw = assay metadata, target-independent by construction
```

These are different objects. The production `kappa` reads the **support's own
labels**, which is coherent for a few-shot prior and incoherent for a
cross-fitted location estimand. `config.py` already anticipates the fix:
`kappa_edges_context_cont` is documented as a "transitional pre-P4 context
quantizer" with one reachable bin so that "continuous assay covariates cannot
alter `B`". Admission should activate exactly that third factor, mapping
`kappa_raw` into a registered finite codomain and giving
`n_context = 3 * 2 * |strata bins|`.

**Z0.11 currently fails.** `m = 7`, so `p in Delta_7` carries 7 free parameters,
while few-shot deployment has `k <= 5` supports. `THEORY_BIOLOGY_INTEGRATION.md`
§4.4 requires `m <= k`. Rather than change `m` — which would discard the frozen
seven-anchor ladder that `audit_l0_operator_contract.py` has already verified —
use the ladder's own structure. Because the six logistic anchors are
**stochastically ordered**, restrict few-shot adaptation to a two-parameter
family:

```text
(ladder position tau in [0,1],  abstention mass p_0 in [0,1])
```

with the full `Delta_7` reachable only by the source-trained population model.
Then `d_adapt = 2 <= k = 5` with margin, and §4.1's monotonicity result
guarantees that moving `tau` moves both endpoints of the emitted affinity-mean
interval in a fixed direction. This makes the sign of the affinity increment a
property of the deployment rather than a parameter to be estimated — which is
precisely the object E-AFF-P0 and H0A failed to estimate.

**Z0.12 has no implementation.** The abstention channel is the seventh proposed
`z_bio` coordinate (support row-space coverage). It is the coordinate that makes
"I don't know" expressible, and it is the one coordinate that is *not* a
chemistry claim and therefore cannot be invalidated by a biological null. It
should be built and tested first.

### 4.2 Gate Z1 — scientific licensing

Admission is **claim-specific**. What a passing stage licenses:

| Passed | `z_bio` may supply | It may **not** |
|---|---|---|
| L1 only | location coordinates that move mass along the ordered ladder | claim within-task ranking, or any interaction semantics |
| X1 only | interaction coordinates driving ligand-dependent band movement within a target | claim absolute level transfer to unseen targets |
| both | both | claim free energy, or any physical thermodynamic reading |
| neither | nothing | `z` stays abstract; `d_z = 28` remains a declared placeholder |

The existing `task.md` criteria are retained unchanged as necessary conditions:
correct protein beats ligand-only by the frozen margin, and beats deranged
protein with a positive lower bound.

### 4.3 Gates Z2–Z4 — the rest of the chain

- **Z2 independent-source replication.** The Z1 result reproduces on BindingDB
  (`D0-B`) under the same governance, with the effect size inside the primary
  study's interval. **[R]**
- **Z3 sealed novel-target transfer.** One evaluation on DAVIS metaval, held-out
  targets, recipient still sealed, `RECIPIENT_LABEL_READS = 0` enforced by the
  F-21 permission mechanism, not by declaration. One shot, preregistered. **[A]**
- **Z4 operator-integrity re-audit and manifest binding.** Re-run Z0 in full
  against the final artifact; bind `frontend_hash`, `state_schema_hash`,
  `view_registry_hash`, `context_registry_hash`, `mechanism_schema_hash`, the
  protein/ligand/pair bank hashes and `B_table_hash` through
  `deployment_manifest`; verify the theory ledger still hash-matches. **[E]**

### 4.4 The admission statement

`z_bio` is admitted only when a single document records: the Z0 conformance
report; the Z1 licensing row that was earned; the Z2 replication; the Z3 sealed
transfer result; the Z4 manifest; and an explicit, signed scope sentence naming
what the admitted statistic does **not** establish. Absent any one of these,
`BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z` stands.

---

## Part V — Honest terminal states

A roadmap that can only succeed is not a roadmap. Three of these are real
possible endings, and each is a legitimate result:

1. **`Var(delta) ≈ 0` in accessible public data.** Then within-task ranking of
   ligands is target-independent to measurement precision, eight historical
   negatives are explained rather than merely recorded, and the correct product
   is a ligand-potency model with a protein-level offset — not a pair model.
   The frozen operator is well suited to exactly that object.
2. **Claim A is unmeasurable on ChEMBL-shaped data.** §2.1 already points here.
   Location requires many proteins measured in shared strata; ChEMBL's document
   structure fights it, in the same way X0-FEAS showed it fights crossing.
3. **Both claims stay indeterminate.** Then the publishable result is the
   governance and negative-evidence chain itself: `92` documented failures with
   preregistration, independent audit and sealed recipients is a stronger
   methodological contribution than most positive DTA papers, and §1.2 gives it
   a single unifying explanation rather than a list.

The roadmap's job is to reach one of these quickly and cheaply, not to defer
them. That is why L1-FEAS and X1-RHO come first: both are cheap, both are
label-light, and either can terminate a whole line before a panel is spent.

---

## Part VI — Recommended immediate order

| # | Stage | Freeze | Cost | Can terminate a line? |
|---|---|---|---|---|
| 1 | Z0.7 repair — dead `kappa` import, two-`kappa` reconciliation | [E] | hours | no, but unblocks everything |
| 2 | L1-FEAS attainability audit | [E] | ~1 day | **yes** — can close Claim A |
| 3 | X1-RHO with abstention | [E] | ~1 day | **yes** — can close Claim B |
| 4 | Z0 conformance harness (12 checks) | [E] | ~2 days | no |
| 5 | Decision on X1-D (DAVIS source) | [A] | user decision | **yes** — can decide Claim B |
| 6 | L1 or X1, whichever survives 2/3/5 | [R] | panel | yes |

Items 1–4 are executable under every current freeze and cost less than a week
combined. Items 2 and 3 between them can close both scientific questions before
any further panel is consumed.

---

## Part VII — What this document does not do

- It creates no information and reverses no verdict. Every historical FAIL,
  CLOSED, NOT-RUN and negative stands exactly as recorded.
- §1.2 is a re-derivation of what the existing results were measuring. It is not
  a reinterpretation of any of them as positive.
- §2.1 and §2.2 are arithmetic on published numbers and should be independently
  reproduced before being relied on; they are the kind of claim that deserves an
  adversarial check.
- The `m <= k` resolution in §4.1 is a proposal about the adaptation family. It
  needs the monotonicity argument of `THEORY_BIOLOGY_INTEGRATION.md` §4.1 written
  out as a proof obligation, not merely cited.
- Nothing here has been executed. No label was read to produce it.
