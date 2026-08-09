# Phase 2B design audit — read-only, pre-execution

Date: 2026-08-10. Repository commit at audit: `0bd1702`.

Subject:
`research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL.md`
SHA-256 `ae6d1a0186bb37af86f3b6eb98c513bce7e67a8745aaf5a3811ce5c9b98ab477`.

```text
STATUS   SUPERSEDED_BEFORE_EXECUTION_DESIGN_DEFECT
RUNS     none — no Phase 2B code was ever written or executed
FILE     kept byte-identical; not edited
```

The **scientific direction** of that document is retained: one frozen
generic-pocket prior plus one small ligand-conditioned residue residual,
supervised by the same-protein ligand differential. The **experimental
contract** contains defects that would have made the result uninterpretable or
unfalsifiable. Each is stated below with the exact clause at fault.

Superseded by `PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md`.

---

## D1 — the prior was ligand-dependent and did not exist on the training split

Clause, §2:

> `b_r(P) := the residue-marginal component alpha_r of the weighted additive
> projection of the sealed B5 logits, as computed in Phase 2A`

Two failures.

**It is not a function of `P` alone.** The Phase 2A additive projection was
fitted *per complex*, i.e. per `(P, L)` pair, on the matrix `G_ra(P,L)`. Its
residue term `alpha_r` therefore absorbs ligand-dependent structure. Writing it
`b_r(P)` asserts a ligand-independence the object does not have, and the
same-protein cancellation the whole differential design relies on would not have
held exactly.

**It does not exist where it is needed.** Phase 2A computed that decomposition
only on held-out A. Phase 2B trains on the training components, where the
quantity was never materialised.

*Repair.* Use the frozen B5 protein-only branch directly:
`b_r^P(P) = b + alpha * w_pi(GELU(W_h(LN(h_r))))`. Every parameter comes from
the frozen B5 checkpoint; the expression contains no atom or ligand term, so
ligand independence is structural rather than asserted. Materialise and hash it
for train and held-out before training.

## D2 — the head was underspecified

Clause, §2: `K <= 8`, and `g(L)` described only as "a permutation-invariant
pooling of the existing 41-dimensional atom features".

`K <= 8` is a family, not a model. "A permutation-invariant pooling" admits
mean, max, sum and attention-free variants that differ materially. A
preregistration that does not pin these leaves a post-hoc choice open.

*Repair.* Freeze `K = 8` exactly, `g(L) = ` the arithmetic mean over the
existing deterministic 41-D atom features, no bias in `U` or `V`, Xavier
initialisation, one fixed seed, `10,568` trainable parameters and no others.

## D3 — the projection span was rank-deficient by construction

Clause, §2.1:

> `delta <- delta - Proj_span{ 1 , b(P) , c(L) * 1 } delta`

`1` and `c(L) * 1` are the same direction in residue space for any scalar
`c(L)`. The span is therefore rank 2, not 3, and the stated projector is
singular. Any implementation would have had to silently resolve the degeneracy,
and the registered orthogonality tolerance would have been checked against an
ill-defined basis.

*Repair.* Build `Q_P` as an explicit orthonormal basis of `{1, b^P(P)}` by
modified Gram–Schmidt in float64, dropping any column whose residual norm falls
below a frozen tolerance — which also handles a constant or near-constant
`b^P(P)`. The projector is protein-fixed, so it commutes with the same-protein
difference and cannot leak ligand information.

## D4 — no consistent score, and no cancellation check

The document defines `logit p_r = b_r(P) + delta_r(P,L)` but never states that
`b` must cancel in the same-protein difference, and never requires that to be
verified numerically. With the D1 prior it would not have cancelled.

*Repair.* Define `s_r(P,L) = b_r^P(P) + delta_r(P,L)` once, derive
`Delta s = s(P,L_a) - s(P,L_b)`, and add a fail-closed numerical check that
`b^P` cancels to machine precision.

## D5 — the primary evaluation was selected by the labels

Clause, §5:

> rank the **symmetric-difference residues** by `Delta delta_r`

The candidate set is defined using the answer. A model is only ever asked to
sort residues that are already known to change, so the metric cannot reward
finding *which* residues change — the hardest and the only deployable part of
the task. It also silently conditions away the majority of the protein.

*Repair.* Score over **all aligned residues**. Report `AP_gain` (score
`Delta s`), `AP_loss` (score `-Delta s`) and `AP_change` (score `|Delta s|`),
each against a per-pair chance level, aggregating
residue → unordered ligand pair → construct → closure component → macro. The
symmetric-difference-only quantity is demoted to a secondary *conditional sign*
diagnostic and may never be quoted as deployment evidence.

## D6 — the non-inferiority arm compared incomparable objects

Clause, §7, gate `D5`:

> full pair AP of `b + delta` must not fall below sealed B5 by more than 0.005

`b + delta` is a residue-indexed vector. Sealed B5 AP is computed on the
complete residue × atom matrix. The two are not the same object, and the
comparison has no defined value.

*Repair.* Define the auxiliary pair score
`G_2B_ra = G_B5_ra + delta_r` with `delta_r` broadcast along the atom axis, so
the frozen B5 atom and pair branches are retained and the gate actually tests
whether the new residue term damages existing pair localisation.

## D7 — the module-participation audit contained a mathematical no-op

Clause, §8:

> gradient blocking: detaching `h_r` must change the outcome measurably

`h_r` is a frozen input with no upstream parameters. Detaching it changes
nothing — no gradient ever flows into it. The requirement can never be
satisfied, so a literal reading would fail every honest run.

*Repair.* Replace it with checks that are actually informative for this
architecture: nonzero `U`/`V` gradient norms, a minimum relative parameter
movement fixed in advance, non-degenerate activation variance, ablation of
`U h_r` and of `V g(L)` each collapsing the differential, residue-context
shuffle and ligand shuffle each degrading it, same-seed determinism, and
synthetic recovery. Only **output-level** `delta` claims are permitted: the
bilinear factorisation has a rotation gauge (`U -> RU`, `V -> RV` leaves
`delta` unchanged for orthogonal `R`), so individual `U`/`V` channels are not
interpretable.

## D8 — the "capacity-matched random head" was untrained

Clause, §7, gate `D3`: "capacity-matched random head (identical parameter count,
frozen random `U`, `V`)".

An untrained head is matched in parameter count but not in *pipeline*. It cannot
distinguish a learned biological signal from anything the optimiser, sampler,
aggregation or evaluation code can manufacture from structure alone.

*Repair.* The capacity control becomes a head with the identical architecture,
optimiser, sampler and budget, **trained** on one frozen within-construct ligand
label permutation. The untrained arm is retained only as a numerical sanity
check and renamed `RANDOM_FEATURE_NULL`.

## D9 — the wrong-ligand control replaced only one side of the pair

Clause, §7, gate `D2`: "`delta` evaluated with a foreign ligand substituted for
`L_a`".

The differential is `delta(L_a) - delta(L_b)`. Substituting only `L_a` leaves the
true `L_b` in the expression, so the control still contains genuine information
about the pair and understates the shortcut it is meant to expose.

*Repair.* Replace **both** ligands with foreign ones, matched by a frozen,
score-blind ligand-size rule, with no fixed points and one map shared by every
compared arm.

## D10 — no sampler, and ordered pairs counted as independent

Clause, §5 speaks of "each held-out **ordered** pair `(L_a, L_b)`", and no
training sampler is specified anywhere.

`(L_a, L_b)` and `(L_b, L_a)` are the same physical comparison; scoring both as
independent observations inflates the effective sample size by two. And with no
sampler, a construct carrying 318 ligands contributes ~50,000 pairs while a
construct with 2 contributes 1 — a single protein family would dominate both the
gradient and the estimate.

*Repair.* Unordered pairs only, with gain and loss scored *inside* the pair.
A deterministic hierarchical sampler enforcing
residue → pair → construct → component balance, with its selected pair IDs
materialised and hashed.

## D11 — "replicate oracle ceiling" was misnamed and over-extrapolated

Clause, §6 calls the replicate agreement a "ceiling" and requires Phase 2B to
report results "as a fraction of that ceiling".

Replicate agreement is not a mathematical ceiling. A model that denoises
annotation error can legitimately exceed the agreement between two noisy
annotations of the same system. Further, the replicate structure exists on only
27 held-out components; extrapolating it to the full held-out panel would be
unsupported.

*Repair.* Rename to `REPLICATE_REPRODUCIBILITY_REFERENCE`, compute it only on
the matched subset where the replicate structure actually exists, state the
subset size, and forbid it from determining any PASS.

---

## What is unchanged

The scientific hypothesis, the frozen-prior-plus-small-residual architecture,
the same-protein differential supervision, the closure-component inference unit,
the prohibition on adding any PLM, attention stack, GNN, geometry branch,
typed-interaction branch, affinity head, PU loss, knowledge graph or parallel
module, and every frozen boundary on affinity, DAVIS, KIBA, recipient data,
few-shot adaptation, production `z`, CSMO, Band, mesh and
`A(F,z) = K(B(z)F(z))`.

## Disposition

`PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL.md` is retained byte-identical and
marked `SUPERSEDED_BEFORE_EXECUTION_DESIGN_DEFECT` in
`report/s7_l2b_r0r/PHASE1_ARTIFACT_SUPERSESSION.json`. It produced no result,
so nothing is withdrawn — only a design is replaced before it was used.
