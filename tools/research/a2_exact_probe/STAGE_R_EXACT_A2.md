# Stage R: the exact episodic A2 operator, tested on real episodes

Authority: `A2_EXACT_meta_val.json`. Frozen A0 (seed 20260815), double-cold
`meta_val` scored **once** after the rank was fixed on `meta_train` component
folds. No trunk parameter trains.

> **Re-run 2026-08-16 after two defects in the first run's episode banks.**
> `SUPERSEDED_A2_EXACT_meta_val.json` retains the earlier numbers as evidence
> of what changed; every figure below comes from the repaired run.
>
> 1. **The banks were not reproducible.** Episode draws were seeded with
>    Python's built-in `hash()` over a tuple containing a string, which is
>    salted per process, so two runs of the "fixed" bank disagreed. Replaced
>    with `scripts.qpsmp_data.stable_seed` (sha256), keyed on the target
>    *name*. `tests/test_bank_reproducibility.py` now launches separate
>    interpreters under hostile `PYTHONHASHSEED` values and requires bitwise
>    equality of episode identities.
> 2. **The banks were not nested.** The query panel was sliced as
>    `order[k : k + query_size]`, so it *moved* with the support size and the
>    k-curve compared different populations. Repaired to the nested contract:
>    one panel per (target, draw) taken after the largest support, shared
>    across k, with each k receiving a prefix of the same support ordering.
>
> **What this changed.** The conclusion is unchanged and the performance gates
> are *stronger*, but two control contrasts lost their resolution and the
> earlier report's claim that corrupting the protein or the labels made the
> operator **resolvedly** better is **withdrawn** — see §"The control gates,
> corrected". Within any single k the first run was still internally paired,
> so its per-k arm ranking stood; its k-curve did not.

## Why this stage was necessary

The v2 representation probe rejected A2 on the strength of a *zero-shot
bilinear delta predictor* — `Δŷ = s·⟨g(P), U(e_i − e_j)⟩`. That shares A2's
feature space but not its structure: it forms no support moment, it reads a
ligand pair rather than a query against an episode, and it has no shrinkage in
k. **Generalising from it to A2 was an over-reach**, and this stage removes the
inference by implementing the operator the plan actually specifies:

```text
z_i  = A_φ(e0(P, L_i))            r_i = stopgrad(y_i − f0(P, L_i))
c_S  = (1/k) Σ_i r_i z_i          η(k) = η_∞ · k/(k + λ)
δ_q  = η(k) · ⟨c_S, z_q⟩          f = f0 + δ
```

Trainable: `A_φ` (one `Linear(D, R, bias=False)`), `η_∞`, `λ`. Nothing else.

## Structural gates — all pass

19 algebraic probes, `tests/test_operator_contract.py`, on random tensors in
float64, at arbitrary parameter values:

| gate | result |
|---|---|
| exact k=0 identity | ✅ `δ ≡ 0` bit-exactly, and survives `η→e⁵⁰`, `λ→e⁻⁵⁰`, `‖A_φ‖×10⁶` |
| shrinkage zero at k=0, increasing in k | ✅ |
| **non-scalar k=1 correction** | ✅ `δ_q = η(1)·r₁·⟨z₁, z_q⟩` varies with `q` |
| support permutation invariance | ✅ |
| query permutation equivariance | ✅ including a single query alone |
| linear and odd in the support residuals | ✅ `δ(−r) = −δ(r)`, `δ(2r) = 2δ(r)` |
| zero residuals give zero correction | ✅ |
| no query-label path | ✅ `forward` takes three arguments, none a query label |
| residuals carry no gradient | ✅ |
| the scalar baseline really is scalar at every k | ✅ range exactly 0 |

**The k=1 gate is the one A2 has that A0 does not.** A0's k=1 transport is
provably a pure level shift (`sar_adaptation ≡ 0`, DATAFLOW_AUDIT F4). The exact
A2 operator is structurally capable of a query-specific k=1 correction. So the
family deserved this test, and the v2 rejection was premature.

## Real-data result — every performance and control gate fails

Population: double-cold `meta_val`, 41 targets / 19 components, 4 episode draws
per target, query panel 16, k ∈ {0,1,2,3,5}. Equal-component, equal-target
weighting. Metrics on the full prediction `f0 + δ` in pK.

### The preregistered gates, as paired component bootstraps

Episodes are identical across arms (same nested bank, same draws, same seeds),
so each contrast is paired per (target, draw). Positive = `a2_embed` better;
the quantity is MSE, so the sign is flipped to make "better" positive.

| gate | k=1 | k=2 | k=3 | k=5 | verdict |
|---|---:|---:|---:|---:|---|
| **beats Tanimoto transport** | **−0.042** | **−0.173** | **−0.217** | **−0.252** | **FAIL, resolved at every k** |
| beats the scalar level baseline | −0.140 | −0.118 | −0.117 | **−0.134** | **FAIL** at k=5; nominally worse everywhere |
| beats the shared-moment control | −0.136 | −0.118 | −0.119 | **−0.136** | **FAIL** at k=5 |
| beats a frozen random projection | −0.066 | +0.192 | **+0.355** | **+0.515** | PASS at k=3,5 |
| degrades under a wrong protein | +0.040 | +0.022 | +0.019 | +0.008 | **unresolved** |
| degrades under shuffled labels | −0.146 | −0.117 | −0.112 | −0.129 | **unresolved** |

### The control gates, corrected

The superseded run reported the wrong-protein and label-shuffle contrasts as
**resolved and inverted** — corrupting the input made the operator measurably
better. Under the repaired nested banks both are **unresolved**. That claim is
withdrawn.

What the repaired evidence supports is weaker and still sufficient:

* there is **no measurable dependence on the correct protein** (the interval
  spans zero at every k), and
* there is **no measurable dependence on correct support labels** (likewise) —
  the point estimates favour the shuffled arm by ~0.12 pK², but the intervals
  include zero, so the honest statement is *absence of evidence for label
  binding*, not *evidence that shuffling helps*.

A mechanism whose entire purpose is to exploit support labels under the correct
protein, and which shows no resolvable dependence on either, has failed its
control gates — but by being inert, not by being inverted.

### MSE (pK²), lower is better

| arm | params | k=0 | k=1 | k=2 | k=3 | k=5 |
|---|---:|---:|---:|---:|---:|---:|
| `tanimoto` (parameter-free) | 0 | 2.1608 | 1.5566 | **1.2974** | **1.1703** | **1.0341** |
| `scalar_level` (2 scalars) | 2 | 2.1608 | **1.4582** | 1.3529 | 1.2695 | 1.1520 |
| `a2_max_state` | 770 | 2.1608 | 1.4710 | 1.3518 | 1.2590 | 1.1412 |
| `a2_label_shuffled` | 770 | 2.1608 | 1.4519 | 1.3535 | 1.2744 | 1.1575 |
| `shared_moment` | 778 | 2.1608 | 1.4624 | 1.3529 | 1.2680 | 1.1503 |
| **`a2_embed`** — the candidate | 770 | 2.1608 | 1.5982 | 1.4704 | 1.3868 | **1.2861** |
| `a2_wrong_protein` | 770 | 2.2076 | 1.6384 | 1.4919 | 1.4054 | 1.2942 |
| `a2_random_projection` | 2 | 2.1608 | 1.5321 | 1.6626 | 1.7419 | 1.8012 |
| `a2_random_feature` | 770 | 2.1608 | 2.5455 | 2.4508 | 2.4373 | 2.4219 |

k=0 is **2.1608** for every correct-protein arm — the exact k=0 identity holds
numerically, and the value sits within retraining noise of A0's recorded 2.1488
on its own bank.

### The verdict, read off that table

1. **A2 loses to a parameter-free kernel at every k, resolved.** 1.2861 vs
   1.0341 at k=5 — 24% worse than fixed Morgan/Tanimoto transport.
2. **A2 loses to two trainable scalars.** 1.2861 vs 1.1520.
3. **Removing the support chemistry improves it.** `shared_moment` uses only
   the mean residual and a learned direction: 1.1503.
4. **No resolvable dependence on the protein or on correct labels.**

### Why it fails — measured, not inferred

The `query_spread_pk` column is the diagnosis. It is the standard deviation of
`δ` across the queries of one episode: the entire query-specific content of the
correction.

| arm | k=1 | k=5 | against a label spread of 0.884 pK |
|---|---:|---:|---|
| `scalar_level` | 0.0000 | 0.0000 | scalar by construction |
| `tanimoto` | 0.0000 | **0.2865** | 32% of the label spread |
| **`a2_embed`** | **0.0021** | **0.0027** | **0.3% of the label spread** |
| `a2_max_state` | 0.0113 | 0.0152 | 1.7% |
| `a2_random_feature` | 0.4613 | 0.4521 | 51% — the operator *can* do this |

The exact A2 operator is structurally query-specific and **numerically almost
constant**: 0.0027 pK of query variation against Tanimoto's 0.2865 pK, a factor
of 106. It has degenerated into a level shift — and a worse-calibrated one than
the plain shrunken mean, which is why it loses to `scalar_level`.

The `a2_random_feature` row rules out the obvious alternative explanation. Given
Gaussian noise features the same operator produces 0.45 pK of query spread and
a catastrophic MSE (2.42 at k=5). The mechanism is not dead, the shrinkage did not collapse,
and the optimiser is not stuck. **Trained on the real representation, gradient
descent drives the operator toward constancy because the moment carries nothing
usable.** That is the strongest available statement that `e0`'s support moment
is uninformative: the operator turned itself off, on its own, on held-out
components.

`max_state` — the representation Phase 3 measured as retaining the most
protein-differential of any stage — behaves the same way: MSE 1.1412 at k=5, above
Tanimoto's 1.0341 and barely distinguishable from the two-scalar baseline.

## Rank selection

Chosen on `meta_train` component folds, before `meta_val` was read. The curve is
flat, so no rank rescues the operator:

| rank | `embed` fold MSE | `max_state` fold MSE |
|---:|---:|---:|
| 4 | 1.2613 | 1.2233 |
| 8 | 1.2371 | 1.2298 |
| 16 | 1.2548 | 1.3119 |

## Conclusion

**A2 is closed, on its own operator and its own preregistered gates.** The
family's S1 admission conditions were: beat the scalar level control, be no
worse than Tanimoto at k≥2, and weaken materially under wrong-protein and
permuted-label controls.

It fails all four. The performance failures are resolved — it loses to a
parameter-free Tanimoto kernel at **every** k with intervals excluding zero.
The two control gates fail by **inertness**: neither the correct protein nor
correct support labels produce any resolvable change in its output. It is not
that corrupting the input helps; it is that the operator does not measurably
notice.

This is a stronger and more legitimate closure than the v2 probe's. The v2
result rejected a different operator and generalised; this one runs the exact
mechanism the plan specifies, gives it a learned coordinate system, a moment
over support residuals, k-dependent shrinkage, and a rank chosen on held-out
components, and finds that it converges to a worse level shift.

### What is *not* concluded

* Not that no protein-conditioned few-shot mechanism can work. This tests one
  operator over one frozen representation.
* Not that the trunk is at fault rather than the objective. The trunk was never
  trained to carry protein-conditioned SAR (DATAFLOW_AUDIT F6/F7). Stage P is
  the test that separates those, and it is unaffected by this result.
* Not a statement about `meta_test`. Its labels were used for no fitting, selection or reported metric; the process-isolation incident remains open.

## Commands

```bash
conda run -n drug python -m tools.research.a2_exact_probe.extract_ligand_features --split meta_train --output tools/research/a2_exact_probe/features
conda run -n drug python -m tools.research.a2_exact_probe.extract_ligand_features --split meta_val   --output tools/research/a2_exact_probe/features
conda run -n drug python -m tools.research.a2_exact_probe.run_probe --features tools/research/a2_exact_probe/features --output tools/research/a2_exact_probe/A2_EXACT_meta_val.json
conda run -n drug python -m pytest tools/research/a2_exact_probe/tests -q
```
