# I-2 — coupling identifiability audit: where the recoverable structure actually is

Date: 2026-08-09. Label-side only. No training, no GPU, no affinity read.
Held-out A: 2,409 complexes with positives, 196 protein closure components;
2,314 entered the coupling test (95 had an active submatrix below 3×3).

## 0. Result in one table

Macro-AP over components, complete residue × heavy-atom matrix:

| Object | macro-AP |
|---|---:|
| ligand-only (measured model) | 0.00450 |
| **Oracle-A** — true *atom* marginal, oracle-exposed | **0.00850** |
| **B4** (measured model, non-PLM residues) | **0.02295** |
| **Oracle-R** — true *residue* marginal, oracle-exposed | **0.21633** |
| **Oracle-RA** — additive projection of both marginals | **0.39075** |
| exact pair (sanity upper bound) | 1.00000 |

## 1. The ANOVA machinery is correct, and it is not double centering

The registered coupling object is

```text
G_ij = mu + alpha_i + beta_j + C_ij ,   C = G − Proj_W(additive space)
```

`Proj_W` is solved as a genuine weighted least-squares ANOVA on the actual mask
by weighted alternating least squares with an identified re-centring step.

Self-test, both directions:

* complete mask + uniform weights → the solver reproduces classical double
  centering to **2.78e-16**, so it is correct;
* non-uniform weights → the solver **differs** from double centering by up to
  **0.423**, so it is genuinely the weighted projection and not double centering
  wearing a different name.

Naive double centering is used only as this oracle in the case where the two
provably coincide, and never on real data.

## 2. Residue localisation dominates atom propensity by ~25×

Oracle-R = **0.216** against Oracle-A = **0.0085**. Knowing which *residues*
contact anything is worth about twenty-five times more average precision than
knowing which *ligand atoms* contact anything.

This is a fresh measurement on our own corpus, split, evaluator and code. It is
**not** a reproduction of the unverified consolidated-report figures
(0.25262 / 0.02494), which remain external claims; the qualitative direction
agrees, the numbers are ours and differ.

## 3. The measured model is nowhere near the marginal ceiling

B4 reaches **0.0229** against an oracle-marginal ceiling of **0.3907**. It
recovers roughly **6 %** of what exact marginals would deliver. Headroom to
Oracle-RA is **+0.368 AP**; headroom to Oracle-R alone is **+0.193 AP**.

The bottleneck is therefore **recovering the residue marginal** — predicting
*which residues bind at all* — not resolving *which residue pairs with which
atom*.

## 4. Coupling beyond the marginals is weak in the labels

Statistic: leading-singular-value share of the marginal-orthogonal residual `C`
on the active submatrix. Null: degree-preserving bipartite rewiring by
checkerboard swaps, 20 independent rewirings per complex, 30 × (positives) swap
attempts, holding every `d_i` and every `e_j` exactly fixed.

| Quantity | Value |
|---|---:|
| true statistic, mean | 0.6124 |
| degree-preserved null, mean | 0.5921 |
| difference | **+0.0203** |
| median z against each complex's own null | **+0.41** |
| complexes above their own null | **63.1 %** (chance = 50 %) |

There is a real but **small** excess. The typical complex's residual structure
sits within half a standard deviation of a degree-matched random matrix.

**Consequence for the registered architecture.** Making the coupling term
marginal-orthogonal and having it *replace* the free pair term is still the right
construction — it stops the pair term from silently absorbing marginal effects,
which is an identifiability improvement regardless of effect size. But this
measurement predicts, in advance and falsifiably, that **the coupling head will
not be the source of a large AP gain on these labels**, because the labels carry
little coupling beyond their margins. Anyone reading a large gain from a coupling
head on this corpus should first suspect marginal leakage into the coupling term.

## 5. Failure localization — now much sharper

| Candidate | Status | Evidence |
|---|---|---|
| **data / label insufficiency for marginals** | **EXCLUDED** | the labels carry Oracle-RA = 0.391; there is a great deal to find |
| **data / label insufficiency for coupling** | **PARTIALLY CONFIRMED** | median z = 0.41 against a degree-preserving null; coupling beyond margins is weak |
| optimization failure | EXCLUDED | trainability control recovers a known function at 0.759 macro-AP |
| closure / inference units | EXCLUDED | 196 components, largest 15.2 %, ligand-disjoint from train |
| atom correspondence | EXCLUDED | I-1: 14,585 admitted, 4 enumerated and quarantined |
| **biological representation failure** | **LEADING EXPLANATION** | B4 recovers 6 % of the marginal ceiling; the gap is 0.368 AP and lies in the residue representation |
| section non-identifiability | NOT REACHED | no adaptation has been attempted |

The previous state record listed representation failure as one of several
competing explanations. This audit removes "there is nothing to find" from that
list: there is a large, oracle-verified marginal signal that the current residue
features do not recover.

## 6. What this authorises, and what it does not

It **supports** the registered B5 design precisely as written: change **only**
the residue features, leave the ligand branch, head, rank, sampler, budget and
evaluation mask fixed. The measurement says the residue side is where the
recoverable structure is, so the one permitted change is the one aimed at it.

It **does not** authorise attention, a geometry branch, a larger PLM than the
registered 650M, an affinity head, or any additional branch. It does not
authorise treating rewiring as a training negative — rewiring is an evaluation
control here and was used only as one.

It makes no affinity, ranking, transfer, few-shot or `z`-admission claim.

## 7. Remaining execution blockers

Discharged by this work: `ATOM_CORRESPONDENCE_NOT_FULLY_VERIFIED` (I-1).

Still open: tie-aware AP and per-pair prediction materialisation; negative and
control manifest completion; publication/time closure; ESM2-650M weight
acquisition. B5 remains operationally blocked.
