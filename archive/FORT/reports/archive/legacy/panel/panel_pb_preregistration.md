# Gate PB preregistration - zero-shot dual-cold interaction on the dense panel

Registered 2026-07-25 after `PANEL_GATE_PA_PASS` and after the panel power audit froze the
threshold, and before any target-conditioned arm was scored.

## Position in the gate order

task.md fixes the order: power audit, then zero-shot Gate Z, then few-shot. PA answered a question
that precedes both -- whether the interaction exists and whether the frozen protein representation
aligns with it -- and passed with adaptive `p=0.000488` against the unchanged `0.01` threshold, a
projected residual SD of `0.6629` pK against the unchanged `0.5` minimum, and top-1% axis energy of
`0.084` (target) and `0.065` (ligand) against the unchanged `0.50` maximum. PB is the predictive
counterpart: can a model *use* that alignment to reorder ligands for a target it has never seen?

## Frozen threshold

`research/panel_power.py` fitted the identical ligand-only base at four seeds under the identical
leave-component-out protocol before any interaction arm existed. Dual-cold development target-macro
Spearman was `0.2786`, `0.3022`, `0.2938`, `0.2913` over 101 components; the per-component retraining
noise SD is `0.0948` and the paired MDE80 is `0.0181`. The Gate PB threshold is therefore
`max(0.03, 0.0181) = 0.03` -- the task's nominal effect size, resolvable on this substrate for the
first time in the program (PLINDER `0.0954`, sparse ChEMBL `0.0586`). Report:
`reports/active/panel_power.json`.

## Protocol

Five deterministic folds over the 101 sequence-homology components, fold map fixed by
`research/panel_power.py` and reused unchanged. For each fold, every arm is fitted on the TRAIN
cells (anchor ligands) of the other folds' components and scored on the DEVELOPMENT cells
(scaffold- and connectivity-disjoint query ligands, max anchor Tanimoto `0.9091`) of this fold's
components. Every arm sees identical rows in identical folds, so all comparisons are paired over
components. Residual arms consume out-of-fold base predictions produced by an inner 5-fold
cross-fit over the outer fold's training components; the provenance audit must report zero
violations. Panel confirmation cells are never read.

Budgets are matched across arms and fixed here: base `4,000` steps, inner cross-fit `5 x 2,000`
steps, every interaction arm `4,000` steps, at most `256` rows per training episode.

## Arms

`B0` ligand-only; `T0` target-only; `A0` additive `b(d)+a(t)` (cannot reorder within a target);
`I0` jointly trained `b+g` without cross-fitting; `R0` plain cross-fitted residual regression;
`CFRI` cross-fitted residual with rank-reversal, target-centering and base-orthogonality;
and the controls `CFRI-Tshuffle` (protein derangement inside the panel target set),
`CFRI-Lshuffle` (ligand derangement in the interaction branch only), `CFRI-Tpool` (amino-acid
composition and length instead of ESM-2) and `CFRI-Trandom` (random frozen protein embeddings).

## Criteria

| id | criterion | threshold |
|---|---|---|
| PB1 | `CFRI - B0` paired component macro Spearman | `>= 0.03` (frozen power audit) |
| PB2 | grouped component bootstrap of that difference | `LCB95 > 0` |
| PB3 | fraction of components with a positive difference | bootstrap `LCB95 > 0.5` |
| PB4 | RMSE safety versus B0 | `<= 1.02 x B0` |
| PB5 | damaged by both derangements | `CFRI - Tshuffle` and `CFRI - Lshuffle` both `LCB95 > 0` |
| PB6 | real protein beats random protein embeddings | `CFRI - Trandom` `LCB95 > 0` |

PB1, PB2, PB4 and PB5 are the unchanged Z1, Z2, Z4 and Z5 criteria. PB3 replaces Z3
("positive in every component"), which cannot be met by any real effect on 101 components and is
therefore not a test; the sign-majority form has no tuned constant. PB6 is new and strictly harder:
it is exactly the criterion `BM1_RR_FAIL_STOP` failed (`+0.0035` against random protein), so the
route is not allowed to survive on a repeat of that failure. Z6 (stratum improvement) is dropped
because the panel has no date axis and its similarity stratum is already implied by the
scaffold-disjoint query design.

A single failed criterion stops the zero-shot panel route and is recorded in `history.md` with its
cause. A pass authorizes the separately preregistered few-shot Bayesian gate PC on this substrate
and nothing else: no multi-seed run, no Hierarchical MoT, no long training, no confirmation access
and no threshold change.

---

## Amendment, 2026-07-25, written after the Gate PB result and before Gate PC was run

Gate PB failed: `PANEL_GATE_PB_FAIL_STOP` (PB1 `+0.0180 < 0.03`, PB2 `LCB95 -0.0145`, PB3 sign LCB
`0.4158`, PB5 derangement damage not significant; PB4 and PB6 passed). Under the clause above, that
result withdraws the authorization this document was going to grant to Gate PC.

Gate PC is nevertheless run once, as a registered amendment rather than as a rescue, on these
grounds:

1. task.md authorizes the Bayesian meta-adapter independently of any zero-shot interaction result.
   Its BM0 contract fixes `mu0(t,d) = b(d)` precisely so that few-shot adaptation does not inherit
   the failed zero-shot interaction. PC is therefore a different mechanism, not a retry of PB.
2. No PC design choice can have been influenced by PB. `reports/active/panel_pc_preregistration.md`,
   `research/panel_episodes.py`, `research/panel_gate_pc.py` and the k=4 power audit
   (`reports/active/panel_power_k4.json`, frozen threshold `0.0367`) were all written and, for the
   power audit, executed before the Gate PB result existed. Nothing about PC's arms, criteria,
   episodes or threshold is tunable after the fact.
3. PC's own criteria are strictly harder than BM0's, because PC6 requires beating both protein
   shuffle and random protein -- the criterion `BM1_RR_FAIL_STOP` failed.

This amendment does not reopen PB, does not authorize a second PB run, and does not change any PB or
PC threshold. If PC fails, both failures stand and the panel route stops for review.

---

## Governance note, 2026-07-25 (audit item 7)

The 2026-07-25 audit fixes the evidentiary status of both runs on this substrate. Gate PB remains a
preregistered failure and is unchanged. Gate PC was executed beyond PB's original sequential stop
under the amendment above and is **exploratory**: it may generate hypotheses, but it cannot serve as
confirmatory evidence and cannot authorise any later gate. Neither result is reinterpreted as
positive, and no threshold in either document is altered.
