# A2S-TRACE Mechanism Exploration Prompt (Q2)

Date: 2026-08-01
Status: **exploration objective + preregistration**. Written after Q1 returned and before any
adaptation parameter was fitted.
Supersedes for the mechanism phase: §6 Q2 of `A2S_TRACE_EXPLORATION_PROMPT_2026-08-01.md`.
Does **not** relax any admissibility condition in `A2S_EXPLORATION_PROMPT.md` or in the TRACE prompt.

Companions:
`reports/active/A2S_TRACE_MECHANISM_ANALYSIS_2026-08-01.md` (history/data analysis and mechanism
derivation), `reports/active/a2s_trace_q1_stratum_2026-08-01.json` (the Q1 measurement),
`research/a2s_trace_stratum.py` (Q1 code), `research/a2s_trace.py` (Q2 code).

---

## 0. Why this prompt exists

The programme's binding blocker was **C9**: two credible measurements disagreed about whether
correctly assigned support labels carry any transferable ranking information. The balanced ChEMBL v2
gate said no; a BindingDB branch said yes and large. The TRACE prompt made the consequence binding:
*no mechanism may be trained before the admissible stratum is measured.*

Q1 has now been measured on one corpus, one frozen base, one analytic estimator and one bootstrap
(`research/a2s_trace_stratum.py`, 12,246 probe episodes, 110 homology components).

**FACT — the Q1 answer.** Support information is not a property of the dataset. It is a property of
the **support→query chemical relation**. Under every one of the three declared support policies the
fixed-Tanimoto-KRR ranking gain over the identical frozen base is:

| nearest support Tanimoto of the query | k=3 ΔCI LCB | k=5 ΔCI LCB | verdict |
|---|---:|---:|---|
| `< 0.20` | −0.003 | −0.007 | null in all policies |
| `0.20 – 0.35` | −0.003 | −0.006 | null in all policies |
| `0.35 – 0.55` | −0.003 | +0.006 / +0.012 | marginal, policy-dependent |
| `≥ 0.55` | **+0.023 / +0.036** | **+0.031 / +0.048** | admitted in both non-degenerate policies |

(LCB = paired component bootstrap 95 % lower bound, 2,000 draws; the two numbers are
`scaffold_disjoint` / `random_within_target`.)

The old v2 policy (`provenance_disjoint`) draws queries with mean nearest support Tanimoto 0.19–0.30
and therefore lives almost entirely in the two null bins. **Its null was correct and is now
explained, not overturned.** The BindingDB positive was correct and is now localised, not imported.

**Binding consequence for this prompt.** The mechanism phase runs inside the **admitted stratum**
and is judged against the **fixed analytic smoother that already exploits it**, not against the
frozen base. Beating the frozen base in the admitted stratum is no longer an achievement; Tanimoto
KRR does it with zero learned parameters.

---

## 1. What is fixed

1. **Task.** Abundant-to-scarce DTA transfer. Meta-train on abundant source targets. At meta-test an
   unseen target supplies `k ∈ {1,3,5}` support affinities under a declared policy; the model ranks
   that target's query compounds.
2. **Paradigm.** Meta-learning. The adaptation rule is amortised across source episodes. No
   recipient-specific SGD and no recipient-specific refit of a free parameter.
3. **Primary endpoint.** Target-macro **CI** and **NDCG@10**, restricted to the admitted stratum and
   also reported over the full query set. RMSE/MAE are secondary and may never be substituted.
4. **Label firewall.** `locked` source role and the A2S recipient roster stay sealed. Only `fit` and
   `probe` roles are opened. `fit` trains; `probe` measures; nothing else exists.
5. **Statistical unit.** Protein-homology component. Aggregation is fixed:
   episode draws → seed/target mean → component mean → paired component bootstrap (≥2,000 draws).
6. **Support policy is part of the estimand.** Every number names its policy and stratum.
7. **The Q1 substrate.** ChEMBL-37 dualcold pKi TRAIN, v2 balanced lock
   (`a2s_source_information_gate_lock_v2_2026-08-01.json`, 222/110/107 fit/probe/locked components,
   zero target/homology/document/assay overlap across roles), frozen component-cross-fitted ridge
   base. The base is **not** retrained for the mechanism phase; it is the same frozen object Q1 used.

## 2. What is open

The adapter's functional form, its features, its parameterisation, its objective, and the episode
construction inside the sealed contract.

---

## 3. Objective

> Produce **one mechanism with one load-bearing claim** that, inside the Q1-admitted stratum, beats
> **fixed Tanimoto KRR at equal support information** on target-macro CI and NDCG@10, with a paired
> component 95 % lower bound above the preregistered MDE, and that fails cleanly when the claimed
> innovation is ablated.

Admissibility — all seven must hold, as in the TRACE prompt:

1. **Learned** across source episodes.
2. **Identifiable at the deployment budget.** State the number of target-specific quantities
   estimated at meta-test. Zero is the preferred answer; then the burden moves to showing the
   amortised object transfers to unseen components.
3. **Query-dependent.** Not an episode constant. Q1 measured the episode-constant channel at
   *exactly* 0.0000 CI in every policy × k × stratum cell — this is now a verified structural fact,
   not an assumption.
4. **Structurally abstaining.** `r_S ≡ 0 ⇒ Δ ≡ 0` to floating-point exactness, provable from the
   functional form.
5. **Bounded** by an observed quantity — no unbounded learned scale on a residual aggregate (C2).
6. **Nested-falsifiable.** Removing the claimed innovation recovers a *named* baseline exactly.
7. **Not shortcut-driven**, with structural controls preferred over empirical ones.

---

## 4. The question this phase answers

> Q1 showed that a *fixed isotropic chemical similarity* already transports residuals usefully when
> the query is close to the support. Is the **reliability of that transport** itself a learnable,
> protein-conditioned function of the support→query relation — i.e. does the model know *which* close
> pairs to trust and which to abstain on — and does that knowledge transfer to unseen homology
> components?

This is the medicinal-chemistry content of the activity-cliff phenomenon stated as a meta-learning
problem. It has zero target-specific parameters, so C3 does not bind at all.

---

## 5. Preregistered decision rules

**MDE and thresholds.** Component-level SD of the paired ranking difference in the admitted stratum
is measured from the Q1 records. With 66–84 admitted components, MDE80 ≈ 0.005 CI. The mechanism
phase preregisters:

| Gate | Quantity | Threshold |
|---|---|---|
| **M1 (headline)** | `ΔCI = TRACE − fixed Tanimoto KRR`, admitted stratum, k∈{3,5}, paired component LCB | `> +0.005` |
| **M1b** | same for NDCG@10 | `> +0.005` |
| **M2 (nesting)** | TRACE restricted to the log-Tanimoto scorer reproduces the NW smoother | exact to 1e−5 |
| **M3 (abstention)** | residual-null arm | bitwise `Δ ≡ 0` |
| **M4 (assignment)** | `TRACE(correct) − TRACE(deranged)` CI, k≥3, paired LCB | `> 0` |
| **M5 (magnitude)** | `TRACE(correct) − TRACE(norm-matched wrong-target)` CI, paired LCB | `> 0` |
| **M6 (no harm)** | TRACE − base RMSE in the null strata (`t < 0.35`) | non-inferior within 0.02 |
| **M7 (protein claim)** | protein-shuffle and protein-zero decrement | reported; if ≈0 the method is renamed ligand-only support-conditioned adaptation and **not** called protein-conditioned DTA |

**Stop rules.**

- If M1 fails, TRACE is not a mechanism. Report the null with its measured upper bound and the
  required component count. Do not run a second architecture on the same probe role.
- If M3 or M2 fails, the harness is wrong and nothing downstream counts.
- If M4 fails while M1 passes, the gain is a label-free reranker and must be reported as such
  (a ligand-similarity prior, not support adaptation).
- `probe` is a development role. A passing result authorises **freezing** the protocol and opening
  the `locked` role once. It does not authorise opening the recipient roster.

**Registered predictions (stated before the run).**

| # | Prediction | Falsifies |
|---|---|---|
| P1 | TRACE ≥ KRR in the admitted stratum at k=5, LCB > 0.005 | the whole idea if it fails |
| P2 | TRACE − KRR is larger at k=5 than at k=1 | the "learned reliability needs contrasts" story |
| P3 | The null-slot ablation loses more in the *null* strata than in the admitted stratum | the abstention claim |
| P4 | Derangement destroys the entire TRACE−base gain, exactly as it destroys the KRR−base gain | the assignment claim |
| P5 | The level channel remains exactly rank-null | the harness |
| P6 | Protein-zero costs less than the null-slot ablation | nothing — it calibrates the honesty of the protein claim |

**Prohibited claims** carry over verbatim from `A2S_TRACE_EXPLORATION_PROMPT_2026-08-01.md` §3, plus:

- Do not call a stratum-restricted gain a global DTA result.
- Do not report the admitted-stratum number without the null-stratum number beside it.
- Do not describe Q1's `t ≥ 0.55` bin as "scaffold-cold"; it is explicitly local SAR.

---

## 6. Required output

Mechanism statement; mathematical definition with the meta-test estimand and its dimension;
identifiability argument at k = 1/3/5; the exact increment over the closest prior art with links;
source meta-training procedure with strict OOF construction; one-pass meta-test procedure; ablation
ladder with named fallbacks; structural vs empirical control inventory; registered predictions; the
decisive falsification experiment with MDE and stop rule; and the maximum scientific risk including
the shape of the honest null.

Label every substantive statement **FACT**, **INFERENCE** or **HYPOTHESIS**. Do not force a winner.
A positively-controlled null with a measured upper bound is an acceptable deliverable.
