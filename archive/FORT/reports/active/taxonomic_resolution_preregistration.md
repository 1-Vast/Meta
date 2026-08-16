# Taxonomic-Resolution-Calibrated Ranking (TR) — preregistration

Date: 2026-07-26. Route code **TR**. Authorized by the estimand-reframe decision taken after the
program-wide re-audit returned verdict ② `SIGNAL_PRESENT_EVIDENCE_INSUFFICIENT`. TR is **not** a
revival of BM0/BM1-RR/CFRI/HQ-GBMA: those required target-specific ligand reordering to be
identifiable from protein features or same-target support, which four independent mechanisms have
falsified (covariance orientation, diagonal precision, unconstrained head `I0`, Grassmann subspace —
HQGBMA Stage D `true−global −0.232 [LCB −0.083]`). TR moves the identifiability requirement to the
level the record shows is real — the coarse **KLIFS group** — and treats the target-specific level as
a calibrated, evidence-gated deviation that abstains to the group prior by default.

## 1. Load-bearing admissibility assumption (user-vetoable)

The dual-cold firewall, as **implemented** in the KirHub strict A1 run and the registries, holds out
for each query target: accession, full-sequence homology cluster (4-mer containment union-find,
threshold 0.40), KLIFS **family**, Bemis–Murcko scaffold and Morgan Tanimoto ≥ 0.50 (ligand axis),
and ChEMBL document. It does **not** hold out the KLIFS **group** (TK, CMGC, AGC, CAMK, CK1, STE, TKL,
Other, Atypical). Evidence: A1 computed a `group_centroid` arm for held-out targets
(`0.0830 [0.0628,0.1028]`, n=308 components), which is only possible because same-group training
targets survive the family/homology holdout.

TR therefore treats KLIFS group as a **legal test-time input**, because: (i) the group label is
derivable from the query target's sequence alone, with no affinity supervision and no knowledge of
which specific targets are in training; (ii) group is strictly coarser than the held-out family;
(iii) the README requires gains from homologous proteins to be *distinguished* from target-specific
transfer (honestly attributed), not forbidden. TR's setting is precisely
**family/homology/accession-cold, group-warm** dual-cold — a realistic, weaker-but-honest transfer
regime: a genuinely novel kinase still belongs to a known group.

This assumption is load-bearing. If the intended contract forbids using group, the **G0 group-cold
falsification control** (Section 6) still returns an honest answer and the route collapses to the
existing ② negative. Flagged for veto before TR-0 runs.

## 2. Scientific hypothesis

**H_TR:** The coarse KLIFS-group taxonomy carries a transferable within-target ligand-ranking
direction (which ligands a group prefers) that (a) beats the ligand-only baseline B0 on
family/homology-cold targets and (b) is specific to the correct group; whereas finer-than-group
(target-specific) ranking resolution is either not identifiable or requires same-target support that
current open panels do not chemically cover (median max query→support Tanimoto 0.185 at k=4). A
predictor that makes the group level its load-bearing, always-available prior and gates the
target-specific deviation by a coherent evidence ratio will (i) beat B0 at all k including k=0, and
(ii) match — but not exceed — the group ceiling unless the target level is genuinely identifiable,
with the resolution decision itself calibrated.

## 3. Estimand

For a query target `t` (accession/homology/family/scaffold-cold; group `G(t)` known from sequence)
and query ligands `d`, predict a calibrated within-target ranking of `{d}`:

```text
y_hat(t,d) = b(d) + f_{G(t)}(d) + omega_t * delta_t(d)
```

* `b(d)`: frozen B0 ligand-only reference (unchanged; Morgan-1024 ++ 10 physchem).
* `f_{G(t)}(d)`: group-level residual ranking direction — an empirical-Bayes partially-pooled
  estimate over **train** targets in group `G(t)`, shrunk toward the global residual direction.
  Available at `k=0`.
* `delta_t(d)`: target-specific ranking deviation inferred from the `k` support residuals
  `r_i = y_i − b(d_i) − f_{G(t)}(d_i)`.
* `omega_t ∈ [0,1]`: evidence-derived resolution weight (the innovation). `omega_t = 0` at `k=0`
  (exact fallback to the group-resolved ranker). Set by a coherent marginal-likelihood ratio between
  the group-only and group+deviation residual models; **not** a free neural gate.

Primary ranking metric: component-macro Spearman; statistical unit = homology component; grouped
bootstrap. Calibration/selective metrics: Section 6 (frozen before running). At `k=0` the predictor is
a pure group-resolved zero-shot ranker — the new zero-shot object whose floor must beat B0.

## 4. Innovation module (exactly one; on the main path; affects inference)

**Calibrated Resolution-Abstention Operator (CRA-Ω).** Given support residuals after group
correction, CRA-Ω computes `omega_t` as a monotone map of
`log[ p(support | group+deviation) / p(support | group-only) ]` under a fixed model prior, with
train-fitted calibration hyperparameters (evidence temperature, model prior) selected only on
training components. It affects the learnable inference process (posterior resolution weight) and
lies on the main prediction path. The group empirical-Bayes prior (partial pooling / James–Stein
shrinkage) and B0 are mature, non-innovative components.

**Difference from BM0's Bayes-factor gate (closed):** BM0 mixed B0 with a target-adapted function,
so full abstention gave zero gain and the route bet entirely on target-specific identifiability.
CRA-Ω mixes the **group prior** with the target deviation, so full abstention (`omega=0`) still beats
B0 via group transfer; the bet is only on whether the target level *adds to* the group level, and the
abstention default is a positive, not a null.

## 5. Why identifiable / difference from failed routes

TR does not require "protein feature → target-specific reordering" (falsified 4×). It requires only
"KLIFS group → group-level reordering," already supported by the frozen A1 arms
(`group_centroid 0.083` vs `ligand 0.043`, non-overlapping CIs, 308 components). The group level is
identifiable because many training targets share each group, giving a well-estimated pooled
direction; the family/homology holdout guarantees the specific query target and its close homologs
are unseen, so the transferred signal is coarse-taxonomy transfer, not memorization. Literature
grounding: hierarchical/partial-pooling shrinkage (Efron & Morris, James–Stein), empirical-Bayes
group means; kinase group taxonomy (Manning et al. 2002, KLIFS Kanev et al. 2021); selective
prediction / abstention (Geifman & El-Yaniv 2017).

## 6. Staged gates (frozen before any execution)

Cheapest-first. **TR-0** is analysis + two cheap controls on the frozen A1 predictions (no model
training, no new selection, no confirmation/sealed access). **TR-1** is the model, contingent on TR-0.

### TR-0 — premise, specificity, ceiling (analysis gate)

Re-analysis of the frozen KirHub strict A1 evaluation (immutable JSON
`90a72159c7ac1626df3f437fe63b6ec94d8163c549dc64b336715b5bf054d4d3`, components registry
`0122780145a5504a61d78a709829d65d093f5636ccb4bca117d5e39cc2c4a7d4`) plus two preregistered new
controls (`random_group`, `shuffled_group`), all on the same 308 homology components, grouped
bootstrap, one seed (1729), no tuning. MDE80 = 0.016 (SD 0.10); `max(0.03, MDE80) = 0.03`.

* **G1 (premise):** paired `group_EB − B0` ≥ 0.03 with grouped LCB95 > 0.
* **G2 (specificity):** paired `group_EB − random_group` and `group_EB − shuffled_group` both have
  grouped LCB95 > 0 (the correct group's direction transfers, not any group-sized effect).
* **G3 (ceiling, non-gating diagnostic):** paired `true_protein − group_EB` and `support_k4 −
  group_EB`, grouped CI. If LCB95 ≤ 0 (expected), the fine level is not identifiable and TR-1's
  contribution is defined as calibration, not raw ranking gain.
* **G0 (group-cold falsification control):** recompute `group_EB` with the query target's GROUP
  additionally held out (pooling over other groups only); the gain over B0 must collapse, proving the
  signal is group-level and protecting the Section 1 admissibility assumption.

A **G1 or G2 failure stops TR** (the premise is false). A G3 pass would be a genuine surprise and
would authorize a separate target-level ablation.

### TR-1 — HTS-CRA calibrated model (contingent on G1 & G2 pass)

Model selection on **train components only** (nested train-component CV for shrinkage strength,
evidence temperature, model prior). One single-seed evaluation on the development substrate with
`k ∈ {0,4,8,16}`.

* **G4 (calibration & selective value — primary for TR-1):** (i) at `k=0` beat B0 with grouped
  LCB95 > 0 (floor = group prior); (ii) never fall below `group_EB` ranking (no harm from the
  abstention machinery); (iii) calibrated resolution decisions — the selective-risk curve (ranking
  quality vs abstention rate) dominates a random-abstention baseline, and `omega_t` is monotone in
  true target-level improvement on train; (iv) exact `k=0` fallback, support-permutation invariance,
  finite positive variance.
* **G5 (mechanism specificity):** the `omega_t` machinery must be damaged by support-label
  permutation and wrong-target support (it should abstain, not adapt, under bad support).

TR-1 passing authorizes only a review for powered independent confirmation, which the program has
established does not currently exist (`NO_OPEN_POWERED_INDEPENDENT_PANEL`). TR-1 is therefore expected
to terminate as a train/dev-validated deliverable, not a confirmed one. No sealed/confirmation
access, no multi-seed, no long training without a separate review.

## 7. Expected source of improvement & honest expected outcome

Expected gain source: coarse-taxonomy (group) ranking transfer, real and ≈2× B0 in the frozen
record. **Honest expectation:** G1 passes, G2 passes, G3 fails (target ≯ group), so TR delivers a
calibrated group-resolved ranker that beats B0 and correctly abstains from unidentifiable
target-level refinement. This moves the program's stated reproducible signal from "B0 only" to "a
validated, calibrated group-resolved ranker over B0, with a sharp identifiability ceiling at the
group level" — meaningful as a better-identified estimand and a narrowed solution space, **not** as
target-specific transfer. TR will not be described as target-specific unless G3 passes.

## 8. Risks and controls

* **Homology/taxonomy shortcut:** mitigated by honest attribution (group-level, not target-specific),
  G2 specificity controls, and the G0 group-cold falsification.
* **Dev exploitation:** TR-0 is one-shot re-analysis of frozen predictions with frozen gates; TR-1
  selects only on train components; no iterative dev tuning.
* **Within-source substrate:** KirHub cannot isolate assay/document, so TR-0/TR-1 on KirHub are
  within-source mechanism/deliverable probes, never cross-assay confirmation. A non-gating
  replication on the frozen Metz/Reinecke group arms may be added if already computed.
* **Pseudoreplication:** statistical unit = homology component; grouped bootstrap; support rows never
  treated as independent pairs.
* **Triviality:** the group prior alone is a mature estimator; the claimed contribution is the CRA-Ω
  calibration/abstention, evaluated by the calibration/selective metrics preregistered above.

## 9. Budget & stop-for-review

TR-0: analysis only, ready to run on approval. TR-1: single seed, train-only selection, contingent on
TR-0. **Stop for review** before TR-1 training, before any multi-seed, and before any
confirmation/sealed access. `sealed_test_consumed=false`; `confirmation_labels_read=true`
(pre-existing; TR reads no confirmation labels).
