# Literature-Informed Solution Menu

Status: **PROPOSAL. Not a registered stage, not an artifact, not hashed.**
Companion to `IDENTIFICATION_ROADMAP_AND_Z_ADMISSION.md` (2026-08-08).
External claims are sourced; internal claims cite the ledger. Nothing executed.

---

## Part A — The measurement problem has already been solved in the literature

### A.1 Interaction Concordance Index (Turku ML, Oct 2025)

`arXiv:2510.14419` introduces the **IC-index**, and it is a direct match to
MetaSieve's situation — arguably the most useful single paper for this project.

**It formalizes the exact 2x2 decomposition derived in Part I of the roadmap.**
For any two drugs and two targets it decomposes the quadruple of outputs into
grand mean, drug main, target main and interaction:

```text
y_DxT = 1/4 (y - y' - y* + y'*)
```

which is the double difference `DD`, up to the factor `1/4`.

**Proposition 1** states that interaction concordance "is invariant to additions
of constant, drug symmetric, target symmetric and additively separable
functions." IC-index therefore measures `delta` **and nothing else** — it is the
interaction-side mirror of what E-AFF-R0 discovered empirically about
within-task concordance. The paper's Table 2 gives the full invariance matrix
across accuracy, global concordance, drugwise concordance, targetwise
concordance and IC. R0 rediscovered one row of that table; the paper supplies
all of it, with proofs.

**Proposition 2 is the result that explains the entire MetaSieve failure chain.**
It partitions off-training-set data four ways — IDIT / IDOT / ODIT / ODOT (in- or
off-training-set drugs x targets) — and proves that a learning algorithm's
**drug and target permutation equivariance implies it is uniformly badly
aligned** with a specific set of these problems, meaning expected utility is
*exactly* `0.5`. The proof is three lines: equivariance makes the predictor
distribution symmetric under swapping two off-training-set targets, so the
expected sign of the predicted interaction is exactly balanced.

The paper's stated remedy is precisely MetaSieve's programme: "In practical
applications, this equivariance is remedied via incorporation of appropriate
side information on drugs and targets."

**What this buys, concretely.**

1. A published, peer-reviewable **readout for Claim B**, with proven invariance,
   an `O(min(n_D^2 n_T log n_T, n_D n_T^2 log n_D))` algorithm, and reference code
   at `github.com/TurkuML/IC-index-experiments`.
2. A **theoretical necessity argument for side information** under unseen
   targets — which is the justification P1B has been missing.
3. A **standard vocabulary** (IDIT/IDOT/ODIT/ODOT) that maps onto the project's
   source/metaval/recipient splits and would make the work legible to reviewers.
4. Confirmation that the roadmap's `delta` reframing is not idiosyncratic —
   an independent group formalized the same object in the same year.

### A.2 The same group supplies the estimation layer

The IC-index authors are the Turku group behind the standard reference on
realistic DTI evaluation (Pahikkala et al., *Briefings in Bioinformatics* 16(2)
325–337, which introduced the four-setting evaluation paradigm) and the
pairwise-kernel toolchain: the generalized vec trick (`arXiv:2009.01054`) and,
as of June 2026, **SPaiK** (`arXiv:2606.16979`), a stochastic generalized vec
trick with a bundle-method optimizer that trains Kronecker-kernel pairwise
models on mini-batches of pairs. That closes the scalability objection to §B.2
below.

---

## Part B — Solutions, ranked by decisiveness per unit cost

### S1. Estimate `Var(delta)` as a variance component, not by counting rectangles **[R]**

**The single most important idea in this document.**

X0-FEAS proved that *one estimator* of interaction — the disjoint-rectangle
variance ratio — is unattainable, because D1 closure unions every target pair
sharing a document, capping effective units at Ki `36` / Kd `12` against a
requirement of `245`. It did **not** prove that `Var(delta)` is unestimable. It
proved that a rectangle-counting estimator is.

Plant breeding solved this exact problem decades ago. The reaction-norm model of
Jarquín et al. (2014) writes

```text
y = mu + g_i + E_j + (gE)_ij ,      Cov(gE)  =  K_G  (x)  K_E
```

and fits it by REML, where `K_G` is a genomic kernel from markers and `K_E` an
environmental kernel from covariates. Because the interaction covariance is
structured by **kernels built from side information rather than from identity**,
the model predicts into *untested environments* — the exact analogue of unseen
targets. Translated:

```text
y[s,t,l] = mu[s] + alpha[t] + beta[l] + delta[t,l],
Cov(delta) = K_ligand (x) K_protein
```

with `K_protein` from ESM plus P1B geometry and `K_ligand` from the frozen GINE.
This is mathematically the Kronecker pairwise kernel — Kron-RLS — so the Turku
toolchain fits it directly, and SPaiK makes it tractable on `152,737` rows.

**Why this changes the power calculation.** A REML variance component uses the
whole matrix at once and borrows strength through the kernel structure. It is not
restricted to cell-disjoint rectangles, so the `36`/`12` ceiling does not apply
in the same form. Information about `Var(delta)` comes from every observed cell,
not only from the sparse set of clean 2x2 sub-designs.

**Caveats to preregister.** (i) The likelihood-ratio test for a variance
component at the boundary has a non-standard null — a 50:50 mixture of `chi^2_0`
and `chi^2_1` — and clustering must enter the null, ideally by a
document/homology-cluster parametric bootstrap rather than an asymptotic
approximation. (ii) `K_protein` built from ESM can smuggle target main effects
into `delta`; the protein-derangement control is what tests this and must be
retained. (iii) `sigma_assay` must come from replicates, not from the residual.

**Verdict: this is the highest-value single change available.** It attacks the
`CROSSED_INTERACTION_EXISTENCE_NOT_YET_TESTED` boundary at the estimator level
rather than requesting more data.

### S2. Get genuinely crossed data without touching the DAVIS seal **[E/R]**

My earlier recommendation — break the DAVIS freeze for the source split — is
**superseded and should be withdrawn.** It was unnecessary.

The kinase field contains several complete, single-platform, crossed profiling
matrices that are *not* the sealed recipient:

| Panel | Scale | Endpoint |
|---|---|---|
| Davis et al. 2011 | 72 inhibitors x 442 kinases = **31,824** Kd | Kd, KinomeScan — **this is the sealed set** |
| Karaman et al. 2008 (*Nat Biotechnol*) | 38 inhibitors x 317 kinases | Kd, KinomeScan |
| Metz et al. 2011 (*Nat Chem Biol*) | 178 inhibitors x ~300 kinases | Ki, functional |
| PKIS / Elkins et al. 2016 (*Nat Biotechnol*) | published kinase inhibitor set | % inhibition |

Karaman, Metz and PKIS are independent of the sealed DAVIS recipient, are
crossed **by construction** rather than by lucky document co-occurrence, and
carry one protocol each — so `mu[s]` is constant and the assay-context
confounding that dominates ChEMBL simply does not arise.

X0-FEAS's structural insight was that crossing and document-disjointness are
produced by opposite kinds of study. The resolution is not to search ChEMBL
harder; it is to use the study design that produces crossing. **Claim B is an
existence question, and one clean crossed panel can answer it.**

Two things I have **not** verified and that must be checked before registering:
the current licence and machine-readable availability of each panel, and whether
their targets can be homology-governed against the sealed DAVIS recipient at the
project's `<40%` threshold. Karaman and Davis are both KinomeScan, so kinase
overlap with the sealed set will be substantial and the exclusion may bite hard.
Metz, being a different assay technology and compound set, is the more likely
candidate.

### S3. Stratify by `kappa` — now with external support **[E]**

The roadmap identified that `run_eaff_l0r.py:219` pools one population band
across all assay strata, leaving `beta_0(z) = b_pop[kappa(z)]` unused. The
literature says this is not a minor inefficiency.

Landrum & Riniker (*JCIM* 64(5) 1560–1567, 2024), "Combining IC50 or Ki Values
from Different Sources Is a Source of Significant Noise", report that for
minimally curated IC50 data **~65% of repeated measurements differ by more than
0.3 log units and 27% by more than one log unit**. Their conclusion is that
pooling across assays is scientifically risky even for the same compound and
target.

L0R pooled across strata and then asked a ligand signal worth `0.034` log units
to be visible against that. The `+0.03421 [-0.03304, 0.10793]` positive-control
failure is the predicted consequence.

### S4. Externally validate the noise ceiling, and use it as a stopping argument **[E]**

Kramer, Kalliokoski, Gedeck & Vulpetti (*J Med Chem* 2012), "The Experimental
Uncertainty of Heterogeneous Public Ki Data", estimate for public ChEMBL Ki a
**mean error of 0.44 pKi units and a standard deviation of 0.54 pKi units**,
after filtering for unit-transcription errors, undifferentiated stereoisomers
and repeated citations of single measurements — which they found accounted for
90% of all pairs.

MetaSieve's independently measured `sigma_assay = 0.47971 [0.47034, 0.48946]`
from `4,261` replicate cells sits squarely inside that published range. This is
worth stating explicitly in any write-up: the project's noise estimate is
externally corroborated, so the ceiling arguments built on it are sound rather
than an artifact of governance choices.

It also supports a *quantitative* stopping argument. Kramer et al.'s framing —
experimental uncertainty "defines a natural upper limit to the predictive
performance possible" — is exactly the roadmap's §2.2 point that a `0.240`
log-unit Gate margin sits at 32% of the total reducible range.

### S5. Decompose the interaction before testing it — AMMI **[R]**

P1R2A already performed the ANOVA half of an AMMI analysis: source candidate
variance is `77.58%` ligand main, `20.76%` protein main, `1.67%` interaction.
That is precisely the additive-main-effects step. What was never done is the
**multiplicative** step.

AMMI (additive main effects and multiplicative interaction) follows the ANOVA
with a PCA of the interaction residual, giving ranked components with
explained-variance shares and per-component significance tests. In published
multi-environment trials the leading interaction PC routinely captures 50–75% of
total interaction sum of squares.

If MetaSieve's `1.67%` interaction concentrates similarly, the first interaction
PC is a far more testable object than the undifferentiated residual: a
concentrated `~1%` signal with a defined direction, rather than `1.67%` spread
over 288 coordinates. The modern replacement — **factor-analytic mixed models,
FA(k), fitted by REML** — handles unbalanced data, which AMMI's complete-grid
requirement does not, and is the natural companion to S1.

### S6. Turn the derangement control into the contribution **[E]**

This reframes F-59 and may be the most publishable thing the project owns.

F-59 evaluated the frozen PSICHIC pair latent (Koh et al., *Nature Machine
Intelligence* 6, 673–687, 2024) on the Core Episode V3 cohort and found:

```text
pair-minus-ligand CI      = +0.05001     (looks like a real gain)
correct-minus-deranged CI = -0.00033     (no protein specificity at all)
```

That is precisely the pattern Proposition 2 of the IC-index paper predicts for
an algorithm evaluated on off-training-set targets. The gain is main effects; the
interaction content is at chance.

Now look at the current state of the art. AdaMBind (*Nature Communications*,
2026, `s41467-026-70554-5`) is a direct peer to MetaSieve — few-shot DTA for
unseen targets, DAVIS/KIBA/BindingDB, support sizes 5 and 40, and a CD-HIT 40%
sequence-identity "novel task split". Its own reported numbers under the hardest
setting are mixed: on the novel split with `k=5`, its CI is `0.6980` on
BindingDB against a best baseline of `0.7240`, and `0.7591` on Davis against
`0.7784` — **worse than baseline on two of three datasets.**

More importantly: **neither AdaMBind nor any of its eight baselines runs a
protein-derangement control.** They report MSE, CI, R², Spearman and Pearson
against baselines, and never test whether the *correct* protein is required. The
IC-index theory says that is exactly the test their reported metrics cannot pass
by construction on off-training-set targets.

**MetaSieve has built the evaluation harness this subfield lacks** — homology
closure, document closure, sealed recipients, one-to-one `<40%` derangement maps
with reuse zero, coupling-null marginals, and preregistration with independent
post-run audit. Applying it to published checkpoints is cheap (F-59 already did
it once), needs no new corpus, and would produce a result of general interest
whether or not MetaSieve's own statistic ever passes.

Suggested framing: *"Derangement-controlled evaluation of drug-target affinity
models: reported gains over ligand-only baselines are main effects."* The IC-index
paper supplies the theory; MetaSieve supplies the governed empirical harness and
92 documented negative controls.

### S7. Adopt IC-index as the Claim B readout **[E]**

Independently of S1, replace or supplement within-task concordance with IC-index
in every future interaction-facing stage, and report results in the
IDIT/IDOT/ODIT/ODOT partition. This is a pure measurement change, is label-free
to implement, has reference code, and removes the R0 blindness problem by
construction rather than by argument.

### S8. Keep the equivalence arm **[R]**

Restating from the roadmap because the literature strengthens it: with an
externally corroborated noise floor (S4) and a properly powered estimator (S1),
an interval that *excludes* the registered effect is a real finding.
`INTERACTION_ABSENT_AT_REGISTERED_SENSITIVITY` should be a permitted terminal
verdict everywhere, or the project can only ever accumulate not-runs.

---

## Part C — Recommended sequence

| # | Action | Freeze | Cost | Decisive? |
|---|---|---|---|---|
| 1 | Implement IC-index; recompute it on **existing** frozen artifacts (P1R2A panel, H0A, H0C, F-59 PSICHIC latents) | [E] | days | reinterprets the whole ledger at near-zero cost |
| 2 | `kappa` repair + Z0.7 dead-import fix (roadmap Part IV) | [E] | hours | unblocks everything |
| 3 | Licence/availability audit of Karaman, Metz, PKIS; homology-governance against the sealed DAVIS set | [E] | days | decides whether S2 is live |
| 4 | Kronecker-kernel REML power study on the governed corpus, label-blind where possible | [R] | ~1 week | decides whether S1 rescues Claim B |
| 5 | Whichever of S1 / S2 survives 3–4 | [R] | panel | **yes — closes Claim B** |
| 6 | Derangement-controlled benchmark write-up (S6) | [E] | weeks | independent of 4–5 succeeding |

Item 1 deserves emphasis: IC-index can be computed on **already-collected,
already-hashed** data. It cannot resurrect a failed Gate, and it must be
registered as a re-analysis with that limitation stated — but it can tell you
whether the interaction signal was present and unmeasured, at the cost of a few
days and no new labels.

---

## Part D — Honest assessment

The literature moves the estimate in two opposite directions, and both should be
recorded.

**More pessimistic.** Kramer et al. and Landrum & Riniker independently confirm
that public affinity data carries `~0.5` log units of irreducible noise, so the
project's ceiling arguments are real and externally grounded. AdaMBind, the
current state of the art, shows CI *regressions* against simple baselines under
cold splits — consistent with the whole subfield operating near a noise floor
with main effects doing most of the work.

**More optimistic.** Three specific things were being done sub-optimally, and
each has a published fix: the readout was blind to interaction (IC-index), the
estimator was inefficient by orders of magnitude (Kronecker-kernel REML), and
the corpus was the wrong shape when correctly shaped ones exist (Karaman, Metz,
PKIS). None of these is a model-capacity change, and none requires weakening a
Gate.

**The reframing that matters most.** Proposition 2 of the IC-index paper says
that for unseen targets, an algorithm without side information gets *exactly*
chance on interaction — not approximately, exactly. Every MetaSieve stage that
found `correct-minus-deranged ≈ 0` was measuring a quantity the theory says is
pinned at chance unless the side information carries genuine partner
specificity. P1B demonstrably does carry it in the geometric sense (AUPRC
`0.43885` vs `0.05149`). The unresolved question is whether that geometric
specificity survives the projection into affinity — and that is now a
well-posed, adequately-tooled question rather than an open-ended search.

---

## Part E — What I did not verify

- Licence, format and current availability of Karaman, Metz and PKIS data.
- Whether those panels survive `<40%` homology governance against the sealed
  DAVIS recipient. Karaman is the same KinomeScan platform as DAVIS, so overlap
  is likely substantial.
- The IC-index paper is an arXiv preprint (Oct 2025); I have not confirmed peer
  review status. Its Propositions 1 and 2 have short, checkable proofs and should
  be verified directly rather than cited on trust.
- Whether SPaiK's public implementation is mature enough for a governed pipeline.
- Exact power of a Kronecker-REML variance-component test on *this* corpus. S1
  argues the ceiling differs from X0-FEAS's; it does not yet quantify by how
  much. That is precisely what step 4 of Part C is for, and the argument should
  not be relied on until that study is done.

---

## Sources

- Interaction Concordance Index — https://arxiv.org/abs/2510.14419
- Generalized vec trick for pairwise kernel models — https://arxiv.org/abs/2009.01054
- SPaiK, scalable pairwise kernel learning — https://arxiv.org/abs/2606.16979
- Pahikkala et al., Toward more realistic drug–target interaction predictions — https://www.semanticscholar.org/paper/Toward-more-realistic-drug%E2%80%93target-interaction-Pahikkala-Airola/5118d7d97271bcae11b56ac36c6375925c8ce8e8
- Kramer et al., The Experimental Uncertainty of Heterogeneous Public Ki Data — https://pubs.acs.org/doi/10.1021/jm300131x
- Landrum & Riniker, Combining IC50 or Ki Values from Different Sources Is a Source of Significant Noise — https://pubs.acs.org/doi/10.1021/acs.jcim.4c00049
- Jarquín et al., reaction norm model for genomic selection — https://pubmed.ncbi.nlm.nih.gov/24337101/
- Penalized factorial regression as a reaction norm model for GxE — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11953130/
- AMMI / GGE interaction decomposition — https://pmc.ncbi.nlm.nih.gov/articles/PMC6483959/
- Davis et al., Comprehensive analysis of kinase inhibitor selectivity — https://www.nature.com/articles/nbt.1990
- Karaman et al., A quantitative analysis of kinase inhibitor selectivity — https://www.researchgate.net/publication/5667440_A_quantitative_analysis_of_kinase_inhibitor_selectivity_Translated_from_Eng
- Metz et al., Comprehensive assay of kinase catalytic activity — https://pmc.ncbi.nlm.nih.gov/articles/PMC3230241/
- Elkins et al., Comprehensive characterization of the Published Kinase Inhibitor Set — https://www.nature.com/articles/nbt.3374
- AdaMBind, meta learning and task adaptive DTA — https://www.nature.com/articles/s41467-026-70554-5
- Rethinking generalization of DTA algorithms via similarity-aware evaluation — https://arxiv.org/pdf/2504.09481
