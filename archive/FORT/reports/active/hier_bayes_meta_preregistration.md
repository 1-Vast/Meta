# Hierarchical protein-conditioned Bayesian meta-learner — preregistration (2026-07-26)

This registers one structural architecture change and its development diagnostic **before any result is
read**. It is written after `TRANSFORMER_BAYES_META_SHORT_FAIL_REVIEW`
(`reports/active/transformer_bayes_meta_decision.md`) and follows that decision's binding instruction:
the next change must *structurally* strengthen protein-conditioned interaction information, not add
epochs and not add a loss that merely makes destruction controls worse (that contrast loss was already
tried and reverted).

## Localized defect being addressed

In the global-prior run, shuffled protein scored 0.2929 versus true protein 0.2983 and random protein
0.2819; `TBM − protein_shuffle` was `+0.0048 [-0.0031, +0.0135]`. The improvement is a
support-conditioned ligand kernel: a wrong protein barely hurts because the prior over the target
function is **global** and the protein enters only the shared interaction basis, which — even under a
wrong protein — is still a serviceable ligand feature space. Protein identity is statistically
redundant given the support residuals. This is the same redundancy diagnosed in `BM1_RR_FAIL_STOP`
and `PANEL_GATE_PC_FAIL_STOP`.

## Structural change (the only change)

The prior over the target function variable `w_t` becomes **hierarchical**:

```text
w_t | t ~ Normal(0, Sigma_t),   Sigma_t = (L0 + Delta_t)(L0 + Delta_t)^T + diag(exp(d))
Delta_t = reshape( MLP(ESM2_pooled(t)) , rank x rank )
```

* `L0`, `d` are shared across targets (the global prior, retained). `Delta_t` is a low-parameter
  protein-conditioned correction to the **covariance factor**, so the frozen protein orients (rotates
  and reshapes) the full prior covariance instead of only rescaling shared coordinates.
* The prior **mean stays structurally zero**. Therefore the exact `k=0`/`k<=1` fallback to B0 is
  preserved, and the "no signed zero-shot correction" contract (a signed prior mean is a separately
  locked route) is not touched.
* `Delta_t` is initialised at ~0, so the model begins indistinguishable from the global prior and must
  *earn* any protein dependence from the episodic query objective. No adversarial/contrast term and no
  extra epochs are added; the meta-objective, optimiser, capacity elsewhere, seeds, folds, episodes and
  evaluation rows are byte-identical to the global-prior run (`--protein-prior` is the only flag).
* Because the same `Sigma_t` enters the marginal-likelihood Bayes factor, a wrong protein also lowers
  the adaptation weight `pi`, coherently, with no neural gate.

Retained core components: Transformer cross-attention interaction basis over frozen ESM-2 segment
tokens (Transformer); exact FP32 Cholesky posterior with a hierarchical prior (Bayesian); episodic
support/query training of shared parameters (meta-learning). This is the "Hierarchical Residual"
element of the long-horizon BHR-MoT-DTA name. Hierarchical MoT remains out of scope and unauthorised.

## Protocol (unchanged from the global-prior diagnostic)

* Substrate: Metz dense kinase panel `panel_metz`, registry sha256
  `94da6bb5a59c2911672fde982530c8dd6a673c194b2b2d7b4638df7768c8173e`; **spent panel development rows**.
  This is a development diagnostic, not a confirmatory gate. Panel confirmation, Davis and sealed
  labels are not read.
* Five frozen homology-component folds; chemically hard label-blind episodes at `k=4`; matched arms
  B0, treatment, wrong-target support, label-permuted support, protein-shuffle, protein-random on
  identical query rows.
* Single seed (1729). Multi-seed is authorised only after a one-seed pass, per the standing rules.

## Frozen pass criteria (identical thresholds; no relaxation)

Primary — `treatment − B0` paired target-component macro Spearman:

1. `effect >= 0.03` (the task's nominal minimum; panel k=4 MDE80 is 0.0367);
2. grouped component-bootstrap `LCB95 > 0`.

Specificity (the defect under test):

3. protein specificity: `treatment − protein_shuffle` **and** `treatment − protein_random` both
   `LCB95 > 0`;
4. wrong-support specificity: `treatment − wrong_support` `LCB95 > 0`;
5. label specificity: `treatment − label_permuted` `LCB95 > 0`.

Safety and invariants:

6. RMSE `<= 1.02 * B0`;
7. exact `k=0` fallback (max abs deviation 0), support-permutation invariance (`< 1e-4`),
   support-offset invariance (`< 1e-4`), finite positive predictive variance.

Verdict is `HIER_BAYES_META_SHORT_PASS` only if all criteria hold, else
`HIER_BAYES_META_SHORT_FAIL_REVIEW`. A pass authorises only a review of a multi-seed repeat; it does
not authorise Hierarchical MoT, threshold changes, long training, or confirmation/Davis/sealed access.
A partial result (e.g. protein specificity now resolved but the effect still under-powered at this
panel size) is reported honestly as evidence about the mechanism, not relabelled as success.
