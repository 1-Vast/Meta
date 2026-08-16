# LPCT one-seed architecture probe — frozen before scoring

Date: 2026-07-26. This is a mechanism-development run on the already-spent Metz panel development
rows. It is not independent confirmation and cannot authorize a sealed-test claim.

## Hypothesis and innovation

The localized failure is query-dependent within-target reordering: the Bayesian adapter learns a
support-conditioned ligand kernel, while protein-conditioned priors/subspaces repeatedly fail their
shuffle and random-protein controls. The proposed **Learned Pairwise Contrast Transport (LPCT)**
therefore makes observed support SAR contrasts, rather than a protein prior, the task-adaptation
primitive.

For every unordered support pair `(i,j)`, LPCT constructs an order-invariant oriented SAR direction

```text
c_ij = EncPair(
    sign(r_i-r_j) * (z_i-z_j),
    |z_i-z_j|,
    (z_i+z_j)/2,
    tanh(|r_i-r_j| / learned_scale)
)
```

where `z` is an end-to-end learned ligand representation weakly FiLM-conditioned on the frozen
protein sequence representation, and `r` is the cross-fitted base residual. Learned query-to-pair
attention and a learned value network transport these contrast tokens to each query ligand and
directly output its affinity residual. Thus LPCT is a parameterized predictive module optimized by
episodic query loss, not a fixed rule, retrieval step, static concatenation, threshold, or closed-form
post-processor.

The structural invariants are exact `k<=1` fallback, support-row permutation invariance and invariance
to an additive shift of all support labels. No target id, split label, query label, confirmation label
or assay/source identifier enters the model.

## Frozen run

- registry SHA-256: `94da6bb5a59c2911672fde982530c8dd6a673c194b2b2d7b4638df7768c8173e`;
- five existing homology-component folds and matched hard episodes, seed `1729`;
- base steps `600`, inner component-cross-fit steps `300`, LPCT meta steps `800`;
- width `64`, AdamW learning rate `1e-3`, weight decay `1e-5`;
- unchanged ranking + centered-magnitude + `0.05` base-orthogonality query objective;
- no CBTU scheduling, protein prior, contrast margin, extra epoch, confirmation/Davis/sealed access;
- output: `reports/active/lpct_short.json`.

Command:

```text
conda run -n drug python research/transformer_bayes_meta.py --pairwise-transport
  --base-steps 600 --crossfit-steps 300 --meta-steps 800 --seed 1729 --device cuda
  --output reports/active/lpct_short.json
```

## Frozen success criteria

LPCT is retained for an independent-panel gate only if all are true:

1. paired target-component macro-Spearman gain over its in-run B0 is at least `0.03`, with grouped
   bootstrap LCB95 above zero;
2. LPCT point Spearman exceeds the retained TBM value `0.2983` on the identical folds/episodes;
3. RMSE is no worse than `1.02 * B0`;
4. correct support beats wrong-target support and label-permuted support with LCB95 above zero;
5. true protein beats shuffled and random protein with LCB95 above zero, or LPCT is explicitly
   reclassified as a protein-free support-adaptation model and cannot claim protein-specific transfer;
6. negative-transfer rate is below retained TBM's `0.356`;
7. all three structural invariants and finite positive predictive variance pass.

The historical zero-shot I0 (`0.3278` Spearman under its different full-query protocol) is reported
as a broader performance reference, not treated as a paired estimand. Failure of any mandatory
criterion records `LPCT_SHORT_FAIL_STOP`; no hyperparameter rescue is permitted on these spent rows.
