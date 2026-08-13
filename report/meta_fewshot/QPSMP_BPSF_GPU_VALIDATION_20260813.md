# QPSMP-BPSF GPU Validation

## Scope

This report records an implementation and development smoke evaluation of the
trainable QPSMP Bipartite Pair-Section Former (BPSF). It is not a confirmatory
Cold Target result and does not authorize G2, G3, or a biological mechanism
claim.

## Integrated model

The active path is:

`protein residue states + ligand atom states -> persistent bipartite pair field -> section latents -> quotient-preserving learned support operator -> scalar affinity`.

The support operator is the primary learnable meta-learning module. It is
permutation invariant, removes constant residuals, returns zero SAR state for
zero centered evidence and one-shot support, and keeps the state in the
centered support row span. The analytic positive ridge remains a comparator.

The geometry head is source-only. It predicts contact and distance labels from
the same pair field during structure pretraining; holo coordinates are never
required by affinity inference.

## GPU and training evidence

- Focused model and interface tests: 42 passed.
- BPSF/SectionFormer endpoint smoke: 20 steps, CUDA AMP, peak memory about
  396 MB, validation MSE 0.8632 at the final step.
- Cached BPSF Stage A: 500 steps, 16 episodes per step, CUDA AMP, about 124
  seconds of optimization, peak memory about 396 MB, GPU utilization observed
  near 78% and power near 76 W.
- Frozen-zero-shot SectionFormer meta-stage: 1,000 steps, best validation
  checkpoint at step 100, about 368 seconds total optimization, peak memory
  about 517 MB.
- Source-only geometry teacher smoke: 32 train and 16 validation complexes,
  one CUDA epoch, validation contact-plus-distance loss 3.0913.

## Development results

On the consumed 18-episode BindingDB development smoke:

- Stage A full MSE: 1.5854; level-only MSE: 1.4516.
- Frozen SectionFormer best-checkpoint full MSE: 1.5536; level-only MSE:
  1.4516; SAR-cut MSE: 1.5924.
- Isolated SAR gain (`SAR-cut - full`): +0.0388 MSE.
- Permuted-support MSE: 1.6625; foreign-state gap: +0.0596 MSE.
- Wrong-protein-state gap: -0.0022 MSE, so target-specific SAR is not
  identified.

The learned section therefore produces a measurable development-only SAR
contrast, but it does not beat the level baseline and does not pass target
specificity controls. The architecture is integrated; the biological and
Cold Target utility gates remain closed.

## Limitations and next gate

The structure corpus has a protein homology split, but the current teacher
manifest does not apply a joint exact-ligand, Murcko-scaffold, and PDB/document
dependency closure. Geometry pretraining must be re-governed before its output
can be used as scientific evidence. The next authorized experiment is a
matched three-seed nested-k comparison of pooled, atom-residue, BPSF, analytic
ridge, and learned SectionFormer arms on untouched confirmation components.
