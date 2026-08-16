# AnchorDelta P0 Interface Audit

## Reuse points

- `model/interaction.py:55-85` provides the frozen label-free feature path:
  `encodeprotein(tokens)` followed by `pairfromprotein(protein, ligand)`.
  A comparator head can consume these features without changing the existing
  posterior API.
- `model/reorder.py:453-477` caches protein and support features in
  `ReorderingState`; `rankfeatures` at `480-495` exposes support-centred
  coordinates, but its centring is not needed for an antisymmetric comparator.
- `scripts/train.py:486-504` (`episodetensors`) is the canonical mapping from
  episode indices to source rows, normalized ligand arrays, labels, and B0.
  Pair training should use the same source-row and normalization convention.
- `scripts/train.py:100-160` (`maketrainroster`) is suitable for episode
  construction, but it intentionally keeps only five support rows and a capped
  query span. P0 pair supervision must instead group all TRAIN rows by target
  and sample a bounded, target-balanced number of ordered pairs.

## Existing prototype and limitation

`research/rankgate.py:95-110` already computes

```text
query_B0 + mean(support_y - support_B0)
          + score(query) - mean(score(support))
```

using `research/pairprior.py`. Its score is an absolute protein-ligand head,
not an explicitly antisymmetric pair operator; PairPrior also has an
independent protein/ligand encoder and only supports k=0 in `predict`. Reusing
it would confound the frozen-feature kill test and should be avoided. Its
`rankgate.v1.json` result is a useful negative control, not a P0 success claim.

## Minimal implementation boundary

Add a separate comparator module (for example `model/anchordelta.py`) and keep
`TargetAdapter`/`ReorderingPosterior` unchanged. The module should:

1. accept frozen pair features `z(p,x)`;
2. compute `h(z_q,z_i)` and return
   `(h(z_q,z_i) - h(z_i,z_q))/2`;
3. expose anchor aggregation `mean_i(y_i + Delta(q,i))` with uniform weights
   as the first control; and
4. train only the small comparator head. The feature extractor must be in
   eval mode with all parameters frozen and all feature tensors detached.

The strict train CLI currently calls `buildreordering(...,
proteinconditioned=args.protein_conditioned)` with the flag defaulting to
false. In that mode `ReorderingModel.interaction` is `None` and no protein
feature exists. A protein-conditioned P0 must explicitly load/build the
conditioned encoder (or use a checkpoint such as `adaptjointfull.v1.pt`); do
not silently treat the global-basis path as protein-conditioned.

## Data and sampling risks

- TRAIN pKi contains approximately 185k rows / 559 targets; all-pairs
  enumeration is quadratic (one target has thousands of rows). Use a fixed,
  target-balanced cap per epoch and include both pair orders with opposite
  labels.
- Pair differences can absorb assay/document effects because the active
  interaction encoder has no source covariates. Record or stratify pair source
  metadata and report same-source versus cross-source bins; do not claim a
  biological relative operator from unstratified edges.
- Keep support/query labels out of feature extraction and support selection.
  At test time labels may enter only as `y_i` in the anchor term. Wrong-target
  and permuted-label controls must preserve the feature tensors and alter only
  target protein or support labels.

## Required tests

- Antisymmetry: `Delta(q,i) == -Delta(i,q)` and `Delta(i,i) == 0`.
- Support permutation invariance for uniform and label-free learned weights.
- Constant support-label shift translates predictions by exactly that constant.
- Correct, wrong-protein, protein-free, wrong-support, and permuted-label
  controls use identical query features and aggregation code.
- Absolute-head and parameter-matched frozen-encoder controls.
- Similarity/source bins and component-bootstrap metrics; reject gains that
  occur only in high Tanimoto or one assay/document bin.

