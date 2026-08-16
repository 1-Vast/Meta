# FACTOR-C F0-C1 preregistration

Date: 2026-07-26  
Route: user-supplied predictive self-supervised carrier audit; no agent candidate slot consumed.  
Role: unlabeled ligand-adapter experiment only. F1-C remains locked.

## Authorization basis

F0-C0 ended as `FACTOR_F0C0_REPRESENTATION_UNIDENTIFIED_STOP`, not as a real chemical-support
failure. Its structural reconstruction, effective rank and true-minus-decoy gap passed, but fixed
features failed inner scaffold-OOD coverage and common neutral atoms made the atom decoy calibration
non-identifying. This is the preregistered condition under which F0-C1 is permissible.

## Strict three-fold corpus firewall

The KIRHub2026, Reinecke2024 and Papyrus-Christmann2016 molecular graphs are the only corpus. For each
leave-one-source-out fold:

1. delete every candidate training molecule whose parent connectivity or Murcko scaffold occurs in
   the held source;
2. split remaining scaffolds 80/20 with the frozen F0-C hash rule;
3. train only on inner-train graphs;
4. use inner-validation solely for masked-structure checkpoint choice, pseudo-OOD bandwidth
   calibration and anti-cheat gates;
5. encode the held source once after the checkpoint and bandwidths are frozen.

No activity/inhibition value, target identity, protein feature, source ID, document ID or held-source
structure is supplied to training.

## Frozen encoder and objective

- four residual GINE message-passing layers;
- hidden width 128, edge width 16, LayerNorm and dropout 0.10;
- categorical atom inputs: element class, formal-charge bucket, degree, aromaticity, hybridization and
  ring state;
- categorical edge inputs: bond order, conjugation and ring state;
- connected-subgraph masks: four SHA-256 masks per molecule, 15% of heavy atoms, cycled by epoch;
- reconstruction heads: masked element/charge/degree/aromaticity/hybridization, directed bond order
  and per-atom BRICS attachment bin;
- loss = mean masked-atom cross entropy + 0.5 bond cross entropy + 0.25 attachment cross entropy;
- AdamW, learning rate 1e-3, weight decay 1e-5, batch size 64, 40 epochs, seed 1729;
- checkpoint = minimum inner-validation structural loss; external coverage cannot choose it.

This is a small predictive masked-graph encoder, not a DTA model. No contrastive augmentation,
Transformer, affinity head or protein path is allowed.

## Frozen multiresolution carrier

The frozen encoder emits L2-normalized atom vectors. Pharmacophore-pair carriers concatenate
role-ordered endpoint vectors with topological distance and same-ring indicators. BRICS motifs pool
mean and maximum atom vectors plus size/ring/attachment summaries. At most 96 pair carriers per
molecule are retained by the same label-blind SHA-256 cap used in C0. This cap is an audit bound only;
if the mechanism is proven, an uncapped high-memory replication is required before final training.

Coverage remains equal-weighted across atom, pair and motif levels and rarity-weighted within level.

## Non-null chemistry-broken decoys

Every evaluated decoy carrier must differ from its true carrier; unchanged records are excluded and
their fraction is reported.

1. role derangement: cyclically replace every role among at least two supported roles;
2. graph chemistry corruption: cyclically change every formal-charge bucket and bond-order category,
   re-encode the corrupted graph, and retain only changed vectors;
3. environment/motif mismatch: cyclically derange atom embeddings within each molecule before pair
   construction and motif pooling; retain only changed vectors.

These decoys preserve graph size and carrier counts. They are invalid chemistry controls, not
additional training examples.

For each level-role, `tau` is fixed on inner-validation exactly as in C0 from the 5th percentile of
changed-decoy nearest-atlas distances, targeting at most 5% similarity >=0.5. Sparse roles inherit
only a train-fold level-wide bandwidth.

## Structural and anti-collapse gates

Every fold must satisfy:

- masked element accuracy >= frequency baseline +0.15;
- directed bond-order accuracy >= frequency baseline +0.10;
- BRICS attachment-bin accuracy >= frequency baseline +0.10;
- 5-NN atom-role macro-F1 from frozen embeddings >= frequency baseline +0.15;
- atom-embedding participation rank >=16 and rank/nonconstant-dimension >=0.10;
- all losses and gradients finite;
- changed-decoy fraction >=0.95 for every decoy family;
- calibration false coverage <=0.05;
- inner scaffold-OOD molecule coverage median >=0.85 and q10 >=0.60.

## External gates

Sources receive equal total weight and molecules equal weight within source. All must pass:

1. median functional coverage >=0.90;
2. q10 functional coverage >=0.70;
3. source-stratified true-minus-mean-decoy coverage LCB95 >0 (10,000 draws, seed 1729);
4. no source weight >0.40;
5. inherited primitive graph spans all sources and grouped MDE80 <=0.03;
6. all structural, rank, changed-decoy and inner pseudo-OOD gates pass;
7. confirmation labels remain unread and sealed test remains unconsumed.

Pass: `FACTOR_F0C1_PASS_AUTHORIZE_F1C`.  
Coverage failure with all representation/decoy gates valid:
`FACTOR_F0C1_REAL_CHEMICAL_SUPPORT_FAIL`.  
Any representation/decoy/calibration failure:
`FACTOR_F0C1_REPRESENTATION_UNIDENTIFIED_STOP`.

No additional architecture, epoch, seed, bandwidth, public pretraining corpus or threshold may be
introduced after this run.

