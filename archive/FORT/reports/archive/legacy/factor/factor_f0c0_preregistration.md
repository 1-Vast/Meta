# FACTOR-C F0-C0 preregistration

Date: 2026-07-26  
Route: user-supplied FACTOR-C; this does not consume an agent candidate slot.  
Role: fixed, label-blind chemical-carrier audit. No affinity/inhibition value is read and no
interaction model is trained.

## Frozen question

Does the discrete F0 failure reflect exact-token aliasing, or is ligand-side chemical support truly
insufficient across the three public development sources?

The original F0 remains failed and its F1 remains locked. F0-C0 changes only the ligand adapter. The
KLIFS anchor adapter, sources, source balancing, 0.90/0.70 external coverage gates, document graph and
grouped-MDE contract are unchanged.

## Data and strict corpus firewall

Sources: KIRHub2026, Reinecke2024 and Papyrus-Christmann2016, using the exact metadata filters from
F0. Papyrus-ChEMBL31 and the quarantined ChEMBL confirmation chain remain excluded.

For each leave-one-source-out fold:

1. remove from the candidate atlas every molecule with parent connectivity or Bemis–Murcko scaffold
   present in the held source;
2. split the remaining training scaffolds 80/20 by SHA-256 (`seed namespace=factor-c0-inner-v1`);
3. fit scalers, distance bandwidths and structural proxies on inner-train only;
4. use inner-validation only for pseudo-OOD calibration;
5. freeze the adapter and evaluate the held source once.

No held-source structure, scaffold frequency or coverage result may choose a distance scale.

## Fixed multiresolution carrier

Every molecule is decomposed into three levels.

### A. Atom-role carriers

One carrier is emitted for each heavy-atom × RDKit role. Roles are Donor, Acceptor, Aromatic,
Hydrophobe, PosIonizable, NegIonizable, HalogenDonor, MetalBinder or Structural. The fixed vector
contains element class, formal charge, degree/valence, aromaticity, hybridization, ring-size flags,
one-hop bond-order counts, and two-/three-hop element counts. Values are standardized by the
inner-train median and IQR.

### B. Pharmacophore-pair carriers

Unordered role pairs at topological distance 1–8 are represented by the two role-ordered atom vectors
plus distance, same-ring and shortest-path bond summaries. Pair identity contains only the two roles,
not scaffold, molecule, source or assay.

Pre-result engineering amendment: the first execution produced no output and approached the host
memory limit because all atom pairs were materialized. To prevent molecule-size weighting and bound
memory, at most 96 pair carriers per molecule are retained by SHA-256 of the label-blind role/vector
record. The cap was frozen before any F0-C0 metric existed; no feature, bandwidth or gate changed.

### C. BRICS-motif carriers

BRICS is used only to define a pooling region. The carrier is the mean and maximum of its atom
vectors plus heavy-atom count, ring count and attachment count. Exact fragment SMILES is not used as
an identity.

Each level contributes one third of molecule coverage. Within a level, carrier weights are
`1/sqrt(inner-train role frequency)` and are normalized per molecule.

## Train-only bandwidth calibration

For every level-role with adequate support, nearest-neighbour distances are measured in standardized
inner-train space. Three deterministic chemistry-broken pseudo-OOD decoys are generated from
inner-validation:

1. cyclic role permutation while keeping carrier vectors fixed;
2. formal-charge and bond-summary permutation across carriers while preserving marginal counts;
3. atom-environment/motif-attachment mismatch by permuting pooled endpoint or attachment blocks.

For a decoy nearest distance `d`, similarity is `exp(-d^2/(2*tau^2))`. `tau` is fixed from the 5th
percentile of the pooled decoy-distance distribution so that at most 5% of calibration decoys reach
similarity 0.5. Sparse roles inherit the level-wide train-only `tau`; no held-source fallback is
allowed.

Inner pseudo-OOD feasibility requires molecule coverage median >=0.85 and q10 >=0.60. Failure stops
the corresponding outer fold and therefore fails F0-C0.

## Anti-collapse and reconstruction proxies

Using only inner-train structural labels and inner-validation scaffold-OOD carriers:

- 5-NN masked atom-role macro-F1 must exceed the frequency baseline by >=0.15;
- 5-NN incident-bond-signature accuracy must exceed the majority baseline by >=0.10;
- 5-NN motif-attachment-bin accuracy must exceed the majority baseline by >=0.10;
- the atom-vector covariance participation ratio divided by nonconstant dimension must be >=0.25.

These are representation checks, not activity prediction.

## External metrics and inference

For every held-source molecule, functional coverage is the rarity-weighted mean nearest-atlas
similarity, equally averaged over A/B/C. Molecules are equally weighted inside a source and sources
are equally weighted. True and chemistry-broken coverage are paired by molecule. The true-minus-decoy
LCB95 uses 10,000 source-stratified bootstrap draws with seed 1729.

Samples are reported as:

- `supported_composition`: coverage >=0.70;
- `low_support`: coverage <0.70;
- the old discrete `novel_primitive` flag is retained only for comparison.

## Frozen F0-C0 gates

All must pass:

1. every outer fold passes inner pseudo-OOD median >=0.85 and q10 >=0.60;
2. external source-balanced median functional coverage >=0.90;
3. external source-balanced q10 functional coverage >=0.70;
4. calibration decoy false-coverage rate <=5% in every fold;
5. source-stratified true-minus-decoy coverage LCB95 >0;
6. all three reconstruction margins pass in every fold;
7. effective-rank ratio >=0.25 in every fold;
8. no source contributes >40% effective weight;
9. the original primitive graph remains connected across all sources and grouped MDE80 remains
   <=0.03;
10. confirmation labels remain unread and the sealed test remains unconsumed.

Pass: `FACTOR_F0C0_PASS_AUTHORIZE_F1C`.  
Fail with valid anti-cheat checks: `FACTOR_F0C0_REAL_CHEMICAL_SUPPORT_FAIL`.  
Fail because calibration/reconstruction/rank is invalid:
`FACTOR_F0C0_REPRESENTATION_UNIDENTIFIED_STOP`.

F0-C1 is not automatically authorized by a C0 failure. It is permissible only when failure is
specifically attributable to semantic aliasing while the corpus-support question remains
identifiable.
