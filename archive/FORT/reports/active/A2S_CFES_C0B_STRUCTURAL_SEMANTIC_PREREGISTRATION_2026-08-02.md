# A2S-CFES C0B Structural Semantic Gate Preregistration

Date: 2026-08-02  
Branch: `research/a2s-conformational-free-energy-state-20260802`  
Status: **frozen before feature materialization or model fitting**

## Decision question

Gate C0B asks whether the locally available, outcome-free PLINDER structure
annotations contain a transferable ligand-pocket interaction signal. The
admitted object is deliberately narrower than the final CFES claim:

> On held protein, pocket, ligand, and PDB-provenance groups, does an explicit
> ligand-by-pocket term predict the observed eight-type contact profile better
> than matched additive and protein-free predictors, and is that incremental
> effect destroyed when the physical ligand-pocket pairing is broken?

A pass authorizes a fit-only affinity representation gate. It does not show
that affinity is a mixture over conformational states, that a state-population
logit is identifiable from k<=5 labels, or that a learned support operator is
better than analytic inference.

## Literature-constrained rationale

The biological premise is supported but not sufficient:

- dynamic conformational ensembles are a general mechanism of biomolecular
  recognition (Boehr, Nussinov and Wright, 2009,
  `doi:10.1038/nchembio.232`);
- conformer populations contribute thermodynamically to ligand recognition,
  so a single best structure need not be the relevant observation (Suruzhon
  et al., 2024, `doi:10.1021/acs.jcim.4c01612`);
- ensemble docking is established and its value depends on selecting useful
  conformations rather than merely adding structures (Amaro et al., 2019,
  `doi:10.1021/acs.jpcb.8b11491`; Gan et al., 2022,
  `doi:10.1007/s10822-021-00433-2`);
- PLINDER supplies interaction annotations and similarity-aware structural
  splits suitable for testing whether protein information is load-bearing
  (`doi:10.1101/2024.07.17.603955`).

These results support the plausibility and breadth of the mechanism, but none
establishes few-shot population adaptation. C0B therefore tests the first
necessary link only: transferable state-specific interaction semantics.

## Outcome firewall

The run may read only these raw PLINDER annotation columns:

- `system_id`;
- `entry_pdb_id`;
- `system_id_no_biounit`;
- `system_pocket_UniProt`;
- `ligand_rdkit_canonical_smiles`;
- `ligand_neighboring_residues`;
- `ligand_interactions`;
- `ligand_is_proper`.

It may read only these split columns:

- `system_id`, `uniqueness`, `split`, `cluster`,
  `cluster_for_val_split`;
- `system_pass_validation_criteria`,
  `system_pass_statistics_criteria`;
- `system_proper_num_ligand_chains`,
  `system_proper_pocket_num_residues`,
  `system_proper_num_interactions`.

No column containing `affinity`, no affinity-presence flag, and no PLINDER
`test` or `removed` row may be loaded. The processed PLINDER registry is
prohibited because its affinity column was exposed during schema exploration.
ChEMBL affinity, source `probe`, source `locked`, and recipient outcomes remain
sealed. Source split metadata may be read only to enumerate all 467 A2S target
accessions, which are excluded from PLINDER before model fitting.

UniProt amino-acid sequences are label-free public covariates. Missing local
sequences may be fetched from the UniProt REST endpoint, cached with a content
hash, and used only to map PLINDER residue indices to amino-acid identities.

## Structural split

Only official PLINDER `train` and `val` rows that pass both quality flags,
contain one proper ligand, a non-empty pocket, non-empty interactions, a valid
canonical molecule, and a resolvable protein sequence are eligible.

The official `val` rows form the untouched C0B audit. Before any feature
normalization or fit, purge from `train` every row sharing any of the following
with an audit row:

- UniProt accession;
- exact canonical SMILES;
- Bemis-Murcko scaffold;
- PDB entry or no-biounit system identifier;
- `cluster`, `cluster_for_val_split`, or `uniqueness`.

This makes the evaluated relation exact-protein-, exact-ligand-, scaffold-,
pocket-cluster-, and provenance-disjoint. Five deterministic audit folds are
formed by hashing `cluster_for_val_split`; they are reporting strata, not
model-selection folds. A deterministic 15% hash split of the purged PLINDER
train clusters is the only inner-validation role.

## Representation and endpoint

The ligand vector contains fixed RDKit physicochemical descriptors and a
Morgan fingerprint compressed on purged train only. The pocket vector contains
the normalized composition of neighboring amino acids, residue chemistry
groups, pocket size, and coarse normalized sequence-position occupancy. No
target identifier, PDB identifier, split cluster, contact count, or interaction
annotation is an input.

For each row, the structural label is

\[
y_j=\log(1+n_j),\quad j\in\{\text{hydrogen bond, hydrophobic,
water bridge, salt bridge, pi stack, pi-cation, metal, halogen}\}.
\]

Each coordinate is standardized with purged-train statistics. The primary
per-row loss is the mean squared error across the eight standardized contact
coordinates. This macro-standardization prevents hydrogen bonds from defining
the result merely because they are common.

## Matched predictors

All fitted arms use identical train rows, normalization, seeds, optimizer
budget, and early-stopping role.

1. `ligand_only`: a compact ligand predictor.
2. `additive`: the sum of compact ligand and pocket predictors.
3. `cross`: the frozen additive predictor plus a rank-16 bilinear residual
   `W_o[(W_l l) * (W_p p)]`.
4. `no_cross_capacity`: the frozen additive predictor plus a parameter-matched
   residual that receives separate ligand and pocket projections but no
   multiplicative term.
5. `frozen_random_cross`: the same bilinear interface with frozen random
   projections and only the output map fitted.

The nesting is load-bearing: the cross residual is trained only after the
additive predictor is frozen. Increasing rank, depth, encoder size, or epoch
budget after seeing C0B is prohibited.

## Physical destructions

Each destruction is deterministic per seed and preserves row count and target
labels:

- `pocket_shuffle`: replace the pocket only inside the cross term with a
  size-matched pocket from another accession;
- `ligand_shuffle`: replace the ligand only inside the cross term with a
  heavy-atom-matched ligand of another scaffold;
- `structure_transplant`: replace the full pocket vector with a size-matched
  pocket from another accession and compare with the correspondingly
  transplanted additive arm;
- `residue_randomization`: permute amino-acid identities over the selected
  pocket residue positions before constructing the pocket vector;
- `state_duplication`: duplicating an identical state must be an exact no-op
  under normalized state pooling;
- `state_order_permutation`: reordering state inputs must be an exact no-op.

The last two are interface invariance tests. They cannot create evidence for an
ensemble when only one state is scored.

## Statistics

For model `m`, let `L_m(g)` be mean row loss in audit cluster `g`. The primary
incremental effect is

\[
G=\operatorname{mean}_g[L_{additive}(g)-L_{cross}(g)].
\]

Confidence intervals use 10,000 paired bootstraps over
`cluster_for_val_split`, seed 1729. Five training seeds
`{1729,1731,1733,1741,1753}` are aggregated within cluster before bootstrap.
The same paired unit is used for all model and destruction contrasts.

For destruction `d`,

\[
R_d = 1-G_d/G,
\]

where `G_d` is the cross-minus-additive effect under the matched destruction.
Negative destroyed gains are retained rather than clipped.

## Positive control

Before real training, run a synthetic generator with the measured dimensions,
sample count, group-size distribution, and noise scale, where the target
contains a known rank-4 ligand-by-pocket term. The harness passes only if:

- the lower 95% bound of cross-minus-additive gain is positive;
- all five audit-fold gains are positive;
- pocket and ligand shuffling each remove at least 70% of the point effect.

Failure stops the run as `C0B_HARNESS_INVALID`; it is not biological evidence.

## Admission criteria

C0B passes only if all conditions hold on untouched official validation rows:

1. `cross` beats `additive` with paired cluster-bootstrap lower 95% bound > 0;
2. `cross` beats `no_cross_capacity`, `frozen_random_cross`, and
   `ligand_only`, each with lower 95% bound > 0;
3. the cross-minus-additive point effect is positive in every one of the five
   deterministic audit folds and in every training seed;
4. at least four of eight contact coordinates improve, including both hydrogen
   bond and hydrophobic coordinates;
5. pocket shuffle, ligand shuffle, structure transplant, and residue
   randomization each remove at least 70% of the cross increment;
6. state duplication and state-order permutation are numerical no-ops;
7. all forbidden-column, split, overlap, and A2S-accession checks pass.

The verdict is either
`CFES_C0B_SEMANTICS_ADMITTED_PROCEED_C1` or
`CFES_C0B_SEMANTICS_NOT_ADMITTED_STOP_CFES`. A failed real-data gate may not be
rescued with extra capacity or threshold changes. Per the branch charter, the
next branch must start from a different general biological principle.

## Generalizability and value decision before execution

CFES remains **biologically generalizable in principle** because ensemble
population shifts are not family-specific and C0A measured broad source
structure coverage. It remains **statistically unverified** because neither a
held structural interaction term nor a few-shot population state has yet been
admitted. Its potential value is **substantial but conditional**: a successful
state mechanism would intervene on chemically distant queries through a shared
physical target state, directly addressing the limit measured by TRACE. No
breakthrough wording or code promotion is permitted before C0-C3 all pass.
