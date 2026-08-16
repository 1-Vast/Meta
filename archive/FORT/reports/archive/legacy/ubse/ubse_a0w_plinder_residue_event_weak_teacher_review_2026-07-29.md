# UBSE-A0W PLINDER Residue-Event Weak-Teacher Review

**Date:** 2026-07-29  
**Status:** source and semantics review only  
**Outcome access:** no affinity value or protected outcome was loaded

## Decision

The local PLINDER `ligand_interactions` annotation is useful only as a
low-resolution residue-by-event weak label:

\[
M_{tlik} \in \{0,1\}.
\]

It adds event type and a small number of protein-side flags relative to the
binary BioLiP contact labels used by UBSE-G1. It is not an independent
residue-by-functional-group-by-event teacher. It is a lossy marginal of the
same holo-coordinate/PLIP evidence needed by A1 and cannot recover:

- ligand atom identity;
- ligand functional-group identity;
- a residue-functional-group pairing;
- most atom-level geometry and direction;
- a unique BioLiP/mmCIF ligand instance.

The binding disposition is therefore:

`STOP_UBSE_A0W_AS_A1_REPLACEMENT`

A separately preregistered, fit-only A0W auxiliary probe remains admissible:

`ALLOW_PREREGISTER_UBSE_A0W_FIT_ONLY_AUXILIARY_PROBE`

Even a passing auxiliary probe may freeze only a future A1
residue-by-event marginal consistency head. It cannot unlock affinity,
replace UBSE-A0/A1, fabricate a functional-group axis, or support a claim
of a complete typed 3D teacher.

## Audited sources

The local table is:

`dataset/public/plinder_2024_06_v2/raw/annotation_table.parquet`

Its footer reports 1,357,906 rows. The audit streamed explicit identity and
interaction columns only. No affinity field was selected. In particular,
`ligand_binding_affinity` was never loaded.

The upstream implementation was checked at PLINDER commit
`85b3f1cb1763530a6cfd934f4263a1777c41afa4`. In
`src/plinder/data/utils/annotations/interaction_utils.py`,
`get_plip_hash` constructs eight residue-event strings: hydrogen bond,
hydrophobic contact, water bridge, salt bridge, pi stacking, pi cation,
halogen bond, and metal interaction. Ligand atom identifiers are absent,
and several atom/group attributes are explicitly omitted from the stored
hash. The ligand-instance and residue formatting logic also shows that
PLINDER assembly instance identifiers and resolved-chain indices are not
BioLiP A0C `mmcif_auth_seq_id` or canonical UniProt positions.

Primary references:

- [PLINDER paper](https://doi.org/10.1101/2024.07.17.603955)
- [PLINDER repository](https://github.com/plinder-org/plinder)
- [PLIP atom-level interaction paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4489249/)
- [LINKER paper](https://pubs.acs.org/doi/10.1021/acs.jcim.6c00527)

The official PLINDER repository also warns that the 2024-06/v2
`ligand_binding_affinity` annotation has a BindingDB parsing bug. This
review did not use that field for any purpose.

## Reproduction contract

The A0 rows were corrected in memory using the accepted A0C
BioLiP-column-20 locator semantics. The PLINDER scan selected only:

- exact PDB entry;
- exact ligand CCD;
- exact ligand author chain;
- exact receptor author chain membership;
- non-stereochemical RDKit connectivity equality;
- non-empty `ligand_interactions`.

For each candidate system, receptor assembly chains were restricted to the
A0 receptor author chain. `ligand_interacting_residues` mapped the stored
resolved-chain residue index, and interactions were collapsed to an
`(index, event_type)` multiset. A row received a consensus weak label only
when every candidate produced exactly the same non-empty multiset.

The stream read only identity, protein-chain, ligand-identity,
`ligand_interactions`, and `ligand_interacting_residues` columns. A split
diagnostic may read only `system_id` and `split`; it must not read
`system_has_binding_affinity`.

## Coverage

| A0 role | rows | any candidate | exactly one candidate | consensus non-empty typed signature |
|---|---:|---:|---:|---:|
| fit | 3,130 | 2,897 | 2,287 | 2,587 (82.65%) |
| validation | 140 | 123 | 84 | 98 (70.00%) |
| audit | 197 | 168 | 121 | 140 (71.07%) |
| total | 3,467 | 3,188 | 2,492 | 2,825 (81.48%) |

Some A0 rows had as many as 93 PLINDER candidates. Examples with
candidate-dependent signatures include PDB entries `5san`, `6mb5`, `6mb7`,
and `1jik`. Consequently, the high "any candidate" coverage is not exact
coordinate-instance coverage.

Blind use of the full PLINDER table is also forbidden:

- it directly contains labels matching A0 validation and audit rows;
- PLINDER's native split is not aligned to A0 roles;
- the annotation table has no PubMed/citation column from which to enforce
  the existing source closure;
- its stored residue index is not automatically a canonical target-sequence
  coordinate.

## Admissible A0W probe

A0W may run during the remote-coordinate wait only after a separate
preregistration. It must use A0 fit rows only and create a new internal held
split closed on PDB, exact target, homology, scaffold, and PubMed. Official
A0 validation and audit event labels remain unread during selection,
training, and threshold setting.

The minimum deployable inputs are protein sequence or a predicted monomer,
the frozen P0A proposal if P0A passes, and the ligand 2D graph. The primary
comparison is against an exact null containing the target-pocket marginal
and a position-free typed pair burden. Required controls are:

- wrong ligand;
- wrong protein;
- event-type shuffle;
- residue-position or monomer-position destruction;
- a parameter-matched position-free model.

Before GPU training, at least 2,000 consensus-labeled fit rows and 50
internally held contrast panels must survive source closure, and the
fixed-measure double-centered typed residual must be nonzero.

All of the following are binding:

1. Cross-model panel macro AP and centered cosine must improve over the exact
   null with target/panel bootstrap lower 95% bounds above zero.
2. Correct-minus-wrong-ligand and correct-minus-wrong-protein margins must
   each be at least 0.05 with lower 95% bounds above zero.
3. The position-free model must not outperform the position-aware model.
4. Three seeds must be finite and stable.
5. No source-closure or outcome firewall may fail.

Failure gives:

`STOP_UBSE_A0W_NO_INCREMENTAL_TYPED_RESIDUAL`

A pass does not change the A0 remote-coordinate wait or authorize Stage-2.

## Novelty boundary

LINKER already maps PLIP-derived residue-by-functional-group-by-event labels
to a sequence-plus-SMILES student and uses the representation downstream.
Neither a PLINDER residue-event weak head nor generic
"structure-supervised, structure-free inference" is a defensible novelty
claim.

The surviving A1 contribution must be evaluated as the complete combination
of source closure, independently certified event reliability, deployment-side
predicted monomer structure, partial/unbalanced OT, an exact typed
pair-burden null, fixed-measure purification, and strict dual-cold
placement-identifiable rectangles. A parameter-matched LINKER-form baseline
remains mandatory.
