# UBSE-G0R cross-publication binding-residue reliability preregistration

Date: 2026-07-29  
Status: frozen before any binding-residue similarity calculation

## Question

The proposed Universal Binding-State Embedding (UBSE) changes the Stage-1
supervision object from an affinity value to a ligand-conditioned interaction
state. Before downloading coordinates or training a teacher/student, G0R asks
whether the interaction label already present in the firewalled BioLiP2
registry is repeatable and ligand-specific:

> For the same exact protein sequence and ligand connectivity observed in
> independent PubMed records, are reindexed binding-residue sets more similar
> than those of a size-matched wrong ligand on the same exact target?

This is a necessary gate for the existing binding-residue-list substrate, not
a sufficient test of fine-grained atom/residue interaction tokens. Passing
G0R may only request a small coordinate reliability pilot. It cannot authorize
affinity access, UBSE training, or a predictive claim.

## Frozen source and firewall

- Source:
  `dataset/public/biolip2/processed/closed_registry.parquet`
- Required SHA-256:
  `7905e4edf88073f564baa4b2d4fb50d496432bc4e15e97cccbfa0766b1b0638d`
- Allowed columns:
  `target_key`, `conn`, `scaffold`, `heavy_atoms`, `pdb_id`, `pubmed`,
  `binding_residues_reindexed`
- Forbidden:
  `affinity_presence`, all affinity values, development/confirmation
  features or labels, sealed outcomes, and structure downloads.

The D1C registry has already removed ligands occurring on more than 50 exact
targets and collapsed each `(PubMed, exact target, ligand connectivity)` to
one deterministic PDB representative.

## Frozen parsing and units

1. Parse `binding_residues_reindexed` as whitespace-delimited amino-acid
   one-letter code plus integer sequence position. Reject malformed tokens;
   compare only integer positions so residue substitutions cannot create an
   artificial mismatch for a sequence-exact target.
2. An observation must have nonempty exact target, ligand, PDB, PubMed, and
   parsed residue set.
3. A repeated unit is `(target_key, conn)` with at least two distinct PubMed
   IDs and two distinct PDB IDs.
4. The correct similarity of a unit is the median Jaccard similarity over
   all observation pairs with different PubMed IDs and different PDB IDs.
5. For every observation in a repeated unit, form wrong-ligand candidates
   from the same exact target, a different ligand, and a different PubMed.
   Prefer candidates with the same scaffold token when any exist; otherwise
   use all candidates. Retain the five candidates with the smallest absolute
   heavy-atom-count difference, breaking ties by ligand, PubMed, and PDB.
   The observation-level hard-negative score is the maximum residue-set
   Jaccard among those candidates. The unit wrong-ligand score is the median
   observation-level hard-negative score.
6. The ligand-specific margin is `correct - wrong_ligand`. Aggregate one
   margin per repeated unit, then one median per exact target. All uncertainty
   resamples exact targets, not rows or structures.

## Frozen retrieval diagnostic

For each observation in a repeated unit, construct an optimistic
within-target candidate pool from observations in other PubMeds and different
PDB entries. Collapse candidate observations by ligand and assign each ligand
the maximum contact Jaccard to the query. Keep queries with at least five
candidate ligands and at least one true ligand. Rank descending by score with
lexical ligand order as the deterministic tie break.

Report Recall@1, MRR, median candidate count, and the deterministic random
chance values `mean(1 / candidate_count)` for Recall@1 and
`mean(H_n / n)` for MRR. This retrieval uses observed contact labels and is
therefore only an optimistic label-discriminability ceiling; it is not a
deployable sequence/SMILES model.

## Frozen uncertainty

Use seed 1729 and 2,000 target-cluster bootstrap replicates. Each replicate
samples exact targets with replacement and reports the median of their
unit-median ligand-specific margins. The two-sided percentile interval is
reported; admission uses its 2.5th percentile.

## Frozen gates

All gates must pass:

1. **R1 scale:** at least 200 repeated units, 80 exact targets, 100 PubMed
   IDs, and 200 distinct PDB entries.
2. **R2 hard-negative coverage:** at least 160 units and 64 exact targets
   have a legal same-target wrong-ligand control.
3. **R3 repeatability:** median correct cross-PubMed Jaccard is at least
   `0.50`.
4. **R4 ligand specificity:** median unit margin is at least `0.10` and the
   target-bootstrap lower 95% bound is greater than `0.05`.
5. **R5 optimistic retrieval:** at least 200 legal queries/50 exact targets
   with at least five ligand candidates, Recall@1 at least `0.35`, and MRR at
   least `0.50`.
6. **R6 firewall:** no affinity field/value, coordinate archive,
   development/confirmation feature/label, or sealed outcome is read.

Pass:
`REQUEST_UBSE_G0C_SMALL_COORDINATE_RELIABILITY_PREREGISTRATION`.

Failure:
`STOP_UBSE_EXISTING_BINDING_RESIDUE_SUPERVISION_NOT_LIGAND_SPECIFIC`.

Failure stops use of the current BioLiP binding-residue lists as the UBSE
teacher target and forbids the proposed 20k-50k complex training run. It does
not prove that a separately specified fine-grained coordinate-event source
cannot work; such a source would require a new small, preregistered
cross-publication reliability gate rather than threshold relaxation.

## Compute boundary

This gate is CPU-only parsing, set algebra, deterministic matching, and
cluster bootstrap. CUDA becomes relevant only after a later coordinate gate
admits model training.
