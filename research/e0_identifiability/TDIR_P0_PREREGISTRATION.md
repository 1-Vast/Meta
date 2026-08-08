# T-DIR-P0 Annotation and Learnability Pilot

Registered: 2026-08-07

## Scope

This is a small, structure-only feasibility pilot. It tests whether open RCSB
coordinates can be deterministically mapped to PLIP interaction labels and
whether fixed frozen-P1B features contain learnable signal for at least one
typed interaction channel. It reads no ChEMBL affinity, DAVIS label, recipient
label or support label. A PASS does not admit a feature to `model/`, `scripts/`
or biological `z`.

## Frozen Inputs

- Governed BioLiP2/RCSB corpus: `pilot20k_homology_split_v2`.
- P1B checkpoint: `p1b_pilot20k_seed17_v1/best.pt`.
- PLIP `3.0.1`; its unavailable Open Babel InChIKey writer is replaced only
  for metadata serialization and is recorded in the manifest.
- CUDA execution in the `drug` environment for frozen P1B inference.

## Label-Blind Selection

Select `24/8/8` train/validation/test complexes before PLIP annotation by the
SHA-256 order of `"TDIR-P0|source_entry_id"`. All 40 complexes must have
different homology groups, PDB IDs and exact protein sequences. Validation and
test Murcko scaffolds must not occur in already selected earlier splits; empty
scaffolds, CCD identifiers longer than three PDB columns and metal-containing
ligands are ineligible because the pilot removes explicit metal environment.
Exact ligand/connectivity overlap is audited. Failed
PLIP records are reported and are not replaced. The frozen P1B checkpoint used
only its governed train split; therefore selected validation and test records
must be P1B-held-out records.

## Pair and Label Contract

Candidate pairs are canonical ligand heavy atom and explicit mapped protein
residue pairs with minimum holo-structure heavy-atom distance `<= 8.0 A`.
Residues retain canonical sequence indices; their P1B slot is
`floor(sequence_index * 128 / sequence_length)`.

PLIP channels are:

1. protein-donor H-bond;
2. ligand-donor H-bond;
3. hydrophobic contact;
4. ligand-negative salt bridge;
5. protein-negative salt bridge;
6. pi stacking;
7. cation-pi;
8. halogen bond.

This is an oracle near-pair candidate rule, so the experiment measures
conditional type learnability rather than deployment end-to-end pair finding.
Multi-atom ligand groups or rings mark each participating canonical heavy atom
and are reported as group-derived weak labels, not direct atom labels.

Hydrophobic contact is the sole primary channel because it has direct atom
semantics and is expected to be sufficiently prevalent in the fixed small
sample. Other channels form a complete preregistered descriptive matrix. A
channel is evaluable only with at least `50` train positives from `8` complexes,
and at least `10` validation and test positives from `3` complexes each, with
both classes present in every split.

## Frozen Arms

- `D0`: frozen P1B contact probability plus five distance probabilities.
- `D1`: D0 plus the 32-dimensional canonical ligand atom feature and a fixed
  six-class residue chemistry one-hot.
- `D2`: D1 plus frozen 128-dimensional ligand-atom and P1B residue-slot states.

Each channel/arm uses a train-only standardizer and scikit-learn logistic
regression with `C=1`, balanced class weights, L2 penalty, fixed seed `1701`
and `max_iter=500`. Every complex has equal total sample weight during fitting
and pooled evaluation. Validation is descriptive only; it performs no model or
threshold selection. Within-complex label shuffle and pair-feature shuffle are
fixed nuisance controls. D0/D1/D2 have different capacities; their differences
are descriptive and are not attributed causally to local-state biology.

## Pilot Success Conditions

The pilot verdict is `PILOT_FEASIBILITY_SIGNAL_OBSERVED` only if:

- at least `36/40` selected complexes annotate and map successfully;
- the preregistered primary hydrophobic channel is evaluable;
- primary D2 AUPRC exceeds prevalence on both validation and test;
- primary D2 AUPRC exceeds D0 on both validation and test;
- primary test D2 AUPRC exceeds its prevalence by at least `0.10`;
- at least `5/8` test complexes individually have D2 AP above prevalence;
- all provenance, selection, mapping and feature artifacts are hash-bound.

Pair observations are correlated; no pair-level bootstrap, confidence interval
or significance claim is permitted. This condition authorizes only a separately
preregistered, adequately powered
T-DIR Gate with corrected derangement and component-bootstrap controls. It does
not establish directionality, affinity semantics, transfer or production use.

## Stop Conditions

No hyperparameter retry, sample replacement, typed-channel deletion based on
test performance, PLIP alternative, real-affinity read, DAVIS evaluation or
production integration is allowed in this stage.
