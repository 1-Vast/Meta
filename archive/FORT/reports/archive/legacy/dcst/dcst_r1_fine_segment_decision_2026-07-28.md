# DCST-R1 fine-segment decision

Date: 2026-07-28  
Decision: `STOP_R1__TARGET_IDENTITY_AND_COORDINATE_MISMATCH`

## Registered result

The seed-1729, 4,000-step source-only run produced:

| Quantity | P0: 8 segments | R1: 32 segments |
| --- | ---: | ---: |
| True centered alignment | 0.0033 | 0.0200 |
| Wrong-target centered alignment | 0.0434 | 0.0290 |
| Wrong-ligand centered alignment | 0.0063 | -0.0147 |
| Privileged certified bands | 2/4 | 2/4 |
| No-privileged certified bands | 1/4 | 1/4 |

R1 failed its source mechanism gate and therefore did not authorize generation
of a ChEMBL 32-segment cache or any new downstream-label load.

## Root-cause audit

The PLINDER upstream source defines a formatted residue as
`<instance.chain>_<residue_number>_<resolved-residue index>`. The third value
is an enumeration over resolved residues, not a UniProt coordinate. The P0/R1
projection divided that value by the length of an arbitrarily selected
cluster-representative UniProt sequence.

On 59,554 firewalled contact residues:

- the structure-chain/representative-sequence length ratio had median 0.710;
- 66.3% were outside the interval `[0.9, 1.1]`;
- 58.3% changed 8-bin assignment when only the denominator was corrected;
- 74.4% changed 32-bin assignment.

The target identity was also corrupted. PLINDER `cluster` is a split/pocket
similarity grouping, not a protein identifier:

- 64.6% of processed rows lie in a cluster containing multiple accessions;
- only 50.1% of rows contain the accession chosen as that cluster's ESM
  representative;
- one cluster contains 146 accessions.

Increasing positional resolution therefore made an invalid target/coordinate
projection more precise. The R1 failure is not evidence against
structure-privileged pretraining after correct entity alignment.

## Authorized successor

Only the separately preregistered R2 entity-aligned route may proceed. R1
artifacts remain diagnostic and must not be used for a predictive claim.

