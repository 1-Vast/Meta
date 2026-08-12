# R0-C Preregistration: Fresh Exact-Distance Confirmation

Status: frozen on 2026-08-11 after coordinate-free metadata selection and
protein closure, before downloading any R0-C coordinate or CCD file and before
reading any R0-C geometry value.

R0-B stopped before fitting because its 53 heldout components gave
`MDE80=0.0074377 > delta*=0.00616076`. R0-C does not lower either threshold.
It supplies a fresh, dependency-closed structural confirmation population for
the exact residue--atom residual registered in
`PREREG_R0_EXACT_DISTANCE_RESIDUAL.md`.

## 1. Frozen source and exposure boundary

The immutable BioLiP2 metadata snapshot is:

- `BioLiP.txt.gz`: `c92229bbc8c55c3bd84a9813c3e278ba62f4cfa44e6315cc98d9bf63ed64b6ec`;
- `ligand.tsv.gz`: `7d107cfaf0cc873d393e9dcee7bbe16fc7dc024a9748495a92adf6e7658f61ec`;
- `protein.fasta.gz`: `f1c1ba7b9838b8b8cde9de3b231440f2cadaca986f85e1c1af9a26e0b9617cd0`.

The exposure registry is the complete 20,000-entry coordinate-free P1B
candidate pool, not only the 14,906 governed P1B records. Any exact source
entry, PDB ID, protein sequence, or ligand CCD ID present in that pool is
excluded. R0-B's 2,845 scored records are a subset of this registry.

This is geometry-blind confirmation within the same frozen BioLiP snapshot. It
is not a temporal holdout and must not be described as one.

## 2. Frozen coordinate-free selection

`prepare_r0c_candidates.py` used namespace `R0C-FRESH-CANDIDATE-v1` and read no
coordinates or affinity values. It hash-selected 3,000 records with mutually
distinct PDB IDs, exact protein sequences, and ligand CCD IDs.

- candidate file SHA256:
  `41e85726746a134c47e40724ef3bede0345e098eaa1852fec15f99634b9e333e`;
- eligible records after old-exposure metadata filters: 91,864;
- selected records/sequences/PDBs/CCD IDs: 3,000/3,000/3,000/3,000.

Protein dependency closure used MMseqs2 candidate edges at identity 0.30,
followed by frozen parasail local-alignment confirmation requiring identity at
least 0.40 and aligned length at least 0.80 of **both** full sequences. The
same all-edge rule was applied candidate-to-exposure and candidate-to-candidate;
representative-only clustering is forbidden.

The closure produced 1,420 internal components. It identified 1,734 exposed
candidate sequences in 600 components, leaving 820 fresh components. Namespace
`R0C-DOWNLOAD-REP-v1` selected one representative from each of 512 components.

- predownload panel SHA256:
  `754c4a9980e535d447a9da541d79b7aac4b8109af8c8b3e40bdbb3fc0d5757ee`;
- MMseqs2 SHA256:
  `ad319fc9fc67500ca2da9b9247d9611d2113e66cb800c15b603949f5652e2a0b`;
- confirmed exposure/internal edges: 36,865/7,551.

No rejected representative may be replaced after coordinate or chemistry QC.

## 3. Download, exact mapping, and chemistry quarantine

Only the frozen 512 representatives may be downloaded from RCSB. Every file
URL, byte count, SHA256, and failure is persisted. Admission uses the existing
P1B structural contract without relaxation:

- X-ray structure at resolution at most 3.0 Angstrom;
- protein length 50--1022 and exact coordinate-to-sequence mapping coverage 1.0;
- noncovalent regular ligand with 6--96 heavy atoms;
- one-to-one equality between canonical CCD heavy-atom names and coordinate
  atom names;
- no duplicate active atom or residue index.

New CCD graphs use the frozen P1B featurizer and are hash-audited. Before any
distance score is read, a record is excluded if its CCD hash, exact canonical
connectivity, or nonempty Murcko scaffold occurs in the complete 15,645-record
P1B ungoverned structural corpus. Within R0-C, records sharing any of those
chemical keys are dependency-connected; namespace `R0C-FINAL-REP-v1` retains
one record per resulting component. Empty scaffolds never form a shared key.

At least 120 final scorable components, zero dependency straddles, and largest
component share below 0.20 are required. Otherwise the terminal verdict is
`R0C_NOT_RUN_FAIL_CLOSED`.

## 4. Pre-fit audit and fitting boundary

The frozen P1B checkpoint and its five ordered distance bins remain the prior.
R0-C first runs the unchanged R0 pre-fit audit on the final confirmation
population:

- `delta*=0.05*S_prior`;
- `S_add_star >= delta*`;
- registered component-bootstrap `MDE80 <= delta*`;
- every bin supported, exact mappings complete, and N3 movable fraction at
  least 0.50.

Any failure stops before GPU fitting. If all preconditions pass, Full, N1 and
N2 are trained only on the already frozen R0-B train split and selected only
on its validation split, using the same three seeds and optimization budget.
R0-C is opened once for the seed-ensemble G1/G2/G3 and NLL guard. No R0-C value
may choose an epoch, hyperparameter, seed, temperature, or representation.

All parent R0 nulls, estimands, Gates, margins, and failure verdicts are
unchanged. A complete PASS records only
`R0C_EXACT_DISTANCE_RESIDUAL_CONFIRMED`. It authorizes a separately
preregistered measured-affinity R1; it does not authorize Meta-Section
integration or production migration.
