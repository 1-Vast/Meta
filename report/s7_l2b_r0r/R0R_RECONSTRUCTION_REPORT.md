# S7_L2B_R0R — reconstruction report

Date: 2026-08-09.
Stage: `S7_L2B_R0R_RECONSTRUCT_DATA_BASELINE_AND_CONTRACT`.

```text
R0R-1  MONN raw provenance reproduction ......... PASS
R0R-2  New ligand identity and closure .......... FAIL CLOSED
R0R-3  Publication/time closure ................. NOT REACHED
R0R-4  B4 definition and materialization ........ NOT REACHED
R0R-5  Unified preregistration freeze ........... NOT REACHED
R0R-6  GPU/runtime preflight .................... NOT REACHED

TERMINAL VERDICT:  NEW_CLOSURE_TOPOLOGY_INSUFFICIENT
```

Execution stopped at the earliest failed boundary. No closure was relaxed to
continue, no B4 was fabricated, no GPU work was performed.

## R0R-1 — PASS, with one serious acquisition defect found and repaired

MONN was cloned and pinned to `f2b62ccf49c18a9502aa0eb0d582c6e0735ef200`
(HEAD verified). Licence recorded verbatim: *"The algorithm and data can be used
only for NON COMMERCIAL purposes."* Use is declared non-commercial research; the
clone lives under `dataset/raw/monn/`, which `.gitignore` excludes, so **no MONN
byte is committed to this repository**.

### The defect

**All six source files failed their SHA-256 on first clone**, each *larger* than
the manifest size. Cause: `core.autocrlf=true` in the inherited Windows Git
configuration, and MONN ships no `.gitattributes` marking its pickles binary. Git
applied LF→CRLF to binary pickle files.

| File | Corrupted bytes | Expected bytes |
|---|---:|---:|
| `out7_final_pairwise_interaction_dict` | 58,935,197 | 53,017,418 |
| `independent_dataset_interaction_dict` | 4,719,169 | 4,234,610 |
| `mol_dict` | 59,934,686 | 59,740,405 |
| `independent_dataset_mol_dict` | 2,267,662 | 2,259,203 |

Repaired with `core.autocrlf=false`, `core.eol=lf`, `git rm --cached -r .`,
`git reset --hard <commit>`. After repair **all six hashes match exactly**.

**Without the supplied manifest hashes this corruption would have been silent**
and would have produced a subtly wrong corpus. Any future acquisition of a
binary-bearing repository on this machine must disable `autocrlf` before
checkout. This is a property of the acquisition environment, not of MONN.

### Reproduction

`rebuild_monn_edge_corpus.py --strict-hashes` reproduced every target exactly:

| Quantity | Development | Additional PDB |
|---|---:|---:|
| Raw dictionaries | 12,987 ✓ | 1,853 ✓ |
| Mapped complexes | 12,738 ✓ | 1,851 ✓ |
| Binary edges | 195,798 ✓ | 9,832 ✓ |
| Typed edges | 202,766 ✓ | 9,832 ✓ |
| Missing atom references | 0 ✓ | 0 ✓ |

The deterministic outputs are **byte-identical to the supplied
`VERIFIED_RAW_AUDIT.json` hashes** (`9489540b…`, `cbb0ed98…`) and **bit-identical
across two independent runs** into separate directories.

Affinity firewall: the script's execution path touches only the two interaction
pickles. The two TSVs were hashed as byte streams and **never parsed**.

## R0R-2 — FAIL CLOSED

### Ligand identities: complete

All **10,972** required CCD codes resolved to a molecule, sanitized, canonical
SMILES and Bemis-Murcko scaffold. **Zero failures, zero quarantined complexes.**

Frozen policies: RDKit 2023.09.6; `pickle.load(..., encoding="bytes")` — note
that `latin1` and the default both fail with *"Bad pickle format: bad endian
ID"*, so this is the only working policy; exact-graph identity is
`sha256(MolToSmiles(isomericSmiles=False, canonical=True))`, which collapses
stereoisomers, the conservative direction for leakage; an **empty** Murcko
scaffold (acyclic ligand) generates **no** closure edge, because treating `""` as
a shared scaffold would merge every acyclic ligand into one artificial blob.

Homology: MMseqs2 was used strictly as a **candidate generator** (5,552 candidate
pairs); the authority is parasail Smith-Waterman, BLOSUM62, gap_open 10,
gap_extend 1, identity = matches / local alignment length, coverage = alignment
length / min(len). At ≥40 % identity and ≥80 % coverage, **1,883 of 5,552
candidates were accepted and 3,669 rejected** — the candidate generator alone
would have over-merged by a factor of three.

### The topology, under the frozen rule set

| Quantity | Value |
|---|---:|
| Total complexes | 14,589 |
| Union closure components | **360** |
| Largest component | **13,595 = 93.19 %** |
| Median component size | 1 |
| Components with ≥5 complexes | 42 |
| Components shared by development and additional-PDB | **18** |

Component-level inference is not supported. A paired whole-component bootstrap
over 360 units where one unit holds 93 % of the data has an effective sample size
near one; an equal-weight component macro-average is simultaneously dominated by
singleton components carrying a handful of edges each. And because 18 components
straddle both cohorts, **the additional-PDB set is not an independent
confirmation cohort under this closure** — before publication/time closure is
even considered.

### Why the giant forms — ablation

The cause is not protein homology. It is the ligand-side relations.

| Closure relations | Components | Largest | Fraction | ≥5 |
|---|---:|---:|---:|---:|
| **protein only** (PDB + sequence + UniProt + 40 % homology) | **1,994** | 453 | **3.11 %** | 534 |
| protein + exact ligand graph | 722 | 11,120 | 76.22 % | 161 |
| protein + scaffold | 556 | 13,185 | 90.38 % | 58 |
| **frozen full rule set** | **360** | **13,595** | **93.19 %** | 42 |
| exact sequence only | 2,847 | 371 | 2.54 % | 622 |
| homology only | 13,736 | 22 | 0.15 % | 48 |

Protein closure alone yields an excellent partition. Adding **exact ligand
graph** alone pushes the giant to 76 %; adding scaffold pushes it to 90 %.

The mechanism is chemical, not statistical: promiscuous ligands and cofactors
(ATP/ADP/NAD/heme-like components) recur across unrelated protein families, so
"same exact ligand graph" transitively bridges the protein universe. Scaffold
sharing compounds it, since common ring systems are shared by thousands of
otherwise unrelated ligands.

**This is a real property of the estimand, not an implementation artifact.** A
closure that simultaneously closes protein identity and ligand identity over a
promiscuous-cofactor corpus does not partition.

### One number that did reproduce

Development **exact sequences = 2,404**, matching the consolidated report's
claimed 2,404. That claim graduates from external claim to **reproduced**. It is
the only derived quantity from that document reproduced here; `4,067/701`,
`8,646`, `524/157` and every AP value remain unverified external claims and were
not used.

*(Neutral note: the protein-only projection below happens to yield 325
components, numerically equal to a "325 satellite components" figure in the
separate, unreproduced U0-U3 summary. Different corpus, different construction —
a coincidence with no evidential weight.)*

## Corpus characteristics (recorded for any future registration)

All seven PLIP channels clear a 5 % development-prevalence bar, so all are
evaluable: Hydrogen Bonds 94.0 %, Hydrophobic 86.5 %, Water Bridges 53.5 %, Salt
Bridges 33.7 %, pi-Stacking 31.0 %, pi-Cation 10.4 %, Halogen Bonds 7.2 %.

The development complete residue×atom matrix is **403,454,851 cells** with a
positive rate of **4.85e-4** (about 1 in 2,060). Complete-matrix evaluation is
computationally feasible but the extreme imbalance is the dominant property any
AP estimator must respect.

## The adjudication this fail-closed requires

The frozen rule set does not partition. There is exactly one decision that could
unblock R0R-2, and it is a **scientific decision about the estimand that I am not
authorized to take**:

> Define the inference partition on the **protein side alone** (exact PDB, exact
> sequence, UniProt, 40 % homology), and control ligand leakage by *evaluation
> design* — held-out ligand graphs and the `BX` wrong-ligand control — rather
> than by closure.

Projected consequences, measured (`PARTITION_FEASIBILITY_PROJECTION.json`):

| | Development | Confirmation candidate |
|---|---:|---:|
| Complexes | 12,738 | 710 |
| Components | 1,669 | 325 |
| Largest component fraction | 3.46 % | 1.97 % |
| Components ≥5 complexes | 458 | **37** |
| Exact sequences | 2,404 | 348 |
| Exact ligand graphs | 9,697 | 346 |
| Positive binary edges | 195,798 | 3,383 |

The confirmation candidate passes **6 of 7** capacity checks. It fails
`components_ge5_complexes` (37 against a required 60).

So even the alternative partition is **not a free pass**. It would require either
accepting fewer well-populated components, or a differently defined evaluable
unit, or more confirmation data — each of which is a registered decision, not a
default.

**What this changes about the claim.** Under a protein-only partition, the
experiment would generalize over *proteins*, not over *protein–ligand pairs*.
That is a weaker and different scientific statement, and it must be written into
the estimand rather than discovered afterwards.

## What must not be concluded

This is outcome **2 — closure cannot support inference**. It is emphatically not
outcome 5. Nothing here says sequence-plus-2D inputs lack the required
information. No model was trained, no AP was computed, no `B4` exists, and no
`B5` Gate was attempted. The PLM hypothesis is untested, and the runtime probe
from the previous session showed the compute path is largely available.

## Frozen surfaces

`theory/`, `model/`, `contracts/`, `scripts/`, `weights/`, `config/` unmodified.
`K(B(z)F(z))`, CSMO, Band, simplex, positive ridge, fixed mesh and production `z`
untouched. No affinity, DAVIS, recipient or few-shot support label read. Nothing
committed, nothing pushed. The pre-existing dirty worktree was preserved; no
unrelated change was reverted or deleted.
