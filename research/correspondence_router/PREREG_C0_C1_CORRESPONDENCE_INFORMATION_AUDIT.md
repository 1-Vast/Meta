# Preregistration — C0/C1

## Untouched correspondence corpus and audit-only exact-coupling information test

Stage identifier: `CORR-C0_C1_EXACT_CORRESPONDENCE_INFORMATION_AUDIT`

Written 2026-08-10, after `PHASE2B_S5D_GATE.json`
(`LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED`, commit `b6d265d`) and before any C0
or C1 code, census, statistic or Gate value exists.

This document registers **audit-only** work. It trains nothing and authorizes no
model. Stage C2 requires its own separate preregistration and may only be
written if C1 passes.

## 1. Question

> Does frozen P1B atom-by-128-slot geometry contain enough information to
> recover ligand-specific **exact** atom-residue correspondence?

C1 answers the prior half of that: is there any exact-coupling structure to
recover at all, beyond additive marginals and a degree-preserving null, on a
corpus no MetaSieve stage has ever touched?

## 2. What P1B actually is — and one correction carried into the design

`scripts/build_structure_supervision.py` defines the slot target as

```text
slot(r)      = min(127, sequence_index(r) * 128 // sequence_length)
distance(i,s)= min over ALL atoms of ALL residues assigned to slot s
contact(i,s) = distance(i,s) <= 6.0 A
```

So `contact_prob(i,s)` is a Bernoulli-like **"any residue in slot s contacts
atom i"** prediction. It is **not** additive contact mass. This stage therefore

- never imposes `sum_r p(i,r) = contact_prob(i,s)`,
- never calls any such equality a conservation law,
- uses P1B only as a frozen gate, offset or admissibility prior,
- and requires every exact-residue output to permit multiple residues in one
  slot to contact the same ligand atom.

A second consequence is structural and is recorded here because it constrains
what C1 can conclude: **P1B is constant across residues within a slot.** Its
contact probability and distance bin cannot, by themselves, discriminate
residues inside a slot. In the C2 estimand its within-slot contribution is
confined to gating which bilinear channels `b` are active. C1 must therefore
measure whether within-slot structure exists at all, and must not be read as a
test of whether P1B can supply it.

## 3. Exposure policy — frozen

Every PDB or system ID consumed by any prior stage is **EXPOSED** and is
excluded from this corpus, including IDs only used for control construction,
acquisition or QC. The registry unions:

```text
pilot20k governed / ungoverned / homology-split corpora  (P1B)
preflight and preflight_qc100 corpora                    (QC)
MONN development + additional-PDB edge corpora           (S7, B5, S3R, S4R, S5D)
the ssl_b2 independent structural set                    (S-programme)
```

Heldout-A is **permanently consumed**. It is not opened, not re-scored, and no
metric, threshold, architecture decision or Gate in C0, C1 or any later stage
may reference it.

## 4. Corpus construction — frozen

Source: raw RCSB mmCIF coordinate files already present under
`dataset/raw/open_structures/pilot20k/mmcif/`, restricted to entries **not** in
the exposure registry. These files are raw coordinates; being on disk is not
exposure, and the audit records that they were never parsed, scored or
governed by any prior stage.

- **BioLiP2** is used **only** as a biological-relevance filter to choose which
  non-polymer entity is a genuine ligand rather than a crystallisation
  additive. No BioLiP2 numeric field is read. Licence role `ANNOTATION_ONLY`,
  not redistributed.
- **PLINDER is not used.** `research/ssl_b2_structural_observability/LICENSE_AND_PROVENANCE_AUDIT.md`
  already records that its annotations carry a separate licence and that
  programme policy forbids pulling them. That standing decision is honoured
  rather than reopened.
- **No affinity field is opened** from BioLiP2, PLINDER, PDBbind, ChEMBL,
  BindingDB, DAVIS, KIBA or any recipient dataset. Affinity reads must be `0`.

### 4.1 Selection rules

A candidate system is one (PDB entry, protein `label_asym_id`, ligand
`label_asym_id`) triple satisfying all of:

```text
S1  the PDB id is absent from the exposure registry
S2  the entry has an X-ray _refine.ls_d_res_high <= 2.5 A
S3  the ligand comp id appears in BioLiP2 for this PDB id and is not in the
    frozen additive blacklist below
S4  the ligand has >= 6 and <= 80 heavy atoms after altloc canonicalisation
S5  the protein entity canonical sequence length L satisfies 150 <= L <= 1200
S6  the ligand has >= 1 heavy atom within 6.0 A of a protein heavy atom
S7  no _struct_conn covalent link joins the ligand to the protein entity
S8  the ligand CCD parses under RDKit and yields a Murcko scaffold
```

Additive blacklist (frozen): `HOH GOL EDO PEG PG4 PGE SO4 PO4 CL NA K MG CA ZN
MN FE CU NI CD HG ACT ACY DMS MPD TRS EPE IMD BME NAG BMA MAN FUC GAL GLC XYL
IOD BR NO3 CO3 FMT OXL TLA CIT MES CAC SCN AZI UNL UNX`.

### 4.2 Mapping rules

```text
M1  altloc canonicalisation: model 1 only; per (asym, seq, atom) prefer blank
    altloc, then highest occupancy, then lowest altloc id
M2  hydrogens are dropped
M3  the protein canonical sequence is _entity_poly.pdbx_seq_one_letter_code_can
    for the entity of the chosen protein asym
M4  a resolved residue's sequence index is int(label_seq_id) - 1, which is the
    mmCIF definition of an index into _entity_poly_seq
M5  slot(r) = min(127, sequence_index(r) * 128 // L), byte-identical to the
    P1B slot policy
M6  residue coordinates are all heavy atoms of that label_seq_id
M7  exact distance d(i,r) = min over heavy atoms of residue r of ||x_i - x_r||
M8  exact contact(i,r) = d(i,r) <= 6.0 A, the frozen P1B contact threshold
```

`M4` replaces P1B's parasail alignment. It is an **exact equivalence** whenever
the record sequence equals the entity canonical sequence, which is how the P1B
corpus was built. The audit must verify `M4` reproduces the parasail mapping on
a sample of already-exposed structures and must fail closed otherwise. This
deviation is declared here, not discovered later.

`M8` reuses P1B's own `CONTACT_THRESHOLD_ANGSTROM = 6.0`. Using any other
threshold would break the correspondence between the slot gate and the exact
edges, which is the entire object of study. A `4.5 A` recomputation is reported
as a **non-gating** sensitivity only.

### 4.3 Exclusion and replacement rules

Exclusions are recorded with a reason and are never silently imputed:
unparsable mmCIF, missing `_entity_poly` sequence, protein mapping coverage
below `0.90`, zero resolved residues, no BioLiP2 ligand, RDKit failure,
covalent link, or dimension violations. **There is no replacement rule**: an
excluded system is dropped, never substituted, and the exclusion census is
published.

## 5. Deconvolution-unit census — frozen definitions

```text
candidate slot     a slot containing >= 2 mapped resolved residues
unit               an (atom i, candidate slot s) pair
positive unit      a unit where >= 1 residue in s has exact contact with i
multi-contact unit a unit where >= 2 residues in s have exact contact with i
multi-contact rate multi-contact units / positive units
valid 2x2          atoms i != j and residues r != r' in one candidate slot with
                   contact(i,r) = contact(j,r') = 1 and
                   contact(i,r') = contact(j,r) = 0, all four read from raw
                   coordinates
replicate system   two systems sharing an exact protein sequence and an exact
                   ligand CCD id, from different PDB entries
```

Ambiguous crossed pairs are excluded by this frozen rule and are never treated
as negatives.

## 6. Closure and splits — frozen

Similarity edges between systems:

```text
E1 protein   global sequence identity >= 0.40, parasail nw_trace_striped_16,
             blosum62, gap open 10, extend 1, identity = matches / min(len)
E2 ligand    identical Murcko scaffold SMILES, or identical CCD id
E3 construct identical exact protein sequence
E4 series    same protein cluster under E1 and same scaffold under E2
```

`E1` is computed with a 3-mer containment prefilter. The prefilter threshold is
frozen at `0.03` and its **recall must be validated by brute-force alignment on
a random 200-sequence subsample**; a prefilter that misses any true `>= 0.40`
pair fails closed.

Reported but **not** part of the primary closure, because they are label-derived
and using them to define the split would be circular: pocket composition
similarity, PLI/contact-pattern Jaccard, deposition date and citation year.

### 6.1 Primary closure and the registered fallback

```text
PRIMARY   connected components of the union of E1, E2, E3
FALLBACK  a DataSAIL-style two-dimensional block partition: cluster proteins by
          E1 and ligands by E2, assign the (protein cluster, ligand cluster)
          grid to train / development / sealed blocks, and DISCARD every system
          whose protein cluster and ligand cluster are assigned to different
          blocks
```

The union closure is reported first together with its giant-component
diagnostic. The fallback is used only if the union fails `C0.7`, and both are
always published so the giant-component assumption is tested rather than
assumed.

## 7. C0 admissibility Gate — frozen

```text
G0a  independent inference components >= 60
G0b  largest component fraction <= 0.25
G0c  minimum detectable effect at 80% power <= 0.05 absolute AP
```

`G0c` power calculation, frozen: with `N` inference components and per-component
statistic dispersion `sigma`, the paired component-bootstrap minimum detectable
effect at 80% power and a one-sided 95% bound is
`MDE = (1.645 + 0.842) * sigma / sqrt(N)`. `sigma` is estimated **from the
degree-preserving null arm only**, never from the empirical-minus-null contrast
that C1 tests. If `MDE > 0.05` the panel is underpowered.

If `G0a`, `G0b` and `G0c` cannot all be satisfied by the primary closure or the
registered fallback, terminate with
`CORRESPONDENCE_DATA_OR_CLOSURE_NOT_IDENTIFIABLE`.

## 8. C1 — audit-only information test

Inference unit is the complex, aggregated complex -> closure component. Rows of
the atom-residue matrix are **never** bootstrapped as IID.

### 8.1 Statistic

Within a positive unit `(i, s)` with `n_s >= 2` candidates, rank the candidate
residues by the **residue additive marginal** `m_r` = the number of ligand atoms
of this complex contacting `r`. The unit's score is the exact tied-AP of that
ranking against the true contact labels of the candidates. `AP_within` is the
mean over positive units, then averaged complex -> component.

The atom marginal is constant inside a unit and therefore cannot affect
within-unit ranking; it is used only for the complete-edge statistic.

### 8.2 Nulls and shuffles

```text
N1 fixed-degree rewire   Curveball / checkerboard swaps on the atom x residue
                         bipartite contact matrix, >= 100 swaps per edge,
                         degrees verified exactly preserved and chain movement
                         verified non-zero
N2 additive marginals    the rank-1 marginal-product model
N3 ligand shuffle        contacts from a different ligand, same protein
N4 protein shuffle       contacts from a different protein, same ligand
N5 atom shuffle          ligand atom identities permuted within the complex
N6 geometry shuffle      slot assignment permuted within the complex
```

`N1` is an **evaluation null**, not a biological non-binder, and is described as
such in every artifact.

### 8.3 C1 Gates — frozen

```text
C1a  AP_within(empirical) - AP_within(fixed-degree rewire) >= +0.05,
     one-sided 95% component-bootstrap lower bound > 0
C1b  positive deconvolution units >= 10,000 and valid 2x2 checkerboards >= 1,000
C1c  replicate agreement is defined and its ceiling exceeds the rewire null,
     so replicate-normalized headroom has a positive, stable denominator
```

Bootstrap: paired, closure components as units, 10,000 resamples, seed
`20260903`, one-sided 95% lower bound. The `+0.05` margin is the project's
standing practical-effect convention, reused unchanged, and is admissible only
if `G0c` shows the panel can detect it.

Also reported, non-gating: complete-edge AP, valid-checkerboard prevalence,
multi-contact rate, the `N2`-`N6` arms, the `4.5 A` sensitivity, and
replicate-normalized recoverable headroom
`(AP_emp - AP_rewire) / (AP_replicate_ceiling - AP_rewire)`, which may be
reported **only** when its denominator is positive and stable and which can
never replace `C1a`.

## 9. Terminal verdicts

Exactly one, by earliest failed boundary:

```text
CORRESPONDENCE_DATA_OR_CLOSURE_NOT_IDENTIFIABLE   C0 fails G0a, G0b or G0c
SLOT_ROUTING_ESTIMAND_INVALID                     the slot/mapping contract fails
                                                  its own equivalence or
                                                  determinism checks
EXACT_EDGE_COUPLING_NOT_SUPPORTED_BY_TEACHER      C1a or C1b or C1c fails
CORRESPONDENCE_INFORMATION_PRESENT_C2_AUTHORIZED  all C0 and C1 Gates pass
```

The remaining verdicts in the programme's terminal set —
`P1B_COARSE_GEOMETRY_INFORMATION_INSUFFICIENT`,
`GEOMETRY_ROUTED_CORRESPONDENCE_NOT_IDENTIFIED` and
`GEOMETRY_ROUTED_CORRESPONDENCE_IDENTIFIED_IN_DEVELOPMENT` — belong to C2 and
are unreachable from this document.

## 10. Stopping rules

One run. No selection rule, mapping rule, threshold, seed, margin, closure
definition or exclusion rule may change after any C0 census or C1 statistic is
read. A failed Gate may not be rescued by widening the corpus, relaxing the
contact threshold, changing the closure or re-running the prefilter.

`CORRESPONDENCE_INFORMATION_PRESENT_C2_AUTHORIZED` authorizes exactly one thing:
writing the C2 preregistration for the single Geometry-Gated Coarse-to-Exact
Correspondence Router described in the parent instruction. It does not authorize
training by itself, does not open affinity, energy, selectivity, few-shot
adaptation or `z`, and does not permit any additional module.

## 11. Boundary

`model/`, `scripts/`, `theory/`, CSMO, Band, the mesh, `z` and the frozen
operator

```text
A(F, z) = K(B(z) F(z))
```

are not modified. Prior-stage code may be imported for equivalence checking but
never edited. This stage identifies at most the presence or absence of
exact-coupling information in an untouched structural corpus. It identifies no
energy, no affinity, no selectivity, no few-shot section, no biological `z` and
no validated end-to-end DTA model.
