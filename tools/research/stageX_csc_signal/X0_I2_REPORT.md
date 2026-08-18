# X0 I2 corrected report — representation-capability instrument

Generated 2026-08-18 (round 2). Preregistration
`STAGE_X0_PREREGISTRATION.md` SHA-256
`03cdc907df3e778f5fe79fb1a238d35ebb6ece5e9e743db181728ba6b25e9683` (verified
frozen). Implementation: `x0_pair_table.py` + `x0_i2.py`. Round-1 artifacts
`X0_INSTRUMENTS.json`, `x0_instruments.py` are preserved as negative evidence
and are NOT the basis of this report.

## 1. What was reproduced and fixed (round-1 defects)

| # | Defect (round 1) | Evidence | Fix |
|---|---|---|---|
| 1 | WT local window used the sequence midpoint; mutant used the mutation coordinate (different positions compared) | round-1 `local_window` ratio exactly 1.0 (= sqrt(2)/sqrt(2): the "pair" distance equalled the inter-protein scale, i.e. midpoint window vs mutation-site window of unrelated positions) | windows extracted at the SAME verified canonical coordinate for WT and mutant; window contract recorded per pair |
| 2 | same coordinate defect in ESM local windows (ratio 1.32) | `X0_INSTRUMENTS.json` | same fix |
| 3 | mutation_token parent-parent denominator = 0 -> ratio 1.41e12 | `X0_INSTRUMENTS.json` mutation_token ratio 1414213538169.86 | mutation_token reported as pair-conditioned edit descriptor with denominator None; excluded from the representation gate |
| 4 | mutation_token counted as an admissible representation | round-1 gate counting | excluded from gate (frozen rule 3 requires admissible representations; an edit descriptor is not an independent protein representation) |
| 5 | >= 4 point mutations did not match the downloaded reference residues: BRAF(V599E) (historical numbering) and the three PDGFRalpha variants (wrong-species accession Q9DE49 = Danio rerio pdgfra) | direct residue check | BRAF: cited historical renumbering V599E -> canonical V600E (Davies 2002 Nature 417:949; UniProt P15056 VAR_018629; P15056 residue 599 = T, 600 = V). PDGFRalpha: S1 GenBank NP_006197 + KLIFS resolve human PDGFRA = P16234 (fetched); D842V/T674I/V561D all pass old-residue checks on P16234 |
| 6 | `hash(name)` seeds random control (process-dependent) | round-1 code | SHA-256-derived seeds only; cross-process test included |
| 7 | KLIFS representation/coverage absent | round-1 code | KLIFS kinase_ID + 85-position pocket fetched per parent; pocket mapped onto the canonical sequence with gapped Needleman-Wunsch (identity >= 0.80); mutation-position coverage census with per-row reasons |
| 8 | no admission rules for multi-mutation/deletion/insertion/truncation/fusion/unknown notation/isoform/long sequences | round-1 code | explicit admission policy in `X0_PAIR_TABLE.json`; every excluded row keeps its reason and stays in the census |
| 9 | distance ratio certifies sensitivity only, not biological direction | — | scope note recorded in `X0_I2.json`; direction is tested by I1 (planted recovery) and I3 (ID-equivalence), not I2 |

## 2. Pair census

97 assay rows (21 WT parents + 76 mutant constructs). 65 admitted
single-point pairs (old residue verified at canonical coordinate, construct
coordinate resolvable and residue-verified). Exclusions:

| Status | Count | Reason |
|---|---|---|
| excluded_multi_point | 3 | KIT(V559D/T670I), KIT(V559D/V654A), EGFR(L858R/T790M) |
| excluded_deletion | 2 | EGFR(d746-750), EGFR(d752-759) |
| excluded_deletion_plus_point | 3 | EGFR(d746-750/T790M), EGFR(d747-749/A750P), EGFR(d747-752/P753S) |
| excluded_insertion | 1 | FLT3(ITD) |
| excluded_construct_unresolved | 2 | PDGFRalpha(D842V) (S1 clone "Cytoplasmic (668-1210)" exceeds canonical length 1089; not silently repaired), FGFR3(K650E) (clone field "cytoplasmic", no numeric range) |

Notation quirk handled with recorded basis: S1 Mutation column "T6741I" ->
T674I (construct name column authoritative). BRAF renumbering recorded with
evidence URLs. Fusion class: no fusion constructs present (count 0, policy
defined). Isoform differences: recorded per row (UniProt SV field); all rows
map to canonical human sequences.

## 3. Corrected metrics (65 admitted pairs; per-representation pair subsets as listed)

r_pair = ||x(WT) - x(mutant)|| at the SAME verified mutation coordinate;
denominator = median inter-parent distance of the same representation type
(local denominators use each parent's windows at its own verified mutation
centers; same-parent pairs excluded). Bootstrap: 2000 parent-cluster
resamples, seed 20260820, numerator and denominator recomputed per draw.

| Representation | Type | n pairs | median pair dist | denominator | ratio [95% CI] | pass (>= 0.05) |
|---|---|---|---|---|---|---|
| pair_centered_local_window (one-hot, radius 6) | local | 65 | 1.414 | 4.899 | 0.2887 [0.2774, 0.2887] | PASS |
| esm_local_window (ESM-2 150M, radius 6) | local | 49 | 0.508 | 4.229 | 0.1201 [0.1094, 0.1335] | PASS |
| klifs_pocket (85-position aligned one-hot) | local | 27 | 1.414 | 9.592 | 0.1474 [0.1393, 0.1562] | PASS |
| global_esm (mean-pooled ESM-2 150M) | global | 49 | 0.0229 | 1.928 | 0.0119 [0.0105, 0.0159] | FAIL |
| composition (amino-acid composition) | global | 65 | 0.00136 | 0.0675 | 0.0201 [0.0169, 0.0241] | FAIL |
| random (SHA-256-seeded Gaussian, dim 128) | control | 65 | 7.929 | 7.855 | 1.0095 [0.9801, 1.0275] | (sensitivity control only) |
| mutation_token (edit descriptor) | edit | 65 | 1.414 | n/a (no legal inter-protein scale) | n/a | excluded from gate |

ESM-excluded pairs (16): mutation position beyond the 1020-token ESM window —
ALK(4: 1156/1174/1196/1275), C-MET(8: 1200-1250), LRRK2(3: 1441/2019/2020),
TIE2(1: 1108). KLIFS-excluded pairs (38): 35 with mutation position not
represented in the 85 aligned pocket positions; 3 with parent pocket
unaligned (LRRK2, identity < 0.80). KLIFS parent census: 21/21 parents have
KLIFS entries; 20/21 pockets align to the canonical sequence with identity
>= 0.90 (KIT 0.976, MAPK14 0.988, PDGFRA 0.906 with 8 mismatches recorded);
LRRK2 unaligned (recorded with reason). Gatekeeper-class mutations map to
KLIFS pocket index 45 (EGFR T790M, C-KIT T670I, RET V804L/M, FGFR1 V561M,
FGFR4 V550E/L, C-SRC T341M, PDGFRalpha T674I, P38alpha T106M, ALK L1196M),
independently validating the pocket mapping.

## 4. I2 gate (frozen)

Rule: at least three admissible representations pass (ratio >= 0.05),
including at least one local representation. Result: 3 passing
(pair_centered_local_window, esm_local_window, klifs_pocket), all three
local. **I2 gate: PASS.**

Global mean-pooled ESM remains insensitive to single-residue edits
(ratio 0.0119, as in round 1) — this is the honest, corrected conclusion:
global pooled PLM representations are insensitive at the I2 threshold, while
mutation-site-local representations express the edit. The one-hot local
window passes by construction (a single substitution changes exactly two
one-hot positions); its CI is degenerate by construction and is reported as
such. Random control passes the ratio metric, showing the metric alone is
sensitive — a pass there is explicitly NOT a biological capability (it is
excluded from the gate and from any biological interpretation).

## 5. Scope and limitations

- I2 certifies representation sensitivity to verified single-residue edits
  relative to an inter-protein scale of the same representation type. It does
  not certify biological direction; direction/transferability are tested by
  I1 and I3.
- WT parent rows have no construct annotation in Duong-Ly S1; WT references
  use the canonical full-length sequences, and construct offsets are verified
  on the mutant rows. Construct-level WT sequence differences remain an
  unresolved assay-level confound, recorded in the pair table.
- KLIFS pockets come from the KLIFS API representative structures; a pocket
  carrying a structure mutation would fail the old-residue check at that
  position and the pair is excluded with that reason (none observed among the
  27 included pairs).
