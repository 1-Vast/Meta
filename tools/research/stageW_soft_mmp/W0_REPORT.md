# Stage W W0 — soft controlled chemical-change surfaces: Davis FAIL, KIBA PASS

Preregistration SHA-256:
`ae96762e319521f30aa09eb1a79fb8bb0e3ea324b21d4b40868aa6826a45dc71`.
Artifact: `W0_SOFT_MMP_CENSUS.json`. No W1 model was designed or trained.

## Redirect

Stage U/V closed the exact-MMP route on BindingDB-Ki double-cold. Stage W
opens the redirected goal on independent open datasets, each separately: first
measure whether a soft but controlled chemical-change family supplies a
protein-component-cold surface, then — only if it passes — build the local
protein × ligand representation.

## Frozen soft family

Single-cut MMP fragments (RDKit Hussain-Rea). Family key:
`sha256(murcko_core | attachment_element | attachment_aromatic |
attachment_in_ring | category(R_a) >> category(R_b))`, with `category` =
quantized heavy-atom / aromatic / ring / HBD / HBA / charge class. This is a
soft family, **not** an exact transformation; the same-core residual is
reported as the chemical-control price.

## W0 results

| dataset | Davis | KIBA |
|---|---:|---:|
| source rows | 25,772 | 117,657 |
| unique targets | 379 | 229 |
| unique ligands | 68 | 2,068 |
| CD-HIT40 components | 206 | 127 |
| same-target MMP observations | 2,653 | 349,100 |
| exact keys / soft families | 7 / 7 | 10,565 / 2,428 |
| rich families (>=3 targets & >=3 components) | **7** | **2,420** |
| cross-component D rows | 498,967 | 5,504,678 |
| top-1 / top-10 family share | 0.143 / 1.000 | 0.031 / 0.137 |
| top-1 / top-5 target share | 0.003 / 0.013 | 0.021 / 0.091 |
| within-family same-target across-core residual median | **not identifiable** | **0.667** (p95 2.951) |
| between-component MS vs permuted-component null | observed > null, p=0 | observed 0.511 vs null 0.304, p=0 |
| **frozen W0 gate** | **FAIL** | **PASS** |

**Davis is insufficient for this surface**: its 68 ligands form only 7 exact
MMP families, below the frozen threshold of 20 rich families, and the
same-core residual is unidentifiable. Davis is recorded as closed for this
soft-MMP surface; no Davis W1 model may be trained.

**KIBA passes all ten frozen checks.** Its 2,428 soft families cover 127
CD-HIT40 components with 5.5M cross-component D rows, concentration is healthy,
and the within-family same-target across-core residual median is 0.667 KIBA
units (below the frozen 1.0 cap). The between-component variation is outside
the component-permuted null (p=0), which is a necessary but not sufficient
condition for protein information.

## Decision

Proceed to W1 on **KIBA only**, under a new W1 preregistration frozen before
any training metric is read. Davis stays closed for this surface. The W1
representation contract is fixed in outline by the redirecting instruction
(multi-pocket residue states × ligand pharmacophore subgraphs, local
cross-attention, independent level/shape heads, all mandatory protein
counterfactuals). Exact W1 hyperparameters and gates are frozen before W1
training.
