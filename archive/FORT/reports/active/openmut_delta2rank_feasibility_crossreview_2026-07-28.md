# OpenMut / Delta2Rank feasibility and adversarial cross-review

**Date:** 2026-07-28

**Decision:** `OPENMUT_DATA_RECOVERY_CONDITIONALLY_FEASIBLE__DELTA2RANK_NO_TRAIN`

**Current result category:** 3 - the present data do not identify amino-acid substitution
geometry as a transferable target-ligand mechanism. A new source-resolved substrate or prospective
measurements are still required.

## Scope and success criteria

Three user-supplied reports were audited:

| report | SHA-256 |
| --- | --- |
| diagnostic report | `e53b2d063b4407c35b754656ab9b85d90a7d93064c89b169c0bf823ec7d0ca4a` |
| Delta2Rank proposal | `b456cf86c25ce1a6df078402fb1f28bdd9dfd6a65e6931f3f247e1e9581412e8` |
| public-data reconstruction | `2f8becab8a7015fdc90132513a039ed3e01e7e03eca1a040dd58bf38fb175418` |

The audit asked whether each proposal adds identifiable information, rather than whether it is
computationally possible. Hardware is not a rejection criterion. Large downloads and expensive
models remain allowed when a preceding information gate justifies them.

No new affinity value was read, no affinity-bearing file was downloaded, and no predictor was
trained. The DAVIS-Complete affinity table, local Davis target-conditioned confirmation gate,
development labels, and sealed labels remain untouched in this round.
Specifically, `davis_complete.tab` was neither downloaded nor read.

The required scientific claim ladder remains:

1. mutation-conditioned ligand reordering exists above measurement noise;
2. a true substitution coordinate beats family, position, composition, and random geometry;
3. the coordinate predicts reordering across independent proteins, families, and provenance;
4. the mechanism improves strict family-cold plus scaffold-cold DTA.

A pass at one level does not establish the next.

## Diagnosis accepted, with one qualification

The diagnostic report is consistent with the accepted LOCK/CLOCK result. KirHub identifies aligned
position and correct-target signal, but fixed BLOSUM semantics are not identified: fixed LOCK was
worse than aligned identity and did not beat either BLOSUM-label permutation or a
correlation-matched random PSD geometry. Pooled ESM-2 is also not a reliable interaction coordinate.

This is not a biological null. It is a coordinate-identification failure on a kinase-only,
ligand-warm, single-source graph. It closes attempts to rescue the same fixed coordinate by changing
rank, kernel width, interaction head, seed count, or training duration. It does not close a new
dataset containing explicit WT-to-mutant operations.

## Historical overlap and genuine novelty boundary

Most of the proposed Delta2Rank mathematics already exists in the project:

| proposed element | prior evidence | decision |
| --- | --- | --- |
| four-cell cross-difference | REWIRE, MISO, DICE/AXIS, KirHub DD, WTPAIR | valid estimand, not new |
| antisymmetric low-rank bilinear operator | DICE/AXIS and WTPAIR | not new by itself |
| additive-nuisance removal | MISO-OR and R1A | conditional tool, not a new module |
| direct centered operator with regular null | R-MAON G0 synthetic calibration | valid carrier, no real substrate |
| pairwise/listwise relative affinity | historical ranking routes and PBCNet audit | auxiliary objective only |
| Hodge/cycle analysis | existing comparison-graph reasoning | quality control, not a predictor |

The relevant historical results are not weak suggestions:

- REWIRE found 471 independent transformations and an independent-group MDE of `0.0645`, then
  stopped before the proposed interaction block.
- R1A measured raw double-difference reliability near `rho=0.336`; high-confidence filtering left
  only 502-513 independent transformations, dominated by one source.
- MISO-OR reduced compound/scaffold nuisance but its selectivity-proxy gain was `-0.001881`
  with interval `[-0.004785,+0.001019]`.
- DICE/AXIS used a bias-free odd low-rank operator; the standard interface screen failed against
  chemistry and global controls over 518 transformation groups.
- WTPAIR directly supervised protein-pair by ligand-pair mixed differences. Its component-macro
  Spearman was `0.0323`, below ligand-only `0.0429`, KLIFS-group centroid `0.0911`, and a matched
  cellwise bilinear model `0.0419`; it did not beat random protein.
- The prior local-edit audit found that an explicit MMP/local-edit representation was worse than a
  whole-molecule difference in the chemistry-cold task, while every arm failed the target-disjoint
  task. An MMP coordinate is therefore a falsifiable arm, not an assumed improvement.

The genuinely untested object proposed here is narrower:

1. an explicitly canonicalized WT-to-single-substitution operation, rather than a difference between
   two pooled target embeddings; and
2. its joint action with a recurring directed ligand change on independent continuous-affinity
   rectangles.

Neither a directed ligand edit nor a bilinear operator is new by itself. The joint substrate has not
yet been tested on compliant continuous-affinity data. Its existence is a data hypothesis, not a model
result or a novelty claim.

## Mathematical correction to the quartet proposal

For one mutation `m` and ligand `l`, define

```text
delta_m(l) = y(mutant, l) - y(WT, l)
```

and remove the mutation-wide shift using train-only weights:

```text
z_m(l) = delta_m(l) - weighted_mean_l(delta_m(l)).
```

The proposed quartet target is then

```text
C_m(a,b) = z_m(a) - z_m(b).
```

Consequently, a complete `L`-ligand panel has at most `L-1` independent reordering degrees of
freedom, not `L choose 2`. Expanding 72 ligands into 2,556 pairs creates correlated training rows; it
does not create 2,556 observations or independent units.

The minimum unlocked model should therefore estimate a centered mutation-response field in frozen,
low-dimensional mutation and ligand bases:

```text
s(m,l) = s0(l) + [u(m) - mean_train(u)]^T Theta [phi(l) - mean_train(phi)]
C_hat(m;a,b) = s(m,a) - s(m,b)
```

`s0` is the same train-only shared ligand/common-mutation baseline in the full and null models.
`Theta` must be optimized directly with a frozen ridge penalty and frozen basis dimensions. It must
not be represented as `UV^T`, `gamma*a*c^T`, or another factorization whose directions exist only under
the alternative. The exact nested null `Theta=0` is then unique, regular, and exactly recovers `s0`.
Full and null must otherwise have identical fitting, weights, and data. If a rank-1 arm is tested, its
two directions must be frozen without affinity labels and only one signed scalar coefficient may be
fitted. This form guarantees ligand-swap antisymmetry and cycle consistency; a separate cycle loss is
identically redundant.

On same-assay complete rectangles, the cross-difference already annihilates additive target,
ligand, mutation-intercept, and assay-offset terms. Functional-ANOVA purification or cross-fitted
nuisance models are useful only for missing or cross-context data, and must be shown not to remove
the interaction. R1A already found a Robinson/R-learner residualization numerically identical to the
raw double difference. These methods cannot manufacture information absent from the comparison graph.

Hodge decomposition has the same restricted role. Pair flows generated from scalar affinities are
gradients by construction. Nonzero curl after stitching records is evidence of incompatible
documents, assay conditions, censoring, or extraction errors; it is not automatically a biological
mechanism.

Continuous centered `delta` is the primary endpoint. Reversal classification, listwise ranking,
uncertainty heads, mixture-of-mechanisms, conformer ensembles, PLM deltas, local structures, docking,
and simulation teachers are all downstream. Outcome-selected "positive rectangles" are prohibited
unless the inclusion probabilities and weighting rule are frozen before fitting.

## Label-free public-source audit

### DAVIS-Complete

The repository `ZhiGroup/DAVIS-complete` is frozen at commit
`799ac5696e7afb9a23d2767cace2352e243b353e`. GitHub exposes no repository license, so its code is not
assumed reusable. The data are separately available from Harvard Dataverse
`doi:10.7910/DVN/RTQGP1`, version 3.0, released 2025-11-27 under CC0 1.0:

| file | bytes | MD5 |
| --- | ---: | --- |
| `davis_complete_full.fasta` | 350,542 | `48dd9b6739c90505ea4e1b3468e12815` |
| `davis_complete_kinase_domain.fasta` | 145,594 | `8d49c499bf5cd2b9272dfff512cde24e` |
| `davis_complete.tab` | 13,477,288 | `153f096a650351a135fdfc55e5875b24` |
| `dict_name_cid_smi.pkl` | 6,586 | `b4df18030939f871164af23de13060fe` |

The two FASTA files were streamed without writing them to disk, and both MD5 values were reproduced
in the `drug` environment. The full-sequence FASTA contains 444 headers. A conservative label-free
parser finds 39 named records that differ by exactly one residue from another named record. Two use
`kit_v559d` as an already-mutated reference; restricting to explicit canonical-WT names leaves
37 clean single-substitution constructs across 10 named WT base proteins. The publication describes
modified entries from 11 proteins overall, so canonicalization rather than filename counting is a
required D1 deliverable.

Even before censoring or missingness is inspected, the independent base-protein count is therefore
about 10-11. At paired SD `0.10`, its optimistic MDE80 is `0.089-0.085`, not `0.03`. The dense
72-ligand panel is useful for within-source mechanism diagnostics, but it cannot independently
establish multi-family or provenance-independent transfer.

Zenodo record `15391611` is open under CC BY 4.0; its main processed archive is about 39.3 GB. It is
not needed for the source/topology gates. Deferral is based on scientific stage, not hardware.

### Davis confirmation firewall

The local asset `dataset/public/chembl_37/processed/panel_davis/manifest.json`, SHA-256
`84035ad45e0a47a8708520acb9caec638e53533a246cf7b8af68465f85199e69`, registers Davis 2011 document
`CHEMBL1908390` as a single-use confirmation panel with `consumed=false`. An arm-blind historical
power audit was permitted, but no target-conditioned confirmation gate has consumed it.

DAVIS-Complete includes the original Davis WT panel plus modified entries. It can never be called an
independent confirmation source for this project. Before any DAVIS-Complete affinity access, F0 must
freeze cell-level overlap keys and a one-way policy:

- preserving `panel_davis` as confirmation excludes every overlapping WT value and generally prevents
  the desired WT-mutant contrast; or
- using those WT values retires the local Davis panel from confirmation before the read and marks all
  Davis-derived records development-only.

No silent reclassification is allowed. A later external claim requires a new source.

### BindingDB

The current download page exposes July 2026 archives and a public TSV schema. BindingDB reports about
3.24 million binding records. Staff-curated records are CC BY 3.0; ChEMBL-derived records retain
CC BY-SA 3.0.

The local file
`dataset/public/open_s/BindingDB_BindingDB_Articles_202607_tsv.zip` matches the published MD5
`f766b06b8d137f7465bcab7ce1f76e56`. It contains one 328,109,536-byte TSV member. Its affinity rows
were not opened. The full July archive is also available with published MD5
`48f8d28a2097ee023946e764ba93de22`.

BindingDB is feasible as a registry source, but mutation strings in target names are not sufficient.
The actual construct sequence, mutation notation, endpoint, assay, DOI/PMID, and ligand parent must
agree. BindingDB and ChEMBL copies of one paper count as one provenance unit.

### ChEMBL variants

ChEMBL documents `variant_sequence_accession`, `variant_sequence_mutation`, mutated residues, and
variant sequence through assays and web services. It also explicitly warns that variant sequences
are not referentially linked to component sequences and that engineered and disease variants require
document review. The local ChEMBL-37 assets are processed pKi/pKd extracts and do not contain the full
variant layer. A full database or bounded API acquisition is therefore required. License is
CC BY-SA 3.0.

### PLATINUM and MdrDB

PLATINUM is reachable and exposes a 687,373-byte flat CSV with more than 1,000 mutation-affinity
records from more than 180 papers, more than 200 ligands, and about 140 UniProt entries. Its current
help/data pages show no explicit data license. It is feasible for schema validation and a manually
reviewed gold subset only after usage and redistribution terms are frozen.

MdrDB is reachable and provides an academic-use license; commercial use requires separate contact.
Its large record count mixes direct binding, cell-line response, resistance annotations, predicted
structures, and computational scores. It may index source papers and direct-binding records, but
GDSC/DepMap response, docking scores, and simulated affinity are permanently excluded from labels.

### Existing KirHub evidence

The historical mutation-anchored KirHub audit already found real mutation-dependent reordering:
222 eligible mutant constructs, 34 genes, 22 kinase families, seven kinase groups, gene-macro
reversal rate `0.1133`, and true-WT minus wrong-WT Spearman advantage `+0.5010`. It uses aggregate
1-uM residual activity, has no raw duplicate measurements, and is one provenance lineage. Its
MDE80 is `0.048` at 34 genes or `0.060` at 22 families for paired SD `0.10`.

KirHub is therefore an ordinal external direction check, not endpoint-compatible pKi/pKd training
data and not broad protein-family evidence.

## Source and observation requirements

`OpenMut-XSource` is feasible to build, but database count must never be confused with provenance
count. Every accepted measurement must retain:

```text
source/version/license/record_id
DOI/PMID/document/assay/institution/provenance lineage
endpoint/relation/value/unit/temperature/pH/replicate
WT accession/construct/sequence
mutant sequence/directed mutation/type
ligand parent/SMILES/scaffold/chemical cluster
evidence location/extraction method/review status
```

Primary labels are direct biochemical exact `Kd` or `Ki`, single amino-acid substitutions, exact
WT-mutant construct pairs, and same-endpoint rectangles. `Kd` and `Ki` are analyzed in separate
strata. IC50, censored values, multi-mutations, insertions/deletions, phosphorylation, and cross-assay
pairs are sensitivity data. Cell response, clinical outcome, qualitative resistance, docking, and
simulated scores are excluded.

Before any model, D1/I0 must report:

- mutations with at least 4, 8, and 16 exact shared ligands;
- independent mutation, canonical base-protein, homology, broad-family, document, institution, and
  provenance-lineage counts;
- repeated directed substitution classes across independent base proteins;
- mutation-position by family and document contingency;
- repeated ligand-edit classes across mutation components and scaffold diversity;
- censoring, replicate covariance, exact-only selection, and reliability;
- comparison-graph connectivity/effective rank and component-level MDE;
- duplicate and source-circularity collapse by original paper, not database.

D1 is label-free only if topology is available in separate metadata or through a preregistered blind
projector that irreversibly drops numeric affinity, relation, and censor fields before any row is
materialized for analysis. The projector must log the source hash, projected schema, dropped columns,
and row count. If this cannot be enforced, the source moves to post-firewall X0/I0 and D1 must report
that its topology is unknown; opening an affinity-bearing table and merely promising not to inspect
one column is not a label-free audit.

DAVIS-Complete alone cannot pass these multi-family/power conditions. The public route remains open
only because BindingDB/ChEMBL/PLATINUM/supplements may supply independent base proteins and documents.
If they yield mainly one-ligand mutation effects, the registry is useful for mutation-main-effect
research but not for ligand reordering.

## Mandatory falsification battery

Every protein coordinate enters the same fixed estimator separately. Coordinates are never
concatenated to claim success.

Protein/mutation arms:

1. position only;
2. WT/mutant identity and aligned one-hot;
3. pocket composition;
4. fixed BLOSUM-derived coordinate;
5. explicit physicochemical mutation delta;
6. pooled ESM-2 difference;
7. BLOSUM-label permutation, directed-substitution-label permutation within matched
   position/family strata, and parameter-matched random PSD;
8. sequence shuffle and position shuffle;
9. matched wrong mutation/target;
10. within-family wrong mutation/target.

Local PLM difference, frozen structure-conditioned CLOCK, and conformer information are deferred
representation arms. A CLOCK arm must use an externally frozen structure map and a
parameter-matched structure-shuffled control; learning extra structure parameters from the same
affinity data is not a coordinate comparison.

Ligand/data arms:

1. ligand potency and ligand identity only;
2. whole-molecule feature difference;
3. physicochemical or Morgan difference;
4. explicit MMP edit where topology supports it;
5. ligand-edit permutation;
6. quartet recombination;
7. additive and parameter-matched random interaction;
8. document/assay/source-only predictor.

Inference is blocked by base-protein, homology, broad family, document, and provenance component.
Pair, quartet, row, seed, and fold counts are descriptive only. Primary performance must improve the
strongest nested non-semantic null by at least the separately frozen material effect, have a
component-blocked lower confidence bound above zero, be detectable at 80% power, survive both
mutation and ligand destruction, and avoid worst-family/source collapse.

## Ordered execution decision

The following stages are conditionally defined and entered into `task.md`. Only metadata-only
`OMUT-D0` is currently active, and it is incomplete. `OMUT-F0` is a current no-go until D0 and the
one-way Davis role policy are frozen; every later stage is blocked.

1. `OMUT-D0`: freeze source versions, rights, hashes, schemas, and safe download plan; no values.
2. `OMUT-F0`: freeze Davis overlap/retirement policy and duplicate-source firewalls.
3. `OMUT-D1`: label-free topology and projected-power audit.
4. `OMUT-X0`: build the evidence-bound BindingDB/ChEMBL/PLATINUM/supplement registry.
5. `OMUT-I0`: audit exact/censored reordering variance, reliability, comparison graph, and empirical
   MDE; no predictor.
6. `OMUT-C0`: run each fixed coordinate in one common low-capacity audit estimator.
7. `OMUT-M0`: only after X0, I0, and C0 pass, fit the directly parameterized centered mutation
   operator with frozen low-dimensional bases and regular `Theta=0` nested null.
8. `OMUT-M1`: multi-seed and held-provenance validation with the full destruction battery.
9. `OMUT-R0`: only after M1, test frozen local PLM or structure-conditioned coordinates separately.
10. `OMUT-T0`: test whether multi-substitution composition is valid; no assumed additivity.
11. `OMUT-T1`: only after T0, test a frozen mutation-derived coordinate in strict dual-cold DTA.
12. `CONFIRM`: use a newly sealed, independent source; DAVIS-Complete is never that source.

The prospective A0 panel remains active in parallel. The existing 384-cell pilot can be made
mutation-aware with six WT-single-mutant pairs from six broad families, 16 shared scaffold-diverse
ligands, two independent sites, and one pKi or pKd endpoint
(`12 constructs x 16 ligands x 2 sites = 384`). Its six independent family pairs give MDE80 about
`0.114` at paired SD `0.10`; it estimates reliability and variance, not a `0.03` predictive gain.
A confirmatory mutation program must expand the number of independently sampled base proteins and
provenance components according to the empirical A0 variance; relabeling the 12 constructs as
12 independent targets is prohibited.

## Cross-review resolution

The new proposals were checked against three independently completed Agent trails rather than treating
the attachment narrative as evidence:

| independent trail | result that survives | challenge to the attachments |
| --- | --- | --- |
| data / LOCK integrity | label-free geometry and the accepted G0 execution chain are reproducible; fixed LOCK still fails the semantic controls | aligned-position signal cannot be relabeled as BLOSUM substitution semantics |
| mathematical / novelty | no pure-model rescue exists; mutation-by-edit algebra overlaps MISO/AXIS/DICE and the regular direct-operator route requires new randomized information | quartet expansion, MMP branding, `UV^T`, cycle loss, and a broad module-novelty claim are rejected |
| paper / adversarial | the direct operator is calibrated only on synthetic outcomes; the real multi-family G0-A substrate is absent and the latest paper audit remains category 3 | structure teachers, relative-pair counts, and database scale do not establish source-independent dual-cold evidence |

The unfinished broader paper-search turns were not counted as evidence. Root cross-review matched the
three attachment claims to the completed independent results and then rechecked the public metadata,
historical ledger, and power arithmetic.

The data/rights conclusion is that public registry construction is feasible but DAVIS-Complete is
neither an independent nor a powered substrate. The mathematical conclusion accepts the estimand and a
directly parameterized centered operator, while rejecting pair expansion, cycle loss, factorized
low-rank nulls, and a broad novelty claim. The adversarial conclusion requires one-way Davis role
freezing, DOI-level source collapse, censor-aware selection, base-protein/family blocking, and matched
wrong-mutation/edit controls.

Cross-checking the three trails produced no surviving disagreement:

- synthetic calibration of the direct operator establishes estimator validity, not real mechanism or
  performance;
- label-free LOCK geometry establishes a reproducible coordinate construction, not fixed BLOSUM
  semantics;
- database availability establishes an acquisition path, not adequate independent components,
  licenses, provenance, or power;
- millions of teacher pairs or database rows do not change the base-protein/provenance inference unit.

The apparent disagreement is therefore resolved by separating data recovery from model authorization:
OpenMut is a feasible data program; Delta2Rank/OMRO, PLM, structure, ensembles, and strict DTA transfer
remain conditional. The current scientific conclusion is therefore category 3.
