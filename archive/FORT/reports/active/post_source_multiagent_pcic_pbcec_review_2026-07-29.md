# Post-source multi-agent review: PCIC, PB-CEC, and the remaining DTA route

Date: 2026-07-29
Status: complete research review; no predictor authorized

## Executive decision

The completed PBCNet2.0, BioLiP2, and PSICHIC audits do not provide the
independent, row-resolved bridge required by the proposed PB-CEC collective
completion model. PB-CEC is therefore stopped before implementation on the
current public substrate. Its useful remainder is a prospective
cycle-closing acquisition and sensitivity-analysis framework.

The proposed PCIC null-space construction is not a new information source in
this repository. `P0-Cycle-A` already projected raw ChEMBL TRAIN exact cells
onto the full orthogonal complement of

```text
Z = [target, ligand, assay, document]
```

using a sparse Frisch-Waugh-Lovell/Hodge projection. That projection covers
all of `ker(Z^T)`, not only four-cell rectangles. It already stopped on
biological amplitude and concentration:

| Endpoint | Exact cells | Projected-dimension LCB95 | Projected affinity SD | Dominant failure |
| --- | ---: | ---: | ---: | --- |
| pKi | 231,090 | 103,877.16 | 0.35596 | below 0.5; top 1% ligands carry 51.45% of residual energy |
| pKd | 28,749 | 14,049.88 | 0.52342 | top 1% ligands/documents carry 64.01%/93.10% |

The recorded verdict is `P0_CYCLE_A_BIOLOGICAL_FAIL_STOP`. Enumerating a
different cycle basis cannot add information or turn its basis vectors into
independent samples.

One narrow increment remains untested:

```text
A = P_perp(D Z) D X
M = A^T A
```

where `D = W^(1/2)` and `X` is one predeclared direct target-ligand operator
coordinate. The remaining question is whether the already-known
nuisance-free edge space contains a stable, provenance-replicated,
strict-dual-cold-estimable operator direction. This question is frozen as
`PCIC-O0`.

## What the three reviews changed

The user-supplied reviews were evaluated independently for mathematical
identifiability, novelty/falsification risk, and local data feasibility. The
data review was completed by the primary agent against the local artifacts
after the dedicated worker exceeded its bounded run.

The reviews agree on four corrections:

1. Nonzero cycle-space dimension is algebraic redundancy, not evidence of a
   transferable biological interaction.
2. Adding provenance fixed effects can only retain or reduce the old
   nuisance-free space on the same rows.
3. Provenance cancellation is not provenance replication.
4. Strict dual-cold prediction requires query estimability in the frozen
   operator coordinates; ridge cannot manufacture it.

## PCIC-to-P0 delta audit

### Already completed by P0

`research/p0_cycle.py` and `reports/active/p0_cycle_a.json` already provide:

- raw ChEMBL TRAIN pKi/pKd exact-cell construction;
- target, ligand, assay, and document nuisance columns;
- the full sparse orthogonal projection rather than rectangle enumeration;
- LSMR/LSQR agreement, KKT residuals, and idempotence checks;
- Hutchinson projected-dimension estimates;
- label residual amplitude and group-energy concentration;
- leave-dominant-ligand-out diagnostics;
- the explicit warning that cycle rank is not independent sample size.

Consequently, these claims are prohibited:

- that prior work only tested rectangles;
- that arbitrary cycles supply a newly discovered information source;
- that circuit expansion increases the independent sample count;
- that nonzero nullity establishes a transferable interaction;
- that PCIC is the first nuisance-null DTA construction in this project.

### Genuinely untested increment

Only the following operator-level quantities are new:

- numerical rank of `P_Z_perp X` for one frozen semantic coordinate;
- fold-common operator directions under lineage-disjoint deletion;
- information concentration by lineage, homology, and scaffold;
- strict-dual-cold query row-span residual;
- semantic-coordinate stability relative to entity-shuffled controls.

These quantities can fail even when `dim ker(Z^T)` is large. Conversely, a
positive operator rank alone is insufficient: arbitrary random coordinates
can also have nonzero rank, so semantic stability and provenance replication
are mandatory.

## Algebraic corrections

### Nested metadata do not create a new projection

If DOI, PMID, institution, patent family, or a transitive lineage identifier
is a deterministic function of document, then

```text
B_lineage = B_document R
```

and its columns are already in the document indicator span. Homology and
scaffold are similarly nested in target and ligand identity. Such fields are
still essential for blocked replication and effective-sample accounting, but
they do not add a new nuisance-null dimension.

If a genuinely nonnested context column is added, then

```text
col(Z_old) is a subset of col(Z_new)
P_new <= P_old
X^T P_new X <= X^T P_old X.
```

The new metadata cannot increase residual label energy or operator
information on the same row population.

### Weighted projection

For predeclared reliability weights, define:

```text
D = W^(1/2)
Z_tilde = D Z
X_tilde = D X
P = I - Z_tilde Z_tilde^+
M = X_tilde^T P X_tilde.
```

An original-scale contrast `c^T y` must obey `Z^T c = 0`. If `Q` is a
whitened null basis, its pseudo-outcome is `Q^T D y`, not `Q^T y`.
Circuits sharing edges are correlated; their covariance is `C Sigma C^T`.
They may not be treated as independent training rows.

The primary implementation should therefore use an implicit projector and a
direct operator. Searching for a sparsest basis in a general augmented
categorical design is not the main algorithm; sparse null-basis problems are
NP-hard in general. Circuits may be shown later only as interpretive
examples.

### Query estimability

For a strict-dual-cold query with operator feature

```text
x_q = v_l tensor u_t,
```

the direct linear functional is estimable only when `x_q` lies in the row
span of `A`. The required diagnostic is:

```text
r_q = ||(I - A^+ A) x_q|| / ||x_q||.
```

Adding `lambda I` makes every ridge leverage finite and therefore hides null
directions. Ridge is prohibited in the O0 rank and query-span gates.

## PB-CEC review

The proposed shared entity factors and relation-specific heads are a form of
collective matrix factorization/completion. Missing-modality marginalization
is a valid likelihood operation, but it does not create evidence in an
unobserved modality. An MNAR observation model is likewise not point
identified without overlap, exclusion restrictions, or an explicit
sensitivity assumption.

The current exact-bridge audit is decisive:

- PBCNet2.0 has no original BindingDB row, assay, document, or source-lineage
  mapping for its generated pairs;
- BioLiP2 provides real complexes and contacts, but its closed topology has a
  79.9767% giant component and does not provide the required independent
  affinity bridge;
- PSICHIC exposes weights but not complete XL training membership or
  row-level source lineage.

PB-CEC therefore cannot currently justify a globally shared physical latent
space or a strict-dual-cold point prediction. It is parked rather than
silently weakened.

If revisited after new bridge measurements, a bounded direct convex operator
may produce a predeclared sensitivity envelope. A deep-network ensemble or
loss-level set must not be called a statistical identified set without a
coverage theorem.

## Local data feasibility

### Present locally

The ChEMBL 37 substrate contains:

- `chembl37_pKi.jsonl.gz` and `chembl37_pKd.jsonl.gz`;
- a TRAIN/development/confirmation registry with target, ligand
  connectivity, endpoint, scaffold, accession, homology cluster, assay and
  document summaries, and split assignment;
- a 44.18 MB frozen target ESM-2 archive;
- a 34.90 MB frozen ligand-feature archive;
- raw row identifiers for activity, molecule, target, assay, and document.

These assets are enough for a label-blind operator design after a dedicated
safe projection is implemented.

### Missing or not yet certified

The local registry does not contain a certified row-level DOI/PMID/patent
family/institution/site/transitive-lineage closure. The raw row builder also
contains numeric `value_nM` and `pK`, so O0 must not reuse its ordinary JSON
decoder. It needs a guard that rejects or byte-skips protected value keys and
emits identifiers only.

Any ChEMBL document/assay metadata enrichment must be source-versioned,
hashed, and bound to the ChEMBL 37 rows. Unversioned live API contents cannot
enter the primary result. If exact closure cannot be recovered, O0 stops at
metadata feasibility rather than substituting document count for independent
provenance.

## Frozen route order

1. Run `PCIC-O0-P`: protected-key-safe exact-cell and provenance closure.
2. If O0-P passes, run `PCIC-O0-I`: CUDA operator information, blocked common
   rank, concentration, shuffled controls, and query-span audit.
3. A failure ends public-data mathematical reparameterization. The next
   deliverable is a provenance-separated cycle-closing acquisition design,
   not another affinity predictor.
4. A pass may authorize a new T1 preregistration. HCRR remains an
   equal-budget reliability baseline; it is not restored as the main route.
   Rectangle-only, full direct operator, entity-shuffled, no-provenance, and
   PB-CEC arms must be included.
5. Only a T1 pass under wrong-target, wrong-ligand, provenance-disjoint, and
   risk-coverage controls can support the candidate framework:
   stable identifiable operator, selective abstention, and cycle-closing
   acquisition.

## Compute boundary and schedule

For pKi, `m = 231,090` and the primary `8 x 8` operator has `p = 64`.
The dense `X` buffer is approximately 59 MB in float32. The Gram calculation
is approximately `m p^2 = 9.5e8` multiply-adds. GPU memory is not the limiting
factor; sparse residualization, metadata closure, blocked repetitions, and
controls dominate wall time.

- CPU only: compressed-stream key-safe parsing, RDKit canonicalization,
  metadata joins, transitive lineage closure, sparse incidence construction,
  and a small number of LSMR/KKT verification probes.
- CUDA in `D:\anaconda\envs\drug`: feature centering/projection where
  applicable, multiway demeaning/scatter residualization, Gram matrices,
  eigendecomposition/SVD, fold-common subspace calculations, and any later
  trainable model.

Expected milestones:

- one working day: safe projection, metadata availability decision, frozen
  block registry, and a reproducible O0-P pass/stop;
- three working days total: pKi/pKd O0-I, five blocked folds, concentration,
  entity-shuffled controls, query span, and one `STOP` or `REQUEST_T1`;
- seven working days total, only if O0 passes: equal-budget exploratory
  direct-operator/HCRR/rectangle/PB-CEC controls and risk-coverage result.

Seven days cannot establish a publishable high-innovation predictor. With an
O0 and T1 pass, independent-source validation and repeated software runs are
estimated at two to six weeks. Prospective cycle-closing measurements are a
month-scale dependency.

## Claim boundary

The strongest defensible innovation candidate is not "a new cycle model."
It is a provenance-audited strict-dual-cold DTA loop that:

1. learns only in a cross-source stable identifiable operator subspace;
2. refuses queries outside that estimable subspace;
3. reports a bounded sensitivity envelope rather than unjustified certainty;
4. converts identification failures into cycle-closing measurement choices.

This remains a conditional framework claim. No high-accuracy predictor,
transferable mechanism, or first-in-literature claim is established.

## Primary literature used in the audit

- Godsil and Guo, cycle-space algebra:
  https://arxiv.org/abs/1609.09118
- Jochmans and Weidner, network fixed effects:
  https://arxiv.org/abs/1608.01532
- Gunasekar et al., collective matrix completion:
  https://proceedings.mlr.press/v38/gunasekar15.html
- Agarwal et al., causal matrix completion:
  https://proceedings.mlr.press/v195/agarwal23c.html
- Durrande et al., ANOVA kernels:
  https://doi.org/10.1016/j.jmva.2012.08.016
- Coleman and Pothen, sparse null bases:
  https://epubs.siam.org/doi/10.1137/0607059
- Chiang et al., dyadic DML:
  https://arxiv.org/abs/2110.04365
