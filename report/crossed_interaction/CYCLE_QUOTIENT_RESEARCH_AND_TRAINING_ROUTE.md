# Cycle-quotient research and minimal training route

Updated: 2026-08-10.

## Decision

The interaction-quotient proposal is mathematically valid and is the strongest
remaining route, but it does not by itself repair the failed ChEMBL dependence
Gate.

```text
CYCLE_QUOTIENT_ALGEBRAICALLY_AVAILABLE_BUT_DEPENDENCY_NOT_REPAIRED
```

The immediate action is a new source census, not GPU training.

## What the frozen ChEMBL graph says

The label-blind census uses the byte-verified X0 cells and computes, for every
bipartite panel graph,

```text
d_cycle = |E| - |P| - |L| + number_of_graph_components.
```

| endpoint | panel cycle dimension | dependency components | largest component share | exact-assay cycle dimension |
|---|---:|---:|---:|---:|
| Ki | 29,677 | 36 | 0.4888 | 0 |
| Kd | 3,279 | 12 | 0.4608 | 0 |

All 597 Ki and 34 Kd document/context panels contain cycles. However, none of
the 1,563 Ki or 162 Kd exact assay graphs spans more than one target, so exact
assays contain no crossed cycle at all. Nearly half of each endpoint's raw
cycle dimension sits in one frozen dependency component. Cycle coordinates
therefore recover algebraic efficiency within panels but do not create new
independent biological units or reverse X1A-R.

Machine evidence:
`report/crossed_interaction/cycle_quotient_feasibility/census.json`.

## Correct mathematical object

For a panel edge vector `y` and protein/ligand incidence matrix `X`, use a
weighted orthonormal basis `U` of the left null space:

```text
U^T W X = 0.
```

The affinity target is `U^T W y`, not a list of enumerated cycles. Given a
frozen measurement covariance `Sigma`, whiten only supported directions:

```text
R = eigentruncated_inverse_sqrt(U^T W Sigma W U)
t = R U^T W y.
```

This makes the loss invariant to the arbitrary cycle basis. It does **not**
make document, assay, homology or scaffold dependencies disappear. Those stay
as inference and split units. Ill-conditioned covariance directions must be
discarded before values are inspected.

## One model, not a module stack

Keep ESM2, the ligand graph encoder and P1B frozen. P1B is a geometric support
measure, not an affinity teacher. Construct one four-dimensional bounded
mechanism vector from the existing atom-local, residue-local and distance-bin
states:

```text
x(P,L) = [polar complementarity,
          signed charge compatibility,
          hydrophobic enclosure,
          aromatic/steric packing] in R^4.
```

Each coordinate must depend jointly on ligand atom and protein residue. Before
aggregation, remove atom- and residue-degree main effects under the actual pair
mask. At panel level, apply the same protein/ligand quotient used for affinity.
No raw PLM embedding, pair map or target identifier may enter the output.

The only trainable affinity object is a small shared response on these four
coordinates. Its primary loss is

```text
L_cycle = mean_panel ||R U^T W (y - q_theta(P,L))||^2,
```

balanced first within panel and then across dependency components. Additive
protein-only or ligand-only functions have exactly zero quotient gradient.
This is the route's main advantage over adversarial shortcut penalties.

## Data route

1. **BindingDB curated articles, metadata census.** Use the release-pinned
   article-only TSV, assay mapping and target sequences. Exclude BindingDB rows
   sourced from ChEMBL, PDSP and PubChem. Freeze DOI/patent, assay, construct,
   endpoint, sequence, connectivity and publication-time provenance before any
   value is opened.
2. **BindingDB patents remain a separate stratum.** They are larger but cannot
   rescue article evidence by pooling after results are seen.
3. **Kinobeads 1,183-compound profiling is an auxiliary positive control.** It
   is dense and useful for selectivity, but two-dose apparent binding values
   and censoring prevent it from serving as clean universal Ki/Kd evidence.
4. **Davis remains frozen.** It may be opened only by a separate explicit
   positive-control authorization and can never be independent confirmation.
5. **PLATINUM/mutation data are attribution validation only.** They can test
   whether ligand-conditioned residue directions agree with mutation effects;
   they do not replace population affinity training.

BindingDB now exposes publication dates, assay identifiers, target sequence and
construct annotations, which makes this governance feasible. It also contains
large ChEMBL-derived and patent strata, so source-level deduplication is
load-bearing rather than optional.

## Execution gates

### CQ-R0: source identifiability, no values

- at least 60 dependency components;
- largest component share at most 0.25;
- conservative effective quotient rank at least 245;
- exact target sequence/construct, ligand connectivity, document/time and assay
  provenance materialized and hashed;
- ChEMBL overlap removed by source plus DOI/measurement closure;
- licence and redistribution status recorded.

Failure means the source cannot support this model claim. It does not authorize
a larger network.

### CQ-R1: interaction existence, CPU

Open uncensored Ki and Kd separately. Estimate measurement noise from
replicates. On frozen quotient coordinates, require noise-corrected interaction
variance above zero with component-level confidence and above the predeclared
biological effect floor. A dense profiling panel can verify implementation but
cannot pass the cross-family source claim.

### CQ-R2: frozen-basis sufficiency, CPU or small GPU probe

Using the same folds, compare the four frozen mechanism coordinates against:

- quotient chance/null;
- protein-only and ligand-only marginals;
- foreign ligand;
- deranged protein;
- degree-preserving pair null.

Only a correct-pair increment with a component-level lower bound above zero
authorizes trainable `q_theta`.

### CQ-R3: minimal GPU training

Train only the four-coordinate response. Do not unfreeze ESM2/P1B, add
cross-attention, introduce a second GNN, or jointly train a few-shot adapter.
Measure end-to-end examples/cycle-coordinates per second and data-loader stall;
GPU percentage is diagnostic, not a reason to enlarge the model.

### CQ-R4/R5: confirmation and few-shot section

After independent source confirmation, fit the target section in
`[1,x]^T in R^5` by closed-form positive-ridge regression. Report support rank,
conditioning, query row-space coverage and abstain off coverage. Only then can
a bounded, source-frozen summary be considered for biological `z`.

## Claim boundary

Cycle-space projection and two-way main-effect removal are established linear
algebra/statistical ideas, not a new theorem. The project-level contribution is
their integration with marginal-audited biochemical coordinates, partner
controls, support-identifiable sections and the unchanged law operator
`A(F,z)=K(B(z)F(z))`.

This route preserves the original goal: few-shot prediction for a new target.
Dense kinase panels are used to establish learnability, not to redefine the
model as kinase-only. Cross-family BindingDB components remain necessary before
the biological statistic can be admitted.

## Primary sources

- BindingDB 2024 FAIR knowledgebase and current downloads:
  https://www.bindingdb.org/rwd/bind/gkae1075.pdf
- BindingDB TSV provenance/construct fields:
  https://www.bindingdb.org/rwd/bind/chemsearch/marvin/BindingDB-TSV-Format.pdf
- ChEMBL 37 and assay provenance:
  https://www.ebi.ac.uk/chembl/beta/
- Kinobeads 1,183-compound selectivity matrix:
  https://www.nature.com/articles/s41589-023-01459-3
- Hodge/cyclic decomposition precedent:
  https://arxiv.org/abs/0811.1067
- PLATINUM mutation-affinity database:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4384026/
