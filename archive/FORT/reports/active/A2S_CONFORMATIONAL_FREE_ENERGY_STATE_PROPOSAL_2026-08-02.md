# A2S Conformational Free-Energy State (CFES)

Date: 2026-08-02  
Branch: `research/a2s-conformational-free-energy-state-20260802`  
Status: **mechanism hypothesis; implementation restricted to Gate C0**

## Exploration objective

Learn a transferable target-adaptation mechanism from abundant source targets
that uses only k={1,3,5} measurements from a strictly unseen target to make a
small, target-specific intervention in query-compound ranking without changing
the frozen support-free DTA path.

The branch asks one narrow question:

> Can source tasks teach how sparse affinity measurements update the population
> of a few physically observed protein conformational states, so that state-
> specific ligand response surfaces improve ranking for an unseen target?

This is not similarity weighting, retrieval, target calibration, uncertainty
estimation, or a larger sequence-ligand encoder. The prospective learned object
is a support-to-conformational-state operator.

## Why the previous mechanism stopped

PIRS learned two query-dependent interaction coordinates from ligand features
and coarse ESM sequence segments. Its synthetic harness passed, but its held
probe full-support oracle had CI gain +0.00460 with 95% interval
[-0.00137,+0.01023], and its low-similarity oracle mean was -0.00077. Correct
support did not beat deranged labels, norm-matched wrong-target residuals,
protein-zero, or protein transplant. Segment conditioning did not beat
ligand-only, pooled protein, or frozen random coordinates.

Therefore PIRS failed at representation admission. No learned support operator
can rescue it, and no capacity or training extension is allowed.

The sequential record now rules out:

- local residual transport outside its measured Tanimoto range;
- a small shared ligand-response subspace;
- arbitrary coordinate sparsity and low-rank heads;
- uncertainty or margin selection as the primary mechanism;
- assay grouping and explicit MMP relations on the current passive coverage;
- coarse protein-segment interaction channels.

It does not test a physically observed conformational ensemble. That is the
new information object.

## General biological principle

Proteins occupy ensembles of conformational substates. Ligands bind with
different preferences to those substates and shift their populations. The
observed equilibrium affinity is therefore not, in general, the score of one
static protein structure. In a finite-state approximation,

\[
G_t(x)=-\beta^{-1}\log\sum_{c=1}^{C}
\pi_{t,c}\exp[-\beta E_\theta(P_{t,c},x)],
\]

where `P_(t,c)` is a physical conformation, `E_theta` is a transferable
state-specific interaction score, and `pi_t` is the target/construct/assay
state population. A support set can act as a set of chemical probes of
`pi_t`; the query compound need not be chemically close to a support compound
if both have discriminable preferences across the same physical states.

This principle is general across protein families, but its deployable value is
conditional on three facts that must be measured rather than assumed:

1. physical state ensembles are available at broad target coverage;
2. state-specific ligand scores differ in a reproducible, ligand-dependent way;
3. one or two centered population coordinates are identifiable at k<=5.

## Current feasibility evidence

Label-blind preliminary coverage on source `fit` metadata found:

- 232/232 targets map to unique UniProt accessions;
- AlphaFold DB v6 returned structures for 230/232 targets (99.1% in the first
  request; one transient request error and one 404);
- PDBe returned at least one experimental PDB structure for 217/232 targets;
- 209 targets have at least two PDB structures, 193 have at least three, and
  the median is 14 distinct PDB entries;
- the local PLINDER preflight records 9,264 leakage-filtered structural train
  systems, 1,375 protein clusters, 2,696 pocket clusters, and 1,828 certified
  experimental apo-linked systems.

These are coverage observations, not mechanism results. The previously
reported PLINDER atom caches are not currently present in the workspace and
must be reconstructed or replaced before a semantic gate. A preliminary schema
inspection also exposed the public processed PLINDER `affinity` column without
using or retaining its values. That processed registry is now treated as
outcome-exposed exploratory data and cannot be an untouched CFES confirmation
set. Gate C0 uses raw structural annotations only.

## Candidate mechanisms

| Candidate | Transferable learned object | k<=5 state | Prior-art boundary | Decision |
|---|---|---|---|---|
| C1. Conformational Free-Energy State (CFES) | state-specific physical response scorer and support-to-population operator | k1 no rank state; k3 one logit; k5 two logits | conformer ensembles and mixture models are known; the candidate claim is a measured few-shot population update that beats analytic inference | **selected** |
| C2. Pocket-choice response state | support selects among a few label-blind pockets | one or two pocket logits | close to mixture-of-experts and allosteric-site selection; many assays have a fixed known site | auxiliary control |
| C3. Reorganization-penalty state | support estimates one target-specific bound/free reorganization penalty | one scalar | closest to IPBind/PReorg-FEP; a scalar often cannot change ranking | nested baseline |
| C4. Contact-topology switch | support activates one of a few residue-ligand contact graphs | one discrete switch | resembles interaction fingerprints and response modes; contact assignment may be unavailable at deployment | defer |
| C5. Conformer-gradient program | support emits one or two learned gradient steps on state logits | one or two directions | closest to Meta-SGD/LEO and can hide an estimator inside optimization | equal-budget meta baseline |
| C6. Kinetic residence-state adapter | support infers slow/fast binding mode | one class | pKi is an equilibrium endpoint and does not identify kinetic state | reject |

C1 is selected because it changes the biological observation, has explicit
state capacity, and makes a distant-query prediction through a shared physical
state rather than support-query similarity.

## Mathematical mechanism

Let `mu(p,x)` be the frozen support-free prediction. For target `t`, construct
at most three label-blind conformational states from public experimental
structures, with AlphaFold as a typed fallback. A shared scorer emits bounded
state energies:

\[
e_{t,c}(x)=E_\theta(P_{t,c},x),\qquad c\in\{1,2,3\}.
\]

The support-free population prior is `pi^0_t`, obtained without affinity labels.
The support residual is `r_i=y_i-mu(p_t,x_i)`. The learned adaptation operator
receives centered residual evidence and state-energy contrasts:

\[
u_t=A_\psi\left(
\{(e_t(x_i)-\bar e_t(x_i),\ r_i-\bar r_S)\}_{i=1}^{k},
p_t,k\right).
\]

The budget mask is structural:

\[
M_1=(0,0),\qquad M_3=(1,0),\qquad M_5=(1,1).
\]

The adapted population is

\[
\pi_t=\operatorname{softmax}(\log\pi^0_t + M_k\odot u_t),
\]

with a reference-state logit fixed to zero. The bounded ranking intervention is

\[
\hat y_q=\mu(p_t,x_q)+
b\tanh\left(
\frac{F_\theta(\pi_t,e_t(x_q))-F_\theta(\pi^0_t,e_t(x_q))}{b}
\right).
\]

Thus support removal, k=1 ranking, or `u_t=0` exactly recovers the frozen base.
The support labels are never averaged into the query prediction and no query
uses support-query chemical similarity in the main path.

## What meta-learning learns

The final trainable object has two load-bearing parts:

1. `E_theta`: a transferable state-specific interaction scorer trained first
   on external structural supervision and then shaped on source fit episodes;
2. `A_psi`: a learned map from sparse state-discriminating residual evidence to
   one or two population-logit interventions.

An empirical-Bayes, ridge, or maximum-likelihood solve for `pi_t` is a required
baseline. CFES is a meta-learning mechanism only if `A_psi` improves over the
identical admitted `E_theta` with those analytic solvers.

## Why k<=5 could identify it

After removing the unknown residual level, k observations provide at most
`k-1` independent contrasts. CFES uses zero rank degrees at k=1, one at k=3,
and two at k=5. The state dimension is fixed in code and cannot expand through
sparsity regularization or a hidden task embedding.

Identifiability additionally requires support ligands to have different energy
contrasts across conformations. Gate C1 must report the singular values of the
support state-contrast design and the fraction of episodes for which its
effective rank is at least the active state dimension. No learned operator is
allowed to operate on rank-deficient episodes except by exact no-op.

## Prior-art boundary

- Conformational selection and population shift are established biology
  (`10.1038/nchembio.232`).
- Protein reorganization thermodynamics can select ligand-compatible
  conformations (`10.1021/acs.jcim.4c01612`).
- Ensemble docking is standard and cannot be claimed as new
  (`10.1021/acs.jpcb.8b11491`).
- PLINDER provides similarity-aware protein-ligand structural data and linked
  apo/predicted structures (`10.1101/2024.07.17.603955`).
- IPBind already learns single-state bound-minus-free interatomic potentials
  (`10.1109/OJEMB.2026.3667030`).
- HyperPCM, DrugBAN, and PSICHIC already cover target conditioning and learned
  interaction representations.

The possible contribution is therefore narrow:

> a budget-matched, support-identified conformational population state for
> unseen-target affinity ranking, with exact frozen-base nesting and evidence
> that a learned support operator adds value beyond analytic population fitting.

Neither structural encoding, conformer pooling, bound/free subtraction,
softmax mixtures, nor meta-learning alone is novel.

## Ordered gates

### C0A - reproducible label-blind coverage

Use source metadata, UniProt mappings, AlphaFold/PDBe public APIs, and raw
PLINDER structural annotations only. Do not load ChEMBL affinity or any source
`probe`, `locked`, or recipient label. Record target, component, accession,
structure count, structure dates, and hashes.

Admission requires broad fit-component coverage, at least two physical states
for a majority of fit targets, and an AlphaFold fallback for at least 90% of
targets. Coverage may justify a restricted applicability domain but cannot by
itself justify training.

### C0B - external structural semantic gate

On protein-, pocket-, ligand-, and provenance-disjoint raw structural splits,
test whether the minimum physical state scorer distinguishes the observed
ligand-compatible pocket/conformation from:

- a single-structure arm;
- uniform state averaging;
- structure transplant from another protein;
- coordinate randomization;
- protein-free ligand scoring;
- ligand permutation;
- state duplication and state-order permutation.

The mechanism must beat the strongest matched control on every held fold and
lose at least 70% of its incremental effect under the corresponding physical
destruction. No affinity label can rescue C0B.

### C1 - fit-only affinity representation gate

Freeze a component-level split of source `fit` into meta-train, inner
validation, and untouched fit-audit roles before reading C1 outcomes. The
historical source `probe` role is not reused. The frozen support-free base is
unchanged.

Require both a full-support population oracle and k=3/5 analytic inference to
improve held-scaffold ranking, including support-query Tanimoto below 0.35.
Correct structure must beat AlphaFold-only, single-state, protein transplant,
state permutation, ligand-only, PIRS, and random-state controls. A gain confined
to local chemistry fails.

### C2 - learned support operator

Only after C0 and C1 pass, freeze `E_theta` and train `A_psi` on meta-train
episodes. Require learned-minus-analytic lower confidence bounds above zero at
k=3 and k=5 on untouched fit-audit components. Required controls are support
removal, wrong target, random support, label permutation, coordinate
permutation, state duplication, frozen/random operator, KRR, ridge/EB, and an
equal-budget gradient meta-learner.

### C3 - complete source breakthrough

Only after C2 passes may a one-time source `locked` gate be preregistered. A
major breakthrough additionally requires positive CI, Spearman, NDCG@10 and
pairwise-proper effects; reduced negative transfer; no support-free
degradation; and exact support-removal no-op. Recipient labels remain sealed.

## Stop rules

Stop CFES without capacity rescue if any of the following holds:

1. C0 lacks deployable structural coverage or the external semantic gate fails;
2. state-specific scores are matched by protein-free, single-state, or random
   structure controls;
3. the full-support physical-state oracle lacks held-component ranking
   headroom;
4. k=3/5 state-contrast designs are usually rank deficient;
5. gains are restricted to chemically local queries;
6. an analytic state solve matches the learned operator;
7. the effect depends on structures containing the scored ligand or on a
   target/ligand/template/provenance overlap;
8. the frozen support-free path changes.

## Generalizability and value assessment

**Biological generalizability: plausible and substantially stronger than
PIRS.** Conformational ensembles and free-energy population shifts are not
protein-family-specific, and current fit targets have broad structural
coverage.

**Statistical generalizability: unverified.** Nothing yet shows that two state
logits explain held-target affinity residuals or that random k<=5 support
ligands discriminate those states.

**Potential value: high if C1 and C2 pass.** It would provide a mechanism for
distant-query intervention grounded in a shared physical target state, address
the central TRACE distance limit, and make the support module auditable through
state destruction. It would also be directly useful for ranking compounds
whose preferences differ across target conformations.

**Maximum risk:** public structural ensembles may be biased toward well-studied
targets and ligand-bound conformations; a state mixture may collapse to a
single static score; approximate pocket/ligand compatibility may be too noisy;
and k<=5 random supports may not span the state contrasts. The branch is not a
breakthrough until those risks are falsified experimentally.

## Promotion rule

No CFES file may enter `model/` or `script/` before C0-C3 all pass. If that
happens, core and utility files are copied once, tested, hashed, and frozen.
Until then all work remains under `research/`, `tests/`, and `reports/active/`.

