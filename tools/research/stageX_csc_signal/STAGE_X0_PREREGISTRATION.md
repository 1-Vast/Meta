# Stage X0 preregistration — pipeline correction, data acquisition and instrument qualification

Frozen **2026-08-17, before any Stage X0 primary result was read.** This is a
new, independent preregistered cycle. Stage S–W artifacts are read-only
evidence and are not modified. Stage X0 only runs statistical audits and small
normally gradient-trained diagnostic models; ridge, pseudoinverse, closed-form
adapters and test-time query gradients are forbidden.

## 0. Scope and ordered gates

X0 is a qualification stage. Its outputs are instruments, not biological
claims. X0 passes only if **all** six instruments below pass their frozen
checks. X0-P and X1 are separate preregistrations and are not entered until
X0 passes.

## 1. Six required instruments

### I1. Synthetic planted-signal control
On the real data graph (targets, ligand pairs, protein components) and the
real censoring pattern, plant an additive protein × ligand interaction of
strength `tau in {0.2, 0.4, 0.8, 1.6}` log units. For each strength, a
normally trained low-capacity bilinear diagnostic must recover:
- sign accuracy vs planted truth >= **0.70** at `tau >= 0.8`;
- Spearman between fitted interaction component and planted interaction
  >= **0.30** at `tau >= 0.8`;
- correct protein dependence: correct-protein feature arm must beat
  ligand-only by >= **0.05** sign accuracy at `tau >= 0.8`.
If the pipeline cannot recover a planted `tau=0.8` signal at the expected
graph size, **no real-data negative may be interpreted**; the pipeline must be
repaired first.

### I2. Representation-capability pre-check
For every admissible protein representation, report
`r_pair = ||x(p) - x(q)|| / median_inter_protein_distance`. A representation
passes expression capability only if the median over labelled WT/mutant or
single-residue-substitution pairs is >= **0.05**. At least five
representations are measured:
global pooled ESM; mutation-position residue token; mutation-centered local
window; KLIFS aligned-pocket token; capacity-matched random representation.
Representations that fail expression capability are **excluded** from
explaining a biological null.

### I3. ID-equivalence test
On the largest admissible single-platform panel, compare real protein
representations against: free target embedding, family id, nearest-pocket
lookup. The real representation must beat family id / nearest-pocket lookup
with a component-cluster bootstrap lower bound >= **0.05** in at least one
preregistered metric (signed-difference sign accuracy or Spearman) on
protein-cold pairs. If it cannot, the representation is recorded as
identity-equivalent and not admitted to X1.

### I4. Censoring instrumentation
All matrices retain raw censoring direction (`>`, `<`, missing) and detection
floor semantics. Four analyses are implemented and reported:
- determinate-only;
- interval-censored likelihood (one-sided Huber/quadratic with explicit mask);
- sign/ranking-only;
- floor-imputed sensitivity (negative control only, never primary).
Davis `Kd=10000 nM`, Metz `pKi=4.0`, Klaeger `pKd=5.0` are never treated as
exact continuous labels.

### I5. Cluster-level inference
Every primary interval resamples protein/pocket components (bootstrap 2,000
draws, seed `20260820`); scaffold, protein-pair and family levels are reported
as sensitivity. Raw rows are never presented as independent samples.

### I6. Live integrity assertions
A test suite fails the whole run if any assertion fails:
label orientation (named anchors: staurosporine broad; imatinib
ABL1/KIT/PDGFR; lapatinib EGFR/ERBB2); Kd↔pKd / Ki↔pKi direction; ligand-pair
order; protein/ligand id mapping; CSC antisymmetry and identity-zero; dead
regularizers; gradient coverage; permutation controls actually destroy the
intended information.

## 2. Stage X0-D external data acquisition contract

Systematically verify and fetch, with DOI/URL/file/licence/SHA-256/date/semantics:
Duong-Ly mutant panel; Anastassiadis panel; Davis/KINOMEscan WT/mutant
constructs; Karaman panel; PKIS/PKIS2; EGFR/ABL1/ALK/KIT/BRAF gatekeeper and
resistance panels; Stanford HIVdb genotype–phenotype; same-assay ortholog
panels; KLIFS and GPCRdb aligned pockets; any matched-variant panels found in
primary sources.

Rules:
- Same vendor does not imply same comparable experiment. Duong-Ly and
  Anastassiadis are linked only after compound identity, construct, label
  direction, ATP conditions and measurement semantics are verified.
- Different label systems are analysed separately; only direction replication
  is allowed, never merged regression.
- CC BY-NC-ND / no-derivative data: only acquisition script + manifest +
  local cache are retained; derived data are not committed.
- A source is `blocked` only after the official URL, mirror and alternative
  repository are all tried and recorded with timestamps and error codes.

## 3. Frozen X0 pass conditions

X0 passes iff:
1. I1 planted-signal recovery passes at `tau=0.8`;
2. I2 at least three admissible representations pass expression capability,
   including at least one local representation;
3. I3 passes for at least one representation;
4. I4 all four analyses execute without integrity failure and produce
   machine-readable outputs;
5. I5 cluster bootstrap is the only primary interval;
6. I6 all live assertions pass.

If X0 fails, the specific failed instrument is fixed in a corrected successor
stage; no real-data biological conclusion is drawn from the failure.

## 4. Artifacts

`X0_PREREGISTRATION_SHA256.txt`, `X0_INSTRUMENTS.json`, `X0D_DATA_AUDIT.json`,
`X0_RESULT.json`, `X0_REPORT.md`, `commands.jsonl`, `tests/`.
