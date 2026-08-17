# Stage W0b preregistration — Core Task 1 data, estimand-hierarchy, censoring and positive-control audit

Frozen **2026-08-17, before any W0b census metric was read.** This is a
corrected successor to `stageW_soft_mmp` W0, which measured only a single soft
MMP layer on Davis/KIBA. Stage W0b does not modify Stage W's frozen
preregistrations; it supersedes their GO/NO-GO by completing the mandatory W0
audit. Stage W W1 is paused at data/model preparation; no W1 training or
training metric has been read.

## 0. Scope

CPU/statistical audit only. No neural model is trained here. Each dataset is
kept separate; labels are never merged across platforms or across datasets.

Datasets audited from local, hash-recorded assets:
* Davis (`dataset/raw/dta/davis.tab`): single-platform Kd panel candidate.
* Metz (`dataset/raw/crossed_panels/kinase_panels/metz_matrix.csv` +
  `metz.xls`): independent pKi panel; original vs cleaned-1421 version is
  reported if distinguishable, otherwise recorded as unverified.
* Klaeger (`.../klaeger_matrix.csv`): Kinobeads apparent pKd; different
  experimental semantics, replication only.
* KIBA (`dataset/raw/dta/kiba.tab`): score units, never combined with pK MSE
  and never a sole positive.
* KLIFS annotations (`protein_annotation/klifs_kinase_information_human.json`):
  aligned 85-residue pocket and kinase groups.
* Anastassiadis: the local acquisition manifest records it as
  `excluded_by_governance` and the files are absent; recorded unavailable.
* HIVdb / resistance / ortholog mutation panels: asset search is performed; if
  absent, W0-P is recorded **NOT RUNNABLE**, which forces W1 GO = NO for any
  biological-null interpretation.

## 1. Frozen audit rules

### A. Asset inventory
Every file listed above: path, SHA-256, row/column shape, label semantics as
parsed, acquisition manifest contents.

### B. Censoring / detection-limit audit
For each numeric matrix: value histogram at the boundary value(s), fraction
censored per protein and per kinase group. Values equal to the known floor are
flagged `censored` (Davis `Kd=10000 nM`; Metz values equal to the matrix floor
after reading the xls sheet; Klaeger `pKd=5.0` apparent floor). Censored rows
are excluded from any Pearson/MSE-like statistic in later stages and reported
separately.

### C. Estimand hierarchy (broad -> strict), per dataset independently
For each target/protein with >=2 measured ligands:
1. `all_pairs`: every unordered within-target ligand pair;
2. `similar_pairs`: pairs with Tanimoto(ECFP4, radius 2, 1024 bits) >= 0.6;
3. `mmp`: pairs that share at least one single-cut MMP core and R group
   (Stage U `mmp.transformation`);
4. `strict_mmp`: Stage U exact key (shared core + full attachment context +
   isomeric R_a/R_b).

For each layer report: targets, protein components (CD-HIT40 where sequences
exist, KLIFS kinase group where pocket annotations exist), within-target
ligand pairs, pair classes repeated across >=3 targets and >=3 components,
cross-component double-difference row counts, effective independent units,
Murcko-scaffold novelty between fit-style and holdout-style strata, censored
pair fraction, and exact/coarse transformation counts.

### D. Interaction variance and noise
Same-layer within-key between-target mean square is computed where repeated
classes exist. Repeated-measure noise is only reported where the raw asset has
repeated rows for the same `(protein, ligand)`; otherwise it is recorded
**not identifiable from this asset**. Cross-platform residual reproducibility
is computed only on classes present in two independent platforms and is
interpreted solely as the cross-platform transfer gate, never as a kill of
single-platform signal.

### E. W0-P positive-control availability
Required: known resistance/gatekeeper mutation panels, HIV resistance variants
or same-experiment ortholog/point-mutant panels. Local filesystem search is
performed and recorded. If no admissible panel exists, W0-P is `NOT_RUNNABLE`
and the overall W0b decision is **NO-GO for W1 biological interpretation**,
independent of any support statistic.

### F. Stage W audit
Stage W files are listed with mtimes and preregistration hashes. Existing W0
and W1 artifacts are retained unmodified. W1 is marked **PAUSED**; no training
metric was read (verified by absence of any `runs/` artifact).

## 2. Frozen GO/NO-GO rules (support-side; not Pearson thresholds)

For a dataset to be `GO` for a W1 identification screen, all must hold:
* at least one estimand layer has >= 500 cross-component repeated-class D
  rows over >= 10 protein components;
* the chosen layer's censored-pair fraction <= 0.25;
* W0-P is runnable **or** the dataset is used only for descriptive
  statistics, never as evidence for biological absence;
* no single class or target carries > 50% of the chosen layer's rows.

Otherwise the dataset is `NO-GO` for that layer. `NO-GO` is a statement about
support/censoring/positive-control availability, not biology.

## 3. Artifacts

`W0B_ASSETS.json`, `W0B_CENSORING.json`, `W0B_HIERARCHY.json`,
`W0B_VARIANCE.json`, `W0B_POSITIVE_CONTROL_AUDIT.json`, `W0B_STAGEW_AUDIT.json`,
`W0B_DECISION.json`, `REPORT.md`, `commands.jsonl`, `tests/`.
Preregistration SHA-256 is recorded in every artifact.
