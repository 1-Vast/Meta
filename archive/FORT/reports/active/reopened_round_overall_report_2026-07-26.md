# Reopened-round overall report

## Terminal category

**② `SIGNAL_PRESENT_EVIDENCE_INSUFFICIENT`.**

All three newly proposed candidates failed at least one frozen gate. Thresholds, documents, folds,
seeds, and endpoints were not changed after results. The round nevertheless found two pieces of
positive evidence that materially refine the root cause.

## Candidate sources and results

### C1 WTPAIR -- changed supervision

Source idea: the user's mutation/double-difference direction plus pairwise relative-effect methods.
WTPAIR directly trained same-KLIFS-group protein-pair x ligand-pair mixed differences using the
existing frozen ESM and Morgan coordinates. It used 256 bilinear coefficients, one seed, and a
matched-capacity cellwise ridge.

Result: 0.0323 component-macro Spearman, below ligand-only 0.0429, matched cellwise 0.0419, and
group centroid 0.0911. It did not beat random protein. This excludes insufficient bilinear
capacity and ordinary loss choice as explanations. The combination of pooled full-sequence ESM
with a whole-molecule ligand basis does not expose transferable within-group SAR.

### C2 CAPIT/ASPIRE P0 -- changed protein information source

Source idea: KLIFS aligned active-site representations and the user's CAPIT synthesis. The
low-cost gate used only the aligned 85-residue pocket and an intentionally favorable shared-panel
oracle before any fragment tensor or neural model.

Result: aligned pocket 0.4539, group centroid 0.4253, pooled ESM 0.3983, shuffled pocket 0.3562,
random target 0.3670. Pocket minus group is +0.0286 [+0.0116,+0.0459], narrowly below the frozen
+0.030 materiality threshold; all destruction controls pass. Aligned active-site composition
therefore carries real target-specific information and is a better coordinate system than pooled
ESM, but the oracle reads held-ligand measurements on training targets and cannot establish strict
ligand-cold prediction. CAPIT training remained unauthorized.

### C3 CROSSDOC -- changed data condition and estimand

Source idea: assemble locally validated panel-relative ranking, endpoint isolation, and document
firewalls rather than require one paper to solve the full task. Exact KIRHub compounds and targets
were compared only with three ChEMBL-37 train documents, keeping pKi and pKd separate.

Result: group-residual reordering correlation +0.4946 [+0.3156,+0.6727]; target permutation
p=0.001499; ligand-permutation null +0.0163 and observed-minus-null +0.4783
[+0.3031,+0.6524]. This is strong independent-source mechanism agreement. Coverage is only 13
target-document units / 11 homology components, below 30/25, so the external replication gate
fails and no model is authorized.

## Explanations excluded

- Ligand-only potency cannot explain the C2 pocket-shuffle destruction or C3 residual agreement.
- Coarse target taxonomy cannot explain C2 pocket over group or C3 leave-group-out residual
  correlation and within-group target permutation.
- Extra parameters cannot explain C1 because the matched-capacity cellwise model is stronger.
- A favorable seed cannot explain the data-only C2/C3 results; no multi-seed model claim is made.
- Endpoint mixing is excluded in C3: pKi and pKd are never pooled on their raw scales.
- Target/ligand identity is not counted as predictive generalization: C3 is explicitly a
  replication estimand, and C2 is explicitly a transductive oracle.

## Most credible root cause

The biological information exists but is poorly matched to the original representation and public
sampling geometry. Kinase selectivity is localized in aligned active-site residue combinations and
their interaction with ligand fragments. Pooled ESM distances and whole-molecule Morgan
neighbourhoods blur that coordinate system. Public independent continuous-affinity documents
contain confirmatory directionality but too few exact factorial overlaps after target homology,
compound, endpoint, document, and saturation controls.

The bottleneck is now specific: **a powered, document-isolated, cross-kinase panel with shared or
compositional-fragment-covered ligands**, not another posterior, support kernel, larger protein
encoder, or longer training run.

## Protocol status and reopening condition

A schema audit displayed five labeled rows from the existing ChEMBL confirmation partition.
Although unused in all candidate calculations, the partition is quarantined and may not serve as
future confirmation. A new confirmation source or newly sealed split is required.

The next exploration round may start only after this report is frozen. It should first search for
public panels that increase independent target-document components and aligned
active-site x ligand-fragment coverage. If no such source exists, the next valuable outcome is a
prospective measurement design; it is not justified to train CAPIT on the current overlap.

Candidate ledger: 3/3. Single-seed trainable runs: 1. Multi-seed runs: 0. Confirmation runs: 0.
Sealed test consumed: false. Existing ChEMBL confirmation partition usable: false.
