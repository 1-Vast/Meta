# Final boundary audit — verification report (no training)

Generated 2026-08-17 by `tools/research/stageN_audit/final_audit.py`. **This file is generated; do not hand-edit it.** It re-derives the load-bearing numbers of report/BOUNDARY_20260817_NIGHT.md from the recorded artifacts, with no training and no meta_test access. Authority: FINAL_BOUNDARY_AUDIT.json.

## Verified numbers

| claim | recomputed | match |
|---|---|---|
| T2 k=0 MSE / level^2 / centered | 2.5961 / 1.7314 / 0.8648 | exact |
| k=0 level share | 66.7% | exact |
| oracle-level k=0 MSE (centered term alone) | 0.8648 | arithmetic |
| within-document level transfer R^2 | +0.4515 (210 targets) | exact |
| K2 pooled k=0 MSE contrast | -0.1118 [-0.1851, -0.0490] | bitwise |
| K2 pooled k=1..5 MSE contrasts | -0.0480 / -0.0273 / -0.0218 / -0.0122 | bitwise |

## meta_test seal — two properties, reported separately

- **Logical exclusion after parsing**: 106 RESULT.json artifacts audited, 0 evaluated, 2 recording `included=True` (the two disclosed legacy R14 artifacts).
- **Physical isolation**: the governed split view is built at `dataset\processed\meta_fewshot\bindingdb_ki_double_cold_v1_views`; its manifest records `meta_test_label_artifact_emitted=False`. 0 recorded artifacts were produced on it — the isolated surface is available to future runs and is not a retroactive relabelling of past ones.
- **meta_test evaluations: 0.**

## Retained trained stages (discovered from the filesystem)

11 retained trained stages, 11 of them preregistered before their results existed. Preregistered but not run: stageR_daviskiba.

Discovery rule: a directory under tools/research/ is a retained trained stage when it carries PREREGISTRATION*.md and at least one evaluation row artifact (*.rows.jsonl).

This supersedes the earlier counts: the 2026-08-18 audit hard-coded 7 stages and the completion inventory hard-coded 8; both lists omitted trained stages that record rows without a *.rows.summary.json sidecar.

| stage | preregistered | prereg before results | row artifacts | reports |
|---|---|---|---:|---|
| stageA_innerloop | yes | yes | 1 | REPORT.md |
| stageB_complementary | yes | yes | 3 | REPORT.md |
| stageD_level_panel | yes | yes | 4 | REPORT.md |
| stageF_pairwise | yes | yes | 2 | REPORT.md |
| stageG_esm650 | yes | yes | 7 | REPORT.md, REPORT_G2.md |
| stageI_lm | yes | NO | 2 | REPORT.md |
| stageJ_assay | yes | yes | 3 | REPORT.md |
| stageK_contrastive | yes | yes | 4 | REPORT.md |
| stageL_gated | yes | yes | 1 | REPORT.md |
| stageP_cpc | yes | yes | 1 | REPORT.md |
| stageQ_frozenhead | yes | yes | 2 | REPORT.md |

### Preregistration-ordering exceptions

**stageI_lm** — at least one evaluation row artifact is older on disk than the stage's preregistration file. This is a disclosed ordering finding, not a corrected one: the artifacts and their mtimes are left exactly as recorded. The check is by file mtime and is weak in both directions (a later edit to the preregistration moves its mtime forward; a restored file loses its original), so it is reported as evidence to inspect rather than as a verdict. The earlier audit asserted only that a preregistration existed *alongside* the results, which is why this was not visible before.

## Completion-inventory consistency

COMPLETION_INVENTORY.json lists 11 stages; this audit discovers 11. Agreement: **True**.

## What this establishes

The final bounded conclusion is reproducible from the raw evaluation rows: the level/shape decomposition, the assay-history transfer measurement, and the three-seed pooled contrast of the strongest mechanism (K-REG) all re-derive exactly. No recorded artifact evaluated meta_test.

What it does **not** establish: this is arithmetic and inventory verification of development measurements. It converts no empirical model failure into an information-theoretic bound, and it makes no claim about untested architectures or other datasets.
