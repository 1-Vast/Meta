# Post-completion review: evidence, limits, and required repairs

Date: 2026-08-18. Status: current interpretation authority.

This review was written after the long-run goal was marked complete on its
second terminal condition. It does not replace leaf `RESULT.json` artifacts or
their numerical reports. It controls how the closing authorities are
interpreted where their language is broader than the recorded experiments.

## 1. Decision

The research cycle is scientifically complete as a **scope-bounded negative
result**, not as a proof that cold-target DTA or every possible model is
incapable of reaching MSE <= 1.00 pK^2.

The supported conclusion is:

> No tested candidate reached MSE <= 1.00 pK^2 at every k under the governed
> BindingDB-Ki double-cold development protocol, sequence and 2D-ligand inputs,
> single-stage differentiable training, and the locally tested legal feature
> families. The dominant measured obstacle at k=0 is target/panel level error,
> much of which tracks assay and document history that the split makes
> unavailable across components.

The completion was conditional at repository level when this review was
written. **The conditions were discharged by the 2026-08-17/18 repair** (see
§6): the physical meta-test seal is implemented and mountable, both audits now
discover the same 11 retained trained stages from the filesystem, the
`GOAL_ACTIVE.md` update is committed, and meta_val checkpoint-selection reuse
is eliminated from the maintained trainer.

## 2. Load-bearing measurements accepted

| Evidence | Recorded result | Supported interpretation |
|---|---:|---|
| T2 k=0 decomposition | MSE 2.5961 = level^2 1.7314 + centered 0.8648 | The current k=0 error is 66.7% target/panel level. |
| Within-document transfer | R^2 +0.4515 over 210 targets | Assay or medicinal-chemistry programme history carries substantial level information inside a document. |
| Tested legal feature probes | best joint measured R^2 +0.259 | The tested transferable inputs explain about 25.9% of level variance; this is not an information-theoretic ceiling. |
| K-REG, three seeds | all-k pooled MSE deltas below zero | A training objective can improve calibration consistently, but its centered gain was not confirmed. |
| Best few-shot band | k=5 MSE 0.944-1.007 | The few-shot target is near the requested boundary; k=0 remains the principal failure. |
| Occupancy representation | held-out ordering r about +0.218 | A small transferable within-target signal exists and is not fully used by the endpoint. |
| Meta-test artifact audit | 104 artifacts, 0 evaluated | Recorded metrics did not evaluate meta-test, although labels were parsed before logical exclusion. |

The bitwise audit verifies stored arithmetic and artifact consistency. It does
not convert empirical model failures into universal causal or impossibility
proofs.

## 3. Claims that must be narrowed

### 3.1 Predictable level fraction

`+0.259` is the largest result among the measured governed probes. Replace
phrases such as "at most 26% is predictable" with "the tested probes explain
up to 25.9%." An unseen representation or newly governed input could explain
more.

### 3.2 The MSE <= 1.00 boundary

The decomposition proves that the target is arithmetically possible: an oracle
level gives k=0 MSE about 0.865. The reported `0.1239` level budget assumes the
current centered error remains fixed. A future model could improve both terms.
The defensible statement is that no tested legal-input candidate reached the
target, not that the target is mathematically impossible.

### 3.3 Level/ranking conflict

Stages E, J, L, and Q repeatedly show a level/ranking trade-off in the tested
shared-trunk, single-stage compositions. This is strong empirical evidence, but
not a theorem about every single-stage end-to-end architecture. Replace
"fundamental" with "reproduced across four tested compositions."

### 3.4 Method-family closure

`tools/research/method_ladder/CLOSURE_MAP.md` maps several named methods to
nearby proxies rather than direct implementations. Orthogonal routing is not a
full OGM or Disentangled Gradient Learning experiment; mean/max panel pooling
is not a Set Transformer; pairwise transport is not a full DrugBAN trunk; and
paired alignment is not FS-CAP. Mark these entries `proxy negative; direct
method not instantiated` instead of `closed by measurement`.

### 3.5 Seal terminology

Meta-test is logically excluded after the all-label artifact is parsed in
**every recorded artifact**; those artifacts are not physically isolated and
must never be relabelled. Zero recorded evaluations make the existing
development metrics usable.

*Discharged 2026-08-17/18:* the goal's physical-seal requirement is now met for
future runs. `QPSMPData(split_view=…)` mounts separately hashed per-split label
artifacts with the meta_test artifact out of tree, and a file-access spy proves
that a train/meta_val workflow opens neither the sealed path nor the all-label
`cells.jsonl.gz`. The two surfaces must always be named separately: say
"logical exclusion after parsing" of the record, and "physically isolated
(governed split view)" only of a process that actually mounted the view.

## 4. Governance and consistency defects (all resolved 2026-08-17/18)

1. `tools/research/GOAL_ACTIVE.md` contained an uncommitted round-19 completion
   update. **Resolved**: preserved verbatim and committed with this repair.
2. `stageN_audit/AUDIT_REPORT.md` reported seven preregistered training stages;
   `COMPLETION_INVENTORY.json` reported eight. **Resolved**: both were stale
   hard-coded lists selecting on `*.rows.summary.json`, which Stages A, B and
   P_cpc never emit. Both now call one filesystem discovery rule and agree on
   **11** retained trained stages, plus `stageR_daviskiba` preregistered and
   not run.
3. Closing documents reported both 147 and 151 research tests. **Resolved**:
   147 = `RUN_SLOW=1 pytest tools/research/stageA_innerloop
   tools/research/stageB_complementary`; 151 = that pair plus
   `tools/tests/test_research_record.py` without `RUN_SLOW` (151 passed / 12
   skipped); 135 = the pair alone without `RUN_SLOW`. All were subsets. The
   authoritative complete research-suite command is now
   `RUN_SLOW=1 pytest tools/research -q` (255 passed / 2 skipped). Neither
   historical number was deleted.
4. Two legacy R14 artifacts record `included=true` but `evaluated=false`. They
   remain disclosed incidents and must not be summarized as a physical seal.
   **Unchanged by design** — `audit_research_record` still exits non-zero while
   they stand.
5. `report/README.md` and the authority date in `CURRENT_MODEL_EVIDENCE.md`
   were stale before this review. **Resolved.**
6. *(Found during the repair.)* The maintained trainer selected checkpoints on
   `meta_val` — the leak Stage B measured at ~0.62 pK^2 at k=0. **Resolved**:
   the leak-free rule is promoted to `scripts/internal_validation.py` and is
   the default.
7. *(Found during the repair.)* 98 Python bytecode caches were tracked in Git
   because `!tools/**` re-admitted them. **Resolved**: re-excluded and removed
   from the index; the files remain on disk.

## 5. Repository authority after repair

The intended compact hierarchy is:

1. This review: final interpretation and repair requirements.
2. `FINAL_STATE_20260818.md`: cycle summary after wording repairs.
3. `BOUNDARY_20260817_NIGHT.md`: numerical boundary after wording repairs.
4. `EVIDENCE_LEDGER.md`: stage decisions and leaf locations.
5. `stageN_audit/`: regenerated arithmetic, seal, preregistration, and inventory
   audit.

`COMPLETION_STATEMENT_20260818.md` records why the goal was closed; it must not
override the narrower scientific scope in this review. Earlier A2 handoffs and
R-series plans are historical inputs, not active instructions.

## 6. Required repair acceptance criteria — verified status

Repair executed 2026-08-17/18. Every criterion below is marked with what was
actually checked, not with an intention.

| # | criterion | status | evidence |
|---|---|---|---|
| 1 | All closing authorities use measured/probe-qualified wording | **DONE** | BOUNDARY, FINAL_STATE, CURRENT_MODEL_EVIDENCE, EVIDENCE_LEDGER, COMPLETION_STATEMENT, task.md, history.md, GOAL_ACTIVE.md rewritten; pinned by `test_no_authority_turns_a_measured_failure_into_a_bound` and `test_closing_authorities_record_the_oracle_level_floor` |
| 2 | Proxy method families distinguished from directly tested methods | **DONE** | `CLOSURE_MAP.md` now classifies each family direct/partial/proxy/not-instantiated: 0 direct, 3 partial, 5 proxy; pinned by `test_the_closure_map_separates_proxies_from_direct_implementations` |
| 3 | meta_test labels in a physically separate artifact; ordinary processes cannot parse them | **DONE** | `QPSMPData(split_view=...)` mounts per-split artifacts and never opens `cells.jsonl.gz`; the sealed artifact is out of tree; proved by the file-access spy in `tools/tests/test_physical_meta_test_seal.py` (29 tests), with a negative control showing the spy fires on the default surface |
| 4 | Final audit covers every retained trained stage and reconciles the counts | **DONE** | both audits now share `stage_inventory.discover_stages`; **11** retained trained stages (the 7 and 8 were stale hard-coded lists that also omitted A, B, P_cpc); `FINAL_BOUNDARY_AUDIT.json`, `AUDIT_REPORT.md`, `COMPLETION_INVENTORY.json` regenerated |
| 5 | One maintained-suite and one research-suite result recorded consistently | **DONE** | `python main.py verify tests` 310 passed / 6 skipped; `RUN_SLOW=1 pytest tools/research -q` 255 passed / 2 skipped; 147/151/135 identified and retained as subset commands |
| 6 | Authorities point to the same hierarchy and scoped conclusion | **DONE** | §5 hierarchy referenced from `report/README.md`, `task.md`, `CURRENT_MODEL_EVIDENCE.md`, `EVIDENCE_LEDGER.md`, `FINAL_STATE`, `COMPLETION_STATEMENT` |
| 7 | Worktree clean after a focused documentation/governance commit | **DONE** | bytecode caches untracked and gitignored (98 `.pyc` removed from the index, files left on disk) so the tree stays clean across test runs |

Added beyond the original list, at the user's instruction:

| # | criterion | status | evidence |
|---|---|---|---|
| 8 | Eliminate meta_val checkpoint-selection reuse | **DONE** | Stage B's leak-free rule promoted to `scripts/internal_validation.py`; `TrainConfig.selection` defaults to `"internal"`; `--selection` CLI on the maintained trainer; every run now records a `checkpoint_selection` provenance block; Stage B imports the promoted module so the two cannot drift |

### Findings raised by the repair itself

1. **stageI_lm preregistration ordering.** The regenerated audit checks
   preregistration *ordering*, not merely presence, and finds that
   `IFROZEN_meta_val.rows.jsonl` (06:45) predates `PREREGISTRATION.md`
   (07:30) on disk while the candidate arm's rows (07:50) postdate it. The
   mtime check is weak in both directions and is reported as evidence to
   inspect, not as a verdict. Nothing was edited to hide it.
2. **All recorded artifacts remain on the weaker surface.** The physical seal
   is available to future runs. It is not retroactive, and no stored
   `RESULT.json` may be relabelled to claim it.
3. **Every recorded figure carries meta_val-selection optimism** (~0.62 pK^2
   at k=0, Stage B). The leak-free rule is now the default, so the recorded
   band and any future band are not directly comparable.

## 7. Next scientific decision

Davis is the highest-information next experiment, but it is a new boundary
replication goal rather than a continuation of the completed BindingDB search.
Before reading labels, its endpoint types and units must be governed separately;
Kd, Ki, IC50, functional response, and composite scores must not be treated as a
single interchangeable pK scale. KIBA acquisition remains out of scope.

An MSA lane or wider structure coverage would also open a new input regime and
therefore a new goal. Restating the objective as centered MSE, ranking, or
per-document calibration changes the scientific task and requires explicit
authorization.
