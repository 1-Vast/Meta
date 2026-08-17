# Governance incident: the meta_test seal was opt-out (2026-08-16)

Status: **repaired and verified**. Numerical impact: **none, demonstrated**, not
assumed. Record impact: **two false claims corrected, one of them mine.**

This is filed as an incident rather than a missing defence because the failure
was not that a guard was absent — it was that the *recorded artifact asserted a
guard that the code did not implement*. An audit trail that states a stronger
property than the code enforces is worse than one that states nothing, because
it terminates further checking.

## 1. What the contract said, and what the code did

`report/meta_fewshot/stageR5_reltransport_20260816/REPORT.md` (row 5) records
the R5 remedy as:

> physical seal: `QPSMPData(include_meta_test=False)`; explicit
> `--include-meta-test`/`--eval-meta-test` opt-in

The class signature implemented the opposite default:

```python
def __init__(self, ..., include_meta_test: bool = True):   # pre-repair
```

So the seal held only for callers that opted *out*. The three trainers
(`train_reltransport.py`, `train_level_shape.py`, `train_grammar_shape.py`)
hardcode `False` and were correct throughout. Analysis scripts written later
took the default.

## 2. Exact scope

### 2a. Scripts that loaded the sealed split

| script | line | wrote a seal claim? |
|---|---|---|
| `scripts/r14_dispersion_audit.py` | 148 | **yes — and it was false** |
| `scripts/stageR6_compare_arms.py` | 323 | no |
| `scripts/stageR9_pair_audit.py` | 113 | no |
| `scripts/stageR2_representation_discriminator.py` | 118 | no |
| `scripts/stageR3_compare_arms.py` | 140 | no |
| `scripts/stageR0_retrieval_falsification.py` | 233 | no (older corpus, no governed split) |

### 2b. Artifacts carrying the literal false claim

The string `"physical: QPSMPData include_meta_test=False"` appears in **37**
recorded JSON artifacts. Classifying by producer rather than by string match:

* **30 are true.** They come from the three trainers, which hardcode the flag.
  The artifact matches the code that wrote it.
* **7 are false.** All are downstream of `r14_dispersion_audit.py`: five it
  wrote directly (`MetaSieve.R14DispersionAudit.v1`) and two stage-level
  `RESULT.json` summaries that copied the string from them.

I wrote `r14_dispersion_audit.py` and its two stage summaries in the previous
cycle. The claim did not match the code I had just written.

## 3. What actually happened to the sealed data

Stated at the precision the evidence supports, and no further.

**Loaded into process, not used in the recorded computation.** For the six
scripts in §2a, `meta_test` cells were decompressed, parsed, retained in
`self.cells`, and indexed into `self.tasks` / `self.components`. No `meta_test`
value entered any recorded metric: every bank is built from
`fixed_nested_episode_banks("meta_val", …)`, every donor is drawn from
`meta_val`, and the label scale is fitted on `meta_train` cells only.

**This is verified, not argued.** `scripts/r14_dispersion_audit.py` was re-run
under the repaired fail-closed seal against the same three frozen A0
checkpoints, the same fixed bank and the same evaluation seed:

```
numeric fields compared : 105
bit-identical           : 105
max absolute difference : 0.000e+00
```

Artifact: `SEAL_REPAIR_REPRODUCTION.json`. If any sealed row had entered the
computation, removing it would have moved a number.

The mechanism behind the exact reproduction is worth recording, because it is
what makes the whole pre-incident `meta_val` record safe. Dropping sealed cells
renumbers `self.cells`, so `EpisodeSpec` support/query integers differ between
the two constructions. But `_unique_ligand_order` permutes an *array of
indices* with a target-seeded generator, and a permutation of an array of
length *n* consumes the same randomness regardless of the values in it. The
same cells are therefore selected in the same order under both constructions —
only their names change. `tools/tests/test_meta_test_seal_contract.py`
asserts this directly on ligand identities and pK values across all five
support sizes.

**What must not be claimed.** The confirmation population was *not* physically
pristine at the process level during those runs, and this document does not
claim it was. The defensible statement is:

> No double-cold `meta_test` label entered any fitting, selection or
> reported metric, and none informed a decision. For six analysis scripts
> the split was parsed **and indexed** in memory. Indexing is now
> impossible; parsing still happens on every load of the legacy loader.
> **Do not write "never opened" or "physically sealed".**

### The repair is not a physical seal, and must not be called one

*(Added 2026-08-16 after review, correcting this document's own first draft.)*

The fix makes the sealed rows **unreachable after construction**. It does not
make them absent from the process. `cells.jsonl.gz` is a single all-label
artifact, so `_governed_cells` decompresses and parses **every** `meta_test`
label on every construction before discarding it. Four properties, routinely
collapsed into one word:

| property | status |
|---|---|
| fail-closed default | ✅ since this repair |
| unreachable after construction (no index, no `materialize`) | ✅ since this repair |
| metric-unconsumed | ✅ verified bit-identical |
| **physically isolated from the read path** | ❌ **false, and was never true** |

The vocabulary is now fixed repository-wide: `seal_record()` emits
`"level": "logical_exclusion_after_parsing"` with
`"physically_isolated": false`, `audit_research_record` prints the isolation
level as its own banner, and
`test_the_seal_does_not_claim_physical_isolation` fails if the stronger word
returns. `SPLIT_ISOLATION_SPEC.md` specifies what a real seal would require —
per-split label artifacts and no `meta_test` shard in the development tree,
the design `scripts/seal_compiled_dataset.py` already implements for the other
corpus lineage. It is **not implemented and not authorized**.

### `violations = 0` does not close this incident

The audit previously reported `0 violation(s)` and exited zero while two
artifacts remained `process_unsealed`. That reads as a clean bill of health for
a record that carries an open incident. `audit_research_record` now prints an
explicit `*** OPEN INCIDENT ***` banner naming the artifacts and **returns a
non-zero status** while any entry remains. The two R14 artifacts keep that
state permanently: the repair prevents recurrence, it does not retroactively
seal a process that already ran.

Two residual exposures remain open by design and are recorded rather than
papered over:

* `_ligand_smiles` and `_protein_sequences` are built from `ligands.jsonl` /
  `proteins.jsonl`, which are split-agnostic. Fingerprints and sequences for
  ligands and targets that appear only in `meta_test` are therefore reachable
  regardless of the seal. **No label crosses**; this is structural identity,
  already public in the corpus manifest.
* The three pre-R5 `SEALED_meta_test_DO_NOT_OPEN.json` sidecars still exist on
  disk. They were computed before the seal was authorised and have never been
  read. The audit tracks them as `sealed_quarantined`, separately from
  `sealed_explicit`, because "metric-unconsumed" and "never computed" are
  different claims and only the first is supported.

## 4. The second defect, found while repairing the first

`scripts/audit_research_record.py` detected an artifact's protocol by looking
at top-level `split_directory` / `split` and `config.split_directory`. The
hand-written stage summaries record it inside a `population` block instead, and
six older stage summaries record no directory at all. Everything that failed
the check was swept into an `older_protocol` bucket whose printed description
is "`bindingdb_ki_main_v0` — a different, consumed population".

**All ten artifacts in that bucket were double-cold stage summaries.** None of
them names `bindingdb_ki_main_v0`. The audit was asserting a consumed protocol
for evidence that is on the sealed one — silence read as a positive claim.

Repaired by scanning the whole document for either corpus name and adding a
`split_undeclared` state for artifacts that name neither.

| audit state | before | after detection | after resolution |
|---|---:|---:|---:|
| `sealed_explicit` | 41 | 42 | **47** |
| `sealed_implicit` (pre-R5) | 15 | 15 | 15 |
| `process_unsealed` (this incident) | — | 2 | 2 |
| `split_undeclared` | — | 7 | **2** |
| `sealed_quarantined` | 3 | 3 | 3 |
| `older_protocol` | 10 | **0** | 0 |
| `violations` | 0 | 0 | 0 |

**The `split_undeclared` artifacts are now resolved from their own evidence.**
Each stage summary aggregates per-seed `RESULT.json` leaves that *do* declare
`config.split_directory` and `split_assignment_sha256`.
`resolve_split_declarations.py` reads the leaves, requires unanimity, and
stamps the parent with a `split_declaration` block naming the leaf files it
inherited from — a derivation, not an assertion. Five stages resolved (R7-R11)
from 3-9 leaf declarations each.

Two are **not** resolvable and are recorded as such rather than guessed:

| artifact | why it stays undeclared |
|---|---|
| `stageR5_reltransport_20260816/RESULT.json` | a contract-and-gates stage: structural gates and a pipeline smoke, no real-data run to inherit a declaration from |
| `stageR13_shape_direct_20260816/RESULT.json` | the family failed its own Stage 1 gates before any real-data run; `population.split` is already explicitly `null` |

Both now carry `population.split_declaration = {"status": "unresolvable", …}`
with the reason, so the audit's remaining count of 2 is a documented terminal
state rather than an unexamined backlog.

`report/BOUNDARY_20260816.md` carried the pre-repair counts and has been
corrected.

## 5. The repair

**Fail-closed default.** `include_meta_test` now defaults to `False`.

**Authorization, not a flag.** `include_meta_test=True` additionally requires a
non-empty `meta_test_authorization` string, recorded in the artifact. A bare
`True` raises; so does an authorization without the flag, rather than the code
silently choosing a meaning. `--meta-test-authorization` is mandatory alongside
`--include-meta-test` in `train_qpsmp.py`, `evaluate_qpsmp.py` and
`evaluate_checkpoint_nested.py`.

**Exclusion precedes every index.** Filtering moved into the same pass that
applies the split assignment, so `self.cells` is never *bound* to a list
containing a sealed cell. Task indices, component maps and `materialize` are
all built from that attribute and therefore cannot reach one. The sealed rows
are still parsed out of the gzip stream — the honest claim is unreachability
after construction, not that the file is never decompressed.

**Artifacts cannot contradict their constructor.** `QPSMPData.seal_record()`
emits the `meta_test` block from the object's actual state. Scripts call it
instead of writing a string. An open run must state `evaluated=` explicitly:
once the split is loaded the object cannot infer what the caller did with it,
and the incident was precisely an unstated claim.

**One legitimate consumer remains.** `build_double_cold_split.py` reads the
whole main_v0 corpus because it *creates* the double-cold assignment; the
double-cold `meta_test` does not exist at that point. It now passes a written
authorization recording that the assignment is label-blind (scaffolds,
fingerprints, panel documents and counts; no `pK` is read).

**Honest vocabulary, enforced.** `seal_record()` emits an `isolation` block
stating `level: logical_exclusion_after_parsing`,
`labels_parsed_in_process: true`, `physically_isolated: false` and why not.
`audit_research_record` prints the isolation level as its own banner and, while
any artifact remains `process_unsealed`, prints an `*** OPEN INCIDENT ***`
banner and **returns a non-zero status** — `violations == 0` no longer reads as
a clean record.

**Contracts.** `tools/tests/test_meta_test_seal_contract.py`, 18 tests:
fail-closed default; a non-zero withheld population; no task index or
materializable cell covers a sealed row; sealed targets absent from every
split; episode draws from `meta_test` raise; bare/whitespace authorization
refused; authorization-without-flag refused; an authorized open really does
expose the split; `seal_record()` truthful in both states; an open run must
state `evaluated`; the reproducibility contract of §3; **and the terminology
guard — `test_the_seal_does_not_claim_physical_isolation` fails if the stronger
word returns to the artifact.** `tools/tests/test_research_record.py` gains
audit-classification contracts.

**Not repaired, by design.** Physical isolation. See §7.

## 6. Verification

```bash
conda run -n drug python -m pytest tools/tests/test_meta_test_seal_contract.py -q
conda run -n drug python -m pytest tools/tests -q
RUN_SLOW=1 RUN_RESEARCH_GATES=1 conda run -n drug python -m pytest tools/tests -q
conda run -n drug python -m scripts.audit_research_record --skip-loading
```

| check | result |
|---|---|
| seal contract | 18 passed |
| maintained suite | 247 passed, 6 skipped |
| record audit | **exit 1** — 0 violations but 2 `process_unsealed` (open incident); 0 older-protocol; 2 split-undeclared with recorded reasons; 12 hashes verified, 0 mismatched |
| R14 reproduction under the repaired seal | 105/105 fields bit-identical |

## 7. What is repaired, and what is deliberately not

| defect | status |
|---|---|
| opt-out default | **repaired** — fail-closed |
| opening by a bare keyword | **repaired** — written authorization required |
| sealed rows reachable via task indices / `materialize` | **repaired** |
| artifacts asserting a seal their code did not implement | **repaired** — derived from the object, guarded by a test |
| `violations = 0` masking an open incident | **repaired** — banner plus non-zero exit |
| ten artifacts misclassified as older-protocol | **repaired** — 0 remain |
| seven artifacts with no declared protocol | **5 resolved from leaf evidence, 2 recorded as unresolvable** |
| future-dated records (2026-08-17) | **repaired** — 28 files corrected to 2026-08-16 |
| **sealed labels parsed on every load** | **NOT repaired.** Requires a new `meta_test`-labelled artifact, which needs authorization. Specified in `SPLIT_ISOLATION_SPEC.md`; recommended trigger is immediately before the confirmation run, not now |

## 8. What changed in the record

* 7 artifacts: `meta_test` block replaced with a truthful state plus a
  `seal_correction` block naming this document. **No numeric field, row or hash
  was modified** — `correct_seal_claims.py` asserts non-seal content is
  unchanged before writing.
* 5 stage summaries: `population.split_declaration` derived from their own leaf
  runs; 2 more recorded as unresolvable with a reason.
* `report/meta_fewshot/stageR14_*/REPORT.md`: correction notes updated from
  "fix planned" to "fixed and verified".
* `report/BOUNDARY_20260816.md`: seal counts corrected; the guarantee restated
  at its true strength for the pre-repair period.
* 28 files: dates corrected from a future 2026-08-17 to the actual 2026-08-16.

## 9. Standing lesson

Twice now the defect has been the *language*, not the mechanism. First a claim
about the seal was written by hand next to code that did not implement it.
Then the repair for that was itself described as a "physical seal" when it is
logical exclusion after parsing — a weaker property that happens to be
sufficient today.

Both are the same failure: a word doing more work than the code behind it. The
structural answer is that the claim is now *derived* (`seal_record()`), *graded*
(`isolation.level`), and *guarded* (a test that fails if the stronger word
returns). `seal_record()` refuses to guess `evaluated` for the same reason.
