# Governance audit of `tools/research/stageU_mmp_interaction/`

**Decision: Stage U is STOPPED and superseded by `stageV_core_mmp`. Its
preregistration is not edited, its thresholds are not retrofitted, and its
frozen text is inherited verbatim where Stage V reuses it.**

## 1. Timeline (filesystem mtimes, `dataset`-local clock)

| time | artifact |
|---|---|
| 16:30:07 | Stage T `PREREGISTRATION.md` frozen |
| 16:39:29 | Stage T `T1_CENSUS.json` (admission statistics read) |
| **17:08:44** | **Stage T `T2_RESULT.json` (all gate metrics read)** |
| 17:10:58 | Stage T `REPORT.md` |
| **17:12:53** | **Stage U `PREREGISTRATION.md` frozen** |
| 17:13 | Stage U `mmp.py`, `provenance.py`, `__init__.py` |
| 17:14 | Stage U `observations.py`, `u0_reliability.py`, `U0_PROVENANCE…gz` |
| 17:15:05 | Stage U `U0_RELIABILITY.json` |
| 17:23 | Stage U `u0_census.py`, `observation_cache.py` |
| 17:25:06 | Stage U `U0_OBSERVATIONS.jsonl.gz` (cache) |
| 17:26:20 | Stage U `u0_census.log` — **0 bytes** |

`PREREGISTRATION.md` SHA-256
`fdc0a830aa92882d07b9aea50f22a4c72fc6d93f92c55a3be6bc15cd6a645c11`.

## 2. What the timeline establishes

* **Stage U is an adaptive correction informed by Stage T, not an independent
  confirmation.** Its preregistration was frozen 4 minutes after Stage T's full
  gate result was read. Its own text says so ("Stage T … its exact
  transformation key omitted the shared core"). This is disclosed, not
  concealed, and it is the correct classification for every claim Stage U or its
  successor makes.
* **No Stage U gate metric has been read.** `U0_CENSUS.json` does not exist,
  `u0_census.log` is empty, `runs/` and `tests/` are empty. The only Stage U
  output is `U0_RELIABILITY.json`, a supervision-reliability audit that gates
  nothing. So Stage U stopped before its own first admission statistic.

## 3. What Stage U got right, and Stage V inherits verbatim

* **Core-inclusive exact key**: `sha256(core_isomeric | repr(attachment
  context) | R_a >> R_b)` (`mmp.py:188`), with hybridization added to the
  attachment environment. Verified present in code, not only in prose.
* **Descriptor consumes the core**: `descriptor()` returns
  `core_counts + R_a + R_b + delta + context`, and `edit_features()` appends
  folded Morgan fingerprints of core, `R_a` and `R_b`.
* **An interaction-variance gate (U1) before any neural training** — exactly the
  "stop before training if interaction variance does not exceed supervision
  noise" requirement.
* **Local protein-region operator** with ordered residue tokens, cross-attention,
  and an explicit prohibition on pooled-protein / target-ID bypass.
* **Structural `D_hat = R(tau,p1) - R(tau,p2)`**, giving identity, protein-pair
  antisymmetry and cycle consistency for free.
* **Paired substitution as the primary correct-vs-shuffled gate** (§4.6 gate 3),
  with the separately trained shuffled arm only as a secondary check (gate 5).
* Domination caps on key/target/component observation share (§2.5.6).
* Three-seed confirmation requirement (§4.6).

Stage V reuses **all** of Stage U's numeric thresholds unchanged, because they
were frozen at 17:12 *before any core-inclusive census number existed*. That is
the property that keeps them non-retrofitted, and it is why they are inherited
rather than restated.

## 4. Load-bearing requirements MISSING from Stage U

These are the reason Stage U is stopped rather than continued. Each was required
by the active instruction and is absent from the frozen text.

1. **Residue-token permutation control.** Absent entirely (`grep` for
   `residue.perm` returns nothing). This is not a nice-to-have for Stage U
   specifically: its whole premise is *ordered* region tokens with sinusoidal
   slot encoding. Without a permutation substitution there is no evidence the
   ordering carries anything, and this repository already carries a standing
   finding that the incumbent protein path was exactly invariant to residue-slot
   permutation.
2. **Capacity-matched random protein representation.** Absent. Without it, a
   gain over `E_mean_tokens` cannot be separated from "any per-target vector of
   this width helps".
3. **A valid target-key shortcut diagnostic.** Stage U's `fit_unsampled` is
   "a frozen 10% sample of fit cross-component `D` rows" (§4.3) and gate 8 rests
   on it. That bank retains the same targets *and* the same transformation keys
   as training, so it cannot detect memorisation of a target key — it is the
   exact construction the instruction rules out. Stage T made the same mistake
   and measured Pearson 0.912 on it against 0.059 out of component.
4. **Identical initialization and identical minibatch order across matched
   arms.** Not stated (§4.5 fixes only the seed, steps, optimizer and schedule).
   Stage T shipped with a per-arm batch seed that gave every arm a different row
   sequence; that defect was caught only by inspection.

Two further defects, recorded for completeness:

5. **Gate 1 is ill-defined as written.** "C_local minus A_zero Pearson >= +0.05"
   is evaluated against a constant-zero predictor, whose Pearson is undefined.
   Stage T hit exactly this and produced `NaN`.
6. **No nested zero predictor in the candidate.** Stage U's arm A is a separate
   constant arm, but `C_local` itself has no way to abstain, so an
   uncontrolled output amplitude can make it lose to the zero predictor on error
   metrics without that meaning anything about protein information. Stage T's
   candidate lost 1.578 vs 0.660 on exactly this axis.

## 5. Rule applied

The instruction is explicit:

> If Stage U's existing preregistration already contains every required
> correction and no relevant gate metric had been read before it was frozen,
> continue it unchanged. If a load-bearing requirement is missing after metrics
> were read, stop Stage U and create a new corrected stage with a new
> preregistration. Do not retrofit thresholds.

Stage T's gate metrics **were** read before Stage U was frozen (17:08 vs 17:12),
and four load-bearing requirements **are** missing. Therefore Stage U is
stopped and `stageV_core_mmp` is created with a new preregistration that
inherits every Stage U threshold verbatim and adds the four missing
requirements plus the two defect repairs.

**No threshold is loosened anywhere.** Every addition is a further control or a
stricter evaluability condition.

## 6. Disposition of Stage U

`tools/research/stageU_mmp_interaction/` is retained read-only as evidence of
the adaptive correction attempt. Its `PREREGISTRATION.md` is **not** edited.
Its `U0_RELIABILITY.json` reproduces Stage T's T0 audit on the same inputs and
is superseded by it for citation purposes. No Stage U number is used as
evidence anywhere.
