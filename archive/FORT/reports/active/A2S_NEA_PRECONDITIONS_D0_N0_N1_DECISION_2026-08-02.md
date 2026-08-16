# A2S NEA preconditions — gates D0, N0, N1

Date: 2026-08-02 · Branch: `research/a2s-transfer-object-20260802`
Runner: `research/a2s_nea_preconditions.py` (`main.py a2s-nea-preconditions`)
Tests: `tests/test_a2s_nea_preconditions.py` (9 passed)
Artifacts: `reports/active/a2s_nea_preconditions_2026-08-02.json`, `..._records_...parquet`
Roles: source **`fit` only**. `probe` deliberately excluded (consumed once by PIRS).
`locked` and the recipient roster never requested. Deterministic basis. No gradient
model trained.

**FACT — decision: `NO_CHEMICAL_ADAPTATION_OBJECT_SURVIVES_SEPARATION_STOP_PROGRAMME`.**

| Gate | Question | Result |
|---|---|---|
| **D0** | does chemical headroom survive simultaneous separation? | **FAIL** |
| N0 | which nuisance group acts? | measured — **affine**, offsets explain 68 % of residual variance |
| N1 | is there a k ≤ 5 deployment path? | **PASS** (85.9 % coverage at k=5) |

---

## 1. D0 — the decisive measurement

The **same 93 fit targets / 92 homology components**, evaluated twice. The only
difference is what the split separates.

| | scaffold-only | **separated** (scaffold + document + assay) |
|---|---:|---:|
| document overlap, support→query | 88.6 % | **0.0 %** |
| assay overlap | 87.0 % | **0.0 %** |
| **document-mean oracle** (chemistry-free) | **+0.0671** [+0.0504, +0.0854] | **+0.0000** [0, 0] |
| own head, all query pairs | **+0.0610** [+0.0448, +0.0784] | **+0.0044** [−0.0161, +0.0242] |
| own head, same-document pairs | **+0.0275** [+0.0125, +0.0435] | **+0.0123** [−0.0041, +0.0281] |

**The harness validated itself.** On the separated split the document-mean oracle
is *structurally* powerless — every evaluation document is unseen, so it predicts a
constant — and it measures **exactly zero**. That is the strongest available
evidence that the split is what it claims to be, and it is a structural check, not
a statistical one.

**The result.** Under scaffold-only separation the chemistry head scores +0.0610 —
and a chemistry-free document oracle scores **more** (+0.0671). Under simultaneous
separation the head's all-pair gain collapses by **93 %**, to +0.0044 with an
interval crossing zero. The same-document gain falls from +0.0275 to +0.0123, also
crossing zero. **D0 fails: the lower bound (−0.0041) is far below the 0.005 MDE.**

**The collapse is not because the separated task is harder.** Base concordance is
essentially unchanged between regimes (0.5499 → 0.5424 all-pair; 0.5037 → 0.5157
same-document), and the separated split actually has **five times more**
within-document pairs (9 290 vs 1 891). The base performs the same; it is the
*head's advantage* that evaporates.

## 2. What this does and does not license

**It does license:** the statement that on open ChEMBL pKi, the per-target
"adaptation object" that nine mechanisms competed to estimate is, under provenance
separation, **not measurable at this power**.

**It does not license:** "the object is zero." The point estimate on same-document
pairs is +0.0123 and consistently positive across both regimes. With SE = 0.0082 at
92 components, resolving an effect of that size to a lower bound above 0.005 would
require **≈ 445 components — about 4.8× the present corpus.**

That is the honest terminal statement, and it is more useful than either
over-claim: *the effect, if real, is ~+0.012 in target-macro CI and needs roughly
five times the current homology-component count to detect.*

## 3. N0 — the nuisance group, measured

Per measurement context (documents and assays give near-identical answers), over
1 421 documents / 1 403 assays with ≥ 5 rows:

| Quantity | Documents | Assays |
|---|---:|---:|
| **variance of residuals explained by per-context offsets** | **68.1 %** | **68.6 %** |
| offset SD | 1.679 pKi | 1.697 pKi |
| within-context dispersion (median) | 0.967 | 0.956 |
| SD of log within-context dispersion | 0.402 | 0.404 |

**FACT.** Per-context offsets explain **more than two-thirds of all residual
variance**, with an offset SD of ~1.7 pKi — larger than most target-specific
chemical effects anyone is trying to measure. The scale term is also material
(log-scale SD 0.40 ≫ 0.15), so the acting nuisance group is the **full affine
group**, not offset-only. NEA would have had to pay the extra degree of freedom.

This is the quantitative explanation for the whole failure sequence: any estimator
consuming absolute residuals spends its budget on a nuisance that is twice the size
of the signal.

## 4. N1 — coverage was never the problem

Passive support draws, fraction containing a context with enough supports to form
an invariant:

| k | ≥ 2 in one document (`G0`) | ≥ 3 (`G` affine) | mean within-context pairs |
|---:|---:|---:|---:|
| 1 | 0.0 % | 0.0 % | 0.00 |
| 3 | 63.9 % | 24.8 % | 1.14 |
| 5 | **85.9 %** | 54.6 % | 3.79 |
| 10 | 97.5 % | 84.7 % | 16.8 |

NEA had a deployment path. k=1 coverage is exactly 0 %, confirming the structural
no-op rather than assuming it. **The mechanism was buildable; there was nothing for
it to learn.**

## 5. Decision

1. **Stop the mechanism track.** Per the registered stop rule in
   `A2S_STAGE1_MECHANISM_SEARCH_V2`, D0 failing ends it. NEA is not implemented.
   Gates N2–N8 are not run.
2. **Nothing promoted** to `model/` or `script/`. No breakthrough.
3. **The terminal deliverable is the negative result**, now with three measured
   components: the document oracle beating the chemical head under conventional
   splits; the 93 % collapse under provenance separation; and the 68 % variance
   attribution to per-context offsets.
4. **`probe` and `locked` remain sealed** — this ran on `fit` only, so a future
   larger-corpus test still has untouched confirmation roles.

## 6. What would reopen this

Not a new architecture. Only more independent components:

- **≈ 445 homology components** with provenance-separated within-target splits, vs
  92 available here. Papyrus 05.7 and BindingDB are the realistic routes; SPD 2023
  is query-depth underpowered and cannot contribute.
- The same protocol is directly portable: the split builder, the structural
  document-oracle validity check, and the same-document estimand.

Until then, the scientifically defensible claim is the measurement, not a mechanism.
