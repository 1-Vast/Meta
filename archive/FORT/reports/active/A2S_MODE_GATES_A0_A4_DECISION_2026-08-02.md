# A2S-MODE Gates A0–A4 — decision

Date: 2026-08-02
Artifacts: `reports/active/a2s_mode_gates_2026-08-02.json`,
`reports/active/a2s_mode_gates_records_2026-08-02.parquet`,
`research/a2s_mode_gates.py`, `tests/test_a2s_mode_gates.py`
Proposal: `A2S_MODE_MECHANISM_PROPOSAL_2026-08-02.md`
Roles opened: `fit`, `probe`. `locked` and the A2S recipient roster never requested.
Device: RTX 4060 Laptop GPU, `D:\anaconda\envs\drug`.

**Decision: `PREMISE_CONFIRMED; K_SHOT_INFERENCE_NOT_YET_IDENTIFIABLE_WITH_AN_UNSHAPED_DICTIONARY`.**

> **SUPERSEDED IN PART — read `A2S_MODE_GENERALIZATION_DECISION_2026-08-02.md` first.**
> Gates G1–G4 subsequently established that (i) §2's stratification axis was wrong — A0 stratified by
> **support→query** similarity, which a per-target head never sees, so the claim "headroom outside the
> local-analogue regime" was **not** established here and is withdrawn; and (ii) the recommendation in
> §8 to implement a separation-shaped dictionary is **retracted**: source-target heads have a nearly
> flat spectrum and low-rank projection onto the dominant source directions destroys the gain, so
> there is no small shared mode structure for shaping to sharpen. §4's A2/A4 failure is explained by
> that absence, not by the estimator. The A0/A1 *numbers* below stand; the interpretation does not.

The gates did what they were built to do: they confirmed the part of A2S-MODE that was speculative
and localised the entire remaining burden onto the one component the proposal named as its innovation
carrier.

---

## 1. Setup

Compact label-free ligand basis `g(x)` — 10 standardised descriptors + 16 Morgan principal components,
all statistics from `fit` rows only, `d = 26`. Per-target heads fitted by ridge on `g` against
frozen-base residuals. Dictionary = k-means over 110 `fit`-target heads (targets with ≥ 40 rows),
sweeping `M ∈ {2,3,4,6}` plus a null mode. Measured `σ = 0.976`, level SD `= 1.898`.
4,029 probe episodes, two policies, `k ∈ {1,3,5}`, 26–76 homology components per cell, paired
component bootstrap, MDE 0.005. Mode ceilings are selected **on a held-out query half**, never with
hindsight on the scored queries.

## 2. Gate A0 — is there anything to select? **PASS**

Per-target head (fitted on that target's abundant labels, evaluated on held-out queries) minus frozen
base, target-macro CI at k=5:

| stratum | components | ΔCI | 95 % interval |
|---|---:|---:|---|
| `t < 0.20` | 31 | **+0.0517** | [+0.0154, +0.0848] |
| `0.20–0.35` | 37 | **+0.0774** | [+0.0451, +0.1131] |
| `0.35–0.55` | 45 | +0.0672 | [+0.0437, +0.0914] |
| `≥ 0.55` | 49 | +0.0832 | [+0.0548, +0.1119] |
| pooled | 54 | +0.0847 | [+0.0624, +0.1080] |

**FACT.** Target-specific ranking structure recoverable by a **query-only** function on a 26-dimensional
label-free basis exists in **every** stratum, including the two where fixed Tanimoto KRR, TRACE and
every other transport operator measure exactly zero.

**INFERENCE.** This is the first positive ranking headroom this programme has measured outside the
local-analogue regime, and it validates the §1.3 reframing: the reachable object is a target-conditioned
function of the query, not a weighted sum of support residuals. It is an *oracle* (it uses the
target's abundant labels), so it is a ceiling, not a method.

## 3. Gate A1 — mode sufficiency **PASS** (M ≥ 3)

Split-half-selected mode minus the **single global head** (the "it is just a better global model"
confound), CI:

| k | stratum | components | ΔCI | 95 % interval |
|---|---|---:|---:|---|
| 3 | `t < 0.20` | 35 | **+0.0192** | [+0.0020, +0.0378] |
| 3 | pooled | 76 | **+0.0271** | [+0.0146, +0.0405] |
| 5 | `≥ 0.55` | 53 | **+0.0227** | [+0.0045, +0.0410] |
| 5 | pooled | 73 | **+0.0301** | [+0.0176, +0.0448] |

**FACT.** `M = 2` fails; `M ∈ {3,4,6}` pass. **INFERENCE.** A small discrete set of response modes is
sufficient to carry a material part of the A0 headroom, and it is not explained by a better global
ligand model — the global head on its own is at or below the frozen base.

## 4. Gate A2 — k-shot identifiability **FAIL**, and Gate A4 explains why

k-shot mode selection minus global head, and selection accuracy against chance `1/(M+1) = 0.20`
(headline `M = 4`):

| k | stratum | accuracy | ΔCI | 95 % interval |
|---|---|---:|---:|---|
| 3 | pooled | 0.259 | +0.0082 | [+0.0006, +0.0160] |
| 5 | `≥ 0.55` | 0.225 | +0.0196 | [+0.0036, +0.0373] |
| 5 | pooled | 0.264 | +0.0136 | [+0.0021, +0.0265] |
| 5 | `t < 0.20` | 0.232 | −0.0042 | [−0.0445, +0.0275] |

Every lower bound sits below the 0.005 MDE, and accuracy is barely above chance everywhere.
**The verdict also flips with the nuisance choice `M`** (`M = 3` passes A2, `M ∈ {2,4,6}` fail), which
is itself a reason not to accept a positive here.

**Gate A4 — synthetic positive control: also FAIL, and this is the decisive measurement.** In a world
where each probe target is *literally generated* from one dictionary mode plus the measured noise,
k-shot selection accuracy reaches only 0.32 at k=5 (chance 0.20) and every gain interval crosses zero.

> **INFERENCE — the localisation.** The A2 failure is **not** evidence that targets lack a discrete
> state. A0 and A1 say they have one. It is a measurement that a **k-means dictionary is not
> separable from `k ≤ 5` noisy residuals**, because k-means minimises within-cluster variance, which
> is unrelated to — and can be antagonistic to — between-mode discriminability
> `D_k = Σ_{i≤k}(ψ_m(x_i) − ψ_{m'}(x_i) − Δ̄)²/σ²` on a random k-subset.
>
> This is exactly the prediction of §3 of the proposal: `ρ_k` is fixed by the data, but `D_k` is an
> object the outer loop must **maximise**, and nothing in this gate maximised it.

## 5. What this means for the mechanism

| proposal component | status after A0–A4 |
|---|---|
| Premise: a query-only target-conditioned correction has headroom outside the transport stratum | **confirmed** (A0) |
| Premise: a *small discrete* state is sufficient, and is not a better global head | **confirmed** (A1, M ≥ 3) |
| The state is inferable from `k ≤ 5` labels **without shaping** | **refuted** (A2 + A4) |
| The separation regulariser `L_separation` is load-bearing | **now the entire claim** |

**INFERENCE.** A2S-MODE's headline claim is therefore sharper than the proposal stated, and it has a
free, exact null baseline:

> **Meta-training a mode dictionary to maximise worst-case pairwise discriminability at the
> deployment budget makes a discrete target state k-shot-identifiable, where an unshaped
> (k-means / likelihood-clustered) dictionary of the same size and capacity is not.**

The unshaped dictionary measured here is the exact nested restriction that the shaped model must
beat — the same role `a2s_bir_global` played for IDA and `R2b_krr` played for TRACE. The claim is
exactly the measured delta, and Gate A4 becomes the mechanism's primary sanity rung rather than an
afterthought: **shaping must first make the synthetic control recoverable.**

## 6. Registered predictions for the implementation phase

| # | Prediction | Falsifies |
|---|---|---|
| Q1 | Shaping raises realised `D_k` (report it per run) and raises synthetic-control selection accuracy well above 0.32 at k=5 | the whole shaping claim |
| Q2 | Shaped real-data k-shot gain over the **unshaped dictionary of the same M** has a paired component LCB > 0.005 at k=3 and k=5 | the mechanism |
| Q3 | The gain does **not** flip sign or significance across `M ∈ {3,4,6}` | robustness; an `M`-dependent result is a nuisance artefact |
| Q4 | k=1 rank action stays ≈ 0 (design decision confirmed by the user) | leakage alarm |
| Q5 | P0 holds: gain flat in support→query Tanimoto, rising in support diversity | the claim is transport in disguise |
| Q6 | Shaped gain survives in `t < 0.35` (Gate A3), where transport is null | the reason to prefer this route over a better KRR |

## 7. Honest limits of this gate

- A0 and A1 ceilings use query labels (per-target oracle; split-half selection). They are ceilings,
  not methods, and are labelled as such throughout.
- Stratum cells hold 26–53 components; MDE80 there is ≈ 0.008–0.010, so single-cell readings at
  `t < 0.20` are weaker evidence than the pooled cells.
- The dictionary comes from 110 `fit` targets with ≥ 40 rows. A different row threshold changes the
  head-estimation noise and therefore the clustering.
- `probe` is a development role and has now been inspected across Q1, Q2 and these gates.
  Confirmation still requires freezing a protocol and opening `locked` once; recipient labels stay
  sealed.

## 8. Recommendation

Proceed to implement A2S-MODE **with the claim narrowed to the shaping objective**, and with these
non-negotiables carried from this gate:

1. The unshaped k-means dictionary at matched `M`, capacity and basis is the **nested baseline**.
2. Realised `D_k` and synthetic-control recovery are **reported every run**; shaping that does not
   move them is a failed run regardless of the real-data number.
3. Sweep `M ∈ {3,4,6}` and require sign stability, since A2's verdict flipped on `M` here.
4. Keep the k=1 rank silence as the built-in leakage detector.
