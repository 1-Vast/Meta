# A2S-MODE Gates G1–G4 — is the adaptation object generalizable, and is it worth building?

Date: 2026-08-02
Artifacts: `reports/active/a2s_mode_generalization_2026-08-02.json`,
`reports/active/a2s_mode_generalization_records_2026-08-02.parquet`,
`research/a2s_mode_generalization.py`, `tests/test_a2s_mode_generalization.py`
Roles opened: `fit`, `probe`. `locked` and the recipient roster never requested.
Device: RTX 4060 Laptop GPU, `D:\anaconda\envs\drug`.

**Decision: `GENERALIZABLE_BUT_NOT_FEW_SHOT_REACHABLE`.**

The target adaptation object is real and survives a scaffold-cold test. It is **not** low-dimensional,
**not** shared across targets, **not** protein-predictable, and its measured label requirement is
**k ≈ 10–20, not k ≤ 5**. A2S-MODE as proposed is refuted, and so is the low-rank-code family it
belonged to. What replaces them is a sharper and better-posed question with a measured baseline.

---

## 0. Correction to the A0/A1 report

**The earlier claim was mis-stratified and is withdrawn.** Gate A0 stratified the per-target head's
gain by **support→query** Tanimoto. That axis is irrelevant to a head fitted on ~130 of the target's
own rows: the head never sees the support set. The statement "the first ranking headroom outside the
local-analogue regime" was therefore not established by A0.

Re-measured on the correct axis — similarity between the evaluation compound and the **head's own
training rows** — the picture is more modest and more honest (§1). The A0 numbers themselves stand;
only the interpretation was wrong.

## 1. Gate G1 — is the object global? **PASS, with a qualified range**

Per-target head minus frozen base, evaluated on a **within-target Murcko-scaffold-disjoint** split
(the head never sees a scaffold it is scored on), CI:

| similarity of eval compound to head-training rows | components | ΔCI | 95 % interval |
|---|---:|---:|---|
| all | 50 | **+0.0524** | [+0.0293, +0.0746] |
| `< 0.30` | 5 | +0.0385 | [+0.0003, +0.0751] |
| `0.30–0.45` | 9 | +0.0093 | [−0.0540, +0.0708] |
| `0.45–0.60` | 22 | +0.0300 | [−0.0146, +0.0794] |
| `≥ 0.60` | 44 | **+0.0706** | [+0.0425, +0.0986] |

**FACT.** The object survives scaffold-disjointness: a head fitted on one set of Murcko scaffolds
transfers to compounds built on scaffolds it never saw, +0.0524 CI on 50 components.

**FACT.** The gain still concentrates where the evaluation compound retains a ≥ 0.60 Tanimoto
neighbour in the head's training rows. The `< 0.30` cell holds **5 components** and is underpowered —
it is uninformative, not negative.

**INFERENCE.** The reach of a per-target head is *longer* than the reach of `k ≤ 5` support transport
(it crosses scaffolds), but it is not unbounded. "Scaffold-cold" and "chemically distant" are not the
same predicate, and this programme should stop treating them as one.

## 2. Gate G2 — is it low-dimensional and shared? **NO. This is the decisive result.**

Spectrum of the 110 source-target heads (variance fraction): `0.154, 0.109, 0.084, 0.072, 0.069, …`
— top three directions hold **34.7 %** of the variance. There is no dominant shared direction.

Projecting each probe target's **own oracle head** onto the top-`r` **source** subspace, then
evaluating (scaffold-disjoint, all):

| projection | ΔCI vs base | 95 % interval |
|---|---:|---|
| rank 1 | −0.0058 | [−0.0153, +0.0038] |
| rank 2 | −0.0033 | [−0.0138, +0.0064] |
| rank 3 | −0.0007 | [−0.0105, +0.0085] |
| rank 5 | +0.0006 | [−0.0096, +0.0109] |
| **rank 26 (full)** | **+0.0524** | [+0.0293, +0.0746] |

**FACT.** A rank-2 projection retains **−6 %** of the full head's gain. Low-rank projection onto the
dominant *source* directions destroys the signal completely.

> **INFERENCE — the core finding.** The directions that matter for a given target are essentially
> orthogonal to the directions along which source targets vary most. Target response heads are
> **high-dimensional and idiosyncratic**, not draws from a small shared basis.
>
> This refutes, on one measurement, three routes at once:
> * **A2S-MODE** — there is no small set of separable shared response modes to cluster, which is
>   *why* Gates A2/A4 failed. The k-means dictionary was not a bad estimator of a real structure; the
>   structure it assumed does not exist.
> * **A2S-IDA** — the rank-`m` code route. Not merely unidentifiable at `k ≤ 5` (as `ρ₅ ≈ 0.147`
>   already said), but **actively harmful**: rank-restricted heads underperform full-rank heads at
>   *every* label budget tested. There is nothing to shape a low-rank basis *toward*.
> * Any future proposal of the form "meta-learn a compact target state" on this substrate.

## 3. Gate G3 — is it predictable from protein alone? **NO**

Head predicted by ridge from the pooled ESM-2 embedding (fitted on `fit` targets, no probe label
read), evaluated zero-shot: **ΔCI −0.0185 [−0.0725, +0.0190]**; absolute CI 0.5362 against a base of
0.5528.

**FACT.** Protein sequence does not predict a target's ligand-response head. **INFERENCE.**
Consistent with this programme's earlier `TR group not resolvable` result, and it means no
protein-conditioned zero-shot shortcut is available. Whatever adaptation happens must come from the
support labels.

## 4. Gate G4 — the label learning curve. **This is the value answer.**

Empirical-Bayes head from `k` labels inside the top-`r` source subspace (prior
`λ_j = σ²/τ_j²` measured on source targets, never tuned on probe), scaffold-disjoint, vs frozen base:

| k | rank 2 | rank 5 | **rank 26 (full)** |
|---:|---:|---:|---|
| 1 | −0.0067 | −0.0067 | −0.0067 [−0.0165, +0.0022] |
| 3 | −0.0065 | −0.0029 | +0.0014 [−0.0100, +0.0122] |
| 5 | −0.0072 | −0.0046 | +0.0111 [−0.0028, +0.0245] |
| 10 | −0.0053 | +0.0017 | **+0.0261 [+0.0122, +0.0410]** |
| 20 | −0.0030 | +0.0058 | **+0.0432 [+0.0252, +0.0613]** |
| 40 | −0.0031 | +0.0191 | **+0.0517 [+0.0332, +0.0735]** |

**FACT — the knee is at k ≈ 10.** The best cell at `k ≤ 5` has a lower bound of **−0.0028**; the first
cell clearing the 0.005 MDE is `k = 10`. At `k = 40` the estimator reaches +0.0517, i.e. essentially
the full oracle head (+0.0524) — the object is fully recovered by 40 labels.

**FACT — low rank never helps.** Rank 1/2/3 are flat-to-negative at every budget including `k = 40`.
Full rank dominates from `k = 3` upward.

**INFERENCE.** The `k ≤ 5` deployment budget sits measurably *below* the knee of this object's
learning curve. That is not a statement about a particular architecture; it is a property of the
object, measured with a prior-regularised estimator that has no free parameters to blame.

## 5. Direct answer to the two questions asked

**Is the mechanism generalizable?**
Partially, and not in the way it needs to be. It generalizes **within a target across scaffolds**
(G1). It does **not** generalize **across targets as a shared low-dimensional structure** (G2), and
it is **not** recoverable from protein (G3). A2S-MODE's premise — a small dictionary of shared,
separable response modes — is refuted by G2, which also explains the A2/A4 failure structurally.

**Does it hold significant potential and substantial value?**
**Not at k ≤ 5.** The measured requirement is k ≈ 10 for a detectable effect and k ≈ 20–40 for most of
the available gain. There is real value at those budgets — but it is delivered by a **closed-form
empirical-Bayes ridge on a 26-dimensional label-free basis**, with no meta-learning of any kind. Under
this programme's own admissibility rules that is a strong baseline, not a mechanism.

## 6. What is worth doing instead

The measurement replaces a vague goal with a sharp one and hands it a baseline curve:

> **Can meta-learning move the label learning curve left — from k ≈ 10 to k ≈ 5 — on an object that
> has no low-dimensional shared structure?**

This is well-posed, decisive and cheap to falsify, because §4 is exactly the curve to beat and it was
produced by an estimator with no tuned parameters. Three concrete levers, in order of measured
promise:

1. **Shape the basis, not the code.** G2 says the *source* head directions are useless, but the basis
   `g` was fixed (10 descriptors + 16 Morgan PCs, chosen for interpretability, not for
   conditioning). A meta-learned `g` optimised so that heads become *low-rank in the learned basis*
   is the one version of the IDA idea that G2 does not refute — it attacks the representation rather
   than assuming a code exists in a given one. Falsifier: the head spectrum in the learned basis must
   concentrate, and the rank-2 retained fraction must rise materially above the measured −6 %.
2. **Raise the budget, honestly.** If a k ≈ 10–20 protocol is admissible for the intended application,
   the closed-form estimator already delivers +0.026 to +0.043 CI today. That should be reported as
   the programme's current best transferable result and used as the baseline for everything after.
3. **Stop treating "scaffold-cold" as "chemically distant."** G1 shows they behave differently. Every
   future stratification should use similarity to the *estimator's own training rows*, not to the
   support set.

## 7. Honest limits

- 50–54 probe targets clear the 40-row threshold; the `< 0.30` similarity cell has 5 components and
  settles nothing on its own.
- The G4 curve uses random support draws within a target. A diversity-aware draw could shift the knee
  left and was not tested; that is a legitimate and cheap follow-up.
- `oracle_rank*` uses the target's own fitted head, so G2's low-rank verdict is about the *source*
  subspace specifically, not about every conceivable low-rank basis — hence lever 1 above.
- `probe` is a development role and has now been inspected across Q1, Q2, A0–A4 and G1–G4.
  Confirmation still requires freezing a protocol and opening `locked` once. Recipient labels stay
  sealed.
