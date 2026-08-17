# Stage L2: the directional-SAR probe, rebuilt fairly

Authority: `STAGE_L2_meta_val.json`. Supersedes `STAGE_L_LIGAND_SAR.md` and its
`STAGE_L_meta_val.json` entirely. Frozen A0 features, one linear directional
head trained on `meta_train` components with weight decay chosen on
`meta_train` component folds, `meta_val` scored once, three fixed pair-sampling
seeds, component-paired bootstrap.

## What was wrong with the first version

Four defects, every one of which flattered the learned arm:

| defect | effect |
|---|---|
| `i` and `j` sampled independently, so orientations were unbalanced | a **symmetric** score picked up an incidental non-zero correlation with the signed gap, and was then compared against a directional arm as a competitor |
| a rank-4/8/16 factorised head `⟨w, U Δe⟩` | algebraically identical to `⟨Uᵀw, Δe⟩` — one linear functional. The "rank selection" was selecting over a reparameterisation of the same hypothesis |
| residualisation slope fitted on `meta_val` using `meta_val` labels | the incremental-value figure was partly fitted to the evaluation split |
| the cliff claim | reported as `embed +0.379` vs `tanimoto −0.370` — the Tanimoto half is meaningless (see below) and the `embed` half **reverses** under fair pairing |

## The construction, and what it guarantees

Every unordered within-target pair contributes **both** orientations `(i,j)` and
`(j,i)`. The signed target is then exactly antisymmetric, so **any symmetric
predictor has identically zero correlation with it — by construction, not by
measurement.** Verified: `tanimoto_signed_structural_zero` measures
**−0.0000 [−0.0000, +0.0000]**.

This is why "Tanimoto scores +0.028 on signed Δy" was a misleading comparison
and why the claim that Tanimoto "points in the wrong direction on cliffs" is
**withdrawn**. A symmetric similarity has no direction. It is scored here only
where it is meaningful — the gap magnitude — and it wins there:

| predicting `|Δy|` | Δ-r | 95% CI |
|---|---:|---|
| **`tanimoto`** | **+0.2989** | [+0.1947, +0.3957] |
| `embed` | +0.2035 | [+0.0851, +0.3353] |

## Result: the directional signal survives, and is real

All arms are the same hypothesis class — one linear functional over a
difference of the same width — so capacity is matched by design.

| arm | Δ-r (seed 20260901) | 95% CI | across 3 sampling seeds |
|---|---:|---|---|
| **`embed`** | **+0.2264** | **[+0.0960, +0.3713]** | **+0.2119 ± 0.0112** |
| `ligand_encoder` (protein-blind) | +0.1404 | [−0.0275, +0.2995] | +0.1191 ± 0.0235 |
| `random_feature` | +0.0538 | [−0.0033, +0.1190] | +0.0541 ± 0.0073 |
| `morgan_difference` | −0.0143 | [−0.1708, +0.1394] | −0.0421 ± 0.0265 |
| `shuffled_labels` | −0.1407 | [−0.2694, −0.0290] | −0.0043 ± 0.1129 |
| `tanimoto` (signed) | −0.0000 | [−0.0000, +0.0000] | structural zero |

`embed` is resolved, stable across sampling seeds (sd 0.011), and beats the
capacity-matched **Morgan fingerprint difference** — which is antisymmetric and
therefore *can* carry direction, and carries none. So the signal is not simply
Morgan bits in another basis.

Incremental value over the protein-blind ligand encoder, with the slope fitted
on `meta_train` and applied frozen: **+0.1880 [+0.0523, +0.3251]** — resolved.
The contact channels add directional information the ligand encoder alone does
not have.

`random_feature` at a consistent +0.054 is a floor worth naming: a random
directional head is not exactly zero on this task. `embed` is about four times
it.

## The two findings that overturn the previous report

### 1. On activity cliffs the signal is *negative*

| stratum | pairs | targets | components | `embed` Δ-r | status |
|---|---:|---:|---:|---:|---|
| **activity cliffs** (Tanimoto ≥ 0.6, gap ≥ 1.0 pK) | 188 | 18 | 10 | **−0.1178** | confirmatory |
| low novelty | 1276 | 28 | 16 | +0.2396 | confirmatory |
| mid novelty | 1288 | 38 | 19 | +0.2448 | confirmatory |
| **high novelty** | 1258 | 30 | 15 | **+0.0523** | confirmatory |

The previous report claimed `embed +0.3793` on cliffs and called it the
direction's best stratum. Under balanced pairing it is **−0.1178**: on the
pairs where highly similar ligands have large potency gaps, the learned
direction points the *wrong way* more often than not.

The stratum clears the preregistered adequacy floor (18 targets, 10 components
against a floor of 8 and 5), so this is a confirmatory reading, not exploratory.

### 2. The signal is concentrated in low-novelty ligands

+0.240 and +0.245 in the two least-novel terciles against **+0.052** in the most
novel one. The double-cold split guarantees **zero shared ligand identities**
with `meta_train` (measured: 0 of 538) and zero shared Murcko scaffolds, so this
is not leakage — it is that the direction generalises by chemical resemblance
and degrades as resemblance falls.

## What this means for the pairwise-operator hypothesis

Stage L's original result motivated a pairwise episodic operator, on the
argument that the signal is pairwise and the A2 moment averages it away. That
argument is unchanged in form and **materially weaker in substance**:

* the directional signal is real, resolved and stable — this survives;
* but it fails exactly where a directional mechanism would earn its keep. Its
  purpose would be to order ligands that similarity cannot separate: activity
  cliffs, where it is now measured at −0.118, and novel chemistry, where it is
  +0.052;
* where it *is* strong is the low-novelty regime, which is also where fixed
  Morgan/Tanimoto transport already performs — and Stage R measured Tanimoto
  transport beating every learned operator at every k.

**A pairwise operator built on this direction would be strongest where the
existing comparator is already strong and weakest where it is needed.** That is
not a promising basis for a model innovation, and Phase 2's conditional gate —
"if the repaired directional-SAR probe remains positive" — should be read
against the strata, not against the headline.

## Claim boundary

1. Ligand-side only. No protein conditioning is tested or claimed.
2. Not a performance result: a frozen-feature pair correlation, not a trained
   DTA model, not comparable to any k=0/k=5 MSE on record.
3. Not meta-learning: no support label enters the predictor.
4. Not an activity-cliff result in any positive sense — the measured cliff
   correlation is negative.
5. `meta_test` labels were used for no fitting, selection or reported metric;
   the process-isolation incident remains open.

## Commands

```bash
conda run -n drug python -m tools.research.a2_exact_probe.stage_l2_directional_sar \
  --features tools/research/a2_exact_probe/features \
  --output tools/research/a2_exact_probe/STAGE_L2_meta_val.json
conda run -n drug python -m pytest tools/research/a2_exact_probe/tests -q
```
