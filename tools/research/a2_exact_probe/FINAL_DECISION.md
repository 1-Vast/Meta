# A2 family: consolidated final decision (2026-08-16)

Supersedes `tools/research/a2_readiness/` (v1) and the A2 verdict in
`tools/research/a2_readiness_v2/` (v2). This is the single high-density record
of what was tested, what was concluded, and which of those conclusions
survived. Machine authorities: `A2_EXACT_meta_val.json`,
`STAGE_L_meta_val.json`, `../a2_readiness_v2/*.json`.

## Decision

**A2 is closed.** Not by inference from an adjacent probe — by running the
exact operator `NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md` §2 specifies, on real
episodes, with a learned coordinate system, a rank chosen on held-out
`meta_train` components, and its own preregistered gates. It fails five of six
with resolved paired intervals, and two of them fail *inverted*.

**One positive result stands, with a material correction.** A
protein-independent directional head on `embed` predicts the **signed**
within-target affinity difference at Δ-r **+0.212 ± 0.011** across three pair
sampling seeds (+0.226 [+0.096, +0.371] on the first), on held-out protein
components. It beats the capacity-matched Morgan-difference control (−0.042),
the protein-blind ligand encoder (+0.119), and a random directional head
(+0.054), and its incremental value over the ligand encoder is resolved
(+0.188 [+0.052, +0.325], slope fitted on `meta_train`).

**But the previous cycle's headline claim about it is withdrawn.** Rebuilt with
balanced `(i,j)`+`(j,i)` pairs — which makes a symmetric score's signed
correlation *identically* zero by construction, verified at −0.0000 — the
activity-cliff stratum reverses from the reported **+0.379** to **−0.118**
(188 pairs, 18 targets, 10 components; adequate, confirmatory). The claim that
Tanimoto "points the wrong way" on cliffs is also withdrawn: a symmetric
similarity has no direction, and the earlier figure was an artifact of
unbalanced sampling. The signal is concentrated in the two least-novel ligand
terciles (+0.240, +0.245) and nearly vanishes in the most novel (+0.052).

## Three cycles, and what each actually established

| cycle | claim | status |
|---|---|---|
| **v1** ordering is interaction-borne, not ligand-borne | +0.186 [+0.057, +0.324] | ✅ **holds**, reconfirmed under five donor strata |
| **v1** the interaction branch's ordering is protein-inert | level 0.215 pK vs centered 0.0007 pK | ✅ **holds**, strengthened |
| **v1** the architecture can express what training removes | one random init shifted 110× more | ❌ **withdrawn** (v2): ten inits move in uncorrelated directions, cosine −0.003 over 1845 pairs |
| **v1** the collapse is in the readout | attention JS 0.241 vs 0.218 | ❌ **relocated** (v2): causal interventions put it at the fusion/pooling inside `ContactGrammar` |
| **v2** A2's premise is falsified | a zero-shot bilinear delta predictor | ❌ **over-reach** — that is not A2's operator. Superseded by Stage R |
| **v2** a protein-independent SAR direction exists in `embed` | Δ-r +0.262 | ✅ **holds**, reproduced at +0.270 and now compared against Tanimoto |
| **R** the exact episodic A2 operator | 6 preregistered gates | ❌ **A2 closed**, 5/6 fail resolved, 2 inverted |
| **L** the `embed` direction beats/complements Tanimoto | orthogonality and strata | ✅ **complements**, does not replace |

Two of v1's four load-bearing claims and v2's central verdict did not survive.
The pattern in both failures was the same: a measurement on one object
generalised to a different one. v1 measured the endpoint scalar and concluded
about the representation; v2 measured a pair predictor and concluded about an
episodic operator. Stage R exists because of that pattern.

## Stage R — the exact operator (`STAGE_R_EXACT_A2.md`)

```text
z_i = A_φ(e0)   c_S = (1/k) Σ r_i z_i   η(k)=η_∞k/(k+λ)   δ_q = η(k)⟨c_S, z_q⟩
```

19 structural gates pass, including the one that made the test worth running:
**the k=1 correction is query-specific**, which A0's provably is not
(`sar_adaptation ≡ 0`).

k=5 MSE (pK²), `meta_val`, equal-component weighting; k=0 is 2.1369 for every
correct-protein arm against A0's recorded 2.1488:

| arm | k=5 MSE | paired contrast vs `a2_embed` |
|---|---:|---|
| `tanimoto` (parameter-free) | **0.9101** | a2 −0.265, **resolved worse** |
| `scalar_level` (2 scalars) | **1.0746** | a2 −0.103, **resolved worse** |
| `a2_label_shuffled` | 1.0820 | a2 −0.092, **resolved worse — inverted** |
| `a2_wrong_protein` | 1.0866 | a2 −0.037, **resolved worse — inverted** |
| `shared_moment` (no support chemistry) | 1.1185 | a2 −0.098, resolved worse |
| **`a2_embed`** | **1.1765** | — |
| `a2_random_projection` | 1.4386 | a2 **+2.011, PASS** |

**Corrupting the protein or the support labels makes A2 better.** The only gate
it passes is against a frozen random projection — so the optimiser works and the
coordinate system is being learned; what it learns is worse than two scalars.

The mechanism is measured, not inferred. `query_spread_pk`, the standard
deviation of `δ` across an episode's queries, is **0.0028 pK** for `a2_embed`
against a label spread of 0.884 pK — 0.3%. Given Gaussian noise features the
same operator produces **0.3497 pK**. The operator is capable of large
query-specific corrections; trained on the real representation, gradient descent
drives it to a constant. It degenerates into a worse-calibrated level shift.

## Stage L2 — the directional result, rebuilt fairly (`STAGE_L2_DIRECTIONAL_SAR.md`)

The first version is deleted, not amended: four defects each flattered the
learned arm — unbalanced pair orientations, an algebraically collapsing
rank-4/8/16 head that was really one linear functional, a residualisation slope
fitted on `meta_val` labels, and a cliff comparison against a symmetric score
that cannot express direction at all.

Rebuilt with balanced pairs, one explicit `Linear(D, 1)` head, weight decay
chosen on `meta_train` folds, `meta_train`-fitted frozen residualisation slopes,
and three pair-sampling seeds:

| arm | Δ-r, 3 sampling seeds |
|---|---|
| **`embed`** | **+0.2119 ± 0.0112** |
| `ligand_encoder` (protein-blind) | +0.1191 ± 0.0235 |
| `random_feature` | +0.0541 ± 0.0073 |
| `morgan_difference` (antisymmetric, capacity-matched) | −0.0421 ± 0.0265 |
| `shuffled_labels` | −0.0043 ± 0.1129 |
| `tanimoto`, signed | **−0.0000** — structural zero, by construction |

`embed` beats the Morgan-difference control, which *can* carry direction and
carries none, so the signal is not Morgan bits in another basis. Incremental
value over the ligand encoder is resolved: **+0.188 [+0.052, +0.325]**. Where
Tanimoto is meaningful — the gap *magnitude* — it wins: +0.299 vs `embed`'s
+0.204.

### The two strata that change the recommendation

| stratum | pairs | targets | components | `embed` Δ-r |
|---|---:|---:|---:|---:|
| **activity cliffs** | 188 | 18 | 10 | **−0.1178** |
| low novelty | 1276 | 28 | 16 | +0.2396 |
| mid novelty | 1288 | 38 | 19 | +0.2448 |
| **high novelty** | 1258 | 30 | 15 | **+0.0523** |

Zero shared ligand identities with `meta_train` (0 of 538) and zero shared
scaffolds, so this is generalisation by chemical resemblance rather than
leakage — and it degrades as resemblance falls.

**The directional signal is strongest where fixed Morgan/Tanimoto transport is
already strong, and fails where a directional mechanism would earn its keep.**
Stage R measured Tanimoto transport beating every learned operator at every k.
A pairwise operator built on this direction would therefore be competing with
the comparator on the comparator's home ground while losing on cliffs and novel
chemistry.

## What is closed, and at what boundary

**Closed:** the protein-conditioned SAR moment family, on the double-cold
protocol, sequence + 2D inputs, frozen A0 representations (`embed`,
`max_state`, `ligand`), ranks 4/8/16 selected on held-out components, k ∈
{0,1,2,3,5}. The operator does not beat two trainable scalars, does not
approach a parameter-free kernel, and improves under both of its falsification
controls.

**Not closed by this evidence:**

1. whether the *architecture* can learn protein-conditioned ordering when an
   objective explicitly demands it. Every wrong-protein control in R0-R14 is
   uncentered and therefore level-only (F7), so no objective has ever asked.
   This is Stage P, and it is untouched by Stage R;
2. whether a **pairwise** learned operator over `(query, support)` could use the
   Stage L direction, where the moment form cannot;
3. anything about `meta_test`, whose labels were used for no fitting, selection or reported metric, and whose process-isolation incident remains open.

## Claims that remain prohibited

* no pocket, contact, binding-site or "biologically localized" language — the
  protein path is exactly invariant to residue-slot permutation (2.4e-08 pK);
* no SOTA or excellence claim — every number is `meta_val` development evidence;
* no protein-specificity claim from an uncentered control;
* no meta-learning claim for the Stage L direction: no support label enters it,
  and the protein contributes +0.017 over a permuted-protein control;
* no performance claim for Stage L at all — it is a frozen-feature pair
  correlation, not a trained DTA model, and Stage R is the only deployment test
  run so far.

## Retained artifacts

| path | why |
|---|---|
| `a2_exact_probe/{operator,run_probe,stage_l_ligand_sar,extract_ligand_features}.py` | the superseding implementation |
| `a2_exact_probe/tests/` | 19 structural gates |
| `a2_exact_probe/{A2_EXACT,STAGE_L}_meta_val.json` | machine authorities |
| `a2_exact_probe/{STAGE_R_EXACT_A2,STAGE_L_LIGAND_SAR}.md` | decision reports |
| `a2_readiness_v2/` | the governance incident, the noise/leakage audit, the causal attention audit, and Stage P's frozen design |
| `a2_readiness/LITERATURE_LEDGER.md`, `tests/` | the literature review and the CPC centering probes, both still valid |
| `*.meta.json`, `FEATURE_CACHE_MANIFEST.json` | provenance for the deleted 50 MB of regenerable caches |
