# Stage R3/R4 preregistration: the two core innovations

Written before any training result exists. Everything runs under the governed
double-cold protocol `bindingdb_ki_double_cold_v1`
(assignment `9d8c7289c1b6162f0e39c0c7ff2222bb45305fe3193bfee7b0b8214c0baf5684`):
`meta_train` 5,643 cells / 346 targets, development `meta_val` 41 targets over
19 components with **zero** exact-ligand, scaffold, component and document
overlap and 81.6% of ligands below Tanimoto 0.40. `meta_test` (22 targets, 10
components) is **untouched** and is not read in this stage.

## Innovation A — level-shape factorized predictor

`model/level_shape.py`:

    f(P, L) = ligand_prior(L) + target_level(P) + centered_interaction(P, L)
    centered_interaction(P, L) = s(P, L) - mean_m s(P, anchor_m)

`anchor_1..M` are learned ligand-side embeddings owned by the model, so the
subtracted term is a per-protein constant computed from parameters, never from
the episode. The branch therefore has exactly zero mean in the anchor basis and
cannot express a target-level constant, while prediction stays inductive.

19 structural gates in `tests/test_level_shape.py` pass, including exact anchor
centering, level-branch constancy across queries, protein-blind ligand prior,
per-query independence of the rest of the batch, exact `k=0` identity, support
permutation invariance, query equivariance, label-locked residuals, and full
gradient coverage.

## Innovation B — counterfactual level-shape gradient-routed training

`scripts/train_level_shape.py`. The squared error decomposes exactly into
`mean(p-y)^2 + var(p-y)`; the method routes each term to the components that own
it while keeping one joint scalar prediction:

    p_level = ligand_prior + target_level + centered.detach()
    p_shape = ligand_prior + target_level.detach() + centered
    L_route = mean(p_level - y)^2 + var(p_shape - y)     ==  MSE(p)  numerically

plus three counterfactual contrasts computed in the same step, each routed so it
cannot be satisfied by the wrong module: protein-shape (ligand prior detached),
protein-level (interaction detached), support-binding (endpoint detached).
Wrong-protein donors are the **most similar** training target from a different
homology component under Stage R2's `esm_whitened` metric. One backward pass,
one optimizer step, single stage.

## Arms — exactly one causal variable per step

| arm | architecture | training | isolates |
|---|---|---|---|
| **A0** | incumbent `similarity_only` | ordinary | the protocol change alone |
| **A1** | level-shape factorized | ordinary shared-gradient | **Innovation A** |
| **A2** | level-shape factorized | + gradient routing | **routing** |
| **A3** | level-shape factorized | + counterfactual = full | **counterfactual** |

Matched seed, matched step budget, matched query-size range. A1 is the
"level-shape decomposition without gradient routing" control; A2 is "routing
without counterfactual training".

Evaluation-time controls on every arm: wrong-protein (similarity-matched),
ligand-only (interaction branch removed), permuted support, level-only
transport.

## Gates

**Zero-shot (k=0)**

| # | gate |
|---|---|
| Z1 | A3 k=0 MSE at least 10% below A0 |
| Z2 | positive component-level bootstrap lower bound against A0 |
| Z3 | CI and Spearman do not regress |
| Z4 | improvement present in every seed |
| Z5 | correct protein beats the similarity-matched wrong protein, component lower bound > 0 |
| Z6 | holds on the `< 0.40` Tanimoto tier (which is 81.6% of this population, so this is nearly the whole test) |

**Few-shot (k >= 1)**

| # | gate |
|---|---|
| F1 | k=1 correction is query-specific, not a scalar level shift |
| F2 | k>=2 full beats the fixed Tanimoto transport's level-only reduction |
| F3 | correct support binding beats permuted support |
| F4 | MSE and ranking improve together |
| F5 | gains are not attributable only to the zero-shot endpoint |

**Training innovation** — the decisive set for Innovation B

| # | gate |
|---|---|
| T1 | A2 beats A1 — routing contributes beyond the factorized architecture |
| T2 | A3 beats A2 — counterfactual contributes beyond routing |
| T3 | both survive three seeds and a component bootstrap |
| T4 | the training method is a **major** source, not a cosmetic regularizer: the A1 -> A3 gap is at least half the A0 -> A3 gap |

## Decision rule

If T1 or T2 fails, the corresponding half of Innovation B is rejected and
reported as rejected. If Z1/Z2 fail, Innovation A is not admitted as a
performance contribution regardless of its structural gates. Nothing is rescued
by changing other modules; a failed arm is reported, not retuned.

`meta_test` is opened exactly once, at the end, only for arms that have already
passed on `meta_val` with three seeds. Whatever it shows is reported.
