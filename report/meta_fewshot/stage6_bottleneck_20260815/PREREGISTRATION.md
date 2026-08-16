# Stage 6 preregistration: chemistry-grounded support weighting

Written before any Stage 6 arm produced data.

## Bottleneck re-diagnosis

Stage 5 concluded "the objective is the blocker". `SIGNAL_meta_val.json` and
`SIGNAL_meta_train.json` **falsify that as stated**. With no learning at all,
Morgan/Tanimoto similarity-weighted support labels beat the support mean:

| split | k | support mean | Tanimoto-softmax | oracle best support | level ceiling |
|---|---:|---:|---:|---:|---:|
| meta_val | 2 | 1.163 | **0.976** | 0.570 | 0.615 |
| meta_val | 3 | 1.096 | **0.887** | 0.357 | 0.615 |
| meta_val | 5 | 1.002 | **0.747** | 0.186 | 0.615 |
| meta_train | 2 | 1.597 | **1.534** | 0.898 | 0.828 |
| meta_train | 3 | 1.317 | **1.202** | 0.488 | 0.828 |
| meta_train | 5 | 1.165 | **0.949** | 0.215 | 0.828 |

Within-target Pearson correlation between Tanimoto similarity and absolute
affinity gap: -0.35 (meta_val), -0.27 (meta_train). SAR continuity is present.

So the MSE optimum is **not** the support mean, and a query-specific weighting
worth 0.19-0.28 MSE at k>=2 exists and is reachable without any protein
information. Every learned transport in this project nonetheless collapsed to
uniform weights. **The bottleneck is the similarity representation used for
support weighting**, not the training objective.

Hard limit of this hypothesis: at k=1 all weightings are degenerate. H1 makes
**no k=1 claim**. Improving k=0 and k=1 needs a different change.

## Hypotheses selected (two, kept separable)

* **H1 (this stage).** Add Tanimoto similarity over Morgan fingerprints as an
  explicit additive term in the support-weighting logits, with a learned scale
  initialised at the value the label-only audit used. Targets k>=2.
* **H2 (only if H1 passes and k=0 remains weak).** Chemistry-grounded features
  in the ligand tower to lift the zero-shot endpoint, whose spread is
  0.087-0.196 pK against a 0.93 pK label spread. Targets k=0.

Not selected, with reasons: Cartesian/equivariant encoders (0/17,717 cells have
a common-frame complex); the signed relative operator (Stage 5, inert);
hypernetworks and support-conditioned trunk modulation (no evidence yet that the
trunk is the binding constraint); noise-aware abstention (does not address a
measured failure); cross-modal pretraining (would confound architecture gains
with extra-data gains, which the contract asks to keep separable).

## Arms

Matched seed 20260812, 800 steps, 4 episodes/step, hidden 384 / embed 192 /
48 contact types / 5 ligand layers, cosine schedule, lr 6e-4, backbone scale
1.0. Identical trunk, encoders and zero-shot endpoint; only the transport
differs.

| arm | `--arch` | transport weighting |
|---|---|---|
| A | `grammar` | learned key + bounded `rho` gate (reused from Stage 5, same seed and budget) |
| E | `similarity` | learned key + `gamma * Tanimoto` |
| F | `similarity_only` | `gamma * Tanimoto` alone |

## Gates, `meta_val`, all eligible episodes

Arm A comparator: `full` minus `level_only` = +0.105 / +0.067 / +0.017 and
permutation gap = +0.163 / +0.110 / +0.081 at k = 2 / 3 / 5.

| id | requirement |
|---|---|
| G1 | `full` beats `level_only` at k in {2,3,5} by more than arm A's margin |
| G2 | permutation gap exceeds arm A's at k in {2,3,5} |
| G3 | **CI and Spearman of `full` are at or above `level_only` at k in {2,3,5}** |
| G4 | k=0 within 0.05 of arm A (the trunk is unchanged; larger drift means the transport is disturbing trunk training) |
| G5 | `full` MSE not worse than arm A at any k in {1,2,3,5} |

G3 is decisive, as in Stage 5: an MSE gain bought by shrinking toward the level
is the pathology already rejected twice.

Advance to a multi-seed run only if G1, G2 and G3 all pass. If arm F matches
arm E, prefer F: the learned key would then be adding nothing and the simpler
transport is the honest one.
