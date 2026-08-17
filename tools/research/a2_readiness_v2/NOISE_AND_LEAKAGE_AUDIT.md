# Phase 1: is the protein-inert ordering real? (2026-08-16)

Authority: `BRANCH_ORDERING_V2_meta_val.json` (+ `.rows.jsonl`). Population:
double-cold `meta_val`, 41 targets / 19 components, k=0, query panel 16, the
frozen R6-R14 bank at evaluation seed 73101. Three trained A0 seeds and ten
independently initialised models. No training, forward passes only.
`meta_test` unreachable by construction (repaired fail-closed default).

Everything in `_frozen.py` — checkpoint hashes, split hash, bank parameters,
donor rule, seeds, metrics, bootstrap, and all thresholds — was fixed before
the first measurement and is written into the output.

## Verdict in one line

**The diagnosis survives every control. My v1 interpretation of the
random-initialisation arm does not.**

## 1. The finding that survives

| quantity | value |
|---|---:|
| `r(full endpoint)` | 0.2128 |
| `r(ligand_value only)` | 0.0272 |
| `r(interaction only)` | 0.2206 |
| within-target sd of `ligand_value` | 0.0147 pK |
| within-target sd of `interaction` | 0.1242 pK |
| within-target sd of the labels | 0.8837 pK |
| **`r(full) − r(ligand_only)`** | **+0.1855 [+0.0566, +0.3236] — RESOLVED** |

Within-target ordering is an interaction-branch property. The protein-blind
ligand branch is nearly flat: 0.0147 pK of spread against a 0.8837 pK label
spread.

## 2. Donor choice is not the explanation

v1 used one donor — the nearest legal cross-component protein. Phase 1 uses
five, at fixed quantiles of whitened protein similarity, so the size of the
perturbation is a measured variable.

| donor stratum | whitened cosine | **level shift** | **centered (ordering) shift** | ratio | shift·truth alignment | `r(int) − r(int │ wrong P)` | verdict |
|---|---:|---:|---:|---:|---:|---|---|
| nearest | +0.335 | 0.2150 pK | 0.0007 pK | 307 | −0.014 | −0.0002 | DECISIVE_NULL |
| q25 | +0.047 | 0.2843 pK | 0.0009 pK | 316 | −0.009 | −0.0011 | DECISIVE_NULL |
| median | −0.002 | 0.2948 pK | 0.0009 pK | 328 | +0.021 | +0.0004 | DECISIVE_NULL |
| q75 | −0.036 | 0.2881 pK | 0.0011 pK | 262 | +0.046 | −0.0002 | DECISIVE_NULL |
| farthest | −0.119 | 0.3421 pK | 0.0011 pK | 311 | +0.065 | **+0.0016** | RESOLVED_NEGLIGIBLE |

The donor rule works — the level shift rises monotonically with distance,
0.215 → 0.342 pK — and the ordering shift does not follow it. The 262:1 to
328:1 ratio is stable across the whole range.

The `farthest` contrast is the one place the interval excludes zero, at
**+0.0016 in correlation**. That is 31× below the preregistered smallest effect
of interest (0.05). Calling it a finding would be a significance test standing
in for a magnitude, which is why the verdict vocabulary separates
`RESOLVED_NEGLIGIBLE` from `RESOLVED`.

**A0 reads the protein for the level and not for the ordering, and this is not
an artifact of which wrong protein was chosen.**

## 3. The measurement floor, which v1 never established

| control | measured |
|---|---:|
| identical-protein substitution, centered | **0.00e+00 pK** |
| identical-protein substitution, level | **0.00e+00 pK** |
| repeated forward pass, centered | **0.00e+00 pK** |

Exactly zero — the forward path is deterministic on fixed inputs, so a 0.0007
pK reading is 100% signal and 0% instrument noise.

## 4. Negative controls

| control | A0 | expected |
|---|---:|---|
| shuffled labels, `r(full)` | −0.0082 | ≈ 0 ✅ |
| shuffled labels, `r(interaction)` | −0.0132 | ≈ 0 ✅ |
| foreign ligand panel, `r(full)` | −0.0097 | ≈ 0 ✅ |
| foreign ligand panel, centered protein shift | 0.0007 pK | same as real ligands ✅ |
| **residue-slot–scrambled protein, level shift** | **2.4e-08 pK** | see below |

The metric collapses when it should. The foreign-panel result matters
separately: the protein-inertness is identical on somebody else's ligands, so
it is not a property of each target's own chemotype.

### The scrambled-protein control is exactly null, and that is architectural

Permuting the residue-slot axis changes the output by 2.4e-08 pK — machine
zero, not a small effect. `ResidueEncoder` pools slots with a **sum** and
`ContactGrammar` reduces the residue axis with a **softmax-weighted sum**; both
are permutation-invariant, so the model cannot distinguish a protein from the
same protein with its slots shuffled. Proved algebraically in
`tests/test_probe_structure.py`.

**Consequence for language, not just for this stage.** No result from this
architecture may be described as pocket-aware, contact-resolved, or
biologically localized, whatever the attention weights look like — there is no
ordered structural information in its input to read. The atom→residue
cross-attention operates on an unordered bag of sequence-window summaries.

## 5. Novelty and scaffold strata

The double-cold split leaves **zero** shared Murcko scaffolds with `meta_train`
(`seen_scaffold_fraction = 0.000` in all three terciles), so the scaffold
stratum is degenerate by design — itself a confirmation of the split.

Terciles of mean max-Tanimoto to `meta_train`:

| tercile | targets | max Tanimoto | `r(interaction)` | centered protein shift |
|---|---:|---:|---:|---:|
| least novel | 14 | 0.396 | 0.289 | 0.00079 pK |
| mid | 14 | 0.325 | 0.295 | 0.00087 pK |
| most novel | 13 | 0.279 | **0.051** | 0.00045 pK |

A0's within-target ordering is concentrated in the less novel two thirds and
collapses on the most novel one (0.289/0.295 → 0.051). The protein-inertness is
uniform across all three. So the ordering A0 does have is **ligand-similarity
transfer**, and it is not protein-conditioned in any stratum.

## 6. Variability, reported at three levels

For `r(interaction)`:

| source | A0 | randinit |
|---|---:|---:|
| between component | **0.3466** | 0.0968 |
| between target within component | 0.1533 | 0.0595 |
| between seed within target | 0.0981 | **0.3800** |

For the trained model the dominant variance is **between components** — which
is why the component-paired bootstrap intervals on `r` are wide (±0.13) even
though the protein-swap contrasts are tight (±0.001). The two are measured on
different quantities: a *paired* difference cancels the target effect, an
*absolute* correlation does not.

For the ten random inits, the dominant variance is **between seeds** (0.3800),
four times the between-component spread. Section 7 is the consequence.

## 7. The correction: my v1 random-init claim was not supported

v1 recorded (E3, now `a2_readiness/SUPERSEDED.md`):

> At initialisation the interaction branch is strongly protein-conditioned in
> its ligand-differential. **Training destroys it.** This is a training
> outcome, not an architectural limit.

That was inferred from **one** random initialisation whose centered shift was
110× A0's. Ten independent initialisations show the shift is real in magnitude
and empty in content. The phase brief asked for four properties to be
distinguished; here they are:

| property | measurement | randinit |
|---|---|---|
| **arbitrary sensitivity** — does the output move? | mean centered shift, nearest donor | 0.0742 pK (vs A0's 0.0007) |
| **aligned** — does it move *with* true affinity differences? | corr(shift, centered truth) | +0.033 (nearest), +0.028 (farthest) |
| **reproducible** — do two inits move the same way? | mean pairwise cosine of the shift vectors, 1845 pairs | **−0.0025** (nearest), **−0.0116** (farthest), sd 0.47 |
| **useful** — does the branch order ligands correctly? | `r(interaction)` | **+0.0225** (vs A0's 0.2206) |

Two independently initialised copies of the same architecture, given the same
protein substitution on the same target, move in **directions that are
uncorrelated** — the cosine is zero to three decimal places over 1845 pairs,
with a spread of 0.47 that is exactly what random directions in this dimension
produce. The movement is also unaligned with truth and produces no usable
ordering.

**A nonzero output shift at random initialisation is undirected propagation of
input variation, not latent representational capacity.** The 110× figure was
real and its interpretation was wrong. E3 is withdrawn.

What replaces it is weaker and more honest: **whether this architecture *could*
be trained to order ligands protein-specifically is untested.** Random
initialisation cannot answer it, and no objective in R0-R14 ever asked for it
(see `DATAFLOW_AUDIT_V2.md` F6/F7).

## 8. Leakage review of this stage

| channel | status |
|---|---|
| `meta_test` | unreachable; `seal_record()` written into the artifact |
| donor pool | `meta_val` only — both proteins equally unseen |
| whitening statistics | `meta_train` only |
| label scale | `meta_train` only |
| query labels | metric targets only; never an input, selector or normaliser |
| episode bank | the frozen R6-R14 bank, unchanged |
| thresholds | frozen in `_frozen.py` before the first run |
| seed aggregation | within target, before the component bootstrap |
| foreign ligand panels | drawn once before any model ran, identical for every arm |

Three checkpoints scored on the same 41 targets and the same query panel are
not three biological samples; `component_bootstrap` averages within target,
then within component, then resamples the 19 components.

## Commands

```bash
conda run -n drug python -m tools.research.a2_readiness_v2.branch_ordering_v2 \
  --output tools/research/a2_readiness_v2/BRANCH_ORDERING_V2_meta_val.json
conda run -n drug python -m pytest tools/research/a2_readiness_v2/tests -q
```
