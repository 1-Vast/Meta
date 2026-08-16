# Stage 2: sequence-derived locality-aware protein refinement — REJECTED

Numerical authority: `seed*/RESULT.json`, `DISCRIMINATOR_meta_val.json`,
`MULTISEED_meta_val.json` (+`.rows.jsonl`), `PAIRED_ANALYSIS.json`.

## Hypothesis

Mac-Diff's transferable idea is **locality-aware alignment**, not diffusion. The
trunk lets every ligand atom attend densely to all 128 protein slots, so nothing
expresses that a binding region is contiguous. The Stage 0 audit established
that the 128 slots are ordered contiguous residue windows
(`floor(i*128/L)`), so a sequence-window prior over slots is legitimate.

Predicted effect: a better zero-shot protein representation, i.e. **lower k=0
MSE**, with the Tanimoto transport held fixed so the two cannot be confounded.

## Exact code change

* `model/interaction_grammar.py`: added a `refine_slots` hook that is an **exact
  identity** by default, so every existing architecture is numerically
  unchanged.
* `model/locality_grammar.py`:
  * `SlotLocalityRefiner` — band-limited self-attention over ordered slots
    (`|i-j| <= 8`) plus a learned global read, behind a **zero-initialised
    residual gate**;
  * `LocalizedContactGrammar` — atom-to-slot cross attention restricted to the
    top-`m` slots per atom, **unioned across atoms** so a slot needed by any
    atom stays available to all;
  * `LocalityGrammarModel` — the accepted `similarity_only` model plus the
    above. Transport inherited unchanged.
* `--arch locality`, opt-in. No default changed.

Naming is binding: this is a **sequence-derived locality prior**. It is not a
contact map, carries no contact supervision, and makes no structural claim.

## Structural tests — all pass (14 gates)

Band mask is a symmetric sequence window; **zero gate reproduces the accepted
baseline exactly**; the refiner is a no-op at gate 0 and active when opened;
padded slots stay zero and perturbing them cannot change live-slot output; the
localizer never selects padding and ignores masked-slot content; top-k >= slot
count reduces exactly to the dense baseline; variable protein (12/17/33) and
ligand (3/5/9) lengths are stable; k=0 returns the endpoint exactly; support
permutation invariance and query equivariance hold with the gate open; protein
slot shuffle changes the zero-shot prediction; no dead trainable tensor; query
labels are not an input. Full suite: 340 passed.

## Result: the gate fails

Three matched seeds, identical budget and banks, complete 44-episode `meta_val`.

### k=0 MSE — the stage's target

| seed | accepted baseline | locality | change |
|---|---:|---:|---:|
| 20260812 | 1.566 | 1.591 | **+0.025 worse** |
| 20260813 | 1.710 | 1.760 | **+0.050 worse** |
| 20260814 | 1.561 | 1.615 | **+0.054 worse** |
| mean | **1.612** | 1.655 | **+0.043 worse** |

**k=0 regressed in 3 of 3 seeds.**

### Paired component-level bootstrap, locality minus baseline (positive = locality better)

| k | MSE reduction | 95% CI | LB>0 |
|---|---:|---|---|
| 0 | **-0.043** | [-0.167, +0.061] | no |
| 1 | -0.011 | [-0.111, +0.086] | no |
| 2 | -0.010 | [-0.105, +0.074] | no |
| 3 | -0.006 | [-0.092, +0.071] | no |
| 5 | -0.013 | [-0.083, +0.049] | no |

Point estimates are negative at **every** k. Nothing is significant in either
direction, but there is no evidence of benefit anywhere and the stated k=0
target moved the wrong way.

### The one-seed signal did not survive

The seed-20260812 discriminator showed k=0 CI 0.575 against 0.556 and lower MSE
at every k>=1, which suggested locality improved endpoint *shape* if not
calibration. Seeds 20260813 and 20260814 reverse it: k=0 CI 0.542 against 0.565
and 0.556 against 0.541, and k>=1 MSE mixed. This is a textbook case of a
one-seed result that is seed variance, and it is why the protocol requires three.

### Consistency check that did hold

The Tanimoto transport still works inside the locality model — its own
within-checkpoint `full` versus `level_only` gives +0.131 / +0.160 / +0.236 MSE
at k=2/3/5 with component lower bounds above zero, matching Stage 6. The
transport was correctly held fixed and is unaffected by the failed refinement.

## Decision

**REJECT.** The preregistered gate was "multi-seed k=0 improves over the accepted
baseline". k=0 regressed in every seed and every cross-arm point estimate is
negative. The module is not promoted to the active path.

Per the staged protocol, Stage 3 (Mac-Diff conformational sidecar) and Stage 4
(support-conditioned conformational routing) proceed **only if** this stage
passes. They are therefore **not authorised**, and no Mac-Diff inference,
weights or ensemble sidecar was run. The Mac-Diff-inspired direction closes here
for this corpus.

## Failure analysis

Three candidate explanations, none tested further because the gate failed:

1. **The prior may be wrong for this representation.** After `_slot_pool`
   mean-pools `L/128` residues per slot, a slot already mixes a whole window;
   an additional 8-slot band may add little beyond what pooling already did.
2. **Dense attention may not have been the binding constraint.** The Stage 0
   audit showed k=0 error is dominated by weak ligand-side and calibration
   signal, not by protein slot mixing.
3. **Capacity cost without payoff.** The refiner adds parameters and a gate that
   must be learned within the same 800-step budget, and the zero-init gate means
   any benefit must be earned late in training.

Resources: 1,256 s per seed, peak 6,527 MB against the baseline's 6,506 MB.

## Retained as a negative result

`model/locality_grammar.py`, `tests/test_locality_grammar.py` and these results
are retained. The `refine_slots` hook stays in `interaction_grammar.py` because
it is an exact identity by default and is covered by tests. `--arch locality` is
opt-in and is **not** part of the active path.
