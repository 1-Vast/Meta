# Phase 3: where the protein's ligand-differential is actually lost

Authority: `ATTENTION_CAUSAL_meta_val.json` (+ `.rows.jsonl`). Frozen A0 seed
20260815 and one random init, double-cold `meta_val`, k=0, no training.

## Why divergence statistics were not enough

The v1 audit reported that trained atom→residue attention is as protein-
sensitive as a random one (Jensen-Shannon 0.241 vs 0.218 nats) while the pooled
features carry 9–20× less protein-differential, and concluded:

> Training has learned an *invariance in the readout* that discards the protein
> signal the attention still delivers.

Two problems. A random model's attention also changes under input substitution,
so the comparison establishes nothing causal. And Phase 2 contradicts the
conclusion directly: the ligand-differential is already protein-invariant at
`occupancy` and `mean_state` (cosine 1.0000), which are the **immediate outputs
of the attention block** — there is no readout between them and the attention.

This audit replaces the correlational statistic with interventions.

## The two channels, driven separately

Inside `ContactGrammar` the protein reaches the ligand-varying path through
exactly two paths:

```text
weight  = softmax( atom_query(atoms) @ residue_key(residues)ᵀ / √d )   # routing
context = weight @ residue_value(residues)                             # content
state   = atom_context( [atoms, context, atoms * context] )            # fusion
```

The audit re-implements this with `residue_key` and `residue_value` fed from
different proteins. `tests/test_attention_intervention.py` proves the
re-implementation is **bit-identical** to the module when both channels agree,
that a full swap equals a plain donor forward, and that a content-only swap
leaves the attention weights bit-identical — otherwise every number below would
be measuring my transcription.

## Result: A0, relative change in the ligand-differential

Level change in brackets. Instrument floor ~2e-07 (float32), four to six orders
of magnitude below every effect reported.

| stage | both channels | routing only | content only |
|---|---:|---:|---:|
| `weight` | 1.4630 [1.3637] | 1.4630 [1.3637] | 0.0000 [0.0000] |
| `context` | 0.4740 [0.0738] | 1.0389 [0.1166] | 0.9621 [0.1144] |
| `occupancy` | **0.0052** [0.0027] | 0.0118 [0.0045] | 0.0100 [0.0034] |
| `mean_state` | **0.0031** [0.0021] | 0.0043 [0.0032] | 0.0044 [0.0027] |
| `max_state` | 0.0503 [0.0016] | 0.0931 [0.0024] | 0.0903 [0.0025] |
| `embed` | 0.0146 [0.0011] | 0.0181 [0.0011] | 0.0091 [0.0002] |
| `section` | 0.0156 [0.0006] | 0.0203 [0.0006] | 0.0106 [0.0001] |
| `interaction` | 0.0044 [0.0219] | 0.0095 [0.0569] | 0.0084 [0.0302] |

**The attention weights change by 146% of their own magnitude and the retrieved
context by 47%. One step later, after fusion and pooling over atoms, the
ligand-differential of `mean_state` has changed by 0.31%.**

That step is a **150× attenuation** (0.4740 → 0.0031), and it happens *inside*
`ContactGrammar`, before `embed`, `section`, `interaction_head`,
`contact_weight` or anything else that could be called a readout.

Both channels contribute comparably (routing-only 1.04, content-only 0.96 at
`context`) and both are attenuated to the same place. Interestingly the
single-channel interventions produce *larger* downstream changes than the full
swap (0.0043 and 0.0044 vs 0.0031 at `mean_state`): driving key and value from
different proteins is a more disruptive, less self-consistent input than a
coherent wrong protein.

`max_state` is the exception — 0.0503, 16× `mean_state`. Max-pooling over ~128
atoms preserves per-atom protein-dependent variation that mean-pooling averages
away. Phase 2 tested `max_state` and found no usable protein-conditioned SAR in
it either (gain vs permuted protein −0.0401), so the survival is real and
uninformative.

## Gradient magnitudes confirm the location

Jacobian norm with respect to the protein residue tokens, normalised by each
quantity's own scale, computed separately for the level and the ligand-
differential:

| stage | ‖∂level/∂P‖ | ‖∂differential/∂P‖ | level/diff |
|---|---:|---:|---:|
| `weight` | 3.02e+01 | 9.13e+01 | 0.33 |
| `context` | 1.20e+01 | **1.13e+02** | 0.11 |
| `occupancy` | 4.70e-02 | **3.36e-02** | 1.40 |
| `mean_state` | 6.70e-01 | 5.48e-01 | 1.22 |
| `max_state` | 2.01e-01 | 1.83e+00 | 0.11 |
| `embed` | 7.62e-05 | 1.01e-02 | 0.01 |
| `section` | 1.64e-05 | 5.61e-03 | 0.00 |
| `interaction` | 3.02e-02 | 1.26e-02 | 2.40 |
| `protein_value` (the level branch) | 1.81e-01 | **0.00e+00** | — |

The protein-token gradient of the ligand-differential is **1.13e+02 at
`context` and 3.36e-02 at `occupancy`** — a factor of ~3,400 across the fusion
and pooling step. `protein_value`'s differential gradient is exactly zero, as
`tests/test_probe_structure.py` proves it must be.

## Comparison with a random initialisation

| stage | A0 | randinit |
|---|---:|---:|
| `context` | 0.4740 | 1.0775 |
| `mean_state` | 0.0031 | 0.3225 |
| **retained fraction** `mean_state`/`context` | **0.0065** | **0.299** |
| `embed` | 0.0146 | 0.3984 |
| `interaction` | 0.0044 | 0.4433 |

An untrained network propagates ~30% of the context perturbation through the
pooling; the trained one propagates 0.65%, a **46× stronger suppression**.

**This does not mean training destroyed something useful.** Phase 1 §7
established that the random network's protein sensitivity is undirected —
uncorrelated in direction across ten initialisations (pairwise cosine −0.003)
and unaligned with true affinity differences. What an untrained pooling
propagates is input variation, not signal. A trained pooling suppressing it is
equally consistent with learning a useful invariance and with discarding
information, and this measurement alone cannot tell them apart. Phase 2 can:
what survives carries no usable protein-conditioned SAR either, at any stage.

## Corrections to the record

| v1 claim | status |
|---|---|
| "the collapse is **downstream of the attention, in the readout**" (v1 E4, now `a2_readiness/SUPERSEDED.md`) | **wrong location.** It is downstream of the attention *weights* and upstream of every readout: at the `atom_context` fusion and the atom pooling, inside `ContactGrammar`. |
| "the trained model's attention is as protein-sensitive as a random one, so the protein is still being read" | the first half is confirmed by intervention (weights change 146%); the second half was over-read — the attention *routing* changes, but nothing downstream of the pooling carries usable protein-conditioned ordering. |

## Language constraint

`tests/test_probe_structure.py` proves the entire protein path is exactly
invariant to residue-slot permutation (measured 2.4e-08 pK). The atom→residue
cross-attention reads an unordered bag of sequence-window summaries. **No
result from this architecture may be described as pocket-aware,
contact-resolved, binding-site-localized, or biologically localized** — there is
no ordered structural information in the input for it to read, and there is no
common-frame complex geometry anywhere in the corpus (0/17,717).

## Commands

```bash
conda run -n drug python -m tools.research.a2_readiness_v2.attention_causal_audit --output tools/research/a2_readiness_v2/ATTENTION_CAUSAL_meta_val.json
conda run -n drug python -m pytest tools/research/a2_readiness_v2/tests/test_attention_intervention.py -q
```
