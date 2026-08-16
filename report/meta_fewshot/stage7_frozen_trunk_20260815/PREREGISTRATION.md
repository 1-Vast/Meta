# Stage 0 audit and Stage 1 preregistration

## Stage 0: verified protein data flow

Read from source, not assumed.

### Artifacts that exist for the active DTA corpus

| artifact | present | shape / form | provenance |
|---|---|---|---|
| raw sequence | yes | `proteins.jsonl`, 499 targets | corpus |
| ESM pooled vector | yes | `[640]` float16 | mean over residue states |
| 128 protein slots | yes | `[128, 640]` float16 | `_slot_pool` in `scripts/cache_structure_proteins.py` |
| slot mask | yes | `[128]` uint8 | `counts > 0` per slot |
| slot residue chemistry | yes | `[128, 4]` | `protein_slot_chemistry`, computed from the sequence |
| attention / contact tensors | **no** | — | none for DTA targets |
| distance bins | **no** for DTA | holo corpus only, disjoint ligands | `pilot20k_structure_supervision_v2` |
| common-frame coordinates | **no** | 0 of 17,717 cells | `GEOMETRY_COVERAGE_AUDIT.json` |

### Are the 128 slots ordered residue intervals or unordered summaries?

**Ordered contiguous residue intervals.** `_slot_pool` assigns residue `i` to
slot `floor(i * 128 / L)` and mean-pools, and `residue_slot_mapping` in
`scripts/build_protein_bank.py` records the matching half-open interval
`[floor(j*L/128), floor((j+1)*L/128))`. Slot `j` and slot `j+1` are adjacent
sequence windows. A slot is masked out only when `L < 128` leaves a bin empty.

Consequence for Stage 2: **sequence adjacency in slot space is real**, so a
band- or window-limited locality prior over slots is defensible. It is a
*sequence-derived locality prior*, not a contact map, because no contact
supervision exists for these targets. That naming is binding.

### Active data flow

`QPSMPData.materialize` -> `EpisodeBatch(protein_pooled[640],
protein_tokens[128,640], protein_mask[128], protein_chemistry[128,4],
support/query atoms[·,A,32], bonds[·,A,A,12], masks, labels, Morgan
fingerprints[·,1024])` -> `ResidueEncoder` (slot MLP + chemistry bias) ->
`ContactGrammar` (atom-to-slot cross attention, **dense over all 128 slots**) ->
shared contact-type dictionary -> interaction embedding `e` -> zero-shot
endpoint `f0` -> support residual `r_k = y_k - f0(P,L_k)` (detached) ->
transport -> `f = f0(q) + s(n) * sum_k w_qk r_k`.

## Stage 1 preregistration

### Hypothesis

The incumbent `grammar` trunk and the validated fixed chemical kernel are
**complementary**: replacing only the support transport on a frozen incumbent
trunk improves k=2/3/5 without degrading ranking. If true, the Stage 6 cross-arm
contradiction was trunk co-adaptation, not a mechanism conflict.

### Exact arm

`scripts/stage1_frozen_trunk_transport.py`. Each incumbent checkpoint is frozen;
its `f0`, its shrinkage `s(n)` and the episode bank are held fixed. Only the
transport varies:

| transport | definition |
|---|---|
| `mean` | `f0(q) + s(n) * mean_k r_k` — the incumbent's own level baseline |
| `tanimoto` | `f0(q) + s(n) * sum_k softmax_k(8 * Tanimoto1024) r_k` |
| `nearest` | `f0(q) + s(n) * r_argmax(Tanimoto)` |
| `incumbent` | the checkpoint's own learned transport |

No training. No parameter is changed. This is an inference-time swap.

### Seeds, banks, budget

Seeds 20260812 / 20260813 / 20260814 (`grammar` checkpoints from Stage 5 and
Stage 6). Development split `meta_val`, complete eligible bank (44 episodes),
`evaluation_seed=73101`, nested k = 1/2/3/5. `meta_test` is run **descriptively
only** and cannot select anything. Zero training cost.

### Metrics

MSE in pK^2, concordance index, Spearman, per k. Paired target-level and
homology-component-level bootstrap, 9,999 draws, pairing within the same
checkpoint and episode. Intervals are conditional on the three trained seeds
because seeds are averaged before components are resampled.

### Stopping rule

**Advance** only if, on `meta_val`, `tanimoto` beats `mean` at k=2, 3 and 5 in
MSE with a component-level lower bound above zero, and CI and Spearman do not
degrade at any of those k.

**Do not** prefer `nearest` over `tanimoto` on a point estimate alone; it is
adopted only if its paired interval against `tanimoto` excludes zero.

**Reject** the complementarity hypothesis if `tanimoto` fails to beat `mean` on
the frozen incumbent trunk. That would mean the Stage 6 gain depended on
co-adaptation and does not transfer.
