# Dataflow audit v2 (2026-08-16)

Supersedes v1's dataflow audit (consolidated into `tools/research/a2_readiness/SUPERSEDED.md`). Every statement is
verified against the working tree or measured; nothing is taken from a report.
Findings that carry over unchanged are listed compactly; the changes are the
point of this document.

## Carried forward unchanged (re-verified)

| id | finding | evidence |
|---|---|---|
| F2 | dead `index` computation in `qpsmp_data.protein_for_target`; raises `StopIteration` rather than a domain error for an absent target | cosmetic |
| F3 | with `use_learned_key=False` (A0), `SimilarityTransport` trains exactly two scalars — `similarity_scale`, `log_shrinkage`; `key.weight` and `log_temperature` are frozen and the residual is `.detach()`ed | `similarity_grammar.py:79-80,197` |
| F4 | at k=1 the support softmax is identically 1.0, so `sar_adaptation ≡ 0` exactly, at every parameter value | `similarity_grammar.py:198-205` |
| F5 | train/eval parity holds; `LabelScale` is fitted on `meta_train` cells only | `train_qpsmp.py`, `compact_episode` on all three paths |
| — | endpoint decomposes exactly as `ligand_value(L) + protein_value(P) + interaction(P,L)`; `protein_value` is constant within a target | now **tested**, not asserted: `tests/test_probe_structure.py` |

## F1 — resolved

The `meta_test` seal defect is repaired, verified and written up as a
governance incident. Default is now fail-closed, opening requires a written
authorization, artifacts derive their seal block from the object, and 16
contract tests hold the property. Numbers verified unchanged (105/105
bit-identical). See `GOVERNANCE_INCIDENT.md`.

A **second** defect surfaced during the repair: `audit_research_record.py`
classified ten double-cold stage summaries as `bindingdb_ki_main_v0`
older-protocol evidence because it looked for the split in three specific
fields and they record it in a fourth. None of the ten names that corpus. Fixed
by scanning the whole document and adding a `split_undeclared` state; the
older-protocol count went 10 → **0**.

## F6 — confirmed in code, with the exact firing conditions

`train_qpsmp.py:774-793`:

```python
if adapt and support_size > 0:
    ...
    wrong_zero = wrong_protein_zero_shot(model, data, episode,
                                         episode.spec.donor_target)
    correct_zero_error = (full.zero_shot - query_y).square().mean()
    wrong_zero_error  = (wrong_zero      - query_y).square().mean()
    loss_protein = binding_contrastive_loss(
        [correct_zero_error, wrong_zero_error], config.binding_temperature)
```

Both errors are **uncentered**. An uncentered protein contrast is satisfied
completely by the level branch: a 0.215 pK level shift already separates
correct from wrong, so the gradient is extinguished by `protein_head` alone and
never reaches the interaction branch's ligand-differential.

Three firing conditions, all newly verified and all relevant to Phase 4:

1. **`support_size > 0`.** `support_size = train_support_sizes[(step-1) % 5]`
   over `(0,1,2,3,5)` (`train_qpsmp.py:734`), so the term fires on **4 of every
   5 steps — 960 of A0's 1200**. The 240 k=0 steps get no protein contrast at
   all, even though the term supervises `full.zero_shot`, a k=0 quantity. The
   objective supervises the zero-shot endpoint but is gated on the episode
   having support: an incidental coupling, not a design.
2. **`adapt`**, i.e. `not (zero_support_only or phase_a)`. A0 has
   `representation_warmup_fraction = 0.0`, so `warmup_steps = 0` and
   `phase_a = step <= 0` is never true for a loop starting at `step = 1`
   (`train_qpsmp.py:728`). **Warmup does not gate this term in A0.** It would
   for any arm that sets the fraction above zero.
3. `protein_contrast_loss_weight` defaults to **0.5** (`train_qpsmp.py:89`).

## F7 — confirmed, and the record correction is now mandatory

`wrong_protein_prediction` (`train_qpsmp.py:538-552`) uses the same uncentered
form throughout the evaluation controls. Every wrong-protein number in R0-R14
therefore measures **level specificity**.

R3R4's "first resolved protein specificity" (+0.4216 at k=2) is a level result.
It is real and it is not evidence of protein-conditioned SAR. Phase 1 measures
the ordering version of the same control and finds
**−0.0002 [−0.0015, +0.0008]**.

## F8 — the contrastive form rewards donor destruction

```python
def binding_contrastive_loss(errors, temperature):
    logits = -torch.stack(errors) / temperature
    return F.cross_entropy(logits.unsqueeze(0), zeros)
```

For two errors this is `softplus((correct − wrong) / T)`, minimised by making
`correct` small **or** `wrong` large. Nothing in the term, and nothing in any
gate recorded in R0-R14, distinguishes "the model got better on the right
protein" from "the model got worse on the wrong one".

This is unguarded in the incumbent and would be unguarded in any CPC variant
built on the same helper. `PREREGISTRATION_V2.md` §4 makes non-degradation of
the donor arm a hard gate rather than a reported statistic.

## F9 — the protein path is exactly invariant to residue-slot order

New, and a constraint on language rather than on numbers.

`ResidueEncoder.forward` pools slots with `residues.sum(1) / gate.sum(1)`, and
`ContactGrammar` reduces the residue axis with a softmax-weighted sum. Both are
permutation-invariant, so the model cannot distinguish a protein from the same
protein with its slots shuffled. Measured at **2.4e-08 pK** — machine zero —
and proved algebraically in `tests/test_probe_structure.py`.

The atom→residue cross-attention operates on an **unordered bag of
sequence-window summaries**. No claim of pocket awareness, contact resolution
or biological localization is available to this architecture regardless of what
its attention weights look like.

## F10 — where the protein's ligand-differential is destroyed

Superseding v1's E4, which named the readout.

The loss is at the `atom_context` fusion and the atom pooling **inside**
`ContactGrammar`: attention weights change 146% under a protein swap, `context`
47%, and `mean_state` 0.31% — a 150× attenuation across one step, with a
~3,400× drop in the protein-token Jacobian of the ligand-differential. Nothing
downstream (`embed`, `section`, `interaction_head`, `contact_weight`) is
responsible. See `ATTENTION_CAUSAL_AUDIT.md`.

## E1-E4 — the measured findings, revised

| id | v1 | v2 |
|---|---|---|
| **E1** ordering is interaction-borne | `r_full` 0.213, `r_ligand_only` 0.027, increment +0.1855 [+0.0566, +0.3236] | **unchanged and confirmed** across 5 donor strata and 10 extra models |
| **E2** the interaction branch's ordering is protein-inert | level 0.2150 pK vs centered 0.0007 pK, contrast −0.0002 [−0.0015, +0.0008] | **strengthened**: holds at all five donor distances (level 0.215→0.342 pK, centered 0.0007→0.0011 pK); measurement floor exactly 0; DECISIVE_NULL at 4/5 strata |
| **E3** "the architecture can express what training removes" | one random init shifted 110× more | **withdrawn.** Ten inits move in uncorrelated directions (pairwise cosine −0.003 over 1845 pairs), unaligned with truth (+0.03), and produce no usable ordering (`r` 0.023). A nonzero shift at init is undirected propagation, not capacity. |
| **E4** the collapse is in the readout | attention JS 0.241 vs randinit 0.218 | **relocated.** Causal interventions put it at the fusion+pooling inside `ContactGrammar`, upstream of every readout. |
| **E5** *(new)* the trunk carries no protein-conditioned SAR at any internal stage | — | differential cosine 0.998–1.000 at every representation; trained probes show a protein-conditioning gain of +0.017 over a capacity-matched permuted protein on `embed`, and negative gains on `mean_state` and `max_state` |
| **E6** *(new)* the trunk carries a **protein-independent** transferable SAR direction | — | `embed`, Δ-affinity `r` +0.2623 [+0.1295, +0.4055] on held-out components from 1,553 parameters; better than the raw ligand encoder's +0.1188 (unresolved) |

## Summary

| id | finding | severity |
|---|---|---|
| F1 | `meta_test` seal was opt-out; a false claim reached 7 artifacts; the audit's protocol classification was wrong for 10 more | **resolved**, incident filed |
| F2 | dead code in `protein_for_target` | cosmetic |
| F3-F5 | transport trains two scalars; k=1 adaptation is identically zero; train/eval parity holds | descriptive |
| F6 | the protein contrast is uncentered, fires on 4/5 steps, and is ungated by warmup in A0 | causal, but see E5 |
| F7 | every wrong-protein control in R0-R14 measures level specificity | **record correction required** |
| F8 | the contrastive form is minimised by degrading the donor | **design constraint on any successor** |
| F9 | the protein path is exactly slot-permutation invariant | **language constraint** |
| F10 | the ligand-differential dies at the fusion/pooling, not the readout | supersedes v1 E4 |
| E5 | no protein-conditioned SAR at any internal stage | **rejects A2** |
| E6 | a protein-independent transferable SAR direction exists in `embed` | the only positive result of this cycle |
