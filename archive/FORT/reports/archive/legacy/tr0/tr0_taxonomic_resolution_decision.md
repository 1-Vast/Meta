# TR-0 taxonomic-resolution gate — decision

Verdict: **`TR0_PREMISE_FAIL_STOP`** / `KLIFS_GROUP_IS_NOT_A_RESOLVABLE_TRANSFER_RESOLUTION`.

Preregistration: `reports/active/taxonomic_resolution_preregistration.md`. Result:
`reports/active/tr0_taxonomic_resolution.json`. One seed (1729), no model trained, no
development/confirmation/sealed affinity label read for selection. The harness reproduces the frozen
KirHub strict A1 arms **bit-identically** (ligand-only `0.0429`, true ESM `0.0719`, group centroid
`0.0830`, over 308 homology components), so the appended control arms and contrasts are trustworthy.

## What the reframe claimed and what the gate tested

The reframe (verdict ② estimand pivot) hypothesised that the coarse **KLIFS group** is a load-bearing,
transferable ranking resolution: group-level transfer beats ligand-only B0 (G1), is group-specific
(G2), and would form the always-available floor of a group→target resolution hierarchy whose
target-level deviation is gated by the CRA-Ω abstention module (TR-1). G0 was the preregistered
falsification control: if the signal is genuinely group-level, removing the query's **own** group from
the pooling must collapse the gain.

## Result

| Contrast | mean | 95% interval | gate |
|---|---:|---|---|
| G1 `group − ligand_only` | +0.0400 | [+0.0181, +0.0616] | PASS (≥0.03, LCB>0) |
| G2 `group − random_group_labels` | +0.0347 | [+0.0119, +0.0578] | PASS (LCB>0) |
| G2 `group − shuffled_group` | +0.0193 | [+0.0022, +0.0361] | PASS (LCB>0) |
| **G0 `group − group_cold`** | **+0.0091** | **[−0.0076, +0.0256]** | **FAIL (LCB≤0)** |
| G0 `group_cold − ligand_only` | +0.0309 | [+0.0107, +0.0515] | — (gain survives without own group) |
| G3 `true − group` (diagnostic) | −0.0110 | [−0.0321, +0.0105] | fine protein does not beat coarse |

Arm ladder (component-macro Spearman, 308 units): random_protein 0.0200 < ligand_only 0.0429 <
random_group_labels 0.0482 < shuffled_group 0.0636 < true_protein 0.0719 ≈ **group_cold 0.0739** <
group_centroid 0.0830.

## Interpretation (honest decomposition)

1. **A real coarse protein signal beats B0.** G1/G2 hold: a coherent protein-kernel predictor beats
   ligand-only by +0.040, and it is specific — it beats a randomised group structure (+0.035) and a
   deranged wrong-group centroid (+0.019). This is genuine protein information, not an artefact.
2. **But the resolution is coarser than the KLIFS group.** G0 fails: a leave-own-group-out centroid
   (`group_cold`, pooling only over *other* kinase groups) retains +0.031 of the +0.040 gain, and the
   own-group increment is +0.009 with an interval crossing zero. Within the kinome the groups are
   similar enough that the query's own group is not a necessary, resolvable transfer level; the
   load-bearing property is "a coherent generic-kinase centroid", not group identity.
3. **G2 and G0 are consistent.** `random_group_labels` (0.048) destroys all group coherence and
   `shuffled_group` (0.064) forces one specific *wrong* group, so both are worse than the own-group
   centroid; but `group_cold` (0.074) — the low-variance average of all other groups — is a good
   generic representation and nearly matches it. Coherence and correctness beat incoherence; own-group
   identity adds nothing significant on top.
4. **The fine level was already the known ceiling.** G3 reproduces the frozen A1 result that per-target
   ESM does not beat the coarse centroid (−0.011).

## Decision

The reframe's identifying control (G0) failed. The charter forbids a mechanism that fails its
identifying control from proceeding to expensive optimisation, and building the CRA-Ω resolution-
abstention module (TR-1) would require a group→target resolution structure that G0 shows does not
exist: there is no significant group resolution to calibrate, only a flat generic-kinase gain over B0
that is already frozen in A1, is not group-specific, and lives on a within-source substrate (KirHub
cannot isolate assay/document, so it is a mechanism probe, never confirmation). This does **not** meet
the charter's bar for meaningful progress (not a new source, not target-specific identifiability, not
a credible improvement over the strongest matched protein arm, no powered independent confirmation).

**TR-1 is not authorised.** No threshold was relaxed, no arm width/rank/seed changed, and no goalpost
moved after seeing the result: G0 was a preregistered falsification whose predicted outcome (gain
collapses) did not occur. The taxonomic-resolution reframe joins the closed-route ledger.

The program verdict returns to **② `SIGNAL_PRESENT_EVIDENCE_INSUFFICIENT`**: a small real protein
signal beats ligand-only B0 (now shown to be a coarse generic-kinase effect, not group- or
target-resolved), but the evidence does not support a stable, attributable, target- or
group-specific improvement, and the remaining plausible routes are reduced to explicit data and power
constraints (`NO_OPEN_POWERED_INDEPENDENT_PANEL`).

Reopening the resolution idea would require a **multi-family** substrate (not within-kinome), where
distinct protein families are far enough apart that a family/group resolution is genuinely
load-bearing and its own-family-cold falsification collapses — i.e. the same open, powered,
document-isolated factorial panel already named as the program's binding data constraint.

`sealed_test_consumed=false`; `confirmation_labels_read=true` (pre-existing; TR read no confirmation
labels).
