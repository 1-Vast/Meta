# UBSE-P0A target-marginal pocket-anchor preregistration

Date: 2026-07-29  
Status: frozen before full-registry contact aggregation or anchor training

## Scope and claim

Post-G1 forensics found that the additive arm retained modest held-domain
pocket localization while joint interaction training destroyed it. P0A asks
only:

> Can a protein-only model trained on the large BioLiP contact registry learn
> a homology/source/chemistry-closed target-marginal pocket anchor that is
> stronger than residue-propensity and shuffled-sequence controls?

P0A marginalizes over ligand. It cannot establish pair-specific information
and cannot solve strict dual-cold affinity by itself. A pass may provide a
frozen pocket proposal for a later three-dimensional event-transport student;
it does not authorize an `A + C + R` pair model, affinity loading, Stage-2,
confirmation, or sealed access.

This is a development-calibrated engineering gate. The G1 validation and
audit results have already informed the route choice; P0A will evaluate only
the existing 64-panel G1 validation role and will not report a new
confirmatory claim.

## Frozen sources and firewall

- BioLiP closed registry:
  `dataset/public/biolip2/processed/closed_registry.parquet`
- Required SHA-256:
  `7905e4edf88073f564baa4b2d4fb50d496432bc4e15e97cccbfa0766b1b0638d`
- G0PB panel manifest:
  `dataset/public/biolip2/processed/ubse_g0pb_panels.parquet`
- Required SHA-256:
  `4fea01e332eb3c60e41d76d5062d33cc95b13bc2e96b01df226532f78fe1b371`

Allowed source columns:

- `target_key`, `sequence`, `pubmed`, `scaffold`;
- `binding_residues_reindexed`.

Forbidden:

- ligand identity or features as model inputs;
- every affinity field or value;
- coordinates or three-dimensional events;
- development/confirmation affinity features or labels;
- sealed outcomes.

## Frozen closure

Reuse the G1 seed-1730 64-panel validation identity and the G0PB 88-panel
audit identity. Before reading contact labels for P0A training:

1. collect all validation and audit target sequences, PubMed identifiers, and
   scaffolds;
2. remove every source target with at least 0.40 containment of unique
   four-mers against any held sequence, using the same conservative project
   homology rule;
3. remove every source row sharing a held PubMed or scaffold;
4. retain the union-complement only.

The completed label-blind count was frozen at 62,849 rows and 38,781 exact
targets after union closure. Any mismatch stops.

The audit contact labels are not evaluated in P0A.

## Frozen target

For each retained target and PubMed, convert every valid residue list to a
binary sequence-length vector and average over its rows. Then average the
PubMed vectors with equal PubMed weight to obtain:

\[
m_{ti}=
\frac{1}{|P_t|}
\sum_{p\in P_t}
\frac{1}{|L_{tp}|}
\sum_{l\in L_{tp}}Y_{tpli}.
\]

This soft target prevents a prolific publication or ligand series from
dominating a target. Reject a target if any contributing row is malformed or
out of range.

For validation, use only the exact ligands in each frozen validation panel.
The binary validation pocket is the union of those ligands' contact residues.

## Frozen model

- Backbone:
  `facebook/esm2_t6_8M_UR50D`
- Revision:
  `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`
- Local checkpoint only.
- Train all backbone parameters plus one dropout-0.10 linear residue head.
- Maximum sequence window 1,000 residues; overlap 100.
- Overlapping validation logits are averaged.
- Length-bucketed batches are constrained by:
  - at most 64 windows;
  - at most 12,000 residue tokens; and
  - at most 4,000,000 summed squared sequence lengths.
- Mixed float16 CUDA training; no CPU model offload.
- Seeds: 1729, 1730, 1731.
- Four fixed epochs; report epochs 1, 2, and 4, with no audit-based selection.
- AdamW:
  - backbone learning rate `1e-5`;
  - head learning rate `3e-4`;
  - weight decay `1e-2`;
  - gradient clip norm `1.0`.

The loss is balanced soft BCE:

\[
\frac12
\frac{\sum m\,\operatorname{softplus}(-z)}{\sum m}
+
\frac12
\frac{\sum(1-m)\,\operatorname{softplus}(z)}{\sum(1-m)}.
\]

Each exact target supplies one target vector. Windows are shuffled
deterministically within each epoch.

## Frozen controls

1. **Residue-propensity null:** fit-only smoothed contact probability for
   amino-acid identity crossed with 20 relative-position bins.
2. **Shuffled sequence:** deterministically permute residues within each
   validation sequence for each seed, retain original labels and length, and
   evaluate the same trained model without retraining.
3. **Constant-position null:** replace a target's predicted probabilities by
   their sequence mean.

No ligand-conditioned control is relevant because P0A makes no pair claim.

## Frozen metrics

Compute each metric per validation target, then macro-average:

- average precision against the union pocket;
- AUROC, when both classes are present;
- oracle-size top-k recall, with `k` equal to the number of union-pocket
  residues;
- soft BCE against the within-panel contact frequency.

Report median across seeds. For correct-minus-control AP deltas, first average
each target across seeds and use 2,000 seed-1729 target-bootstrap replicates.

## Frozen gates

All must pass:

1. **P0A-1 substrate:** exactly 62,849 closed rows and 38,781 retained
   targets before contact parsing; at least 35,000 valid training targets and
   all 64 validation panels after parsing.
2. **P0A-2 absolute anchor:** validation AP at least 0.25, AUROC at least
   0.75, and oracle-size top-k recall at least 0.25.
3. **P0A-3 propensity increment:** AP at least 0.05 above the fit-only
   residue-propensity null, with target-bootstrap lower 95% bound greater than
   zero.
4. **P0A-4 sequence destruction:** AP at least 0.05 above shuffled sequence,
   with lower 95% bound greater than zero and all three seeds positive.
5. **P0A-5 stability:** correct-sequence AP range across seeds at most 0.05,
   all runs finite.
6. **P0A-6 execution/firewall:** CUDA training, no silent CPU offload, input
   hashes match, and no forbidden field or outcome is loaded.

Pass:
`FREEZE_UBSE_P0A_FOR_A1_POCKET_PROPOSAL_ONLY`.

Failure:
`STOP_UBSE_P0A_TARGET_MARGINAL_ANCHOR_INADEQUATE`.

Neither decision changes the UBSE-A0 coordinate-fetch wait or the G1 affinity
lock.

## Hardware record

Sample GPU utilization, framebuffer memory, power, and temperature during
training. Report mean, median, p95, and maximum where applicable. The goal is
to increase useful batching subject to the frozen memory limits, not to
overclock hardware or hide CPU bottlenecks.

Expected wall time on the RTX 4060 Laptop GPU: 30-120 minutes, including
three seeds and preprocessing.

