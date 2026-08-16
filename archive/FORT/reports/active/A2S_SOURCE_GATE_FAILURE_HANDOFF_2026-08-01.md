# English Research Handoff: A2S Source Information-Gate Failure

This handoff is intended for an independent research agent who receives only
partial attachments in a chat window. Do not assume CMAL, CAMP, ECA, CSRIO,
MDK, or any previous proposal is correct. Do not write code until the audit is
complete.

## Research objective

The desired contribution is a genuinely learned, transferable adaptation
mechanism for abundant-to-scarce drug-target affinity prediction. A mechanism
must learn from abundant source targets how to use k={1,3,5} measured support
labels from a strictly unseen recipient target to produce target-specific,
query-dependent compound-ranking improvement over an identical support-free
base. Calibration, retrieval, kernel regression, interpolation, closed-form
Bayesian updates, and simple fine-tuning are baselines only.

Every statement must be tagged `FACT`, `INFERENCE`, or `HYPOTHESIS`.

## Attachments to request first

1. `reports/active/A2S_SOURCE_INFORMATION_GATE_DECISION_2026-08-01.md`
2. `reports/active/a2s_source_information_gate_lock_2026-08-01.json`
3. `reports/active/a2s_source_information_gate_2026-08-01.json`
4. `research/a2s_information_gate.py`
5. `research/a2s_source_lock.py`
6. `reports/active/A2S_META_ADAPTATION_RESEARCH_HANDOFF_2026-08-01.md`

If the chat accepts only three attachments, use items 1-3. The decision report
contains the complete numerical result and the label firewall. Do not attach
recipient outcome files.

## Established facts

- The new metadata-only source lock uses only ChEMBL-37 `dual_cold_split=train`,
  `endpoint=pKi` metadata.
- Components are closed over target homology and all pipe-delimited document and
  assay tokens. Roles are fit=484 targets, probe=41 targets, locked=34 targets.
- All target, homology, document-token, and assay-token role overlaps are zero.
- One provenance component contains 380/559 source targets, so component power is
  a serious limitation.
- A label-free comparison gives 517 homology-only components (maximum size 4),
  but homology plus exact document cells gives 163 components with a
  347-target giant component. Assay cells do not connect targets. The current
  provenance field itself is the bottleneck.
- Only fit and probe labels were opened. Locked and recipient labels were not
  requested.
- The full fit/probe supervised base is component-cross-fitted: every fit row
  is predicted by a model excluding its component; probe rows are predicted by
  a model trained only on fit components. The base is a target-balanced ridge
  over the complete ligand feature vector plus a fixed 16-D protein projection.
- The OOF fold is highly unbalanced: the 380-target provenance component is
  held out as one fold (176,193 held-out rows versus 5,382 training rows), so
  that fold contains 97.0% of all fit-role OOF holdout rows. The exclusion
  contract is valid, but the real G0/G1 result is non-confirmatory until a
  provenance design with adequate fold balance exists.
- A nested linear probe exposes the same query design and support chemistry to
  G0 and G1; only the residual-label channel is masked for G0. G1 residuals are
  correct for the main arm. For k>=3, residuals are rolled among the same
  support compounds for the assignment control. k=1 assignment permutation is
  undefined by construction.
- A synthetic support-label signal is recovered at all k, proving that the
  diagnostic can detect an injected signal.
- Real probe `Delta_label` and `Delta_assign` have no positive component-bootstrap
  lower bound at k=3 or k=5. k=1 has no assignment test and no positive ranking
  label gap. The machine-readable decision is `NO_GO_INFORMATION_NOT_ADMITTED`,
  not an information-theoretic impossibility theorem.

## Questions for the independent agent

1. Is the new provenance closure scientifically justified, or does the giant
   document component create an avoidable power collapse? Propose a predeclared
   alternative using only metadata available at inference time. Do not use
   recipient labels to choose it.
2. Does the G0/G1 feature construction actually isolate support-label information,
   or could chemistry, base residual scale, assay identity, or target frequency
   still create a shortcut? Audit the code line by line.
3. Is the component-level OOF contract sufficient for the linear base, and what
   would be required to claim the same contract for a neural encoder/head?
4. How should high-data oracle headroom be defined so that query labels are not
   accidentally included in oracle fitting? Is the current target-specific
   leave-query-out oracle an admissible upper-bound diagnostic?
5. Given the null real gate and positive synthetic control, distinguish among:
   - no incremental information in this episode distribution;
   - insufficient probe capacity;
   - insufficient independent components;
   - chemically/assay-incoherent support/query construction;
   - a weak frozen base that obscures a real signal.
6. Design one and only one next source-only falsification experiment. It must
   specify the split, support/query closure, target/component unit, controls,
   primary ranking loss, bootstrap, minimum detectable effect, and a stop rule.
7. Decide whether any candidate mechanism remains paper-worthy under the current
   data. Explicitly reject renamed CNP/ANP, kernel, posterior, calibration, or
   generic listwise reranker proposals.
8. Decide whether the OOF fold imbalance alone invalidates the probe or merely
   downgrades it. If it invalidates the probe, specify a label-free component
   construction that avoids a giant provenance component without weakening the
   no-overlap contract.

## Required response format

### Part 1: Audit

List findings in descending severity. For each finding provide `FACT`,
`INFERENCE`, or `HYPOTHESIS`, the file/function, and why it changes the
scientific interpretation.

### Part 2: Identifiability

State what k=1, k=3, and k=5 can plausibly identify under the measured
support/query chemistry and assay metadata. Do not claim a universal
information-theoretic impossibility from this one source probe.

### Part 3: Candidate routes

Give at most three routes. For each route define the learned object, why it is
transferable, why it is not merely calibration/kernel/posterior, and the single
experiment that could falsify it. If none is defensible, say so.

### Part 4: Decision

Choose `REOPEN_WITH_NEW_STRATUM`, `REQUIRE_NEW_DATASET`, or
`STOP_CURRENT_A2S_EPISODE`. Explain why the other choices are rejected. Do not
recommend training on the locked role unless the decision is first changed by a
pre-registered protocol.

## Hard constraints

- No recipient labels.
- No locked-role labels during this review.
- No broad hyperparameter sweep.
- No claim based only on correct-vs-wrong arm separation.
- No model promotion to the `model/` directory.
- No GitHub commit or push.
- Preserve the distinction between a bounded diagnostic null and a theorem.
