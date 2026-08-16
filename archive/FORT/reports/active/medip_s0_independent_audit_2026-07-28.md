# MEDIP S0 Independent Failure Audit

**Date:** 2026-07-28

**Scope:** independent review of the frozen MEDIP-S0 preregistration, model,
simulator, tests, accepted JSON, and decision. No file from the formal run was
edited for this audit.

**Result:** the stored result and STOP verdict are reproducible, but the
experiment supports a narrower interpretation than a fully confirmatory
mechanism calibration.

```text
MEDIP_S0_ENGINEERING_CALIBRATION_STOP
```

## 1. Reproduction

In `D:\anaconda\envs\drug\python.exe`:

- all 10 MEDIP tests pass;
- a fresh formal rerun exactly reproduces the stored JSON and verdict;
- no likelihood-sign, censor-direction, ordinal-index, train/test, or
  optimization mismatch explains the two failed gates;
- no real dataset or affinity value is read.

The STOP decision therefore stands.

## 2. Preregistration Limitation

The prose preregistration freezes dimensions, seeds, optimizer, steps, loss
weights, variants, metrics, and thresholds. It does not numerically freeze the
complete data-generating process:

- interaction and main-effect scales;
- endpoint/source assignment coefficients;
- observation offsets and noise scales;
- censor limits;
- selectivity comparison count;
- the latent-difference selection threshold.

The result binds the prose preregistration hash but does not bind generator and
model source hashes. The run is deterministic engineering evidence, not a
fully confirmatory preregistered mechanism experiment.

## 3. Metadata Destruction Interpretation

The correct and metadata-shuffle arms use the same decoupled architecture.
They vary observation IDs, not whether metadata is permitted inside the
interaction encoder. The arm therefore cannot establish that architectural
decoupling itself is or is not load-bearing.

Correct metadata is demonstrably important for the observation process:

| Metric | Correct | Metadata shuffle | Registered contrast |
| --- | ---: | ---: | ---: |
| Exact-mean RMSE | 0.133927 | 1.793285 | shuffle minus correct = 1.677035 |
| Cross-target ordering | 0.984848 | 0.913826 | descriptive loss = 0.071023 |
| Mixed-difference correlation | 0.998710 | 0.993709 | correct minus shuffle = 0.005862 |

The gate failed because mixed-difference correlation preserved direction and
relative geometry under metadata corruption. The defensible conclusion is:

> Correct endpoint/source IDs are load-bearing for calibration and affect
> ordering, but this redundant simulator did not show that they are necessary
> for preserving scale-insensitive mixed-difference correlation.

It is not defensible to conclude that metadata decoupling is generally
irrelevant.

## 4. Selectivity Destruction Interpretation

The simulator rejects candidate target pairs whose true latent difference is
below `0.20`, although the preregistration says only that pairs are sampled.
This is outcome-conditioned synthetic inclusion. The shuffle preserves the
truth-selected pair identities and permutes only their signs.

The ordering metric also uses the full score, including target main effects.
It is therefore not specific to ligand-conditioned reordering: the separable
null reaches ordering `0.580492` while its doubly centered interaction is
numerically zero.

Every training target-ligand pair already supplies continuous, binary, and
ordinal evidence from the same latent matrix. Selectivity is highly redundant,
and the frozen design has no nested selectivity-omission arm. Selectivity-sign
shuffle changes ordering by `-0.000947` and leaves mixed-difference
correlation essentially unchanged (`0.998710` versus `0.998464`).

The defensible conclusion is:

> This S0 design does not demonstrate incremental interaction information from
> the selectivity loss.

It does not establish that matched selectivity is biologically useless under a
different, nonredundant observation topology.

## 5. Minor Diagnostic Defect

If optimization stopped for a non-finite loss, the runner would report the
number of logged loss fields as `steps_completed`, not the iteration count.
All formal losses were finite and all 500 steps completed, so this defect does
not affect the accepted result. It must be corrected and tested before any
newly preregistered successor run; it is not a reason to alter or rerun S0.

## 6. Claim Boundary and Reopening

S0 establishes:

- the endpoint-separated likelihood implementation is numerically stable;
- correct metadata greatly improves observation calibration;
- the direct interaction carrier can recover this synthetic truth;
- the exact separable comparator has zero doubly centered interaction;
- the current mixed objective does not pass its frozen destruction battery.

S0 does not establish:

- type-I calibration under a generated true interaction null;
- incremental value of metadata decoupling versus metadata injection;
- incremental value of selectivity versus a nested omission arm;
- ligand-conditioned selectivity rather than target main effects;
- biological, open-data, or strict dual-cold predictive validity.

Do not rescue S0 by changing a weight, rank, width, seed, sample size, or
threshold. A valid successor is a scientifically different preregistered
question and must:

1. bind an executable DGP plus generator/model hashes before execution;
2. use outcome-independent or known-probability comparison inclusion;
3. model a realistic nonredundant modality-support topology;
4. include nested component-omission and matched-corruption controls;
5. separate doubly centered interaction or genuine ligand-order reversals
   from observation calibration and target main effects;
6. include a generated true-null regime if type-I calibration is claimed.

Real-label work remains blocked by the existing information, provenance,
topology, firewall, and power gates.
