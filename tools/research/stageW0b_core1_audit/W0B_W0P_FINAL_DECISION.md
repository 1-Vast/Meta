# Core Task 1 W0/W0-P final GO/NO-GO decision

Date: 2026-08-17.
Machine-readable: `W0B_W0P_FINAL_DECISION.json`.

## Decision

**NO-GO: Core Task 1 cannot be resolved with current local assets.**

This is a data/positive-control limitation statement, not a biological claim.
No neural W1 model was trained and none is authorized on current assets.

## Evidence

1. **W0-P positive control FAILED.** Leave-one-pair-out, 3 seeds,
   low-capacity bilinear model. Correct mutation positions sign accuracy
   0.240; random positions 0.156; BLOSUM-approximate unrelated positions
   0.125. Global ESM pooled difference reached 0.760, but with 6 pairs this
   is recorded as unexplained and under-powered, not as a pass.
2. **The local W0-P panel cannot be enlarged.** CD-HIT98 on BindingDB
   sequences yields 7 candidate pairs with >=1 shared ligand (34 rows max;
   6 pairs/32 rows with >=3 shared ligands). No standard resistance/mutation
   or ortholog panel exists locally.
3. **Censoring is severe.** Davis 71.2%, Metz 60.4%, Klaeger 93.5% of labels
   sit at detection floors. Raw Pearson/MSE/cliff statistics on these panels
   would be floor-dominated.
4. **Strict-MMP support is small where labels are interpretable.** Davis 7
   classes; Klaeger 38 classes; Metz 480 classes but 60% censored.
5. **KIBA is uncensored but score-semantic**, and cannot be the sole positive
   evidence for Core Task 1.

## Consequences

- Stage W W1 remains **PAUSED / NO-GO**.
- Any future cold-protein null is not interpretable as biological absence
  until a W0-P panel passes.
- Stage W, W0b and W0P preregistrations are frozen and unmodified.

## Authorized next actions

1. Acquire a licensed standard W0-P panel (gatekeeper/resistance mutations or
   ortholog/point-mutant panels with matched ligand panels) and record source,
   version, license, checksum and label semantics.
2. Re-census the single-platform panels after censored exclusion or
   interval-censoring.
3. Re-preregister W0-P and W1 thresholds from the effective sample size of the
   enlarged panel, then rerun.

## Censored re-census update (round 25)

`W0B_CENSORED_RECENSUS.json`: after excluding detection-floor rows, the broad
all-ligand-pair layer still has ample support on all three single-platform
panels:
- Davis all_pairs: 87,371 pairs / 2,192 classes / 1,988 rich classes /
  3,935,561 cross-component D rows / EIU 205;
- Metz all_pairs: 851,357 pairs / 24,500 classes / 13,552,707 cross-component
  D rows / EIU 10;
- Klaeger all_pairs: 63,823 pairs / 11,910 classes / 334,822 cross-component
  D rows / EIU 10.

Strict MMP layers remain small after censoring (Davis 7 classes / EIU 7;
Metz 480 classes / EIU 9; Klaeger 37 classes / EIU 9). Therefore, when a
passing W0-P panel exists, the W1 screen should be preregistered on the
all-pairs / similarity layer, with strict MMP retained as the high-explainability
confirmation layer — not as the default primary surface.

The NO-GO verdict is unchanged: W0-P has not passed and no W1 training is
authorized.
