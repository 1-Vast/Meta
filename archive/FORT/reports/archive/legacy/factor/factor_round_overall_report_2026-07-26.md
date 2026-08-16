# FACTOR route overall report before reopening exploration

Date: 2026-07-26

## What was tested

The user's FACTOR proposal was treated as one unified route, not as an agent-proposed candidate. Its
first mandatory gate was a label-blind F0 audit across three public development sources:
KIRHub 2026, Reinecke 2024 and the Christmann2016 slice of Papyrus/Kinase200. Protein anchors were
KLIFS-aligned residue identity/physicochemical states; ligand carriers were count-Morgan local
environments, BRICS fragments and localized pharmacophore tokens. Family, scaffold and source were
held out jointly for primitive-coverage accounting.

## Overall outcome

The original **discrete ligand-adapter gate** failed before model training. This is not evidence that
the complete FACTOR mechanism has failed. The precise failure type is
`F0_DISCRETE_LIGAND_CARRIER_COVERAGE_FAIL`.

The positive result is specific: the protein side is not the current bottleneck. Corrected KLIFS
mapping is 98.73%; weighted anchor coverage has median 1.000 and q10 0.956. All 30 endpoint-specific
environments form one primitive graph spanning all three sources. Equal-source weighting limits each
source to 33.33%, and the optimistic family×scaffold×source design envelope is MDE80 0.0226.

The blocking result is also specific: coverage under **exact count-Morgan/BRICS token identity** is
insufficient. Its weighted median is 0.801 against the frozen 0.90 gate, and q10 is 0.604 against
0.70. There are 10,922 discrete-novel-primitive observations (20.11%). The audit cannot determine
whether these are truly unsupported chemical functions or functionally similar environments split by
an ontology that is too exact.

## Why it failed

1. The three releases share many kinase-pocket states, but exact fragment identities differ strongly
   across their ligand libraries.
2. A fully connected document graph is too weak a sufficiency condition: a few common anchor and
   carrier tokens connect all environments while rare, observation-important carrier tokens remain
   uncovered.
3. KIRHub and Reinecke each contribute only one dense experimental environment. The analytic grouped
   MDE is therefore optimistic and does not by itself establish out-of-source uncertainty.
4. Adding interaction-model capacity or private residuals cannot decide whether the gap is real or a
   discrete-token aliasing artifact.

## Frozen decision and next exploration

`FACTOR_F0_FAIL_STOP_BEFORE_F1`, failure type
`F0_DISCRETE_LIGAND_CARRIER_COVERAGE_FAIL`. No original F1 model, Transformer, posterior, extra seed
or threshold change is authorized.

Exploration reopens as a new preregistered `F0-C0` representation audit: fixed continuous atom
environments, pharmacophore roles and pooled motifs, calibrated only on train-source pseudo-OOD
scaffolds. The external gates remain 0.90/0.70 and must be accompanied by chemistry-broken decoys,
structural reconstruction proxies and an anti-collapse rank check. Only if continuous F0-C still
fails should the result be attributed to real public-data chemical support and trigger a search for a
new paired panel or the existing prospective factorial design. Existing ChEMBL confirmation remains
permanently quarantined, and Papyrus-ChEMBL31 remains excluded.

## FACTOR-C completion: fixed and predictive continuous carriers

The reopened continuous-carrier branch is now complete. It does not change the precise original F0
classification above.

### F0-C0 fixed continuous carrier

F0-C0 replaced exact tokens with fixed atom, pharmacophore-pair and pooled BRICS vectors. Structural
reconstruction, effective rank and the true-minus-decoy gap passed, but common neutral atoms were
often unchanged by the nominal charge/bond permutation. Atom false coverage was 0.294 to 0.334, so
the anti-cheat control was non-identifying. Inner scaffold-OOD coverage also failed. Frozen verdict:
`FACTOR_F0C0_REPRESENTATION_UNIDENTIFIED_STOP`; this was not a real-support verdict.

### F0-C1 predictive self-supervised carrier

The only preregistered repair was a four-layer, width-128 masked GINE trained on strict
source-disjoint, connectivity-disjoint and scaffold-disjoint molecular graphs. It read no activity,
target, protein or source feature. All three fixed 40-epoch folds completed on CUDA.

The learned representation contains real local-chemistry information:

- 5-NN pharmacophore-role macro-F1 0.753--0.804, margins +0.724--+0.776;
- bond-order accuracy 1.000 and attachment margins +0.347--+0.390;
- changed-decoy fractions 0.9945--1.000;
- source-stratified true-minus-decoy coverage +0.5091, LCB95 +0.5004.

It nevertheless failed its representation-identifiability gates:

- participation rank 8.15--10.05 versus >=16, ratio 0.064--0.078 versus >=0.10;
- atom false-decoy coverage 0.0660--0.0748 versus <=0.05;
- inner scaffold-OOD median 0.7735--0.8203 versus >=0.85;
- source-balanced external median/q10 0.6710/0.4789 versus 0.90/0.70.

One Reinecke-held element-reconstruction margin was 0.14969 versus 0.15, but this near-boundary
failure is not decision-dominant: rank, calibration and inner coverage fail in every fold.
Frozen verdict: `FACTOR_F0C1_REPRESENTATION_UNIDENTIFIED_STOP`. F1-C is not authorized.

## Whole-route judgment before new exploration

The route has localized, not solved, the bottleneck:

1. protein anchor coverage is already adequate;
2. exact discrete ligand identity is too brittle;
3. fixed continuous carriers have invalid atom decoys;
4. the small predictive encoder learns chemistry and separates broken controls, but collapses the
   atlas into roughly 8--10 effective directions and still does not transfer across scaffolds.

Therefore the evidence does **not** authorize either “FACTOR works” or “public chemical support is
truly absent.” The correct category is **② `SIGNAL_PRESENT_EVIDENCE_INSUFFICIENT`**. The true-decoy
and reconstruction results are credible signal; the representation and calibration gates prevent a
support claim.

No more width, epoch, seed, threshold or pair-cap changes are allowed on this FACTOR-C instance.
Memory is not a scientific failure condition: the 96-pair cap was only an audit bound, and any
future passed mechanism must be reproduced uncapped on a larger machine before full training.

Exploration may reopen only after this report, as a substantively new mechanism rather than another
carrier-encoder tuning pass. The existing ChEMBL confirmation partition remains permanently
quarantined (`project_historical_confirmation_labels_read=true`), the current FACTOR runs did not
read it, and `sealed_test_consumed=false`.
