# FACTOR-U U1-L decision

Date: 2026-07-26  
Verdict: `FACTOR_U1L_REPRESENTATION_UNIDENTIFIED_STOP`

## Outcome

The local strict-corpus proof rebuilt 138,805 globally firewalled molecules, selected exactly 50,000
unique scaffold keys by the frozen SHA rule and produced 50,000 valid graphs. The split was exactly
45,000 train / 5,000 validation. ChEMBL confirmation rows materialized = 0.

Public-unlabeled structural prediction passed:

- masked element margin +0.1757;
- directed bond margin +0.5112, accuracy 1.000;
- BRICS attachment margin +0.3967;
- all losses and gradients finite.

Every evaluation fold retained strong atom-role information (5-NN macro-F1 0.723--0.821, margin
+0.694--+0.792), all role and chemistry decoys changed, and environment mismatch changed more than
99.4%. The source-balanced true-minus-decoy gap remained large at +0.4845, LCB95 +0.4744.

The decisive defects persisted:

- participation rank 8.77--8.91 versus >=16; ratio 0.0685--0.0696 versus >=0.10;
- atom false-decoy coverage 0.0670--0.0742 versus <=0.05, while pair/motif levels passed;
- inner scaffold-OOD medians 0.7345--0.7947 versus >=0.85;
- source-balanced external median/q10 0.6231/0.4078 versus 0.90/0.70.

Compared with F0-C1, increasing the structurally balanced unlabeled corpus from roughly two thousand
to fifty thousand molecules did not raise effective rank or repair atom calibration; external
coverage decreased from 0.6710/0.4789 to 0.6231/0.4078. This falsifies the local proof hypothesis
that corpus scale alone is the missing factor for this ligand-only predictive carrier.

## Decision boundary

Do not run the 138,805-molecule/full-pair replication on the larger machine. The local model did not
prove effective, so the user's resource escalation condition was not reached. Do not add a variance
penalty, larger encoder, more epochs or another seed to U1-L; those would be new mechanisms selected
after observing the failure.

F1-C remains locked. Low external coverage is still not a real-support verdict because the
representation gate failed first.

Current-run activity/affinity/protein fields read = false; current-run confirmation labels read =
false; project historical confirmation labels read = true; existing confirmation remains
quarantined; sealed test consumed = false.
