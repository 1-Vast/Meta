# P1R2B-D0/D1 Release-Pinned Affinity Corpus Report

Updated: 2026-08-06

Decision: `D0-C_PASS; D1_PASS; E0_DATA_READY_AWAITING_AUTHORIZATION;
T/P2-P4_FROZEN; RECIPIENT_LABEL_READS=0`.

## FACT

The official ChEMBL37 SQLite archive contains `5,764,252,857` bytes and
matches the registered official SHA-256
`33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`.
The extracted SQLite database contains `30,480,314,368` bytes, has SHA-256
`4be13df3b68e25dcd0bff44bf094033b5aebe98f415acdc8c1cdf380e0c15142`,
and schema SHA-256
`fb9d838db717e011cc031f6acb24665c309c973d17f854fcff86382af1f8f972`.
The official license file is recorded as CC BY-SA 3.0.

The frozen whole-release SQL selected exact Ki/Kd point measurements from
binding assays with confidence 9, single-protein targets, curated standard
values, no validity warning or potential duplicate, an exact component or
variant sequence, canonical parent structure, and document provenance.
`EnergyPilot.v1` contains `343,562` canonical rows and zero RDKit rejections.
It contains `41,619` target x endpoint x assay x context tasks; `4,092` meet
the preregistered minimum of 20 exact compounds and 12 non-tied comparisons.

Independent verification reconstructed every sequence/context/task hash,
pAffinity conversion, row count and task manifest. Canonical rows SHA-256 is
`c8b2762683a6ce8abbf7c0722a11eefaca4b5d042270f9167963ec308efe9a29`;
task manifest SHA-256 is
`2d804e066bc36855ef6799ddbb8ad9546375e0d3281b80c86d1d153377430f33`.

D1 compared 5,188 candidate source sequences against 159 protected DAVIS
metaval/recipient sequences using MMseqs2 candidate search and the registered
parasail 40% local-identity confirmation. It excluded 401 candidate sequences,
of which 131 occur in E0-Core tasks. The governed corpus contains 3,817 tasks,
697 targets, 3,504 Ki tasks, 313 Kd tasks, 334 retained homology components,
378 document components, and 253 homology-union-document closure components.

The largest closure component contains 1,467 tasks (38.43%). The deterministic
five-fold task counts are `1467/588/588/587/587`. Homology, document and union
closure straddling counts are all zero. Split SHA-256 is
`1e0bf9d6bbd287782d861825278968260980fb615764bc07d57450b309cc6689`.

The historical crosswalk found all old 200 assay IDs in static ChEMBL37, but
only 58 contain candidates under the new Ki/Kd E0-Core contract. This audit is
non-gating and does not alter the historical F0R failure.

## INFERENCE

The live API was the wrong immutable provenance boundary. Static-release raw
hashing plus frozen SQL, normalizer, row schema, governance and split hashing
provides scientific reproducibility without claiming byte recovery of the old
API JSON serialization. D0-C and D1 therefore pass their registered data and
governance contracts.

The 253 closure components support component-closed source OOF, but the 38.43%
largest component makes task-level independence and balanced-fold assumptions
invalid. E0 must report closure-component uncertainty and the fixed fold
imbalance; it may not resplit by target or task to improve metrics.

## UNTESTED HYPOTHESIS

Frozen P1B contact/distance geometry followed by a local-before-pooling energy
map may add correct-protein affinity direction beyond an OOF ligand prior on
this corpus. No energy model was trained and no DAVIS affinity Gate was read in
D0/D1. E0 requires a new explicit authorization. Typed interactions and P2-P4
remain frozen.

## Frozen E0 Data Gate

Before any E0 model result, the data sufficiency floor is frozen at
2,500 governed tasks, 500 retained targets, 200 closure components, and 500
tasks in every outer OOF fold. The realized corpus passes at
3,817 / 697 / 253 / minimum-fold 587. The largest-component fraction is a
mandatory estimand limitation, not a tunable exclusion criterion.

Key machine reports:

- release manifest: `be53dd63672d80e15e4f80790f55dd00aee1d5d78d5f3e1fd52bcd4bb1499043`
- corpus manifest: `5badb7b7ebb66591d696a486420539bd182dfe6758b37193474661e620d479bb`
- corpus verification: `a3717580cca65556a733ef60bd0bdbec0ad0fc5f08b1a14c1a88301de5a5a025`
- governance manifest: `56e2ace6e4c1e1956a8b3a5c2ee8908c255d316c0139b364c2a98f86c4b27c01`
- D1 audit: `ed89d0d8cb55bb9a4148f0eb0e8da1a7d58a8b22ffa8b6a7eaa8a206d339ae7f`
- legacy crosswalk: `1f9a8af23507ea5efe59cd57c6394881f179f82141545fd973faf3ecf740118d`

All D0/D1 machine artifacts record `training_authorized=false` and
`recipient_labels_read=false`.
