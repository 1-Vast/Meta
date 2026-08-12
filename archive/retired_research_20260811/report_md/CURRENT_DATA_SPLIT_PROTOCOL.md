# Current Data Split Protocol

Version: `CORE_FEWSHOT_SPLIT_V2`  
Effective: 2026-08-05

This file is the authoritative split and exclusion policy for new MetaSieve
training, metaval and recipient evaluation. Earlier Gate Z1/Z2 manifests,
homology-component runs and published benchmark reproductions remain immutable
historical evidence under their registered protocols.

## Core few-shot legality

The core estimand is few-shot adaptation to an exact unseen target within a
single assay context. Every new run must satisfy all of the following:

1. **Exact target separation.** An exact target identifier or exact normalized
   protein sequence may occur in only one of source/train, metaval or recipient.
   All assay contexts for that exact target follow the same partition.
2. **Target x assay isolation.** A task is `target x assay/context x endpoint`.
   Its rows, labels and fitted transforms cannot straddle partitions. Ki, Kd,
   IC50 and EC50 remain separate endpoint families.
3. **Support/query disjointness.** Support and query are sampled without
   replacement. They share no measurement row, normalized compound identity or
   exact protein-ligand pair. Query labels are unavailable to adaptation.
4. **Exact duplicate exclusion.** Duplicate measurements and exact normalized
   protein-ligand complexes are assigned as one group; they cannot cross a
   partition or support/query boundary.
5. **Recipient-complex exclusion.** Source and metaval remove every exact
   protein-ligand complex found in the sealed recipient identity manifest.
   This comparison is label-blind and does not authorize recipient-label reads.

These are universal legality requirements. Failure of any item invalidates the
core run.

## Separately reported stress-test strata

The following are not universal exclusions from core training:

- protein homology novelty at the registered 40% identity threshold;
- ligand scaffold novelty;
- pocket novelty;
- protein-ligand interaction novelty.

Each must be computed without labels and reported as a separate stress-test
stratum when the required metadata exist. A sample may belong to multiple
strata. Reports must include the similarity definition, threshold, reference
partition, row/task counts and confidence interval per stratum. Missing pocket
or interaction annotations are reported as `not measurable`, not silently
excluded.

The stress strata answer harder generalization questions; they do not redefine
the core unseen-target estimand, reduce the core training set, or select models
using recipient labels.

## Published reproductions and historical gates

Published reproductions keep the paper's exact split, including CD-HIT 40%
novel-task or other cold-start variants, and are labeled external benchmark
evidence. They do not replace this protocol.

Gate Z1/Z2 and completed K1 runs keep their original homology/document closure
and aggregation rules. Their artifacts and failure conclusions are not
recomputed under V2. Any new training or dataset assembly after the effective
date uses V2 and records the protocol version in its manifest.

## Structural Teacher Governance

Mechanism pretraining uses a separate structural corpus and never reads DTA
affinity labels. Before the structural corpus is frozen:

1. Use sealed DAVIS/KIBA metaval and recipient protein sequences as a
   label-free blacklist.
2. Remove every structural receptor with confirmed sequence identity >=40% to
   a blacklisted protein.
3. Group the remaining structures by protein homology for structural
   train/validation/test splitting; random complex-row splitting is forbidden.
4. Record exact ligand connectivity and Murcko scaffold overlap with DTA
   partitions. The primary Pilot-10K reports these overlaps; a separate
   ultra-strict ablation may enforce scaffold-cold structure pretraining.
5. Record the similarity tool, version, threshold, reference sequences and
   manifest hashes. For corpus scale, use MMseqs2 candidate discovery followed
   by threshold confirmation instead of Python all-pairs comparison.

These exclusions govern leakage claims for the P1 structural teacher. They do
not expose recipient affinity values and do not change historical split
assignments.
