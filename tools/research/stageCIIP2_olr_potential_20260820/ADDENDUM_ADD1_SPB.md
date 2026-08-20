# Addendum ADD-1 to PREREGISTRATION.md (issued before any training)

Date: 2026-08-20. Prereg base SHA-256: a7b17e8a3a6300d1e02bad44233a07d70a91b0b46a7c2a5d4ccd5cdb89489912
Reason: data-contract correction discovered during implementation census; no
result has been produced in this stage yet, so this is a construction-rule
fix, not a results-driven amendment. The frozen base file is never edited.

## Correction 1: SPB is built on the ESM-covered subset

The base prereg listed 12 multi-pair parents including ALK, FLT3, LRRK2,
CMET/MET. Those parents have >= 2 pairs in the full 65-pair table but their
pairs are not ESM-covered (DATA2X2 covered subset = 49 pairs over 18
parents). The parent-disjoint split must be formed on the same covered
subset used by every other analysis (leakage parity with CIIP-1A).

Covered multi-pair parents (>= 2 covered pairs): KIT 9, ABL1 8, RET 8,
EGFR 5, FGFR4 3, FGFR3 2, PDGFRA 2, TEK 2 (8 parents, 39 pairs).

## Correction 2: SPB greedy direction

Base rule "rank by covered-pair count desc, greedy to TEST until >= 5
parents and >= 20% of pairs" would move KIT+ABL1+RET+EGFR+FGFR4 (33/49 =
67% of pairs) to TEST, starving train. Amended rule (frozen):

1. Consider only the 8 multi-pair covered parents.
2. Order them by covered-pair count ASCENDING, ties by name ascending.
3. Move parents to TEST in that order until TEST has >= 5 parents AND
   >= 10 covered pairs (20.4% of 49).
4. All remaining covered pairs stay available to TRAIN/VAL: VAL takes one
   pair from each single-pair TRAIN parent where parent has >1 pair else
   the largest TRAIN parent's pairs by name order until ~15% of remaining
   pairs; simplest frozen form: VAL = every 7th pair (by stable sort of
   (parent, mutation)) among TRAIN-side pairs; the rest are TRAIN.

Resulting frozen allocation (deterministic, verifiable):
TEST parents {FGFR3, PDGFRA, TEK, FGFR4, EGFR} = 14 pairs (28.6%);
TRAIN side {KIT 9, ABL1 8, RET 8} + 10 single-pair parents = 35 pairs;
VAL drawn from TRAIN side per rule 4.

All other rules of the base prereg are unchanged. This addendum is frozen
at issue.
