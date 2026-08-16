# Reinecke 2024 development-panel registration and power decision

Date: 2026-07-26  
Status: frozen before scoring any new target-conditioned candidate

## Decision

The Reinecke et al. kinase panel is admitted as a **development panel only**. It
cannot support an independent-confirmation claim because aggregate label shape
and resolution were inspected before the registry split. Confirmation labels and
the sealed test remain unread.

The registered data are pinned by:

- registry SHA-256:
  `92ca4981efb87ee42ae03fd4f6d2b2a0fc11ee90a678921a1c67d6b0c2768405`
- ligand-feature SHA-256:
  `30253e858eed32bafb318fdc22aec87f7e1253f57022401822ce85a6e1e39dc2`
- endpoint: `pKd_app = 9 - log10(Kd_app,nM)`

## Frozen leakage firewalls

1. Evaluation targets belong to combined sequence-homology components absent
   from the historical ChEMBL-37 and Metz substrates.
2. Development ligand parent connectivity and Bemis–Murcko scaffolds are absent
   from historical affinity substrates and registered anchors.
3. Exact Morgan Tanimoto is below 0.95 versus historical ligands and anchors.
4. The complete held target homology component is excluded from model fitting
   in each evaluation fold.
5. Query rows and the five-fold component map are fixed before candidate
   training.

The registered development set has 826 cells, 109 targets, 104 homology
components, 171 ligands, and 121 scaffolds. There is zero train/development
scaffold or parent-connectivity overlap. Seventy-seven targets have at least
four query ligands.

## Arm-blind power result

The power audit trained and scored only the ligand-only B0 arm for four seeds
and therefore revealed no result for a target-conditioned candidate.

| Seed | B0 target-macro Spearman | B0 RMSE |
|---:|---:|---:|
| 1729 | 0.1028 | 0.4854 |
| 2027 | 0.0575 | 0.4902 |
| 4241 | 0.0434 | 0.4877 |
| 5501 | 0.0623 | 0.4926 |

Across 80 scored homology components, the empirical retraining-noise standard
deviation was 0.2930 Spearman and the paired 80% minimum detectable effect was
+0.0668. The median number of query ligands per target was only five.

## Frozen performance gate

The primary endpoint is paired improvement over B0 in target-macro Spearman on
the fixed development queries, with the sequence-homology component as the
statistical unit.

The candidate passes only if all of the following hold:

1. mean paired improvement is at least **+0.0668 Spearman**;
2. the component-paired uncertainty interval excludes zero;
3. RMSE does not materially regress;
4. correct support labels outperform support-label permutation and wrong-target
   support;
5. the protein-conditioned mechanism outperforms shuffled/random-protein and a
   matched-capacity protein-independent control;
6. provenance checks report zero firewall violations.

The `+0.0668` threshold is `max(+0.03, empirical MDE80)` and cannot be loosened
after observing a candidate. Results below this threshold may be reported as
reproducible failures or weak signals, but not as a performance success.

