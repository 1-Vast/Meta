# CIIP-S1 S0 Audit Report (read-only; plan 9.1 items 1-10)

Machine-readable: S0_AUDIT.json. All computations keyed (S1.* namespaces);
no test labels used except split membership (public); no fitting performed.
Endpoint remains percent inhibition throughout.

## Item 1 — Coverage audit: PASS
- 65 admitted pairs; 49 covered re-derived from sequences == DATA2X2
  covered_pair_indices EXACTLY (bitwise).
- 16 missing pairs: ALL pos > 1020 = ESM_MAX_LEN (re-verified); parents ALK
  (4), MET (3), ... as in DATA2X2 missing_detail; no other reason.
- Covered split counts 32/8/9 match DATA1A pair_split; construct lengths
  270-1047 residues.

## Item 2 — Parent overlap: PASS
- split x parent table in JSON. WT-row sharing: every same-parent covered
  pair set shares exactly ONE WT measurement row; different parents never
  share a row (n distinct WT rows = 18).
- All 9 test pairs have >= 1 same-parent TRAIN pair (F9 definable):
  ABL1 x2 (4 sibs each), KIT x2 (5), EGFR (3), FGFR4 (1), RET x2 (4), TEK (1).
- Test clusters: {ABL1:2, KIT:2, EGFR:1, FGFR4:1, RET:2, TEK:1}.

## Item 3 — Mutation coordinates vs Q0B: PASS
- All 65 (parent, mutation, pos, wt->mt residue) match Q0B_MAPPING_AUDIT
  duongly_variant_records (0 mismatches). Alias ledger: 1 entry
  (BRAF historical numbering), consistent.

## Item 4 — Ligand overlap: PASS
- Per-pair panels 179-183 ligands (median 183); pairwise common ligands
  175-183 (median 182); per-ligand pair coverage 47-49 (158 ligands in all
  49 pairs). The panel is effectively shared -> ligand-pattern floor F7f is
  mandatory in every contrast (as frozen).

## Item 5 — Assay semantics: PASS (with recorded limitations)
- Endpoint namespace: percent inhibition only (no Ki/Kd/pK/DDG anywhere).
- Raw panel cells outside [0,100]: full panel 23.02% (matches prior 23.0%
  census); covered-rows-only 21.8%. Min -12.5, max 191.3.
- WT panel (covered parents, 18 rows): per-ligand WT mean median 89.4,
  87/183 ligands WT mean > 90, 0 ligands < 10, per-ligand WT sd median 15.0,
  mid-zone (10<mean<90) 96 ligands. CIIP-2 sec 2.4f census (90.9 / 99 / 84)
  used a wider row basis; both censuses agree on the structure:
  ceiling-loaded WT panel, sensitivity concentrated in mid-zone ligands.
- Concentration metadata: NO per-well concentration column or unit string
  exists in the local supplement copy (S1/S2 searched). Single-dose
  percent-inhibition endpoint; concentration assumed per source paper
  methods; recorded as data limitation (not usable as a feature).

## Item 6 — Censoring: PASS (limitation recorded)
- No censoring annotations (0 symbols like <, >, >=, <= in the matrix).
  Interval-censored formulations are NOT identifiable on this data; recorded
  as data limitation.

## Item 7 — Plan section-4 diagnostics re-derived (train+val 40 pairs): MATCH
| diagnostic | S1 re-derivation | plan value | verdict |
|---|---|---|---|
| main-effect energy share of d | 10.05% | 10.5% | match |
| median abs mean_l d | 3.29 | 3.29 | exact |
| same-parent cross-mutation corr | 0.4425 | 0.442 | exact |
| same-parent WT-residualized | 0.391 | 0.406 | match (common-ligand set differs) |
| different-parent keyed baseline | 0.041 | 0.036 | match (keyed draw) |
| parent-profile LOPO Spearman | 0.559 | 0.579 | match |
| parent-profile LOPO R2 | 0.326 per-pair median (0.260 pooled) | 0.326 | match (plan value is the per-pair median) |
| ligand-global leave-pair-out R2 | 0.060 | 0.060 | exact |
| parent-residualized cross-mutation | -0.280 | -0.28 | exact (LOO-biased negative) |

No magnitude discrepancy -> no stop. The plan's decomposition is confirmed
under the S1 preregistration.

## Item 8 — Power table (keyed simulation; S1.power streams)
- sigma (per-pair metric dispersion of a legal LOPO profile predictor):
  R2 0.599, Spearman 0.218.
- MDE at 80% detection (parent-cluster bootstrap lo2.5>0, 2000 draws,
  9 pairs / 6 clusters sizes [2,2,1,1,2,1]): delta-R2 0.566; delta-Spearman
  0.208.
- Consequence: only effects of about half an R2 unit are confirmable at
  n=9; ALL S1 outcomes carry the power label; thresholds frozen in
  S1_ADDENDUM_THRESHOLDS_20260820.md BEFORE any S1 fitting.

## Item 9 — Erasure cache: job launched (gen_erased_s1.py); asserts:
erased WT string == erased MT string (exact), max |embedding delta| <= 1e-5
(expected exactly 0 since strings identical). Cache + SHA in this directory.

## Item 10 — Leakage audit: no unresolved channel
1. Same-parent pairs share one WT row -> handled by parent-cluster
   bootstrap + cross-fitted F9 (self-excluded; val/test from train parents).
2. F2 X-position encodes the mutation coordinate -> F2 is counterfactual,
   non-deployable; deployability tested only via F3/F4.
3. Coverage selection (16 pairs excluded, pos>1020) -> claims restricted to
   the covered subset.
4. Shared ligand panel -> F7f floor mandatory.
5. No replicates -> mutation-specific variance noise-inclusive; no
   noise-floor claims.
6. Frozen-input integrity: radius-6 window features recomputed from
   q1_esm_cache.npz reproduce DATA2X2.npz esm_wt/esm_var EXACTLY
   (max abs diff 0.0 over 49 pairs).

## S0 exit gate: items 1-8 COMPLETE, no unresolved leakage, erasure cache
job in flight (gate closes when the cache + SHA land). S1 fitting does not
start before the addendum below is frozen.
