# RB-DR-QMAPD — Stage O1 oracle preregistration (2026-07-26)

Frozen before any O1 number is read. Decisive gate: does additional same-target evidence contain a
reproducible, scaffold-generalising SAR ranking signal beyond the same teacher at `k` shots? Design:
`rb_dr_qmapd_mathematical_design.md`; feasibility: `rb_dr_qmapd_power_audit.md`. Train + spent
development rows only; no confirmation/Davis/sealed label. A predictive number here is
`ARCHITECTURE_MECHANISM_RESULT_ONLY`.

## Frozen teacher (no post-hoc maximum)

**Teacher A**, the existing validated exact-Cholesky Bayesian residual posterior
(`model/descriptor_baseline.py` + `model/bayes_posterior.py`): base `b(d)` (heteroscedastic, ligand-
only), posterior design `[1, z(rep_d)]` (intercept = calibration, linear = ranking). Meta-trained per
CV fold on TRAIN rows with support sizes drawn uniformly in `{2,…,K}` (so it is competent at both `k`
and `K`), then frozen; the closed-form posterior is the identical architecture at every support size.
Teachers B (REML measurement-error) and C (deep-kernel GP) are registered alternatives, implemented
only if A is inconclusive. Features and all hyperparameters are frozen before O1 is scored.

## Protocol

* Primary pair **(k=4, K=16)**, 76 independent components; secondary robustness (reported, non-gating)
  **(k=4, K=32)**, 40 components.
* Five frozen leave-homology-component-out folds (`component_folds`, seed 1729). Base + teacher fit on
  TRAIN rows of the fold's train components; evaluated on the held fold's development targets.
* One label-blind nested episode per held target: greedy scaffold-distinct support `S_K`, fixed
  dual-cold query (scaffold-disjoint, Tanimoto<0.95). `S_k = S_K[:k]` (nested). A second independent
  support permutation (drawB) provides the support-resampling null for empirical MDE80.
* Arms on identical query rows: `B0` (no support), `T_k`, `T_K`. `B_k = B0`.
* `rho` = within-component (tunit) macro Spearman. Grouped component bootstrap (10,000 draws).
* `empirical_MDE80_info` from the paired null `T_k(drawA) − T_k(drawB)` (same size, no added info);
  `empirical_MDE80_total` from `T_K(drawA) − T_K(drawB)`. MDE80 via the frozen `mde_from_spread` rule.

## Frozen pass criteria (Stage O1) — a single failure ⇒ `RB_DR_QMAPD_ORACLE_INFORMATION_FAIL_STOP`

1. `Delta_info ≥ max(0.03, empirical_MDE80_info)`.
2. grouped paired-bootstrap `LCB95(Delta_info) > 0`.
3. `Delta_total ≥ max(0.05, empirical_MDE80_total)`.
4. grouped paired-bootstrap `LCB95(Delta_total) > 0`.
5. `Delta_info` positive in every preregistered major component stratum (low-similarity,
   scaffold-diverse).
6. `RMSE(T_K) ≤ 1.02 · RMSE(strongest baseline)`.
7. `Delta_info > 0` under held-target + held-scaffold evaluation (the scaffold-diverse stratum, which
   is already held-target by CV and held-scaffold by construction).
8. `Delta_info` not confined to high support–query similarity: `Delta_info > 0` in the low-similarity
   stratum.
9. the gain is not solely a teacher-architecture advantage at `k`: `Delta_info` (criteria 1–2) holds
   independently of `Delta_arch`.

If `Delta_arch > 0` but `Delta_info` fails, report/improve the k-shot baseline; do NOT implement
privileged completion distillation. Passing O1 authorises only Stage O2 (completion specificity). No
threshold, teacher, feature, or hyperparameter is changed after the result; three seeds only after a
one-seed pass of the full downstream gate (far in the future). Do not add data, a larger teacher, a
flow, or a protein controller as a rescue.
