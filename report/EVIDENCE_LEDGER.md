# Evidence ledger

Compact authority after the 2026-08-16 cleanup. Deleted pre-R0 evidence is
summarized in `meta_fewshot/LEGACY_PRE_R0_SUMMARY.md` and recoverable from Git.

| Cycle | Retained evidence | Decision |
|---|---|---|
| Legacy through 2026-08-14 | `meta_fewshot/LEGACY_PRE_R0_SUMMARY.md` | Solver, HyperSAR, D-MEMT, CIPF/TERM and K3/ELMT families closed |
| R0 | `meta_fewshot/stageR0_retrieval_falsification_20260815/` | Earlier k0 retrieval claim does not generalize to exact-free cold ligands; selection/transduction caveats established |
| R1 | `meta_fewshot/stageR1_double_cold_split_20260815/` | Governed double-cold split and meta_test seal established (logical exclusion; see the standing corrections) |
| R2 | `meta_fewshot/stageR2_representation_discriminator_20260815/` | Representation identified as the next causal axis |
| R3/R4 | `meta_fewshot/stageR3R4_level_shape_20260815/` | A0/B3 level-shape ladder; real k0 MSE/CI trade-off |
| R5 | `meta_fewshot/stageR5_reltransport_20260816/` | Donor, whitening, gradient and seal contracts repaired |
| R6/R7 | `meta_fewshot/stageR6_reltransport_screening_20260816/`, `stageR7_reltransport_3seed_20260816/` | Learned query-specific transport gates rejected as inert/harmful |
| R8 | `meta_fewshot/stageR8_stronger_shape_20260816/` | Shape-first training yields real cliff signal but fails joint promotion |
| R9 | `meta_fewshot/stageR9_cliffweight_20260816/` | C2 joins Pareto set; no dominant model |
| R10 | `meta_fewshot/stageR10_variance_20260816/` | Variance objective rejected |
| R11 | `meta_fewshot/stageR11_grammar_shape_20260816/` | Grammar-shape family rejected |
| R12 | `meta_fewshot/stageR12_margin_20260816/` | Margin ranking rejected |
| R13 | `meta_fewshot/stageR13_shape_direct_20260816/` | Direct shape family rejected; 18 gates: 16 pass, 2 expected failures |
| R14 | `meta_fewshot/stageR14_diagnostics_20260816/`, `stageR14_screening_20260816/` | A0 ordering floor best; aligned ListCE structurally valid but inert; loss-form axis closed |
| Boundary | `BOUNDARY_20260816.md`, `CURRENT_MODEL_EVIDENCE.md` | No excellent candidate; A0/B3/C2 Pareto, 0.782 cliff record is development-only |
| A2v1 | `tools/research/a2_readiness/` | consolidated into `SUPERSEDED.md`; two of four conclusions withdrawn. Retains the literature ledger and the CPC centering probes |
| A2v2 | `tools/research/a2_readiness_v2/` | governance incident repaired; noise/leakage audit; causal attention audit; Stage P's frozen design. **Its A2 verdict was an over-reach and is superseded by Stage R** |
| A2 exact (R, L2) | `tools/research/a2_exact_probe/` | **A2 closed on its own operator and gates**: loses to parameter-free Tanimoto at every k with resolved intervals; its protein and label controls are **inert**, not inverted (the earlier "inverted" reading came from non-nested banks and is withdrawn). `embed` carries a signed-SAR direction, +0.212 +/- 0.011, but it is **negative on activity cliffs** (-0.118) and near zero on novel ligands |
| Stage P | `tools/research/stageP_cpc/` | **P1 fails** (-0.0066 [-0.0545,+0.0417]). Correct and wrong protein give identical ordering at every k in both arms. The centered objective excluded the level branch as designed and made the protein response reproducible across seeds (+0.316) but unaligned with truth (+0.022). Objective-only protein conditioning is closed for this architecture and budget |
| M0 | `meta_fewshot/stageM0_msa_probe_20260816/PREREGISTRATION.md` | Independent diagnostic, not run |
| Next | `tools/research/a2_readiness_v2/PREREGISTRATION_V2.md` | Stage P: frozen, costed (~3.75 h), **not run**. The only untested protein-conditioning question. `NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md` is closed |
| Stage A/B | `tools/research/stageA_innerloop/`, `stageB_complementary/` | AdaMBind-inspired inner/outer-loop meta-adaptation NOT PROMISING then REJECTED; meta_val checkpoint selection measured at ~0.62 pK^2 (k=0); ligand representation collapse quantified |
| Stage C | `tools/research/stageC_level_shape/` | level/shape decomposition: k=0 is 68% level; boundary measured, not trained |
| Stage D0/D/E | `tools/research/stageD_level_panel/` | Five governing questions answered (D0_REPORT.md); panel-set level head + orthogonal routing REJECTED (G1/G2); attribution ablations complete; document-transfer R^2 +0.451 measured |
| Stage F | `tools/research/stageF_pairwise/` | pairwise learned transport REJECTED; fifth learned-kernel family to lose to fixed Tanimoto |
| Stage G/G2 | `tools/research/stageG_esm650/` | ESM-650M residue-input lane NOT CONFIRMED across 3 seeds |
| Stage H0 | `tools/research/stageH_pocket/` | pocket-prior lane rejected at identifiability gate (2.44 vs 2.62 constant; shuffle 2.49) |
| Stage I | `tools/research/stageI_lm/` | live ESM-150M LoRA lane REJECTED (no resolved MSE gain; two resolved ranking observations) |
| Stage J | `tools/research/stageJ_assay/` | assay-aware level head (journal/panel/protein) REJECTED: best k=0 level on record (1.30) but resolved k=2/3 ranking degradation; D0c journal probe 1.619 |
| Stage K/K2 | `tools/research/stageK_contrastive/` | K-REG: first all-k resolved MSE improvement across 3 seeds; k=0 centered gain does not survive pooling -> NOT CONFIRMED, nothing promoted |
| Stage L | `tools/research/stageL_gated/` | support-gated assay-aware level head: best k=0 calibration on record (2.0997) but resolved k=2/3/5 ranking degradation -> REJECTED; closes the level-head composition axis |
| Stage M0 | `tools/research/stageM_chemberta/` | ChemBERTa-77M ligand embeddings: ordering r +0.147 (below occupancy), level probe = grand mean -> ligand-LM family falsified; external-representation ledger complete |
| Stage N (audit) | `tools/research/stageN_audit/` | final boundary audit: all load-bearing numbers re-derive bitwise from raw rows; 104 seal artifacts audited, none evaluated; 7/7 stages preregistered |
| Stage P0 | `tools/research/stageP_go/` | ProteinKG25 GO annotations (313/387 matched): level probe 2.27 vs 1.43 constant -> protein-function-annotation family falsified; external-representation ledger complete |
| Stage Q | `tools/research/stageQ_frozenhead/` | decoupled frozen-feature level head: Q0 probe 1.3416 (best frozen predictor) but trained composition degrades k=0/2/3 ranking with resolved intervals -> the level/ranking conflict on one shared trunk is fundamental (4th composition to fail) |
| Boundary | `report/BOUNDARY_20260817_NIGHT.md` | k=0 <= 1.00 protocol-conditioned with the measured legal-input families; k=5 at 0.944-1.007 across seeds |

Only `RESULT.json`, formal reports, necessary prediction tables and loadable
admitted checkpoints are authoritative. Progress logs and short smokes are not
performance evidence.

## Standing corrections (2026-08-16)

These override any earlier phrasing wherever it survives:

1. **Every wrong-protein control in R0-R14 is uncentered and therefore measures
   *level* specificity.** R3R4's "first resolved protein specificity" (+0.4216
   at k=2) is a level result. The ordering version of the same control measures
   **−0.0002 [−0.0015, +0.0008]**.
2. **A0 has no protein-conditioned within-target ordering at k=0**, at any
   internal stage of the trunk, and the exact episodic A2 operator over those
   representations fails every performance and control gate — corrupting the
   protein or the support labels *improves* it. A random initialisation's
   protein sensitivity is undirected across seeds, so it is no evidence of
   latent capacity either.
2b. **A protein-independent signed-SAR direction does exist in `embed`**
   (Δ-r +0.270 on held-out components), is orthogonal to Tanimoto and strongest
   on activity cliffs. It is ligand-side transfer, not meta-learning, and the
   moment form cannot exploit it — the signal is pairwise.
3. **The protein path is exactly invariant to residue-slot permutation.** No
   result from this architecture may be described as pocket-aware,
   contact-resolved or biologically localized.
4. The `meta_test` seal was **opt-out** until 2026-08-16. No recorded number is
   affected (verified bit-identical), but the pre-repair claim of process-level
   purity is not supported.
5. The seal is **logical exclusion after parsing**, never physical isolation:
   the corpus is one all-label artifact, so every sealed label is decompressed
   and parsed on every load. `violations = 0` does not close the incident — two
   artifacts remain `process_unsealed` and `audit_research_record` exits
   non-zero while they do. See `GOVERNANCE_INCIDENT.md` and
   `SPLIT_ISOLATION_SPEC.md`.
