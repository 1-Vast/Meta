# R-MAON G0 multi-agent synthesis

**Date:** 2026-07-28  
**Accepted formal artifact:** `reports/active/rmaon_g0.json`  
**SHA-256:** `3276a854339cf9ebb7684df2295bb526681004f78847a8dda34dc6a03dacd4cc`  
**Environment:** `D:\anaconda\envs\drug\python.exe`, CUDA  
**Final category:** **2 — credible positive signal, evidence insufficient**

## 1. Agent judgments

The independent passes converged on the same bottleneck from different directions.

| role | core judgment |
| --- | --- |
| history and contradiction audit | Earlier positive pretraining, denoising, residue-teacher, physical-field, Bayesian/meta, support-adaptation, and protein-prior signals either failed transfer, were not load-bearing, were underpowered, or were produced by an invalid protocol. A model-only retry without new information is closed. |
| data and measurement design | No local or staged public source simultaneously provides multi-family independence, at least about 423 powered components, at least 40 scaffold-diverse query ligands per target, retained inactives, and independent provenance. A prospective randomized panel is the missing object. |
| biology and medicinal chemistry | Transfer requires a target-side mechanism variable whose states recur within and across protein families and change ligand ordering. Sequence taxonomy, pooled embeddings, pocket composition, contacts, and ordinary structure priors have not established this. Mechanism-orthogonal roster design is more valuable than another encoder. |
| mathematics and statistics | K-LBP R3 used `Theta = gamma a c^T`; at `gamma=0`, `a,c` are unidentified and its Fisher information is singular. Direct `Theta` has a unique null and admits an operator-norm score without an optimizer or `V_t` inversion. |
| architecture and training | The smallest defensible model is a standard compact sequence encoder, standard ligand encoder, strong ligand-only/shared-global head, and a directly trainable centered operator. Larger Transformer/Mamba/GNN/MoE capacity introduces no identifying information. |
| adversarial review | A synthetic rank-one pass is conditional estimator evidence, not affinity performance. The coefficient-level gate has no ligand-scaffold score axis, the target coordinate is kinase-only, low-fidelity assay monotonicity is untested, and G0-A independently blocks any real-label model stage. |

## 2. Evidence map

### Signal absent versus signal not transferable

* The Metz panel contains structured target-ligand residual interaction: the historical PA4 held-out
  completion correlation was about `0.408`. This rejects a global "no interaction exists" claim.
* AMOB's staged low-fidelity diagnostic produced target-specific residual order correlation
  `rho=+0.434`, but the files have no license or assay/document IDs and overlap protected accessions.
  RECRO L0 later showed that nominal cross-document residual signal collapses under provenance-family
  isolation. AMOB is therefore evidence of a candidate signal, not certified independent biology.
* Pooled protein features, aligned-pocket composition, mechanism priors, structure/contact models,
  posterior adaptation, and pretraining repeatedly failed to carry a material interaction gain to
  unseen targets under matched destructive controls.

### Power and protocol corrections

* The old `+0.093` few-shot claim was retracted after matched-query and split auditing.
* K-LBP R1's first artifact was invalidated because it materialized labels before dropping them.
* K-LBP R3 stopped in its first null replicate because the factorized null was non-regular; that is an
  estimator no-decision, not evidence that the mechanism effect is zero.
* R-MAON's first two concurrent pre-A1 writes were invalidated. The first materialized Metz development
  affinity through `DualCold.panel()`; the second ran before the repaired protocol was explicitly
  amended. Neither is interpreted.

## 3. Public data audit

| source | license / provenance | useful structure | binding failure |
| --- | --- | --- | --- |
| Metz TRAIN / ChEMBL 37 | ChEMBL CC BY-SA; one panel document | 101 kinase components, dense pKi TRAIN geometry | kinase-only; no independent predictive confirmation lineage |
| Reinecke 2024 | CC BY 4.0, hashed supplements | dense continuous Kinobeads pKd; 77--104 kinase components depending on firewall | within-kinome and previously spent for development; not a 423-component multi-family panel |
| Novartis SPD 2023 | CC BY 4.0 data, MIT code | 101 genes across multiple top-level families; systematic, 91.3% tested negatives retained | median 14 compounds per gene-assay; only 30/144 assays have at least 40; floor/censor dominated |
| Papyrus 05.7++ | public aggregate with source metadata | broad targets and chemistry | exactly one aggregate row per parent-target; zero independent document-replicated cells |
| native BindingDB / ChEMBL raw | open licenses, provenance recoverable in principle | broad observational graph | no demonstrated powered provenance-independent factorial panel after all firewalls |
| Harmonic activity-integration staging | empty LICENSE; no assay/document IDs | single-dose plus exact activity; strong engineering diagnostic | not admissible for training; overlaps development and sealed accessions |
| prospective R-MAON panel | absent | would provide randomized ordinal and independent dose-response views | `manifests/rmaon_prospective_panel.v1.json` does not exist |

The conclusion is not that open databases are small. They are large but observationally aggregated:
target-specific reordering, selection, assay, and provenance cannot be separated at the required
independent unit and query depth.

## 4. Candidate generation and screening

| rank | candidate | genuinely new condition | cheapest result | disposition |
| ---: | --- | --- | --- | --- |
| 1 | R-MAON direct operator plus randomized assay-monotone supervision | prospective randomized, inactive-retaining low-fidelity order plus independent dose-response calibration | G0-A topology + G0-B regular-null score | I2 passes; G0-A stops; I1 untested |
| 2 | MC-ORTHO roster with SENTRY profile blocks and TWIN-SITE replication | mechanism states crossed within/across families, immutable binding-profile blocks, two independent sites | label-blind rank/power audit, then 12-target two-site pilot | best measurement design; no manifest/material yet |
| 3 | CPRV-DTA replicated-view bottleneck | the same cell measured by genuinely independent laboratories | cross-site target-centered residual-rank reliability | conditional runner-up; roughly doubles measurement cost |
| 4 | LEXOR-MC fixed qualitative channels | source-bound non-affinity mechanism evidence | MC0--MC2 extraction/redundancy/compatibility gates | separate Task B only; cannot create labels or components |
| reject | larger attention/Mamba/GNN/MoE, pose/pocket rescue, free RAG/LLM, support meta-learning, causal completion without an instrument | none | historical destruction and information audits | closed unless a new information condition appears |

The candidates are not stacked. Randomization, profile blocking, and site replication are measurement
controls, not extra model innovations.

## 5. Accepted G0 experiment

### Protocol

The accepted A1 run used:

* seed 1729;
* five regimes, 200 replicates each;
* 999 homology-component Rademacher multipliers per replicate;
* the valid, label-guarded R1 deterministic coordinate;
* reader-filtered Metz TRAIN coefficient/covariance geometry;
* exact PSD eigensquare-roots of the empirical, sometimes rank-deficient `V_t`;
* no covariance inversion, no Cholesky jitter, and no optimizer at the null.

A concurrent post-run writer serialized the same deterministic A1 result again. A fieldwise comparison
found every scientific field and all 1,000 replicate records identical; only `wall_clock_s` differed.
The stable hash above is for the final A1-bound serialization.

Focused amended verification returned `38 passed`. It checks the TRAIN-only loader boundary, exact PSD
reconstruction, component balancing/whitening, unique zero operator, bit-identical shared-global
fallback, attached/detached gradients, regular-null score, and debiased recovery.

### Results

| regime | reject / 200 | rate | exact 95% binomial interval | median recovery ratio | frozen gate |
| --- | ---: | ---: | --- | ---: | --- |
| S1 null | 11 | 0.055 | [0.0278, 0.0963] | — | pass (`<=0.10`) |
| S1 active | 200 | 1.000 | [0.9817, 1.0000] | 1.0027 | pass |
| S2 heteroscedastic | 200 | 1.000 | [0.9817, 1.0000] | 0.9993 | pass |
| S3 degenerate | 199 | 0.995 | [0.9725, 0.9999] | 0.9610 | pass |
| S5 wrong coordinate | 9 | 0.045 | [0.0208, 0.0837] | — | pass (`<=0.10`) |

All six G0-B gates pass:

```text
RMAON_G0_NULL_SCORE_AND_RECOVERY_PASS__MODULE_ONLY
```

The null and wrong-coordinate exact upper bounds are both below 0.10. Active-regime recovery is close
to one; even S3's 2.5th--97.5th recovery-ratio quantiles are about `[0.885, 1.037]`.

G0-A fails because the prospective manifest is absent:

```text
RMAON_G0_TOPOLOGY_OR_POWER_STOP
```

## 6. What improved relative to the strongest baseline

The strongest reproducible deployable model remains the ligand-only/shared-global baseline. No
real-affinity R-MAON comparison was authorized, so its Spearman, RMSE, calibration, and paired gain
relative to that baseline are **not estimated**.

The demonstrated improvement is narrower:

* old factorized K-LBP estimator: could not produce a valid null decision because 2/5 cross-fit
  direction fits failed in the first formal null replicate;
* direct R-MAON operator: finite non-iterative null score, 0.055 false rejection, at least 0.995 power
  in every active stress regime, and calibrated magnitude recovery.

This is an estimator-module advance, not a predictive-performance advance.

## 7. Innovation accounting

### I1 — randomized assay-monotone multi-fidelity supervision

Biological hypothesis: within a fixed target and assay, low-fidelity single-dose order preserves
target-specific ligand reordering outside saturation; independent dose-response measurements calibrate
the latent order to pKi or pKd.

Status: **not tested**. No compliant randomized, licensed, provenance-independent substrate exists.

### I2 — regular-at-null direct target-chemical operator

Statistical hypothesis: a directly parameterized centered `Theta` removes the alternative-only
directions that made `gamma a c^T` singular at zero and can detect/recover a true operator under
empirical heteroscedastic noise.

Status: **passes G0-B conditionally**. Zero-operator fallback and gradient dependence are verified in
the trainable module; the score/recovery experiment passes every frozen synthetic gate.

## 8. Alternative explanations not removed

1. The synthetic truth is rank one and aligned with the supplied coordinate. Passing proves conditional
   sensitivity and calibration, not that real affinity uses that coordinate.
2. The empirical `V_t` is a TRAIN-derived plug-in covariance. G0-B does not validate its calibration as
   a measurement model.
3. Coefficient compression removes ligand-scaffold score contributions. G0-B uses an honest one-way
   component bootstrap; future edge-level inference still requires component-by-scaffold dependence.
4. The surviving R1 coordinate is audited only within kinases. It does not identify cross-family
   mechanism states.
5. No deep encoder was trained and no assay-monotone loss was evaluated on real measurements.
6. No independent site, binding-profile sentinel, or prospective randomized tranche exists.

These are reasons to stop at G0, not reasons to weaken the gate.

## 9. Root cause and next admissible condition

The most credible root cause is **measurement topology and target-side information identifiability**,
not optimizer capacity. A valid next stage requires a frozen, label-blind prospective manifest with:

* a 70--155-component multi-family mechanism pilot and a credible path to at least 423 independent
  predictive components;
* at least 40 randomized scaffold-diverse query ligands per target;
* inactive and censored outcomes retained with known inclusion probabilities;
* binding-profile blocks isolated from model features;
* at least two genuinely independent provenance lineages;
* raw curves, exact WT construct and endpoint separation;
* `PA2 >= 0.5 pK` and at least 80% power at the 0.03 paired ranking floor.

The first paid gate remains a small blinded 12-target, at least six-family, 16-ligand, two-site
reliability pilot. It must be preregistered separately and cannot be inferred from G0-B.

## 10. Final conclusion

**Category 2: credible positive signal, evidence insufficient.**

The direct operator is a load-bearing, regular estimator module under the intended synthetic mechanism.
There is still no leakage-safe evidence that R-MAON improves real strict dual-cold affinity prediction
over the shared-global baseline, and current public data cannot support that comparison. Category 1 is
therefore declined. Category 3 is also too strong because structured interaction and a calibrated
operator are present; what is missing is independent multi-family measurement support and real transfer
evidence.

`sealed_test_consumed=false`; no Davis, confirmation, or sealed label was read in the accepted run.
