# Stage W0b — Core Task 1 data / censoring / hierarchy / positive-control audit

Preregistration SHA-256:
`ff23c408d20cc79b1bd5fcd20854a0443280d32d6fc3dbb8abf0733a9a70631f`.
Machine-readable: `W0B_AUDIT.json`. No neural model was trained.

## Verdict: W1 GO/NO-GO = NO-GO on current local assets

**The mandatory biological positive control (W0-P) is NOT RUNNABLE**: no
HIV-resistance, gatekeeper/mutation or ortholog point-mutant panel with
matching ligand affinities exists in the local workspace. Under the frozen W0b
rules, support statistics alone cannot authorize a biological signal or
biological-null interpretation. In addition, the single-platform panels have
severe censoring that the earlier Stage W W0 did not report.

## 1. Asset inventory

Local, hash-recorded assets:
- Davis `dataset/raw/dta/davis.tab`;
- Metz `metz_matrix.csv` + `metz.xls` (3,858 rows in xls; 1,497 with SMILES;
  704×172 matrix);
- Klaeger `klaeger_matrix.csv` + `klaeger_smiles.json` (222×343);
- KIBA `dataset/raw/dta/kiba.tab`;
- KLIFS kinase groups and pocket annotations.
- Anastassiadis is `excluded_by_governance` in the acquisition manifest and
  absent.

## 2. Censoring audit — decisive

| dataset | cells | censored | fraction | per-target median censored |
|---|---:|---:|---:|---:|
| Davis Kd | 25,772 | 18,343 | **71.2%** | 73.5% |
| Metz pKi | 39,216 | 23,690 | **60.4%** | 60.5% |
| Klaeger apparent pKd | 75,117 | 70,243 | **93.5%** | 95.0% |
| KIBA score | 117,657 | 0 | 0.0% | 0.0% |

Davis/Metz/Klaeger labels are dominated by their detection floors (Kd=10000,
pKi=4.0, apparent pKd=5.0). Any Pearson/MSE/cliff statistic on raw values
would be dominated by censored floors. Censored rows must be excluded or
modelled as bounds in any later stage.

## 3. Estimand hierarchy (broad -> strict)

| dataset / layer | within-target pairs | classes | classes >=3 targets & >=3 comps | cross-component D rows | EIU | top-1 class share |
|---|---:|---:|---:|---:|---:|---:|
| Davis all_pairs | 863,362 | 2,278 | 2,278 | 162.4M | 206 | 0.166 |
| Davis similar_pairs | 863,362 | 2,278 | 2,278 | 162.4M | 206 | 0.166 |
| Davis MMP | 2,653 | 7 | 7 | 498,967 | 7 | **54.1** |
| Davis strict MMP | 2,653 | 7 | 7 | 498,967 | 7 | **54.1** |
| Metz all_pairs | 4,219,332 | 24,531 | 24,531 | 304.6M | 10 | 0.007 |
| Metz MMP | 82,560 | 480 | 480 | 5.96M | 10 | 0.358 |
| Klaeger all_pairs | 8,187,753 | 23,871 | 23,871 | 1.14B | 10 | 0.014 |
| Klaeger MMP | 13,034 | 38 | 38 | 1.81M | 10 | **9.03** |

Broad ligand-pair layers have massive support but the censored fraction makes
their signed differences floor-dominated. Strict MMP layers have interpretable
chemical context but very few classes on Davis (7) and high class
concentration on Davis and Klaeger.

## 4. Stage W audit

Stage W (`stageW_soft_mmp`) artifacts retained unmodified. Its W0 prereg
SHA-256 `ae96762e…45dc71` was frozen before W0 statistics; its W1 prereg
SHA-256 `038f4d97…49082` was frozen before W1 split statistics. W1 has **no
trained artifacts and no training metric**; it is marked **PAUSED**. The Stage
W W0 census did not include the censoring audit and did not test W0-P; those
omissions are corrected here without retrofitting Stage W.

## 5. GO/NO-GO

- **W1 biological interpretation: NO-GO.**
- **Davis:** NO-GO for W1 — 71.2% censored, only 7 MMP classes, W0-P absent.
- **Metz/Klaeger:** NO-GO for W1 until censored labels are excluded/modelled
  and a positive control is available; they remain valuable replication
  panels.
- **KIBA:** descriptive only — score semantics differ from pK, and a KIBA-only
  positive cannot satisfy the Core Task 1 SOLVED definition.

## 6. What would change NO-GO to GO

1. acquire or construct an admissible W0-P panel (resistance mutations,
   orthologs, or point mutants with matched ligand panels);
2. re-census with censored rows excluded or interval-censored, then
   preregister thresholds from the post-censoring effective sample size;
3. for Davis, prefer the broad or similarity layer over the 7-class MMP layer
   unless W0-P and a strict layer can both be satisfied.
