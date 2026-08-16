# RAMCI-DTA Stage S0 — source and license audit

Date: 2026-07-26. RAMCI-DTA asks whether a small, unbiased, known-probability measurement **anchor**
can make target-specific dual-cold interaction identifiable despite MNAR historical bioactivity. S0
audits whether an admissible, complete/auditable, dual-cold-capable, open anchor source exists, before
any modeling. Binding prior: RECRO v2 L0 showed the ChEMBL observational graph is MNAR + provenance-
duplicated (`RECRO_SIGNAL_EXPLAINED_BY_PROVENANCE`); the reproducible baseline remains B0.

## Candidate-by-candidate audit

| source | license | families | endpoint | inactives retained? | uniform/complete? | dual-cold-capable at power? | local? |
|---|---|---|---|---|---|---|---|
| **Davis** `CHEMBL1908390` | ChEMBL CC BY-SA | kinase only | pKd (dose-response) | yes (no floor spike) | **yes, complete** | **NO** — 68 cpds → 12 query/target → MDE80 **0.16** | yes (sealed) |
| Metz `CHEMBL1201862` | CC BY-SA | kinase only | pKi | partial | dense, single-doc | spent (dev); coarse ties 0.817 | yes |
| KirHub 2026 | CC BY | kinase only | %inh (single-dose) | yes | dense, single-source | spent; within-kinome; saturation | yes |
| Reinecke 2024 | CC BY 4.0 | kinase only | pKd (Kinobeads) | partial | **sparse ~2.8% fill** (not complete) | spent; within-kinome | yes |
| PKIS/PKIS2 | CC BY 4.0 (Zenodo) | kinase only | %inh (single-dose) | yes | dense profiling | within-kinome; low-fidelity | no |
| KCGS | CC BY 4.0 | kinase only | %inh | yes | matrix not local (only compound/kinase lists) | within-kinome | metadata only |
| KIBA | mixed | kinase only | **composite Ki/Kd/IC50 score** | no | aggregated from public MNAR | **BARRED** (composite + provenance) | no |
| Binder2030 | — | — | — | — | **not a real/indexed dataset** | — | no |
| CACHE | open (challenge data) | single-target per round | prospective hit-finding | n/a | one target per challenge | **not a multi-target panel** | no |
| activity-integration (Harmonic) | **empty LICENSE** | kinase | single-dose+IC50 | partial | ChEMBL/PubChem (MNAR) | not admissible as-is; MNAR | staged |
| kinase-modelling (Papyrus) | empty LICENSE | kinase | pIC50 (barred) | no | aggregated ChEMBL (MNAR) | leakage + MNAR | staged |
| **Novartis SPD** | **CC BY-4.0 (data) + MIT (code)** | **multi-family** (GPCR/kinase/ion-channel/enzyme/transporter) | **AC50 (8-pt dose-response)** | **yes (key property)** | **systematic safety panel** | **UNVERIFIED — needs download** | **no** |

## Key finding

The local ChEMBL extract contains no dual-cold-capable complete panel: the one complete/clean panel
(Davis) is underpowered (MDE80 0.16, 12 query ligands/target) and within-kinome, exactly as the Davis
registration already concluded ("an admissible source needs ~100 components AND ~40 query ligands/target
— a second Metz-scale panel; the local extract does not contain one"). All other local/staged sources
are within-kinome, knowledge-selected, sparse, MNAR-aggregated, or barred by endpoint/license.

**Novartis Secondary Pharmacology Database (SPD)** is the one genuinely-new admissible candidate and the
first that could break both program limitations at once:
- **Open + licensed:** Zenodo `10.5281/zenodo.8103950` (CC BY-4.0, 150 MB), GitHub `Novartis/SPD` (MIT),
  Nature Comms 2023 (Brennan et al., doi:10.1038/s41467-023-40064-9).
- **Multi-family:** 200 secondary-pharmacology assays across GPCRs, kinases, ion channels, enzymes,
  transporters — not within-kinome, so cross-family selectivity (own-family-cold falsifiable) is
  testable, unlike the within-kinome panels where TR-0's group control failed.
- **Non-MNAR:** systematically tested safety panel that **retains inactive/low-activity results**; the
  paper reports 95% of results are unique vs public resources which are biased toward higher activity —
  i.e. SPD supplies exactly the unbiased/complete measurement structure the RAMCI anchor requires.
- **Quantitative:** 8-point concentration-response AC50 (not single-dose), released per-compound-per-
  assay (`final_summarized_activity_data_pub.txt`, `Dataset_S1.xlsx`) with an assay→gene map.
- **Structures recoverable:** 1,958 marketed/prescribable drugs (public structures via the released
  DrugCentral mapping); not proprietary-blocked.

## Unverified items (require the download) and considerations

1. **Density / dual-cold capability at power:** # biochemical single-target assays (vs cell-based, which
   must be excluded), independent homology-component count, and query-ligand depth after a
   scaffold-disjoint split. S0 forbids inferring density from unique counts alone.
2. **Estimand shift:** SPD is *secondary/off-target* pharmacology of marketed drugs. Most drugs' primary
   therapeutic target is not in the panel, so on-target interactions may be sparse and the signal
   skews toward off-target-liability prediction — a legitimate but different estimand from primary
   binding-affinity dual-cold. To flag before committing.
3. **Leakage vs FORT:** marketed drugs overlap ChEMBL heavily; SPD can only serve a **self-contained**
   S1 (dual-cold enforced within SPD), never as a FORT confirmation source without heavy firewalling.
4. **Endpoint:** AC50 is assay-condition dependent; admissible only as within-assay ordinal for the
   RAMCI interaction estimand (not as exact pK), consistent with RAMCI's own contract.

## S0 status and decision

Local sources: **`STANDARDIZED_PANEL_NOT_DUAL_COLD_CAPABLE`** (Davis and all local panels are
within-kinome and/or underpowered). One admissible non-local candidate (SPD) is identified but its
dual-cold capability is unverified and it requires a 150 MB permission-gated download plus a substantial
build (assay→UniProt→ESM mapping, drug→SMILES→scaffold/homology firewalling, density audit) before the
S1 pseudo-prospective MNAR simulation can run. S0 is therefore **not yet terminal**: the decision to
download and verify SPD is escalated for review. `sealed_test_consumed=false`;
`confirmation_labels_read=true` (pre-existing; S0 read no FORT labels).

## SPD download + structure verification (approved 2026-07-26)

SPD activity + mapping files downloaded from Zenodo (CC BY-4.0), md5- and sha256-verified, recorded in
`dataset/public/spd_2023/manifest.json`. Structure audit of `final_summarized_activity_data_pub.txt` +
`assay_group_vs_gene_map.txt`:

* **Admissible and multi-family, as hoped:** 1,948 marketed drugs (DrugCentral IDs + InChIKeys →
  structures recoverable from open DrugCentral); 144 gene-mapped assay groups / **101 genes across
  GPCR, kinase, ion-channel, enzyme, nuclear-receptor, transporter families**; endpoint is 8-point
  concentration-response summarized IC50/AC50 with censoring prefix.
* **Systematic with retained inactives (the non-MNAR property, confirmed):** every gene-mapped tested
  cell has `N CRC total > 0`; 87,757 tested cells (~31% of the 1,948x144 space); **8.7% active, 91.3%
  genuine tested-negatives** (censored `>`). This is exactly the unbiased/inactive-retaining structure
  RAMCI needs, and it is real.

**But it is not dual-cold-capable at power.** The scarce resource is query depth, and SPD is thin on
it: **median 14 compounds per tested gene-assay; only 30 of 144 assays have >=40 compounds.** This is
comparable to the underpowered Davis panel (12 query ligands/target -> MDE80 0.16) and below the ~40
query ligands/target the Davis registration established as necessary for a resolvable dual-cold
interaction test. Component count (~80-100 multi-family genes) is adequate, but — exactly as the Davis
audit proved — component count and query depth are separate resources and only query depth determines
per-component ranking noise. Additional handicaps: 91% floor-domination and an off-target shape (only
751 on-target `is_moa` cells), which further thin the target-specific interaction signal; and the
quantitative unbiased Novartis-only AC50 is ~4% fill (the systematic-negative bulk is censored, not a
dense continuous anchor).

## Terminal S0 verdict

**`STANDARDIZED_PANEL_NOT_DUAL_COLD_CAPABLE`.** No local source is dual-cold-capable, and the one
genuinely-new admissible open candidate (Novartis SPD) — despite correctly supplying a multi-family,
systematic, inactive-retaining panel that addresses the MNAR limitation — is underpowered on query
depth (median 14, like the underpowered Davis) and off-target/floor-dominated. A powered, clean,
dual-cold RAMCI S1 cannot be run on it; building the full SMILES/ESM mapping + 12-arm S1 on an
already-underpowered substrate is inadmissible (a model must not be built on a substrate that cannot
resolve the planned effect).

RAMCI's contribution is nonetheless real and convergent: it identified the closest existing template to
its ideal source and sharpened the requirements for the **prospective S2 panel** (the program's
standing recommendation): a multi-family systematic panel that additionally retains inactives with a
known/randomized inclusion probability AND provides >=~40 scaffold-diverse query ligands per target
across >=~100 independent homology components. No existing open source — within-kinome or multi-family —
satisfies all of these simultaneously. `sealed_test_consumed=false`; SPD `confirmation_labels_read` n/a
(SPD read for source audit only; no FORT dev/confirmation/sealed label read).
