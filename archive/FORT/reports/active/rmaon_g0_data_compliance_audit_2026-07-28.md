# R-MAON G0 data and compliance audit

Date: 2026-07-28

## Decision

The local data do not authorize a real-label R-MAON performance experiment.

```text
RMAON_G0_TOPOLOGY_OR_POWER_STOP
```

This is a data-identifiability stop, not an estimator-effect failure. The separately registered synthetic
module gate passed as `RMAON_G0_NULL_SCORE_AND_RECOVERY_PASS__MODULE_ONLY`. A0, M1, strict dual-cold
training and prediction remain blocked.

No development, Davis, confirmation or sealed label was read for this audit.
`sealed_test_consumed=false`.

## Binding requirements

The binding predictive substrate remains:

* approximately 423 independent multi-family target components;
* `PA2 >= 0.5 pK`;
* at least 40 scaffold-diverse query ligands per target;
* target accession, homology/pocket, binding-profile correlation, ligand parent, scaffold,
  chemical-neighbour, assay/document provenance and duplicate isolation;
* randomized inactive-retaining inclusion with known nonzero probabilities;
* at least two genuinely independent provenance lineages.

No frozen local source satisfies all of these conditions.

## Reader-boundary correction

`src/data/dualcold.py::DualCold.__init__` reads the complete Parquet registry and materializes
`affinity` before callers select a split. The previous
`research/klbp_r3_synthetic.py::load_real_coefficients` invoked `DualCold.panel()` only to construct
component folds, making that path unsuitable for G0.

The corrected loader obtains affinity-bearing rows only through
`research.panel_gate_pa.load_panel_train`, which projects named columns and applies
`dual_cold_split == train` in `pd.read_parquet`. Component folds are now reconstructed from that
already-filtered TRAIN frame. The G0 report records both
`only_train_split_read=true` and `parquet_train_filter_applied=true`.

## Safe TRAIN-only geometry

| source | TRAIN rows | targets / clusters | ligands | components | admissible G0 role |
| --- | ---: | ---: | ---: | ---: | --- |
| ChEMBL-37 dual-cold | 201,827 | 559 targets | 121,401 | 517 | engineering only; fails the registered interaction-identifiability floor and lacks the prospective provenance/profile design |
| Metz | 12,574 | 112 targets | 619 | 101 | empirical coefficient/noise template and synthetic calibration only; kinase-only, one document |
| PLINDER | 5,804 | 2,221 clusters / 3,828 accessions | 2,717 | 1,081 | engineering diagnostics only; only 217 TRAIN clusters have at least four ligands and median depth there is 7 |

The R1 coordinate aligns to 111 eligible Metz targets. Component-balanced whitening retained 11 of 14
coordinate dimensions. Thirty eligible target covariance matrices were rank deficient (rank range 9--64);
G0 used exact positive-semidefinite eigensquare-roots and added no jitter.

## Source-role audit

| source | licence / provenance | admissible role | blocking fact |
| --- | --- | --- | --- |
| ChEMBL-37 | CC BY-SA 3.0; aggregated ChEMBL documents | TRAIN-only engineering and noise geometry | the multi-family graph has `PA2=0.356 < 0.5`; no prospective randomization or independent lineage |
| Metz | ChEMBL CC BY-SA 3.0; document `CHEMBL1201862` | TRAIN-only empirical `V_t` and synthetic calibration | 101 kinase-only components, one document; 64 accessions and 618 ligands overlap ChEMBL TRAIN |
| Davis | ChEMBL CC BY-SA 3.0; document `CHEMBL1908390` | single-use confirmation only | sealed and unconsumed; median 12 query ligands and underpowered |
| PLINDER | CC BY 4.0; BindingMOAD-derived affinity | engineering diagnostics | sparse per target and no independent provenance/profile design |
| BindingDB native articles 202607 | CC BY 3.0 US; publisher archive and native-curation subset | candidate independent source, rejected for prediction power | only 38 pKi and 6 pKd targets have an assay-controlled article block with at least 40 ligands before further firewalls |
| OpenBind EV-A71/CVA16 2A | data CC0-1.0; benchmark code Apache-2.0 | single-target local-SAR mechanism evidence | one protein; cannot test target conditioning, protein destruction or unseen-target transfer |
| KirHub 2026 | article/supplement CC BY-NC-ND 4.0; one release | within-source kinase mechanism evidence only | one-concentration aggregate residual activity; assay/document/publication isolation impossible; median strict query depth 8 |
| Reinecke 2024 | CC BY 4.0; Kinobeads `pKd_app` | completed development source | role was fixed as development after label-shape inspection; kinase-only and no independent-confirmation claim |
| Novartis SPD 2023 | data CC BY 4.0; code MIT | systematic inactive-retaining engineering source | median 14 compounds per gene-assay, 91.3% censored inactive, ChEMBL/DrugCentral-heavy provenance |
| EPA ToxCast v4.3 | CC0; DeepChem file is a projection of the same EPA source | mechanism-proof projection only | not independent affinity supervision; full-scale work would require the official release |
| Papyrus 05.7++ | CC BY-SA; aggregated release | closed provenance audit only | one aggregated row per parent-target; no document-resolved repeated cell for the required firewall |
| KLIFS | structural/mechanism metadata, no affinity outcome | non-affinity mechanism infrastructure | cannot supply randomized affinity supervision or prediction power |

Exact local evidence is recorded in:

* `dataset/public/chembl_37/processed/dualcold/manifest.json`;
* `dataset/public/chembl_37/processed/panel_metz/manifest.json`;
* `dataset/public/chembl_37/processed/panel_davis/manifest.json`;
* `dataset/public/plinder_2024_06_v2/processed/dualcold/manifest.json`;
* `manifests/open_sources.json` and `reports/active/open_data_only_amendment.md`;
* `reports/active/openbind_o0_decision.md`;
* `reports/active/kirhub_spkop_a1_strict_decision.md`;
* `dataset/public/reinecke_2024/processed/panel_reinecke/manifest.json`;
* `dataset/public/spd_2023/manifest.json`;
* `dataset/public/toxcast_invitrodb_v4_3/manifest.json`;
* `reports/active/papyrus_f0_decision.md`.

## Only valid next real stage

The first real-label action is a separately preregistered A0 reliability pilot, not model training:
at least 12 targets across six families, 16 randomized ligands and two independent provenance lineages.
It must establish correct-target residual-order signal, cross-site reliability and the dispersion needed
for the later 0.03 paired component-macro claim. Until those measurements exist, additional architecture
search cannot resolve the stop.
