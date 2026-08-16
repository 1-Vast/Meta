# Open dense continuous-pKd kinase panel acquired — Reinecke 2024 (2026-07-26)

Authorised open-dataset download under the `OPEN_DATA_ONLY` amendment. This overturns the specific
blocker recorded in the July OPEN-S audit, which found the *2017* Klaeger Kd(app) matrix "not openly
recoverable" (publisher 403, PRIDE raw spectra only, Europe PMC "not open access"). The **2024
follow-up is open access (CC BY 4.0)** with the processed affinity matrix in supplementary.

## Source and provenance (recorded in `manifests/open_sources.json`)

* Reinecke et al., *Nat Chem Biol* 2024, 20(5):577–585, "Chemical proteomics reveals the target
  landscape of 1,000 kinase inhibitors." DOI `10.1038/s41589-023-01459-3`; PMC `PMC11062922`;
  licence **CC BY 4.0**. Compounds: PKIS/PKIS2/KCGS/Roche kinase chemical-probe libraries.
* Downloaded 2026-07-26 from the Springer open-access static-content host (PMC `bin/` and NCBI OA FTP
  were bot-blocked / unreachable from this environment; recorded honestly). SHA-256 of all five
  supplementary Excel tables recorded in the manifest and in `raw/sha256.txt`. Path
  `dataset/public/reinecke_2024/raw/`.
* Key tables: **MOESM3 / Table S2** "Drug matrix of target affinity values" (`Kinobeads Drugmatrix -
  all`: 318 genes × 1183 compounds, apparent Kd in nM) + `Target annotation`; **MOESM2 / Table S1**
  compound annotation with **SMILES** (1182/1183 covered); MOESM4 / Table S3 a 10,554-row melted
  `(compound, gene, apparent Kd)` list.

## Cheapest gates run this round (model-free) — both PASS

**Shape feasibility.** From the Kd matrix (10,553 measured cells, 22 non-positive artifacts dropped):

| threshold | targets (all) | protein/lipid-kinase targets |
|---|---:|---:|
| ≥ 10 binders | 245 | 215 |
| ≥ 20 binders | 180 | 162 |
| **≥ 40 binders** | **105** | **96** |
| ≥ 60 binders | 50 | 47 |

Per-target median 25 binders (max 178). This is the first **open** source to meet the amendment's
*shape* requirement (~100 targets with ~40 within-target ligands), and it is **dense-both** (target-
broad and ligand-deep) rather than the target-shallow BindingDB-native or ligand-shallow Davis shapes.

**Label resolution.** Clean pKd (= 9 − log10 Kd_nM) range 5.45–9.95, mean 6.28, SD 0.61, with
**10,513 of 10,531 values distinct (essentially tie-free)**. This directly addresses the program's
label-resolution ceiling: Metz pKi is rounded to 0.1 pK with 81.7% within-target duplicate values,
which caps within-target ranking resolution; Kinobeads apparent Kd is effectively continuous.

## Honest status — NOT yet a confirmed powered independent panel

```
OPEN_DENSE_CONTINUOUS_KINASE_PANEL_ACQUIRED__ADMISSIBILITY_AUDIT_PENDING
```

Shape and label resolution are necessary, not sufficient. Admissibility under `OPEN_DATA_ONLY` still
requires, and this round does **not** yet establish, all of:

1. **Homology-component count.** 96 kinases must be mapped gene-symbol → UniProt → sequence and
   k-mer-containment clustered; kinases cluster by family, so the independent-unit count (the real
   statistical n) will be **below** 96 and must be measured, not assumed.
2. **Independence firewall.** Overlap vs the ChEMBL-37 train extract **and** the spent Metz panel at
   homology-component, UniProt accession, ligand parent-connectivity, Bemis–Murcko scaffold and
   Tanimoto<0.95 levels. PKIS/PKIS2/KCGS are public and partly in ChEMBL, so compound overlap is a
   real risk; kinase targets overlap Metz/ChEMBL by construction, so the *independent-confirmation*
   component count is whatever survives homology-disjointness from Metz.
3. **Power.** A grouped MDE80 audit on the surviving components (as done for Metz/Davis) — only then
   can it be called powered.
4. **Role fixed before any label is scored** (development-pool vs single-use confirmation), per the
   amendment; and pKd here is a different endpoint from Metz pKi, so any pooling keeps panel-specific
   intercepts/scales/endpoints.

Two plausible uses, to be decided by the audit: (a) a **second independent panel** whose Metz-disjoint
components are pooled to raise total independent components toward the ~150–400 the power decomposition
says is needed; (b) an **independent confirmation** source for a model developed on Metz. Neither is
claimed yet.

No confirmation/Davis/sealed label was read. No model was trained. No predictive claim is made. The
next round is the registered admissibility audit (registry build with UniProt sequences + k-mer
clustering + the firewalls + MDE80); its result decides whether the powered-independent-panel blocker
(`NO_OPEN_POWERED_INDEPENDENT_PANEL`) is finally cleared.
