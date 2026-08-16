# BioLiP-2026 LLM/rule strict-matrix metadata decision

Date: 2026-07-29  
Decision: `STOP_BL26_D0_STRICT_MATRIX_ADOPTION`

## Safe source audit

The official BioLiP homepage reports a July-2026 hybrid LLM/rule workflow and
23,502 new records. The independent official
[`all_affinity.tsv`](https://zhanggroup.org/BioLiP/data/all_affinity.tsv) was
audited with the four affinity spans skipped at byte level:

- 44,925,681 bytes;
- SHA-256
  `eea56333b650bb900d175e673cafe938cb04a647d2fea20c52ab7346175320a`;
- 68,890 rows and 23 columns;
- zero affinity-field bytes decoded.

The file is a real new target-ligand affinity source, but its schema does not
make it an experimental exact-complex contact source. In particular, the
LLM-AFDB rows use predicted receptor structures, and the LLM-RCSB/rule-BioLiP2
rows map external ligand SDFs to receptor sites without an observed
ligand-chain/serial atom instance.

## Frozen strict metadata gate

Safe identity required one exact UniProt, one full 27-character InChIKey, one
numeric PMID, and a curation stratum. Duplicate identity triples were
collapsed. This retained 17,679 triples:

| Safe scale | Result |
| --- | ---: |
| Targets / ligands / PubMed IDs | 4,190 / 5,850 / 9,903 |
| LLM / rule triples | 14,316 / 3,363 |
| Exact ChEMBL-TRAIN accessions | 291 / 559 (52.1%) |
| ChEMBL-TRAIN homology components | 279 / 517 (54.0%) |
| Exact ChEMBL-TRAIN ligands | 768 / 121,401 (0.633%) |

The source is well located in the target domain but poorly located in the
joint target-ligand matrix. Requiring both exact TRAIN identities while the
pair itself is new leaves:

| Strict new-edge topology | Result |
| --- | ---: |
| Joint exact target-ligand rows | 557 |
| Already present TRAIN pairs | 299 |
| New matrix edges | 258 |
| Targets / homology components | 94 / 91 |
| Ligands / scaffolds / PubMed IDs | 157 / 103 / 141 |
| Conflict components | 18 |
| Largest component | 88.76% |
| Conflict-free resource ceiling / greedy packing | 91 / 59 |

The largest PubMed, target, homology, and scaffold shares are respectively
10.47%, 10.47%, 10.47%, and 9.69%. The frozen edge-scale, topology, packing,
and concentration gates all fail.

## Evidence and rights boundary

The BioLiP2 paper gives a database-level BSD availability basis, and the
[official download page](https://zhanggroup.org/BioLiP/download.html) makes
the source publicly downloadable. However the 2026 file has no dedicated
file-level license/readme and exposes no row-level evidence sentence/span,
extraction confidence, LLM/model/prompt version, rule-set version, or
validation status. The file hash and official 2026-06-29 release date are the
only exact extraction-version anchors currently available.

Those omissions prevent an auditable literature-extraction teacher even
independently of the failed strict topology. They also require clarification
before assuming that every upstream source and commercial redistribution use
inherits the code repository's BSD terms.

## Nonbinding upper bound

An exploratory relaxation retained exact TRAIN targets with novel ligands and
optimistically assigned every unknown scaffold its own singleton token. It
produced 1,979 rows, 268 homology components, and a greedy packing of 259.
That can exceed the 88-unit mechanism floor but not the 423-unit predictive
floor; its largest component is still 86.56%. Because it changes the frozen
identity question, it does not override the STOP.

## Consequence

BioLiP-2026 remains catalogued as a high-target-coverage, low-strict-matrix
candidate affinity source. It does not:

- supply observed ligand-conditioned contact events;
- rescue RDIB/PD-MVR;
- replace the G0PB same-scaffold contact substrate;
- authorize decoding the four affinity fields; or
- authorize model training.

Reopening requires row-level evidence/version metadata plus a new frozen
precision audit, or a materially larger provenance-separated set of strict
new target-ligand matrix edges.

Machine record: `reports/active/biolip26_d0.json`.
