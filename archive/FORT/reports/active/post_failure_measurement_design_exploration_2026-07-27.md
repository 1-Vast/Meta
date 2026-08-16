# Post-failure raw-source and measurement-design exploration

Date: 2026-07-27

## Scope and status

This is a read-only exploration after the required three-candidate failure
report. It is not a fourth candidate, a model experiment, or an authorization
to contact providers, request compounds, download a new corpus, or train a
model. F1--F4 remain closed.

## Raw-source landscape

| source | useful raw fields | why it cannot enter F0 directly |
| --- | --- | --- |
| PubChem BioAssay | AID/version, source IDs, optional NCBI protein accession, PMID/DOI, SID-level results, endpoint/unit | Target and organism metadata are optional; NCBI accessions are not exact UniProt identifiers; AID is not an independent document or repeat; endpoints are often percent effect, AC50, or IC50; provenance can include ChEMBL mirrors. |
| Guide to Pharmacology 2026.2 | human UniProt, structures, pKi/pKd, relation, PMID/patent, assay description | The public interaction table has no stable assay ID and includes affinity-range summaries. Strict human exact-pKi/pKd filtering leaves 3,543 rows and only 36 cross-PMID repeated parent--target--endpoint cells across 30 targets. |
| Drug Target Commons | reported PMID retrieval, UniProt links, and assay annotations | The current release could not be schema- and checksum-verified from this environment. Historical terms are CC BY-NC-SA 3.0 and its integrated BindingDB/multidatabase content cannot be presumed to be independent raw provenance. |
| ChEMBL raw records | activity, assay, document, and exact accession fields | Already failed the provenance-family audit: apparent document replication was overwhelmingly duplicated measurement lineage. |
| BindingDB native articles | article-level continuous affinity records | Already audited: only 38 pKi and 6 pKd targets have an article block with at least 40 ligands before strict dual-cold firewalls. |

PubChem remains a possible upstream reconstruction resource only. Any future use
would require a source/protocol/PMID-or-DOI provenance-family model, exact WT
UniProt mapping, endpoint normalization without Ki/Kd/IC50 pooling, license
verification per contributing record, and a new preregistered topology audit.
It has not passed any of those checks here.

## Exploration outcome

No currently available open dataset is eligible to enter F0. Renaming or
re-aggregating an existing source would not address the failure modes already
observed in MMP-X, TCOPA, Papyrus, or the earlier raw-ChEMBL provenance audit.

## PRISM-PANEL compatibility delta

The archived PRISM-PANEL is a prospective measurement-design program, not a
fourth model candidate. Its historical `PRISM_P1_READY` package has real
strengths: continuous biochemical Ki dose-response curves, accession and
construct fields, raw curve/QC return requirements, blinded compound codes,
and planned bridge repeats. It is not yet compatible with the current F0
contract without the following amendments.

1. Add a preassigned `provenance_family_id` to every raw measurement, based on
   genuinely independent laboratory/site/campaign lineage. Plate, batch, or
   assay ID alone must not be treated as an independent provenance family.
2. Design repeated parent--target cells across at least two such independent
   provenance families. Technical replicates and bridge compounds estimate
   measurement error but never count as independent replication components.
3. Before external work, run a label-blind topology and power audit after exact
   accession, WT/construct, homology/pocket, parent, scaffold, high-similarity,
   and provenance firewalls. A target count or 500 transformations is not a
   substitute for independent component count.
4. Strengthen the future held-out tranche beyond target/transformation/scaffold
   separation to include parent connectivity, chemical-neighbour, homology and
   pocket-family separation, plus an explicit provenance-family firewall.
5. Keep P1 as a reliability pilot only. Its 80 transformations and one planned
   operational setting cannot establish a powered, provenance-disjoint
   target-conditioned transfer claim. A future P2 must pass the amended F0
   topology/power gate before any encoder comparison.

## External authorization boundary

The historical P1 package deliberately left RFQs unsent, orders unplaced, and
provider/collaborator contact untouched. Executing it requires explicit user
authorization, compound-access arrangements, an assay provider or laboratory,
and funding. This exploration has not taken any of those actions.

## Current conclusion

The only scientifically coherent continuation is a prospective factorial panel
with raw per-measurement records, exact WT target identity, independent
provenance families, and a pre-registered strict dual-cold topology audit. The
next safe local task is to freeze the above PRISM compatibility amendment; it
does not authorize wet-lab execution or Mamba/DTA training.

