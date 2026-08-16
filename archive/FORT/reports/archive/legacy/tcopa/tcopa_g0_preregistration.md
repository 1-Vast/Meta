# TCOPA-G0 preregistration: target-contrast pretraining data gate

Date: 2026-07-27  
Candidate count: agent-proposed candidate 2 of at most 3.  
Role: dataset and split audit before sequence acquisition or model training.

## Candidate

**Target-Contrastive Orthogonalized Pretraining Anchor (TCOPA)** uses systematic ToxCast hit calls
to learn a protein-conditioned interaction anchor from *within-ligand target contrasts*:

`same ligand: active on protein A, inactive on protein B`.

This changes the supervision unit from absolute ligand activity or a broad family average to an
exact-protein response reversal. A future model must decompose

`logit P(hit) = ligand tendency + assay nuisance + g(ligand, protein)`

and train the anchor `g` with both hit-call BCE and a pairwise target-contrast loss. Batch
double-centering and an assay-source/technology nuisance head prevent `g` from becoming another
ligand-promiscuity or protocol carrier.

This is distinct from failed SAFSA (six-class family means) and from MMP-X (directional chemical
edits requiring cross-source repeated transformation--family units).

## Sources

- Official EPA ToxCast invitroDB v4.3 target mappings and assay annotations, CC0,
  DOI `10.23645/epacomptox.6062623.v14`.
- The public DeepChem/MoleculeNet ToxCast CSV is used only as a lightweight local label projection
  for mechanism proof. It is not an independent source and cannot establish cross-source
  replication. Full-scale training, if authorized, must use the official invitroDB release.

Hashes and URLs are frozen in
`dataset/public/toxcast_invitrodb_v4_3/manifest.json`.

## Disclosed feasibility inspection

Before this registration, schema and aggregate counts showed:

- 8,597 structures and 617 projected endpoints;
- 1,538,395 observed binary labels, including 1,411,325 retained negatives;
- exact endpoint-name intersection with the current official target map;
- after a provisional single-human-gene/protein filter, about 57 genes, 10 intended target
  families, 67 endpoints and 88,668 labels meet `n>=100`, `positive>=20`, `negative>=20`.

No within-ligand discordant-pair count, scaffold allocation, target fold, model score or downstream
label was inspected. The thresholds below may not be changed after the audit.

## Frozen endpoint registry

1. Join a projected endpoint exactly to
   `assay_component_endpoint_name.x` in the official target map and to
   `assay_component_endpoint_name` in official assay annotations.
2. Keep human (`ncbi_taxon_id=9606`) `entrez_gene_id` mappings with a nonempty official symbol.
3. Keep an endpoint only if it maps to exactly one human gene and
   `intended_target_type == protein`.
4. Require at least 100 observed labels, 20 positives and 20 negatives.
5. Select one endpoint per exact gene without using outcome values: maximum observed-label count,
   then lexical endpoint name. Record the discarded same-gene endpoints as protocol metadata, not
   independent targets.
6. Canonicalize structures with RDKit, collapse exact canonical duplicates, and mark a gene label
   missing if duplicate rows disagree. Missing values never become negatives.

## Frozen contrast and split units

A discordant target-contrast unit is a unique
`(canonical ligand, active exact gene, inactive exact gene)` triple where both labels are observed.
The two genes must differ. Count each triple once regardless of endpoint aliases.

Ligand holdout key is Bemis--Murcko scaffold; acyclic compounds use their exact canonical SMILES.
Target holdout key is exact human gene. Build five deterministic seed-1729 target folds by greedy
balancing family counts, gene cell counts and positive counts. No exact gene may cross folds.

This G0 audits whether both axes can support a later Cartesian dual-cold validation. It does not
read FORT development, confirmation or sealed labels and does not claim ToxCast hit calls are
affinity labels.

## Frozen gates

All must pass:

1. at least 50 representative exact genes from at least eight intended target families;
2. both NVS and TOX21 assay sources contribute representative protein endpoints;
3. at least 50,000 observed gene-level cells, with tested-negative fraction at least 0.70;
4. at least 1,000 canonical ligands have five or more observed representative genes;
5. at least 500 ligands have both an active and inactive exact gene;
6. at least 10,000 unique discordant target-contrast triples spanning at least 200 distinct
   active-gene--inactive-gene pairs;
7. at least 1,000 ligand scaffolds and no scaffold contains more than 5% of canonical ligands;
8. every target fold contains at least eight genes, at least five intended families, 5,000 observed
   cells and 500 positives;
9. assay source, detection technology and assay design metadata are nonmissing for at least 95% of
   representative endpoints;
10. current-run FORT confirmation labels are unread, historical confirmation remains quarantined,
    and the sealed test is unconsumed.

Pass verdict:
`TCOPA_G0_PASS_AUTHORIZE_SEQUENCE_AND_FIREWALLED_DATASET_BUILD`.

Fail verdict:
`TCOPA_G0_INSUFFICIENT_TARGET_CONTRAST_SUPPORT_STOP`.

On pass, the next registration must define the downstream-entity firewall, UniProt sequence
mapping, exact loss, baselines and mechanism-proof acceptance thresholds before any training. On
failure, do not relax endpoint depth, merge multi-gene endpoints or count missing cells as
negatives.

