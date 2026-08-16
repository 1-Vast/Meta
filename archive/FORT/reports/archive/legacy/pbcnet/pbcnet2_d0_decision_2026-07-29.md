# PBCNet2.0-D0 decision

Date: 2026-07-29  
Decision: `STOP_PBCNET2_D0_IDENTITY_OR_PROVENANCE_INADEQUATE`

## Result

The public release is now materially more accessible than the previously
restricted record:

- Zenodo record `18299525` is open under CC BY 4.0;
- six files total `5,917,812,656` bytes;
- the pair manifest is `38,413,018` bytes with MD5
  `6407a4a64969e8d8d9166c08b1a4e509`;
- the BindingDB-derived source archive is `5,441,875,666` bytes with MD5
  `b1f52cb60c9326e99fff3be7a2834845`.

The first-line-only projection found:

`lig1, lig2, smile1, smile2, similarity, Label1, Label2, Label, abs_L,
0, 相似度区间, 活性差值区间, dir_1, dir_2`.

`lig1/lig2` and `dir_1/dir_2` are sufficient to identify the two pair members,
and the label columns are separable. However, the manifest exposes no
BindingDB row ID, assay, document, DOI/PMID, source, patent, or other lineage
key. The 8.6 million rows are generated pair combinations, so without lineage
they cannot be collapsed into independent measurements or screened for
ChEMBL/BindingDB overlap and provenance duplication.

The GitHub repository was pinned at commit
`3d46e6e594531c5692376e242b606641979e8550`. Its README says MIT, but the root
contains no LICENSE file and the GitHub API detects no license. Code reuse is
therefore not authorized by this audit.

## Firewall

Exactly 110 bytes through the first newline were decoded. Zero CSV data rows
and zero outcomes were read. The 38.4 MB manifest, 5.44 GB archive, model
weights, development/confirmation data, and sealed test were not downloaded
or loaded.

## Consequence

Do not download or train on the PBCNet2.0 package for RBSDD. Reopening requires
an author-provided mapping from each `dir_1/dir_2` member to original
BindingDB row/assay/document/source identifiers and an explicit repository
license file for code reuse.

Authoritative machine result: `reports/active/pbcnet2_d0.json`.

