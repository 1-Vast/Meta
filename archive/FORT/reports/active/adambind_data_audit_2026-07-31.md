# AdaMBind official data provenance audit (2026-07-31)

## Conclusion

The isolated download directory `D:\FORT\tmp\adambind-data` contains the three
training CSVs from the AdaMBind GitHub repository at immutable commit
`01a169a6d62fba0d6c003f47bfba539e55f5b344` and the paper's Figshare supplement.
The CSV byte lengths and Git blob SHA-1 values match GitHub's API metadata;
the SHA-256 values below are the local reproducibility fingerprints. The
checkout under `tmp\adambind-source` is the same logical data with CRLF line
endings introduced by Git's `* text=auto` checkout rule, so its byte hashes are
expected to differ from the downloaded LF files.

The GitHub repository has no declared license and contains no `LICENSE` file.
The Figshare supplement is explicitly CC BY 4.0. The underlying Davis, KIBA,
and BindingDB records should therefore retain their upstream dataset
attribution/terms; this audit does not infer a license for them from GitHub.

## Provenance

- Code repository: <https://github.com/Moohyun-w/AdaMBind>
- Pinned commit: `01a169a6d62fba0d6c003f47bfba539e55f5b344`
- Pinned commit URL:
  <https://github.com/Moohyun-w/AdaMBind/commit/01a169a6d62fba0d6c003f47bfba539e55f5b344>
- Dataset paths in that commit: `data/{bindingdb,davis,kiba}-full-data.csv`.
  The repository also stores byte-identical logical copies under `data/raw/`.
- Downloaded training-data directory: `D:\FORT\tmp\adambind-data\data\raw`
- Paper supplement DOI: `10.6084/m9.figshare.30963823.v1`
- Figshare API record: <https://api.figshare.com/v2/articles/30963823>
- Supplement download URL:
  <https://ndownloader.figshare.com/files/61860844>

## File fingerprints

| file | source URL | bytes | SHA-256 | Git blob SHA-1 / Figshare MD5 |
|---|---|---:|---|---|
| `bindingdb-full-data.csv` | `https://raw.githubusercontent.com/Moohyun-w/AdaMBind/01a169a6d62fba0d6c003f47bfba539e55f5b344/data/bindingdb-full-data.csv` | 32,865,741 | `3ebd8dfabd2a20c0dbceba35cc59ba8e6dd44a90798667fef2c9059bab63fbba` | Git blob `d53dfe26191102ec37d1c35d61f4dfcc24625b85` |
| `davis-full-data.csv` | `https://raw.githubusercontent.com/Moohyun-w/AdaMBind/01a169a6d62fba0d6c003f47bfba539e55f5b344/data/davis-full-data.csv` | 25,810,493 | `dc9331894d5eafa46787632cc0d9754406e5a96eb87980b27d4abe22308a6994e` | Git blob `03c759494a980ab70ab69a9703ea4a953cc11534` |
| `kiba-full-data.csv` | `https://raw.githubusercontent.com/Moohyun-w/AdaMBind/01a169a6d62fba0d6c003f47bfba539e55f5b344/data/kiba-full-data.csv` | 94,281,374 | `7b1e306a2344e38c5d5bbcda6f6112201440bbaa92d5081a4fc054ed83edca24` | Git blob `e02e153fd325d0bc8fdd9feacbebc2d08ab8b2ef` |
| `source_data.xlsx` | `https://ndownloader.figshare.com/files/61860844` | 3,981,306 | `1b73ef01a34578d543070caf0a724c65dcd4d397022c3741864efeaed52cb0ac` | Figshare MD5 `a5e2b3f5d754d169c063f6ecd61b3108` |

## CSV integrity summary

All three CSVs parse with UTF-8 CSV headers and finite numeric parsing except
for two `Infinity` affinity values in BindingDB. Counts are over data rows;
the Davis/KIBA leading empty column is an original row index.

| file | schema | rows | unique targets | unique compounds | unique compound-target pairs | repeated pair rows | non-finite affinity |
|---|---|---:|---:|---:|---:|---:|---:|
| BindingDB | `compound_iso_smiles,target_sequence,affinity` | 42,203 | 1,088 | 9,862 | 42,203 | 0 | 2 (`Infinity`) |
| Davis | `,compound_iso_smiles,target_sequence,affinity` | 30,056 | 379 | 68 | 25,772 | 4,284 | 0 |
| KIBA | `,compound_iso_smiles,target_sequence,affinity` | 118,254 | 229 | 2,068 | 117,657 | 597 | 0 |

The duplicate-pair counts reflect repeated rows in the released CSV and must
not be silently treated as independent observations in a strict benchmark.
The two BindingDB `Infinity` values require an explicit finite-value policy
before loss/evaluation; otherwise runs can produce NaN/Inf metrics. They occur
at CSV rows 28,088 and 28,550 (the two compounds share a target prefix
`MGAGALALGASEPCNLSSAAPLPDG`).

## Re-download and verify (PowerShell)

```powershell
$sha = '01a169a6d62fba0d6c003f47bfba539e55f5b344'
$dest = 'D:\FORT\tmp\adambind-data\data\raw'
New-Item -ItemType Directory -Force -Path $dest | Out-Null
foreach ($name in 'bindingdb-full-data.csv','davis-full-data.csv','kiba-full-data.csv') {
  Invoke-WebRequest "https://raw.githubusercontent.com/Moohyun-w/AdaMBind/$sha/data/$name" -OutFile (Join-Path $dest $name)
}
$supp = 'D:\FORT\tmp\adambind-data\supplement'
New-Item -ItemType Directory -Force -Path $supp | Out-Null
Invoke-WebRequest 'https://ndownloader.figshare.com/files/61860844' -OutFile (Join-Path $supp 'source_data.xlsx')
Get-FileHash (Join-Path $dest 'bindingdb-full-data.csv') -Algorithm SHA256
Get-FileHash (Join-Path $dest 'davis-full-data.csv') -Algorithm SHA256
Get-FileHash (Join-Path $dest 'kiba-full-data.csv') -Algorithm SHA256
Get-FileHash (Join-Path $supp 'source_data.xlsx') -Algorithm SHA256
```

## Failure log

1. **Initial GitHub tree API lookup failed.** Command used the commit SHA in
   `/git/trees/{sha}` and returned HTTP 404 (`Not Found`). The endpoint expects
   a tree SHA, not a commit SHA. The failure was corrected by querying the
   commit API and then `/contents/data?ref=<commit>`, which returned the three
   file sizes, download URLs, and Git blob SHAs above. No downloaded file was
   changed by this failed read-only request; the corrected metadata and local
   Git blob checks passed.

2. **License check found no repository license.** GitHub API returned
   `license: null` and the checkout has no `LICENSE` file. This is not a data
   download failure; it remains an attribution/redistribution limitation and
   is explicitly preserved in the conclusion above.
