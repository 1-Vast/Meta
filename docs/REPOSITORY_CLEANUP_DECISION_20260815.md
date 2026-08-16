# Repository Cleanup Decision - 2026-08-15

This table records the repository-hygiene decisions made before the next model
redesign. `KEEP` means scientific or runtime authority. `DEFER` means the item
may be removable, but current references or user-owned work prevent deletion.

| Path | Status | References / role | Decision | Reason |
|---|---|---|---|---|
| `model/bpsf.py` | tracked, modified | active model and tests | KEEP | Active interaction trunk |
| `model/encoders.py` | tracked | active model and geometry trainer | KEEP | Active encoders |
| `model/qpsmp_meta.py` | tracked, modified | train/evaluate scripts and tests | KEEP | Active zero/few-shot model |
| `model/cartesian.py` | tracked, modified | package export and tests | KEEP | Optional tested extension; inactive on current BindingDB bank |
| `model/bpsf.py` geometry head | tracked | geometry commands and tests | KEEP | Source-only supervision belongs with its pair-field producer |
| `model/bands.py`, `mathematical.py` | tracked | maintained contract tests | DEFER | Not active DTA runtime, but still tested |
| `model/config.py`, `runtime.py` | tracked | data/runtime scripts and tests | KEEP | Live non-model dependencies |
| `model/meta_operator.py` | tracked | maintained contract tests | DEFER | Compatibility API remains tested |
| `research/` | mixed current prototypes | many direct test imports | KEEP | Bulk movement would break research verification |
| `report/` | result lineage | current and failed experiment evidence | KEEP | Numerical provenance; failure is not redundancy |
| `archive/theory/` | untracked user migration | historical theory references | DEFER | Existing user-owned migration; no destructive dedupe authorized in this pass |
| `archive/legacy/retired_qpsmp/` | tracked, moved | historical architecture report | KEEP | Retired implementation needed to interpret old results |
| `archive/FORT/model`, `research`, `reports` | untracked predecessor snapshot | historical evidence | KEEP | Unique predecessor source and reports |
| `archive/FORT/.env` | untracked secret | no valid evidence role | DELETE, DONE | Plaintext API configuration |
| `archive/FORT/tmp/` | ignored temporary payload | no active references | DELETE, DONE | Third-party clones, duplicate downloads, smoke products |
| `archive/FORT/.git/` | nested repository | no runtime references | DELETE, DONE | Removed 733,869,294-byte VCS database after exact-path relocation and long-path cleanup |
| `archive/retired_research_20260811/ephemeral_quarantine/` | ignored cache | archive manifest marks reproducible | DELETE, DONE | Reproducible quarantine payload |
| root `.pytest_cache/`, `__pycache__/` | ignored caches | none | DELETE, DONE | Rebuildable runtime state |
| `.tmp_pytest_final_cleanup_20260815/` | ignored test temp | none | DELETE, BLOCKED | Original path was retired; child ACL denies the current execution account and no Python/pytest process owns it |
| `LLM/llm.json` | untracked local configuration | K3 gateway scripts | KEEP LOCAL | Contains local credentials; covered by `.gitignore`, never publish |
| `LLM/run_k3_framework_*.ps1` | untracked helper scripts | manual K3 review | KEEP LOCAL | Streaming and non-streaming gateway fallbacks |

No active dataset, split manifest, `RESULT.json`, retained checkpoint, model
source, or training/evaluation protocol was deleted in this cleanup.
