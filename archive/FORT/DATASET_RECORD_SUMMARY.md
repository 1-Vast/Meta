# 数据集记录与预处理审计

更新时间：2026-08-01

## 结论

当前数据已经按用途分为两类，并保留原始层：

| 用途 | 路径 | 状态 |
|---|---|---|
| 创新机制/模块测试 | `dataset/innovation_tests/a2s_validation_small.v1/` | `DEVELOPMENT_ONLY_H0_BLOCKED` |
| 正式 pKi 训练语料 | `dataset/formal_training/chembl37_pki_formal.v4/` | `FORMAL_PKI_CORPUS_READY_NATURAL_TAIL_BLOCKED` |
| 原始数据 | `dataset/public/chembl_historical/snapshots/chembl_37/chembl_37.db` | 只读保留 |

正式包可以用于 exact pKi 的 source training；它目前不能支持正式的 prospective natural-tail 评估，因为现有 ChEMBL 37 的 publication/patent year 不是完整 measurement timestamp，且候选 recipient 数量不满足预注册功效门槛。

## 数据来源

主数据源只有 ChEMBL 37 原始 SQLite 关系数据库：

- 原始文件：`dataset/public/chembl_historical/snapshots/chembl_37/chembl_37.db`
- 大小：30,480,314,368 bytes
- SHA-256：`4be13df3b68e25dcd0bff44bf094033b5aebe98f415acdc8c1cdf380e0c15142`
- 许可：ChEMBL CC BY-SA 3.0
- raw 清单：`dataset/raw/chembl_37/manifest.json`
- 原始 SQLite 只读打开，所有清洗和聚合都写入新的派生文件，不修改原始库。

BindingDB、Papyrus、KIBA、Davis、Metz 等没有进入正式 pKi 主监督集。pKd 只作为独立 auxiliary endpoint 保存，不能与 pKi 混合。创新测试包来自已经冻结的 `dataset/ready/a2s_validation_small.v1/`，只用于机制和模块开发验证。

旧候选包、重复构建产生的中间包、过期分类清单和 6 个失败 processing-run 记录已清理；当前仅保留原始 SQLite、封存的正式训练包和创新测试包。旧的 dual-cold registry 仍属于历史输入资产，不是当前正式模型输入。成功的 processing-run 记录继续保留作 provenance。

## 统一代码入口

所有数据预处理代码已整合到：

```text
D:\FORT\scripts\preprocess.py
```

该文件同时提供：

- ChEMBL 37 SQLite 只读抽取；
- exact Ki、censored Ki、pKi/pKd 标准化；
- target、parent compound、document、assay、assay-context 表构建；
- assay-context 重复测量聚合和噪声隔离；
- ECFP4/分子描述符/蛋白序列特征构建；
- metadata-only natural-tail D0 诊断、source/recipient closure 和报告；
- 创新测试包构建、A2S split/leakage 检查和 source-only normalization；原有 `preparetable`、`preparerows`、`preparevectors`、`normalizeligands` 接口。
- `verify_formal` 封存包完整性和泄露契约验证；`organize_datasets.py` 只调用该统一验证函数并写分类索引。

分类和封存校验入口为：

```text
D:\FORT\scripts\organize_datasets.py
```

它不会重新复制受保护的创新包，也不会从旧 registry 派生正式数据；会 fail-closed 检查 raw/manifest/产物哈希、标签约束、censored 隔离、context 唯一性、特征对齐、label-blind roster 和硬 overlap。

新版本重建命令必须使用新的版本目录，不能覆盖已封存 v4：

```powershell
python scripts\preprocess.py --output dataset\formal_training\chembl37_pki_formal.v5
python scripts\organize_datasets.py
```

A2S 创新验证包也使用同一入口：

```powershell
python scripts\preprocess.py --mode innovation --output dataset\processed\a2s_validation_small.v2
```

## 主集过滤规则

正式 pKi 主集来自以下严格条件：

- `standard_type = 'Ki'`；
- `standard_relation = '='`；
- `standard_units = 'nM'` 且 `standard_value > 0`；
- `standard_flag = 1`；
- `pchembl_value IS NOT NULL`；
- `data_validity_comment IS NULL`；
- `potential_duplicate = 0` 或 NULL；
- `confidence_score = 9`；
- `assay_type = 'B'`；
- `variant_id IS NULL`；
- `target_type = 'SINGLE PROTEIN'` 且 tax ID 为 9606；
- target 只有一个 component；
- 序列长度 50--5000，`X` 比例不超过 0.01；
- activity 和 parent 均为 `Small molecule`，parent canonical SMILES 存在；
- 排除 `Targeted Protein Degradation`。

标签自行计算：

```text
pKi = 9 - log10(standard_value_nM)
```

ChEMBL `pchembl_value` 只作一致性检查，误差阈值为 0.05。censored relation `<`、`<=`、`>`、`>=`、`~` 不进入 exact 回归集，单独保存。不同 assay context 不跨 context 平均；同 context 内保留 median、MAD、range、replicate_count 和 activity IDs。

## 数据漏斗与规模

| 阶段 | 行数 |
|---|---:|
| ChEMBL all activities | 24,527,044 |
| `standard_type = Ki` | 887,151 |
| exact relation | 607,253 |
| nM 且正值 | 605,786 |
| `standard_flag = 1` | 605,770 |
| pChEMBL 存在 | 588,370 |
| 无 validity comment | 588,032 |
| 无 potential duplicate | 531,732 |
| 严格 target/assay/compound/sequence QC 后 | 157,613 |
| pChEMBL mismatch quarantine | 0 |
| censored Ki audit | 24,558 |

正式包统计：

- exact Ki activity measurements：157,613；
- exact Ki assay-context rows：152,427；
- 主 context（`pKi_range <= 0.5`）：151,411；
- high-noise contexts：1,016，其中 range > 1.0 的有 467；
- parent compounds：82,646；targets：1,098；documents：7,028；assays：17,704；contexts：15,807；
- pKd auxiliary measurements：22,056；
- pKi 分布：min 2.0000，q01 4.3188，median 7.1024，mean 7.1005，q99 10.0000，max 11.0000，std 1.3047。

## 输出文件

正式包的权威入口是：

```text
dataset/formal_training/chembl37_pki_formal.v4/manifest.json
```

主要内容：

```text
canonical/
  pki_measurements_exact.parquet
  pki_measurements_censored.parquet
  pki_measurements_context_aggregated.parquet
  pki_measurements_context_main.parquet
  pki_measurements_context_high_noise.parquet
  quarantine_label_mismatch.parquet
components/
  targets.parquet
  compounds.parquet
  documents.parquet
  assays.parquet
  assay_contexts.parquet
auxiliary/
  pkd_measurements_exact.parquet
  pkd_measurements_context_aggregated.parquet
features/
  ligand_features.npz
  target_sequences.json
reports/
  preprocessing_report.md
  leakage_audit.md
  natural_tail_roster.md
  overlap_matrix.json
  funnel.json
  assay_noise_audit.json
```

配体特征为 2048-bit ECFP4 加 10 项描述符：molecular weight、logP、HBD、HBA、TPSA、rotatable bonds、ring count、formal charge、heavy atoms、fraction CSP3。formal package 只保存未归一化特征；归一化必须在已封存的 `source_meta_train` 上拟合。

## Natural-tail 与 closure

roster 只使用 target/parent/document/context 数量和 `document_year`，代码明确记录 `labels_used_for_roster_selection = false`。候选 cutoff 统计为：

| cutoff | source candidates | recipient candidates | closed source candidates |
|---:|---:|---:|---:|
| 2018 | 262 | 5 | 260 |
| 2019 | 263 | 4 | 261 |
| 2020 | 266 | 2 | 265 |
| 2021 | 267 | 1 | 267 |
| 2022 | 267 | 1 | 267 |
| 2023 | 269 | 0 | 269 |

选定诊断 cutoff 为 2022，但只作为 metadata-only diagnostic，不是正式 natural-tail test。source/recipient 闭合后的硬交集均为 0：target UID、component、document、parent。support/query parent、support/query scaffold 和 dev/test family 尚未构成，因为 recipient 功效不足；这些轴不能宣称已经完成。

## 数据泄露审计

### 已通过的保护

1. raw SQLite 只读，模型输入策略明确禁止 raw database。
2. natural-tail roster 的 target 选择只读取 metadata counts 和年份，不按 affinity、label variance 或模型结果选择。
3. pKi 自行计算，pChEMBL 只做固定阈值一致性检查；mismatch 不静默删除。
4. assay-context 只在相同 context 内聚合，不把不同实验直接按 target-compound 平均。
5. `normalizeligands` 必须显式接收 `trainrows`，统计量只从这些训练行拟合。
6. formal manifest 和分类入口会检查 exact relation/unit/value、activity 唯一性、censored 隔离、context key、特征 UID 对齐、label-blind roster、source/recipient hard overlap 和 raw 哈希。
7. source/recipient closure 审计结果：target UID、component、document、parent overlap 均为 0。

### 必须遵守的边界

formal exact corpus 是完整监督语料，不是自动完成的 train/dev/test split。它包含标签是正常的，但下游训练必须先使用独立、冻结的 source/recipient split；不能把 `canonical/pki_measurements_context_main.parquet` 同时用于训练、调参和最终评估，不能用全 corpus 的 pKi 分布拟合 scaler，也不能把 report 中的 label summary 当模型输入。

当前 natural-tail roster 未通过功效和完整时间/lineage 门槛，所以不能把现有 recipient candidate 当正式 test。创新包也只允许机制/模块开发，不用于 natural-tail 科学结论。

因此，结论是：**预处理代码在既定 manifest/model-input contract 下已通过数据泄露审计；下游若绕过 contract 直接把完整 canonical 表用于评估，仍然可能人为造成泄露，代码不能替下游调用者承担这一误用。**

## 验证记录

- `python -m py_compile scripts/preprocess.py scripts/organize_datasets.py`：通过；
- `python scripts/organize_datasets.py`：`verification.status = PASS`；
- v4 manifest 中 29 个派生文件 hash 全部通过；raw size/hash 全部通过；
- exact pKi `activity_id` 全唯一，157,613 行均为 `Ki`/`=`/`nM`/正值，pKi 公式逐行一致；
- censored 24,558 行全部隔离；主 context 无重复 key；ECFP4 与 compound UID 完全对齐；
- 统一入口相关测试集合：15 passed；其中 `tests/core/test_htl_topology.py`：2 passed；
- 全量 pytest 仍有既有环境问题：3 个收集错误缺少 `parasail`，7 个 historical projection 测试在现有 SQLite authorizer 环境下 `PRAGMA integrity_check` 返回 `not authorized`。这些失败不涉及本次统一预处理逻辑。

详细机器可读记录：

- `dataset/formal_training/chembl37_pki_formal.v4/manifest.json`
- `dataset/formal_training/chembl37_pki_formal.v4/reports/leakage_audit.md`
- `dataset/formal_training/chembl37_pki_formal.v4/reports/preprocessing_report.md`
- `dataset/registry/DATASET_CLASSIFICATION_REPORT.v2.json`

## 2026-08-01 target-ID-unseen formal-track update

The preceding natural-tail warning remains valid, but it is no longer the only
experimental track. A separate, explicitly scoped target-ID-unseen chain is now
sealed:

- `dataset/formal_training/a2s_d0r_roster.v3`: PASS; 206 source targets, 63
  recipient targets, 55 independent homology components, five document-ordered
  support draws at k=1/3/5;
- source/recipient target UID, accession, document, parent and assay overlaps
  are zero; support/query parent and document overlaps are zero;
- nine recipient targets are homology-warm. The admissible full-cohort claim is
  strict target-ID unseen; the 54-target homology-cold subgroup is reported
  separately;
- every support/query row freezes a deterministic `measurement_uid` selected
  from `pki_measurements_context_main.parquet`; label values were not read;
- `dataset/formal_training/a2s_cmal_episodes.v2`: 30,123 label-blind episodes,
  source meta-splits intact by homology component, three target-mismatched
  support mappings, content SHA-256
  `d8e9a259f594db293c4e46779cb716b8e52a828a989ae20cdbc5571805877f9b`;
- all package audits report zero parent/document/time/nesting/component-split/
  measurement/same-target-negative violations.

This admits external A2S-CMAL training for the target-ID-unseen estimand only.
It does not admit a natural-tail, all-recipient homology-cold, or multi-release
historical claim. No formal recipient training has run locally.
