# CIIP-S1 执行提示词（信号分解与估计量阶梯）

你是 `D:\MetaSieve` 项目的独立研究、审计和实现 agent。你的任务是执行
`report/research_ideas/CIIP_SUCCESSOR_STAGE_RESEARCH_PLAN_20260820.md`
第 9 节定义的 **CIIP-S1（Signal Decomposition & Estimand Ladder）** 阶段。
该计划文档是本任务的唯一科学权威；本提示词只是执行纪律。两者冲突时以
计划文档为准；计划文档未覆盖的实现自由度必须先以日期化 addendum 冻结，
再实现。

你的目标不是训练一个大模型，而是把 CIIP-1A 的失败分解为可分别裁决的
三个命题：

```text
A: mutation-site 表示能否预测突变整体效应（T0/T0m）
B: mutation-erased 蛋白上下文能否提供蛋白状态信息（T1，Form-2）
   —— 所有 B 分支必须标注 "mutation-free but not necessarily biologically causal"
C: 蛋白表示能否改变 ligand ranking / ligand-specific response
   （T1 = C_total，T2 = parent 残差化后的 C_sharp，T3 = rank 形式）
```

端点永远是 functional percent inhibition；禁止改名 Ki/Kd/pK/ΔΔG。

## 0. 开始前的强制核对（不通过则停止并报告）

1. `git status`、进行中进程、`.dsh-jobs`：确认没有并发进程正在写
   `tools/research/stageCIIP_signal_decomp_*`。已知并发工作
   `tools/research/stageCIIP2_olr_potential_20260820/` 与
   `report/research_ideas/ciip/CIIP2_RESEARCH_REPORT_20260820.md`
   属于另一阶段：**只读参考，禁止修改、覆盖或在其目录写入**。
2. 完整阅读：计划文档全文；`PREREGISTRATION_STAGE1_CONTROLS.md`、
   `CONTROL_ADJUDICATION.json`、`CONTEXT_PROPAGATION_REPORT.md`、
   `STAGE1_COLLAPSE_AUDIT.md`、`REPORT_2X2_DIAG.md`、
   `report/research_ideas/CIIP_CONTROL_FAILURE_IMPLEMENTATION_AUDIT.md`、
   `DATA2X2.json` 的 coverage_bias 块；以及 CIIP-2 报告的 §2.4
   （独立复算诊断，与你的 S0 审计互为交叉验证）。
3. 若 CIIP-2 阶段已经产生训练结果：停止，写出重叠分析报告，请求
   治理裁决（S1 先行 vs CIIP-2 增补），不得自行并行训练。

## 1. 阶段目录与冻结程序（任何计算之前）

1. 新建 `tools/research/stageCIIP_signal_decomp_<YYYYMMDD>/`。只允许写
   入该新目录。`model/`、生产 `scripts/`、`dataset/`、既有
   `stageCIIP_*` 目录一律只读。
2. 把计划文档 §9.0–9.4 逐字抽取为 `PREREGISTRATION.md`，补全：
   阶段目录名、日期、环境、冻结输入 SHA-256 清单（见下）、随机流命名
   空间、以及 **S0 功效表出来之前的临时阈值占位符**。
3. 计算 `PREREGISTRATION.md` 的 SHA-256，写入 `PREREGISTRATION_SHA256.txt`
   并记入 `commands.jsonl`。**此后任何阈值/规则变更只能以新的日期化
   addendum 文件追加，且必须在看到对应结果之前冻结。**
4. 冻结输入（沿用既有 SHA，逐项重新校验后记录）：
   - `stageCIIP_potential_bridge/DATA1A.json`（1c2b92df…）
   - `stageCIIP_potential_bridge/DATA1A.npz`（40f69509…）
   - `stageCIIP_potential_bridge/DATA2X2.json` / `DATA2X2.npz`
   - `stageX_csc_signal/stageX0c_measurement_qualification_20260818/q1_esm_cache.npz`（c8b59e33…）
5. meta_test（BindingDB 封存集）**永远不读**；Duong-Ly test pair 标签
   只允许用于冻结 arm 的最终评估，禁止用于拟合、归一化、特征构造、
   模型选择、检索或 checkpoint 选择以外的任何用途（checkpoint 选择
   只能用 val）。

## 2. S0：只读审计（全部通过前禁止任何拟合）

按顺序完成并落盘 `S0_AUDIT.json` + `S0_AUDIT.md`：

1. **覆盖审计**：65 admitted pairs 的 parent/position/construct 长度表；
   16 个未覆盖 pair 的排除原因（pos > 1020）复核；covered 49 的
   split 归属（32/8/9）与 DATA1A `pair_split` 逐位一致。
2. **parent 重叠**：split × parent 关联表；WT-row 共享图（哪些 pair
   共用同一 WT 行）；确认 9 个 test pair 的 parent 全部有同 parent
   train pairs（F9 可定义性）。
3. **突变坐标**：与 Q0B_MAPPING_AUDIT 逐对复核 65 个映射；alias ledger
   检查；任何不一致立即停止。
4. **配侧重叠**：per-pair 公共配体数分布；per-ligand 行覆盖。
5. **assay 语义**：端点命名空间检查（只出现 % inhibition）；越界细胞
   普查（复算 ≈23.0%）；WT panel 饱和度（复核 CIIP-2 §2.4f 的
   ~90.9% 中位与 ~84 个中区配体）；浓度元数据从源补充表只读摘录。
6. **censoring**：确认无删失标注 → 记录"interval-censored 公式不可
   识别"为数据限制。
7. **复算计划文档 §4 的探索性诊断**（train+val 40 pairs，keyed rng）：
   能量分解（主效应 ≈10.5%）、同 parent 跨突变一致性（≈0.442 /
   WT 残差化 ≈0.406）、parent-profile LOPO 上限（≈0.579/0.326）、
   配体全局效应（≈0.060）、parent 残差化一致性（≈-0.28，注明 LOO
   机械偏负）。与计划文档数值的偏差必须解释；数量级不符立即停止。
8. **功效表**：用 train 标签做 keyed 模拟，给出每个冻结对比在
   n=9 test pairs / 6 clusters 下的最小可分辨效应（MDE）；据此填出
   各命题的点估计阈值，写进**日期化 addendum（S1 训练前冻结）**。
9. **erasure 推理**（CPU，唯一允许的新计算）：对 49 covered pairs，
   把 verified 位点替换为 X（WT、variant 两侧），**断言替换后字符串
   完全相等**，用本地 ESM-2-150M 推理，断言每对最大 embedding 差
   ≤1e-5；新 cache 写入本阶段目录并 SHA-256 钉住。
10. **泄漏审计输出**：发现的任何通道必须在 S1 前中和或写成结论限制。

S0 出口门：1–8 完成、无未解决泄漏、erasure cache 有效。任一失败：
停止，写 `FAILURES.md`，不得转入 S1。

## 3. S1：估计量阶梯与表示臂

### 3.1 目标（全部来自同一合法标签）

```text
T0  = d_vl = y_var - y_wt                    # 命题 A（含主效应）
T0m = mean_l(d)                              # 命题 A 纯量（突变严重度）
T1  = c_vl = d - mean_l(d)                   # C_total（CIIP-1A 估计量）
T2  = c - parent_profile^{(-pair)}           # C_sharp（突变特异残差）
T3  = 对 c 的 within-pair 排序目标            # C_total 的尺度自由形式
```

T2 的 parent-profile 必须 cross-fitted：train 内 leave-pair-out；
val/test 的残差只用 train parents 拟合的 profile。评估中的 pair 的标签
永远不得进入它自己的 nuisance。

### 3.2 表示臂与模型形式

```text
F1 mutation-site local ESM（oracle，radius-6）      Form-1/Form-2
F2 mutation-erased（双侧 X 替换；Form-1 下结构恒零） Form-2
F3 full-sequence ESM mean（cache 池化）              Form-2（可部署）
F4 KLIFS pocket one-hot                              Form-2（可部署，结构受限）
F5 family-preserving shuffle（ keyed，parent 内置换） 控制
F6 random local window（重解释：上下文传播值测量）    控制
F7 ligand-only（蛋白置零）                            零地板
F8 protein-invariant constant shift                   零地板
F9 parent-profile predictor（cross-fitted 均值）      上限参考
F10 外部特征（conservation/MSA/PSS）—— 预期不可用，
    只有数据与许可核实后才允许；UniRef 快照当前缺失
```

- Form-1：冻结的 rank-8 / hidden-64 potential，g = s(Pv,L) − s(Pw,L)，
  只做 T1/T3；AdamW 1e-3、wd 1e-4、200 epochs、batch 512、grad clip 10。
- Form-2：小 MLP（hidden 64）回归 (protein-features, ECFP4) → 目标；
  同优化器/预算；只做信息上限/诊断，**永远不是生产机制**。
- 禁止 ridge/闭式解/pseudoinverse、test-time 梯度、query 标签。

### 3.3 必备负控（与所有 arm 同一指标管线）

ligand-label permutation（pair 内 keyed）；same-parent wrong-mutation；
family-preserving shuffle；random window；protein-invariant shift；
erasure null（Form-1 下必须严格输出零）。

### 3.4 执行阶梯（不可跳级）

1. **CPU smoke**：2 个 pair、5 epochs、≤30 分钟；结构测试全绿
   （identity=0、反对称、erasure null 严格零、F7/F8 在地板上）。
2. **单 seed（seed 1）**：全部 arm + 全部控制。
3. **多 seed（{1,2,3}）**：仅当单 seed 的结构与全部负控行为符合冻结
   预期（地板在地板、erasure 在零、permutation 无系统增益）。
4. 任何负控失效：停止该分支，记录，不得"修到通过"。

### 3.5 指标合同（每 pair 每 arm）

nonconstant flag 与 N_nonconstant/N_total、var_true、var_pred、
scale_ratio = sqrt(var_pred/var_true)、centered MSE、centered R2、
OLS slope、dead-zone（10 单位）sign accuracy、Spearman（常数时
**undefined，绝不记 0**）、N_rank_evaluable/N_total、per-parent 聚合。
所有 arm 间对比：paired、parent-cluster bootstrap（2000 draws、keyed）、
leave-one-parent-out 符号稳定性；bootstrap 均值永远不作点估计。

### 3.6 冻结裁决（lo2.5 > 0 且 LOPO 符号稳定；点阈值来自 S0 addendum）

```text
S1-PASS-B        F2(erased, Form-2) 在 T1 上超过 F7/F8 地板
S1-PASS-A        F1 在 T0/T0m 上超过 F2-erased
S1-PASS-Csharp   F1 在 T2 上超过 F2-erased（配对对比）
S1-NULL-ALL      所有格子在地板上
S1-UNRESOLVED-*  CI 跨零且不满足 NULL-ALL（合法终局，禁止继续加模型抢救）
```

## 4. 授权逻辑（裁决后唯一允许的动作）

| 结果 | 授权 | 关闭 |
|---|---|---|
| S1-PASS-B | 起草方向 2（mutation-free family conditioner，parent-disjoint）预注册 | — |
| S1-PASS-A | 报告突变整体效应结果（禁止配体条件化措辞） | — |
| S1-PASS-Csharp | 方向 3（erasure 反事实通道）范围评审 | — |
| S1-NULL-ALL | 更新边界文档 | Duong-Ly 交互路线（两级）最终关闭 |
| 任何结果 | CIIP-1B / BindingDB bridge / 生产集成 / 可部署表示声称 **保持 NOT AUTHORIZED** | — |

## 5. 绝对禁令

1. 不修改 `model/`、生产 `scripts/`、任何冻结工件、CIIP-2 目录。
2. 不把 oracle mutation-coordinate ESM 接入任何部署路径。
3. 不以 correct-vs-random-window 差异声称生物机制。
4. 不以单 seed / 单 parent / 少数 pair 声称成功。
5. 不用更大 backbone、更多层数或更高预算掩盖不可识别——本阶段只许
   缩小问题，不许扩大模型。
6. 不把 % inhibition 写成 Ki/Kd/pK；不把上下文传播幅度当预测值；
   不把失败解释为"生物学上没有蛋白条件化信号"。
7. 不用 test 标签做拟合/归一化/选择/检索/特征构造；不读 meta_test。
8. 不移动已冻结阈值；不加无预注册的新 arm；不删除失败的控制结果。

## 6. 交付物

`PREREGISTRATION.md` + `PREREGISTRATION_SHA256.txt`；S0 addendum
（阈值，S1 前冻结）；`S0_AUDIT.json/md`；erasure cache + SHA；
`RESULT.json`（机器可读全指标）；`REPORT.md`（含裁决表与授权
状态块）；`commands.jsonl`；结构与数据契约测试；`SHA256SUMS`；
`FAILURES.md`（如无失败则显式说明）；最后同步 `history.md`、
`task.md`、`report/EVIDENCE_LEDGER.md`（只追加，不改写旧结论）。

最终报告必须逐条回答：

```text
命题 A（突变整体效应）:        SUPPORTED / NOT SUPPORTED / UNRESOLVED
命题 B（mutation-free 上下文）: SUPPORTED / NOT SUPPORTED / UNRESOLVED
   —— 若 SUPPORTED 必须带 "mutation-free but not necessarily
      biologically causal" 标注
命题 C_total / C_sharp:        SUPPORTED / NOT SUPPORTED / UNRESOLVED
parent-profile 上限可达性:     各臂与 F9 的差距表
授权状态:                      按 §4 表逐项列出
```

目标是对"是否存在不依赖 mutation-coordinate oracle、能改变未见
protein/family 上 ligand ranking、且经 matched controls 证明不是
protein main effect / mutation severity / assay batch / 随机上下文传播
的可迁移 protein-conditioned interaction"给出**可裁决的分解答案**，
而不是制造一个通过数字。
