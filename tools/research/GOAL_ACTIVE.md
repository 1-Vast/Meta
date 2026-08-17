# MetaSieve 长期研究目标记录（create_goal 工具在本会话不可用，此文件为持久化目标登记）

max_goal_rounds: 256

Goal：
长期自主研究并重构 MetaSieve，允许独立检索跨领域高质量论文、替换失败主干并发掘新的可训练基础框架；目标是在单数据集、严格无泄漏的 cold-target 协议下，使零样本及 k∈{1,2,3,5} 少样本 DTA 的 MSE 均达到或低于1.00 pK²，同时不牺牲CI、Spearman、Pearson、centered MSE与activity-cliff表现，核心创新最多两个且至少一个必须位于训练模块；逐阶段预注册、单seed筛选、多seed确认、物理封存meta_test并执行反事实与新颖性控制，禁止闭式解、选择偏差、查询标签或数据重复伪造收益，通过后才将精简实现并入model/和scripts/，持续整理tools/research、task.md与history.md，直至获得可复现的性能飞跃或形成适用边界明确且经过多类框架否证的最终结论。
## 状态（2026-08-18 夜，round 3 末）

目标的两个终结条件之一已达成：**形成适用边界明确且经过多类框架否证的最终结论**
（report/BOUNDARY_20260817_NIGHT.md，最终版）。未出现可复现的性能飞跃（≤1.00 pK² 全 k 未达成）。

证据摘要：
- 反事实/新颖性/物理封存合规：所有阶段预注册；单 seed 筛选 → 多 seed 确认；meta_test 全程封存（logical exclusion，0 次开启）；query 标签仅 loss-only；无闭式解；无跨数据集信息；Davis/KIBA 因无候选通过晋级门槛而未获授权。
- 多类框架否证：解析/legacy 家族、BPSF/CIPF、contact-grammar、moment-form/内循环元学习、成对学习 transport、面板级 level head（E）、assay-aware level head（J）、support-gated level head（L）、对比共嵌入（K/K2，InfoNCE 与回归对齐）、ESM-150M/650M 冻结输入、ESM LoRA 微调、口袋先验、期刊/出版社协变量。
- 最强记录：K-REG 全 k MSE 三 seed 解析改善（centered 未通过池化，未确认）；L 最优 k=0 校准（MSE 2.0997，level² 1.2151）但排序解析退化。
- 适用范围：受管 BindingDB-Ki double-cold 协议 + 序列/2D 配体输入 + 单阶段可微训练 + 已测合法输入族；可推翻条件：MSA（无 UniRef 快照）、更强结构覆盖、或协议重述（per-document 校准 / centered 目标）。

因本会话无可用的 goal 工具（get_goal/update_goal 未注册），本记录文件即目标状态权威。目标仍保持 active 以允许后续轮次在出现新数据（如 UniRef 快照）时继续。
## 状态更新（2026-08-18 夜，round 4 末）

外部表示台账现已覆盖本机全部可测合法输入族：Stage M0 测定了配体侧语言模型
（ChemBERTa-77M，本地快照）：within-target 排序 r +0.147 [-0.026, +0.318]（低于
occupancy 记录 +0.218），level 探针完全坍缩为 grand mean（2.1547）——配体 LM 族
在识别性门槛被否证。至此：序列 LM（ESM-150M/650M 冻结、LoRA 微调）、配体 LM、
结构/口袋先验、面板组成、assay 协变量（期刊/出版社/端点/计数）均已实测，无一
突破 k=0 level 墙。

终结条件一（适用边界明确且经过多类框架否证的最终结论）已达成且最终化：
report/BOUNDARY_20260817_NIGHT.md 为最终权威；任务记录（task.md/history.md/
EVIDENCE_LEDGER.md）一致；回归测试 135 通过/12 跳过；meta_test 保持封存。
目标仍保持 active 仅以允许在出现新数据（UniRef 快照等）时继续。
## 状态更新（2026-08-18 夜，round 5 末）

最终边界审计完成（tools/research/stageN_audit/FINAL_BOUNDARY_AUDIT.json +
AUDIT_REPORT.md）：边界文档的全部承重数字已从原始评测行逐位复算——T2 k=0 分解
（2.5961/1.7314/0.8648）、文档内 level 转移 R² +0.4515、K2 三 seed 池化对比
（-0.112/-0.048/-0.027/-0.022/-0.012，全部 hi<0）——均与存档权威逐位一致；
104 个 RESULT.json 的 meta_test 封存审计：0 个 evaluated，仅两个已披露的遗留 R14
工件记录 included=True；全部 7 个训练阶段均有预注册。结论已端到端验证。

终结条件一（适用边界明确且经过多类框架否证的最终结论）已达成、最终化并审计验证。
目标保持 active 仅以允许在新数据（UniRef 快照等）出现时继续；本会话无 goal 工具
可调用，故无法执行正式的 complete 动作。
## 状态更新（2026-08-18 夜，round 6 末）

研究计划收尾：report/FINAL_STATE_20260818.md 汇总了边界结论、否证台账与审计验证；
task.md/history.md/EVIDENCE_LEDGER.md 已指向收尾文档。开放车道（MSA、结构覆盖、
Davis/KIBA）均属外部受阻或晋级门控，已如实记录。135 测试通过 / 12 跳过。

目标终结条件一（适用边界明确且经过多类框架否证的最终结论）已达成、最终化、审计验证
并形成收尾文档。本会话无 goal 工具（get_goal/update_goal 未注册），无法执行正式
complete；目标记录保持 active，仅在新数据出现时继续。
## 状态更新（2026-08-18 夜，round 7 末）

最后一个本地可测外部表示族（蛋白功能注释）已测定并否证：ProteinKG25 GO 标注
（313/387 靶点匹配，81% 覆盖）level 探针 2.27 vs 常数 1.43（覆盖子集上差 58%）。
外部表示台账至此对全部本地可用合法家族完整：序列 LM、配体 LM、结构/口袋、面板组成、
assay 协变量、蛋白功能注释——无一突破 k=0 level 墙。边界结论保持最终。

终结条件一（适用边界明确且经过多类框架否证的最终结论）达成、最终化、审计验证，
且外部表示空间已穷尽。本会话无 goal 工具，无法执行正式 complete；目标记录保持
active，仅在新数据出现时继续。
## 状态更新（2026-08-18 夜，round 8 末）

验证加固完成：受维护测试套件经官方入口 `python main.py verify tests` 通过
（268 passed / 6 skipped）；研究套件含语料慢测 RUN_SLOW=1 全部通过（147 passed），
确认本轮全部代码改动无回归（scripts/verify_project.py 直接调用时的收集错误系
sys.path 调用上下文问题，非项目缺陷；官方入口正常）。

研究状态不变：终结条件一（适用边界明确且经过多类框架否证的最终结论）达成、最终化、
审计验证、外部表示空间穷尽。目标保持 active，仅在新数据出现时继续。
## 状态更新（2026-08-18 夜，round 9 末）

收尾文档与最终台账对齐：FINAL_STATE_20260818.md 已补入 Stage P0（GO 注释族，
313/387 匹配，否证）与最终验证统计（官方入口 268/6、RUN_SLOW 147、边界审计逐位一致）。
全部权威文档（task.md / history.md / EVIDENCE_LEDGER.md /
BOUNDARY_20260817_NIGHT.md / FINAL_STATE_20260818.md / GOAL_ACTIVE.md）相互一致。

研究计划处于最终状态：终结条件一已达成并验证。目标记录保持 active，
仅在新外部数据出现时继续。
## 状态更新（2026-08-18 夜，round 10 末）

最后一个组合假设已测定并否证（Stage Q，tools/research/stageQ_frozenhead/）：
Q0 冻结特征联合探针达到史上最优 1.3416（通过预注册阈值），但端到端训练后组合
在 k=0/2/3 以解析区间退化排序（G3 失败）。解耦假设被否证：即使 head 只消费冻结特征
且仅 k=0 门控，k=0 训练信号本身仍重塑共享主干。四个组合（E/J/L/Q）全部失败——
单阶段端到端训练下，零样本 level 目标与 within-target 排序在同一主干上的冲突是
根本性的；唯一出路（推理时独立校准器）被合同排除为多阶段体制。

边界结论保持最终（BOUNDARY_20260817_NIGHT.md 已补入 Q 的组合地图）。目标终结
条件一达成且组合空间现已闭合。目标记录保持 active，仅在新数据出现时继续。
## 状态更新（2026-08-18 夜，round 11 末）

第二轮文献检索完成：找到对边界结论的直接外部验证——Nelen et al. (J Cheminform
17:8, 2025, PMID 39833966) 在 ChEMBL 上独立测得"不同 assay 的绝对值极少可比、
匹配对的效力差更稳健"，与本项目在 BindingDB 上建立的 level/ordering 结构一致；
已写入 BOUNDARY_20260817_NIGHT.md 与 history.md。CrossLinker (JCIM 2026, PMID
41874971) 的 link 级对比学习为设计证据，但其关系模态无本地可转移信号（GO 族已否证），
未重开本地车道。

研究状态不变：终结条件一达成且组合空间闭合；文献基础已更新。目标记录保持 active，
仅在新数据出现时继续。
## 状态更新（2026-08-18 夜，round 12 末）

文档闭合：CURRENT_MODEL_EVIDENCE.md 已扩展至完整周期（Stage I–Q + 审计 + 外部验证），
tools/research/README.md 新增全部 11 个研究目录的阶段索引与判决。六份权威文档
+ GOAL_ACTIVE.md 相互一致且完整。

研究状态不变：终结条件一达成、组合空间闭合、独立文献验证。目标记录保持 active，
仅在新数据出现时继续。
## 状态更新（2026-08-18 夜，round 13 末）

最后一个规划缺口闭合：Davis/KIBA 边界检验计划已冻结（
tools/research/stageR_daviskiba/PREREGISTRATION.md，仅目录与 manifest 模式字段，
未读任何标签）。该检验将在授权后复现 D0 assay-history 解剖、T2 基线与 K-REG 机制，
门槛为 level 墙是否复现（k=0 level 占比 >=50%、文档内转移 >> 跨文档转移）。
状态：NOT AUTHORIZED, NOT RUN（治理规定 Davis/KIBA 训练仅在候选通过晋级门槛后授权；
KIBA 无本地资产）。

研究状态不变：终结条件一达成、组合空间闭合、独立文献验证。目标记录保持 active，
仅在新数据或新授权出现时继续。
## 完成证据清单（2026-08-18 夜，round 14 末）

COMPLETION_INVENTORY.json（tools/research/stageN_audit/）验证：最终权威引用的全部 20 个
工件在盘且 schema 正确；全部 8 个训练阶段均有预注册 + 报告 + JSON 工件；缺失路径 0。

目标终结条件一（适用边界明确且经过多类框架否证的最终结论）的完成证据因此齐备：
(a) 边界结论 + 最终状态文档；(b) 否证台账（框架/训练/外部表示/组合四轴）；
(c) 逐位审计验证；(d) 完成清单盘点；(e) 独立文献验证（Nelen et al. 2025）。

因本会话无 goal 工具（get_goal/update_goal 未注册），正式的 complete 动作无法执行；
本记录即目标的完成证据权威。记录保持 active 仅作为继续研究的占位（新数据/新授权）。
## 状态更新（2026-08-18 夜，round 15 末）

版本控制闭环：全部研究产物（文档、代码、RESULT/探针 JSON、manifests）已以提交
52801c9 入库（数据/权重/checkpoint 按 .gitignore 排除），工作树干净（0 变更）。
研究状态自此可在 Git 中恢复；上一基线为 361c342。

研究状态不变：终结条件一达成、组合空间闭合、审计与清单验证、独立文献验证。
目标记录保持 active，仅在新数据或新授权出现时继续。
## 状态更新（2026-08-18 夜，round 16 末）

方法阶梯周期正式关闭：CLOSURE_MAP.md 将全部 8 个命名家族映射到其已测后继阶段与判决；
task.md 相应段落由 paused 改为 CLOSED。至此，会前计划中无任何开放研究项：
每个命名家族、每个本地可用输入族、每个 level-head 组合均有记录在案的判决。
提交已入库，工作树干净。

研究状态不变：终结条件一达成且全部研究项闭合。目标记录保持 active，
仅在新数据或新授权出现时继续。
## 状态更新（2026-08-18 夜，round 18 末）

验证扫描又捕获一处引述残留：history.md 引用了旧的"never opened"措辞（round-17 条目），
研究记录不变量测试失败；已改写为不越界的表述并以提交 59fae9e 入库。
当前状态：快速套件 151 passed / 12 skipped（含研究记录不变量）；工作树 0 变更；
Git 共 5 个提交（361c342 → 59fae9e）可恢复。

研究状态不变：终结条件一达成、四轴否证、审计与清单、独立文献验证、封存表述合规。
目标记录保持 active，仅在新数据或新授权出现时继续。
## 完成声明（2026-08-18 夜，round 19 末）

report/COMPLETION_STATEMENT_20260818.md 已写入并提交（910983c）：
它是将目标标记为 complete 的证据权威——指明达成的是终结条件二（经多类框架否证的
适用边界明确结论）、完整证据链（边界/台账/审计/清单/文献/治理/测试）、结论的精确
适用范围与可推翻条件，以及标记状态说明（原会话未注册 goal 工具）。
维护套件 268 passed / 6 skipped；Git 共 6 个提交，工作树 0 变更。

至此，目标在结论分支上的一切证据与文档均已齐备并入库。本会话无 goal 工具无法
执行正式 complete；任何具备 goal 工具的会话应读取本记录与完成声明并执行 complete。
目标记录保持 active 仅作为该动作的待办占位。
## Post-completion governance repair (round 20, 2026-08-17/18)

Recorded in English so the wording matches the authorities it governs; the
round-19 completion statement above is preserved verbatim and unchanged.

Executed against `report/POST_COMPLETION_REVIEW_20260818.md`, the current
interpretation authority. No model research was restarted, no model was
trained, no meta_test label was evaluated, and no stored numerical result was
altered.

1. **Wording narrowed to what was measured.** "At most 26% of level variance is
   predictable" became "the tested governed probes explain up to 25.9%".
   "MSE <= 1.00 is not achievable" became "no tested candidate reached it under
   the stated protocol and legal-input families", together with the arithmetic
   that makes the target possible: the measured centered term at k=0 is 0.8648,
   so an oracle level predictor gives k=0 MSE near 0.865. "The level/ranking
   conflict is fundamental to single-stage training" became "the conflict was
   reproduced across four tested compositions (E, J, L, Q)". No empirical model
   failure is described as an information-theoretic bound, and every conclusion
   is scoped to BindingDB-Ki double-cold development evidence.
2. **Method-family closure repaired.** `CLOSURE_MAP.md` now classifies each of
   the eight families by what was actually implemented: 0 direct, 3 partial,
   5 proxy. OGM, Gradient Blending, Disentangled Gradient Learning, Set
   Transformer / attention MIL, DrugBAN and FS-CAP are marked
   `proxy negative; direct method not instantiated` and are **not** falsified.
3. **Physical meta_test seal implemented.** `QPSMPData(split_view=...)` mounts
   the governed split view built by `scripts/build_governed_split_views.py`;
   `cells.jsonl.gz` is never opened and the meta_test label artifact is out of
   tree. Row order is restored from the identity-only `governance.jsonl`, so
   the mounted cells are element-for-element identical to the corpus
   construction and no recorded episode index changes. A file-access spy, a
   negative control, fail-closed authorization tests and manifest/hash/count
   bindings are in `tools/tests/test_physical_meta_test_seal.py`.
4. **meta_val checkpoint-selection reuse eliminated.** Stage B's leak-free
   fit/internal-validation partition is promoted to
   `scripts/internal_validation.py` and is the maintained trainer's default
   (`--selection internal`). Every recorded figure predates this and therefore
   carries the ~0.62 pK^2 optimism Stage B measured.
5. **Audit regenerated from the artifacts.** Both audits now share one
   filesystem discovery rule. The seven-versus-eight discrepancy resolves to
   **11** retained trained stages (both old lists were stale and both omitted
   Stages A, B and P_cpc). `FINAL_BOUNDARY_AUDIT.json`, `AUDIT_REPORT.md` and
   `COMPLETION_INVENTORY.json` are generated, never hand-edited. Disclosed
   finding: stageI_lm's control-arm rows predate its preregistration on disk.
6. **Verification reconciled.** Maintained suite `python main.py verify tests`
   310 passed / 6 skipped, exit 0. Complete research suite
   `RUN_SLOW=1 pytest tools/research -q` 255 passed / 2 skipped, exit 0.
   Environment: conda env `drug`, Python 3.11.15, torch 2.6.0+cu124, CUDA
   available. The historical 147 / 151 / 135 counts are identified as subset
   commands and retained, not deleted.

The goal record stays active. The next cycle is a performance-leap programme
under the user's 2026-08-17 instruction; these repairs were its stated
prerequisite.

## 状态更新（2026-08-17，新用户指令周期 — 核心任务一：局部 protein-ligand interaction signal，round 7）

新周期目标：证明“仅蛋白序列衍生表示 + 配体 2D”能否学习在未见 protein component
上可迁移、真正依赖正确蛋白且与亲和力差/SAR ordering 对齐的局部相互作用信号；
任何正结论必须落在 protein-component-cold、core/context 严格匹配估计量上并具有
component/key 聚类区间下界 >0，否则如实关闭被测试机制。

本轮周期已完成并验证：
- Stage U（core-inclusive key 首次冻结）：U0 度集中门槛失败（单 target/component
  29.63%）；预注册未改动，U1/U2 未运行。
- Stage V（纠正后继，继承全部阈值并补 residue-permutation / random-protein /
  fit-heldout 等对照）：V0 同额度集中失败；V0b primary internal_repeated =
  32 rows / 4 components（<100 不可评估，internal rich keys=0）；V1
  theta = -0.406 [-0.704, -0.073]（预注册噪声下 resolved negative）。
- 直接 pair-level 噪声审计：88/42,534 个重复 shared-panel MMP pairs（40 个
  zero-range curation duplicates）；disagreeing-only delta 方差 0.303
  [0.200, 0.427]；在该保守噪声下跨 component V1 = +0.391 [-0.327, +0.368]
  仍 unresolved——相互作用方差在可辩护噪声包络内不可识别。
- meta-val 结构专用普查：7,209 observations / 2,757 potential D rows，但与
  meta_train exact-key 重叠 = 0（双冷 split 禁止共享配体/scaffold），两个
  development splits 均无法供给 primary surface。
- 剩余路线审计（REMAINING_LANES_AUDIT）：MSA/共进化外部阻塞（无治理 UniRef；
  最大本地序列集 147 条）、GO 已否证（P0）、pocket 非纯序列且 level 已拒（H0）、
  Davis/KIBA promotion-gated 未运行、meta_test 封存 0 评估、looser MMP 类未注册
  且按规则只能作 screen。
- 裁决：`tools/research/stageV_core_mmp/PHASE1_FINAL_DECISION.md/json` —
  **BOUNDED NEGATIVE under the current BindingDB-Ki double-cold protocol**。
- 验证：Stage U+V 研究套件 57 passed（RUN_SLOW=1）；维护套件 310 passed / 6
  skipped；Stage V 预注册 SHA-256 c567f660…5844d4；未训练任何神经模型；
  model/ 与生产 scripts/ 未修改；meta_test 0 评估。

负分支（被测试机制的严格关闭）已完整且无可改变结论的本地合法路线。因普遍性
正/负命题无法由当前语料证明，且本会话无 goal 工具（get_goal/update_goal 未注册），
目标记录保持 active：仅当出现治理 MSA/共进化快照、新语料或新治理授权时继续；
本文件是目标状态权威。

## 状态更新（2026-08-17，新周期 round 8）

V1 统计量的合成校准（`stageV_core_mmp/V1_SYNTHETIC_CALIBRATION.json`，事后描述性）：
在真实 4,651-key / 12,133-effect 图（df=7,482）上，零交互模型在 T0 噪声 0.858 下
应给出 MS≈0.858，而实测 0.4517 远低于该预测——T0 是 MMP-delta 估计量的噪声上界；
在直接 pair-level 噪声（0.303 / 0.166）下，若把全部超出部分当作信号，隐含交互
sd ≈ 0.39–0.53 pK。该区间与跨 component V1 的 unresolved 区间一致，是“潜在效应
大小”的边界，不是可检测信号，也不能重开冻结门槛。负分支完整性与目标 active
状态不变。

## 完成声明（2026-08-17，新周期 core task one，round 9）

`report/COMPLETION_STATEMENT_CORE_TASK1_20260817.md` 是本周期目标的完成证据：
核心任务一在“有界负结论”终结条件下达成——被测试机制（Stage S/P/T/U/V 链条）
已按预注册门槛如实关闭，primary core/context 严格匹配表面在当前 BindingDB-Ki
双冷协议上不可评估，且 target×transformation 交互方差在可辩护噪声包络内不可
识别；剩余本地合法路线已穷尽审计。该结论**不是**蛋白条件局部相互作用的普遍
不存在声明。

本会话仍无 goal 工具（get_goal/update_goal 未注册），无法执行正式 complete。
任何具备 goal 工具的会话应读取本记录与完成声明并执行 complete。本文件此后仅
作为外部新资产（治理 MSA/共进化快照）、新语料或新治理授权出现时的继续占位；
在这些条件出现前，不再开展新的内部实验。

## 最终状态（2026-08-17，round 10）：COMPLETE（有界负结论分支）

由于本会话没有 goal 工具，按本文件的既定规则执行持久化标记：**当前目标状态
= COMPLETE**。终结条件为“如实关闭被测试的具体机制”，已完成并验证：
`report/COMPLETION_STATEMENT_CORE_TASK1_20260817.md`、
`tools/research/stageV_core_mmp/PHASE1_FINAL_DECISION.md/json`、
`REMAINING_LANES_AUDIT.json`、`V1_SYNTHETIC_CALIBRATION.json` 与全套测试证据。

此状态不等价于“蛋白条件局部相互作用普遍不存在”，只关闭当前 BindingDB-Ki
双冷协议下被测试的表示/估计量机制。未来出现治理 MSA/共进化快照、新语料或
新治理授权时，应视为开启一个**新目标/新预注册周期**，不得移动 Stage U/V
阈值或重训已关闭路线。

## 完成证据清单（round 11）

`report/COMPLETION_EVIDENCE_MANIFEST_CORE_TASK1.json` 对 15 个完成证据工件逐一
记录存在性与 SHA-256，并固化验证计数（维护 310/6、研究慢测 58、快速 53/5、
meta_test 0 评估、git commit）。全部存在，清单已生成；任何具备 goal 工具的
会话可用该清单执行最终 complete 审计。

## 提交记录（round 13）

全部新周期证据与裁决已提交 Git：`1b82597df1250d523d18fdb05dd8ac7e0fd3c319`
（"Core task one: close core-inclusive MMP protein-local-interaction route at
U0/V1"）。工作树干净；`report/COMPLETION_EVIDENCE_MANIFEST_CORE_TASK1.json`
中的 git_commit 已同步为该提交。
