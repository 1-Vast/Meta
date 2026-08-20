# Core Task 1 — 数据锚定完成导向报告 / Data-Anchored Completion Report (v3, global-history synthesis)
### 证明"模型真的利用了蛋白—配体条件交互" (sole purpose)

- Date: 2026-08-20 (v3：全局本地历史综合 + 新一轮文献；supersedes prior drafts of this filename)
- **唯一目的 / Sole purpose**：完成**核心任务一**——证明一个**可部署**模型**真的利用了**蛋白—配体**条件交互**（within-target, ligand-dependent ordering）做冷靶点 DTA，而非靠 level/target-key/family key/配体全局模式/attention 幻觉。
- **最终核心任务（北极星）**：把可信交互 + **训练机制创新** → 冷靶点零样本(k=0)/少样本(k∈{1,2,3,5}) DTA 显著增益。
- **本版相对 v2 的新增**：(i) 真正**全局**扫描本地历史（R0–R14 + Stage A–X + CIIP-1A/2 整条弧线）；(ii) 关键定量发现 **R_g=1081** 与 **effective rank 18/1700** 提为一阶证据；(iii) 新文献锚点 **TAPB（target prior bias 干预去偏，Nat Commun 2025）** 把本地"target key"失败与文献正式对接；(iv) 方案围绕"证明交互的三个必要条件"重构，并显式映射任务契约要求的**双创新**（训练机制 + 交互表征）。
- Governance：Core Task 1 = **UNRESOLVED**（report/CORE_TASK1_UNRESOLVED_TERMINAL_20260817.md）。本报告是**完成程序**，不是推翻。UNRESOLVED / FALSIFIED-AS-TESTED 均为合法终局。
- **本版实际审阅**：history.md 全弧线、task.md、report/EVIDENCE_LEDGER.md、report/LITERATURE_R14/R15、dataset/processed/meta_fewshot/* 全 manifest、dataset/processed/open_structures/pilot20k_*、model/（encoders/interaction_grammar/qpsmp_meta）、tools/research/stage{A_innerloop,B,C,D,F,G,H_pocket,I,J,K,L,M,N,P,Q,R,S,T,U,V,W0P,X,CIIP,CIIP2,CIIP_context}。

---

## 0. 一页结论 / Executive verdict（全局综合）

| 维度 | 结论（含本地一阶证据） |
|---|---|
| 核心任务一 | 证明可部署模型真的用了蛋白—配体**条件交互（within-target ordering）** |
| **核心障碍（三重，全部本地量化）** | ① **level/ordering 不对称**：蛋白通路携带 level（target-key），ordering≈0（排序对照 −0.0002）；② **目标竞争 R_g=‖g_abs‖/‖g_ctr‖=1081**（level 梯度是交互的 1081 倍）；③ **表征退化**：KLIFS one-hot effective rank **18/1700**、10/13 对结构性零对比、蛋白通路对残基槽置换**精确不变**（非 pocket-aware） |
| 文献对接 | 本地"target key" = 文献 **target prior bias**；TAPB（Nat Commun 2025）已用**干预去偏**对症，但只在 DTI 分类、未做冷靶点 DTA 的交互归因 → **可迁移创新空间** |
| 唯一已解析可迁移信号 | **配体侧成对 SAR**（Δ-r +0.270 held-out，正交 Tanimoto，cliff 最强）——蛋白无关、成对形式 |
| 唯一通过的蛋白表征 | Q1 **pair_centered_local_esm +0.189 [0.033,0.363]**——但依赖 oracle 坐标（非可部署） |
| 两条未关闭通道（history 明确） | (a) **pairwise learned operator over (query,support)** 用 Stage L 方向（须胜 Tanimoto transport）；(b) **表征**（R2：表征是决定性问题） |
| 功率天花板 | CIIP-2 仪器：当前尺度连种植交互场都恢复不出（>+0.03 R² 不可达）；meta_test 仅 10 components |
| 推荐主线 | **A（干预去偏堵 level）+ B（pocket 落地的成对配体门控）** 为核心；**D（梯度/目标重平衡，训练机制创新）**；C 可部署表征、E few-shot、F 可识别参照下限 |
| 完成标准 | CT1-1∧CT1-2∧CT1-3∧CT1-4∧CT1-5 全过；功率不足记 UNRESOLVED-by-power（合法终局） |

---

## 1. 核心任务一定义：五柱判据（可证伪化）

全程预注册 + SHA-256；keyed rng；family-cluster bootstrap 2000 次；leave-one-family-out(LOFO) 符号稳定。

| 柱 | 名称 | 内容 |
|---|---|---|
| CT1-1 | Deployability | 仅序列/合法先验 + 配体 + 少量 support；**无突变坐标、推理期无复合物结构、无 MSA** |
| CT1-2 | Cold-target transfer | 未见蛋白/家族/化学系列/独立数据集；k=0 与 k∈{1,2,3,5} 稳定优于强基线 |
| CT1-3 | Interaction load-bearing | 消融交互通路 → 冷靶点按预注册幅度掉点（bootstrap lo2.5>0） |
| CT1-4 | Interaction-specificity | 反事实方向 + 负控摧毁；非 family/assay/**target-key** 捷径 |
| CT1-5 | Metric superiority | MSE/RMSE、CI/Spearman 稳定增益，**不可归因于重拟合 level** |

> **本地教训直接改写判据**：CT1-4 必须以 **ordering（排序）形式**为主判据——LEDGER 证明 level/uncentered 对照会**虚假显得蛋白特异**（R3R4 level +0.4216 但 ordering −0.0002）。任何"蛋白特异性"若只在 level 形式成立，**一律不计入**。

---

## 2. 本地数据资产清单（实际读到的 manifest）

### 2.1 主动冷靶点语料 processed/meta_fewshot/bindingdb_ki_main_v0
| 字段 | 值 |
|---|---|
| pair observations | **18,331**（cells 17,717；exact Ki rows 21,695） |
| proteins | **500** 序列（CD-HIT40 → **296 clusters**）；median len 356 / max 1450 |
| ligands | **9,880**（max_atoms 128，max MW 1000） |
| endpoint | exact positive uncensored Ki；**pKi = 9 − log10(Ki[nM])** |
| cross-panel pairs | 618（panel_id = BindingDB assay proxy） |

### 2.2 双冷切分 bindingdb_ki_double_cold_v1（label_blind；轴 ligand=Bemis-Murcko scaffold / protein=CD-HIT40）
| split | cells | ligands | targets | components | scaffold clusters | 备注 |
|---|---|---|---|---|---|---|
| meta_train | 5,643 | 3,825 | 346 | 258 | 1,628 | median 9 lig/target |
| meta_val | 1,411 | — | — | 19 | — | 模型选择用 |
| **meta_test** | **768** | 518 | **22** | **10** | 223 | median 22 lig/target；**密封** |
- **closure 0 泄漏**（component/document/exact_ligand/scaffold overlap 全 0）；meta_test 相似度 lt40=456/t40_60=57/t60_80=5。
- ⚠ **功率约束**：meta_test 仅 **10 components** → family-cluster bootstrap 有效簇≈10 → 小交互效应（~+0.01 R²）在 sealed test 几乎不可检出 → 模型选择落 meta_val，test 只读一次。

### 2.3 冻结特征/编码银行
| 资产 | 内容 |
|---|---|
| ..._protein_bank | **ESM-2 t30 150M**，hidden 640，**128 residue slots**，499 蛋白，1022-chunk 池化 |
| ..._ligand_bank_compact | RDKit 图，atom_feat 32 / bond_feat 12，9,880 配体 |
| ..._tbasis_features.npz | **288 维**，arms={correct, foreign_ligand, deranged_protein} |
| open_structures/pilot20k_* | **14,906 holo 复合物**；split 11,926/1,490/1,490；ESM-2 t30 128-slot 银行 + **128-slot 结构监督标签**（方案 B 的 pocket 监督源） |

### 2.4 生产模型面（model/，研究期**不改**）
interaction_grammar.py（394 行）= 当前交互模型：atom-resolved cross-attention 到 128 residue slots + 全局共享 contact-type 字典 + query 依赖 few-shot transport。Stage 0 审计的旧 BPSF 病灶：k=0 塌缩（spread 0.065 pK）、换蛋白只动 0.009 pK（protein-blind）、k=1 残差核退化为 mean(r)。

### 2.5 当前定量基线（task.md authoritative）
- **T2 leak-free double-cold 基线 MSE = 2.5961 / 1.7712 / 1.3245 / 1.2197 / 0.9859** @ k=0/1/2/3/5。
- **k=0 MSE = level² 1.7314 + centered 0.8648**；oracle target level → ≈0.865 MSE。level 主导但受 level wall（≤25.9% level 方差可蛋白预测）。
- 最强可复现 k≥2 query-specific 比较器 = **固定 Morgan/Tanimoto 残差加权（配体-only）** → 任何蛋白条件方案必须超越的下限。
- 任务契约要求**双创新**：① 一个**训练机制**（效果可单独归因）；② 一个**蛋白-配体交互表征**。目标 MSE≤1.00 across k + 竞争性 RMSE/CI/Spearman。

---

## 3. 全局本地史：已证伪机制 + 两条未关闭通道

### 3.1 全弧线证伪史（R0–R14 + Stage A–X + CIIP）
- **Legacy（pre-R0）**：Analytic QPSMP/LIRMS、HyperSAR、D-MEMT/DORM、CIPF/TERM、K3/ELMT 全证伪。反复病灶：k=1 结构性死亡、support 路由均匀/标签不敏感、**蛋白特异性弱**、residual/level 语义含糊、primitive slots 不可识别。
- **R0–R4**：R0 否决配体检索先验；R1 立双冷切分+meta_test 密封；**R2 确立"表征是决定性问题"**；R3/R4 立 level-shape 族（A0/B3/C2 Pareto）。
- **R5–R8 transport/shape**：query-specific gates 部署惰性或伤校准；R8 shape-first 成对/cliff 训练改善 within-target shape 但与 CI 权衡、未晋级。
- **R9–R14 objective 族关闭**：R13 Tanimoto 残差 k=5 RMSE 0.9988（vs 目标 1.0000，但门要全 k）；R14 关闭单支撑解析 transport。
- **Stage A/B（meta-learning）**：inner/outer-loop  screened、complementary 被拒；发现 ligand 表征在 target 内塌缩、评估泄漏。
- **Stage P**：objective-only protein conditioning **CLOSED**。
- **Stage H**：pocket 先验**未过可识别性门**（H0_POCKET_IDENTIFIABILITY）。
- **Stage S（FiLM 势场）**：B_protein 1/6、C_protein_cf 0/6 门，protein=**target key**，FAIL。
- **Stage T→V（MMP×protein）**：estimand 不可识别（40.4% fit 行 core 集不相交），STOPPED before training。
- **W0-P（点突变 bilinear）**：FAIL（correct sign 0.240 < global-pooled 0.760）。
- **Q2d（低秩 bilinear 合成族）**：CLOSED（SNR 1 下优化/估计失败）。
- **CIIP-1A/2（mutation 级）**：oracle 局部 ESM 非可部署；A1 bilinear R²0.105<先验 0.536、A2 router R²−0.327；仪器证明当前尺度恢复不出种植场（>+0.03 R² 不可达）。

### 3.2 一阶机制证据（collapse audit，2026-08-19）
| 量 | 值 | 含义 |
|---|---|---|
| **R_g = ‖g_abs‖/‖g_ctr‖** | **1081** | level 目标梯度是交互（对比）目标的 **1081 倍** → 朴素联合训练永远先学 level |
| C_g（梯度相关） | −0.016 | 两目标近正交，level 不捎带交互 |
| KLIFS one-hot effective rank | **18/1700** | 蛋白表征退化、近零有效秩 |
| 结构性零对比对 | **10/13** | WT/variant one-hot 相同 → 反对称对比恒 0，优化器无解 |
| 蛋白通路残基槽置换 | **精确不变** | 非 pocket-aware，只是更花哨的 target key |

### 3.3 两条未关闭通道（history 明确授权方向）
1. **Pairwise learned operator over (query,support)**，用 Stage L 的配体侧方向——moment 形式抓不住的成对信号；**须胜 Tanimoto transport 才有意义**。
2. **表征**（R2 决定性；collapse audit 指表征为主导瓶颈）——需 pocket/结构落地才可迁移。

---

## 4. 核心障碍与"证明交互"的三个必要条件（全局综合）

**障碍**：核心任务一要的蛋白条件 within-target ordering 在本地 ≈0，且被 level/target-key（R_g=1081）掩盖；唯一 resolved 可迁移信号是配体侧成对 SAR（蛋白无关）。文献把这同一现象命名为 **target prior bias**（TAPB）。

**因此要"证明交互被利用"，方案必须同时满足三条件：**
- **条件①（堵 level）**：结构性排除 level/target-key——within-target 中心化 / 成对差分 / **干预去偏**，让蛋白通路构造上无法成为 target key。（对症 R_g=1081、Stage A/S）
- **条件②（可迁移落地）**：交互必须 grounding 在**可迁移的蛋白特征**（结合位点化学/几何），而非 target 身份——否则跨 family 不迁移（Stage A out-of-component +0.0065）。需要 pocket/结构 + 破置换不变。（对症 LEDGER 3、rank 18/1700）
- **条件③（成对原生）**：交互读出为**成对**形式，对接 resolved 的成对 SAR 信号（Δ-r+0.270），而非 mean/moment。（对症 LEDGER 2b）

功率诚实：先算 power table（基于真实 component 数）；小效应在 sealed test 不可达则记 UNRESOLVED-by-power，CIIP-2 successor 数据规模为解锁条件。

---

## 5. 冷靶点 DTA 文献"怎么做交互"（新一轮；含证明范式与去偏）

调研 2023–2026 冷靶点/零样本 DTA + 去偏文献。**结论：多数"声称"用交互，但证明范式有系统性缺口；而去偏文献（TAPB/DebiasedDTA）只治偏差、未做冷靶点 DTA 的交互归因——两者交汇即创新空间。**

| 论文 [ref] | 刊物/年 | 交互机制 | "证明交互"方式 | 缺口 / 可借鉴 |
|---|---|---|---|---|
| ZeroBind [1] | Nat Commun 2023 | 子图+图注意力池化 | attention 可视化+冷 benchmark | attention≠load-bearing |
| DCI-SiteDTA [2] | BMC Bioinform 2026 | 结合位点检测+位点感知双交叉互作 | 交叉注意力+位点命中 | 位点检测≠因果 |
| EBD-DTI [3] | bioRxiv 2026 | episodic bridge diffusion | 零样本 episodic benchmark | bridge 未因果隔离 |
| CompBind [4] | JCIM 2025 | 复合物引导预训练（推理 structure-free） | 预训练复合物监督 | **可部署范式（grounding B/C）** |
| LABind [16] | Nat Commun 2025 | 配体感知结合位点识别（学配体×蛋白交互） | 位点识别+交互图 | pocket-grounding 借鉴（grounding B） |
| LaPro-DTA [17] | arXiv 2026 | 潜在双视角药物表征+显著蛋白特征 | 冷切分 benchmark | 显著性≠因果 |
| CS-DTA [12] | Front Chem 2026 | LM 驱动，**比较单向注意力变体 Lp→L/L→P** | **交互设计消融** | 少数做交互消融者；但仍非 load-bearing CI |
| **TAPB [18]** | **Nat Commun 2025** | **干预去偏治 target prior bias（do-算子/后门调整）** | 因果干预+去偏 benchmark | **直接对症本地 target-key；但仅 DTI 分类、未做冷靶点 DTA 交互归因（grounding A）** |
| DebiasedDTA [19] | 2021 (arXiv 2107.05556) | 去偏提升泛化 | 偏差建模+泛化 benchmark | 去偏思路（grounding A/F） |
| Task-conditioned DTI [10] | NeurIPS 2022 | 任务条件化 meta-learning | few-shot benchmark | 条件化未消融 |
| AdaMBind [8] | Nat Commun 2026 | hypernetwork/adapter + 任务自适应 + 课程 | few-shot benchmark（R15 已审） | 训练机制借鉴（grounding D/E） |

**系统性缺口**：冷切分 load-bearing 消融（带 family-cluster CI）+ 反事实方向 + nuisance 不变 + **ordering 形式判据**——文献无一同时做到；本地数据证明这正是区分真交互与 target-key 的唯一试金石。

---

## 6. 六条可行方案（每条：机制 / 三重 grounding（本地+生物文献+跨域）/ 满足哪条必要条件 / falsification）

> 已证伪的朴素 bilinear / naive router / FiLM 势场**不再作候选**（§3），仅作 CIVS 负控下限。下列方案围绕"三必要条件"重构，并覆盖任务契约的**双创新**（训练机制 + 交互表征）。

### 方案 A — TP-ID：Target-Prior Interventional Debiasing（干预去偏堵 level）★条件①
- **机制**：建因果图 Target→(level 捷径)→Affinity、Protein×Ligand→(交互)→Affinity；用**后门调整 / do(Target)** 阻断 level 捷径，强迫预测走交互路径；配 within-target 成对训练（level 自动抵消）。
- **grounding**：本地 R_g=1081 + Stage A/S target-key；生物 **TAPB [18]**（target prior bias 干预去偏）+ DebiasedDTA [19]；跨域 Pearl 因果推断 / 后门调整。
- **创新点**：**首次把干预去偏用于"暴露并归因蛋白-配体交互"**（而非仅提精度），直接把 target-key 从竞争者变成被控制的混杂。
- **Falsification**：去偏后交互通路在冷靶点不优于 ligand-only（lo2.5≤0）→ FALSIFIED。

### 方案 B — PG-PLG：Pocket-Grounded Pairwise Ligand Gating（pocket 落地的成对配体门控）★条件②③，主线
- **机制**：(i) 蛋白编码经 **pocket 感知预训练**（用 pilot20k 14,906 holo + 128-slot 结构监督，**破置换不变**；推理 structure-free，CompBind/LABind 范式）；(ii) 以 resolved 的**配体侧成对 SAR**（Δ-r+0.270）为底座，蛋白只学**门控**——within-target 中心化（无法携带 level），决定调用哪条配体-pair 方向。
- **grounding**：本地 LEDGER 2b 成对 SAR + LEDGER 3 置换不变 + Stage A/S target-key；生物 LABind [16]/CompBind [4]；跨域 FiLM/gating（CV）[20]。
- **创新点**：交互 = **蛋白条件的路由**，作用于**可迁移的配体-pair SAR 原语**——把"生成信号"降为"选择信号"，低维、抗 level 捷径、可消融。
- **Falsification**：冷靶点上蛋白门控的成对 SAR 未超 ligand-only 成对 SAR（如 Stage A protein-side flat）→ FALSIFIED；或 pocket 注意力不优于 slot-置换负控 → FALSIFIED。

### 方案 C — DPC：Deployable Pair-Centered distillation（无坐标恢复 +0.189）★条件②
- **机制**：以 Q1 唯一通过的 **pair_centered_local_esm +0.189** 为**教师**，把 oracle 局部信号**蒸馏**进只看序列、在 128 residue slots 上的**可部署 router**；训练目标=匹配 pair_centered 残差，且显式破置换不变。
- **grounding**：本地 Q1（LEDGER 217/463）+ 修 CIIP-2 A2 naive router R²−0.327；跨域 knowledge distillation（CV/NLP）。
- **Falsification**：可部署 router 不超 shuffled-residue/置换负控 → FALSIFIED（=target key）。

### 方案 D — GORT：Gradient/Objective-Rebalanced Training（梯度重平衡，**训练机制创新**）★条件①，训练创新
- **机制**：显式重平衡 level 与交互目标梯度——PCGrad 投影冲突梯度 / GradNorm 自适应权重 / 双时间尺度，把 R_g=1081 压到可学习区间；或 level 由独立（闭合式/support 提供）项吸收，交互成为**唯一可学残差路径**。
- **grounding**：本地 collapse audit **R_g=1081 + C_g=−0.016**（目标竞争已证）；跨域多任务学习（PCGrad NeurIPS 2020、GradNorm ICML 2018、uncertainty weighting）。
- **创新点**：**直接对症已量化的梯度竞争**——这是任务契约要求的"效果可单独归因的训练机制"。
- **Falsification**：重平衡后冷靶点 ordering 归因仍不显著 / 交互通路仍被 level 吸收 → FALSIFIED。

### 方案 E — CINP：Support-Conditioned Interaction Neural Process（few-shot 桥）★条件①③
- **机制**：k≥1 support 提供 level（绕开 k=0 level 墙），模型只学 support 内排序；用**蛋白条件 support 加权**超越最强 ligand-only Tanimoto 残差加权。
- **grounding**：本地 QPSMP meta 契约 + Tanimoto 残差加权是最强 k≥2 下限；生物 AdaMBind [8]/Task-conditioned DTI [10]；跨域 CNP/MAML [21]。
- **Falsification**：support 操纵不改变冷靶点预测，或不超 Tanimoto 残差加权 → FALSIFIED。

### 方案 F — CF-Ref：Cross-Field 因子分解参照下限（可识别性地板）★可识别性
- **机制**：field-aware 双线性交互项 + 交叉拟合残差（DML）剥离主效应，作**精确可识别的参照下限**（所有方案的强制对照）。
- **grounding**：生物 Cross-Field Fusion [5]；跨域 FM/FFM [22]、DML/R-learner [23]；**注**：CIIP-2 A1 已在 mutation 面板证伪朴素 bilinear，故此仅作受控参照，非候选。
- **Falsification**：（作为下限，无独立 falsification；其存在即定义"交互增益须显著>0"）。

> **主线 = A（堵 level）+ B（pocket 落地 + 成对门控）** 为架构核心；**D（梯度重平衡）为训练机制创新层**叠加；C 攻可部署表征、E 服务 k≥1、F 为可识别性地板。
> **双创新映射**：交互表征创新 = 方案 B（PG-PLG）/C（DPC）；训练机制创新 = 方案 D（GORT）/A（干预去偏）。
> 全部 ≤2M 参数、冻结 ESM/ligand 特征、分钟级/seed、keyed rng；复用 interaction_grammar.py 的 128-slot cross-attention 作可部署蛋白接口。

---

## 7. CIVS — 反事实交互验证套件（7 项；负控含本地已证伪机制）

架构无关；热图/位点重叠仅描述性（Jain & Wallace [24]）。**本地负控**：朴素 bilinear(A1)、naive router(A2)、FiLM 势场(Stage S)、slot-置换(LEDGER 3) 必须复现为"不优于配体-only/被摧毁"。

1. **Cold-target transfer**（CT1-2）：双冷 split + DTI-DG [11]；vs {Tanimoto、ligand-only、additive、frozen-ESM}；k∈{0,1,2,3,5}；**ordering 与 level 分开报告**。
2. **Load-bearing ablation**（CT1-3）：置零交互通路；冷靶点 gap 的 family-cluster bootstrap lo2.5>0。
3. **Counterfactual ligand-swap**（CT1-4）：固定冷蛋白换配体；within-protein 排序朝标签方向。
4. **Counterfactual protein-edit**（CT1-4）：固定配体编辑蛋白；within-ligand 排序朝标签方向。
5. **Nuisance-invariance**（CT1-4）：对 assay/family/batch + **残基槽置换** + 随机上下文不变；**level 形式不得显得特异**（LEDGER 1）。
6. **Interaction-scrambling negative**（CT1-4）：训练期打乱交互通路须塌缩到基线。
7. **Identifiability gate**（CT1-4/5）：family-cluster bootstrap + LOFO；不被 pocket/family-key 探针吸收；**附 power table**。

---

## 8. 功率与运行可行性（具体）

- **CIIP-2 仪器天花板**：49-pair 尺度连种植交互场都恢复不出（>+0.03 R² 不可达，标准 0.25）→ mutation 级交互当前功率不足，R3 UNRESOLVED-by-power。
- **双冷 meta_test 功率**：仅 10 components → family-cluster bootstrap 有效簇≈10 → 小交互效应（~+0.01 R²）几乎不可检出。**主读出放 meta_val（19 comp）+ 内部受控面板；sealed test 只终验一次。**
- **结构可用性（方案 B）**：需 CT1-0 审计双冷 meta_test 靶是否有可用结构/AlphaFold；若无，B 降级为"预训练-only 先验"。
- **解锁路径（CIIP-2 successor，冻结可复用）**：① 同终点面板汇到 ≥100 独立 mutation 条件 / ≥30 parents（单剂量 % inhibition），或 Platinum 级 Ki/Kd ΔΔG；② replicate/噪声表征；③ OLR-Potential/SPB 已冻结，先在更大面板重新检定仪器。
- **诚实终局**：合法数据规模+可部署约束下，A–F 的 CIVS-2 全失败 → 记 FALSIFIED-AS-TESTED / UNRESOLVED-by-power（合法终局），**不**外推为"生物上无蛋白条件化信号"。

---

## 9. 分级、晋级门与脚本锚定

- **CT1-0** 治理+数据审计：用 scripts/build_double_cold_split.py / build_governed_split_views.py 核对 §2.2；family 键=CD-HIT40；**power table + MDE**；结构可用性审计（方案 B）。
- **CT1-1** 方案 A（TP-ID，最小干预去偏）+ 端到端验证 CIVS；**CT1-2** D-on-A（GORT 梯度重平衡）。
- **CT1-3** 方案 B（PG-PLG pocket 门控）；**CT1-4** 方案 C（DPC 蒸馏）；**CT1-5** 方案 E（few-shot）。
- **复用而非重造**：interaction_grammar.py 128-slot cross-attention（可部署蛋白接口）；冻结 OLR-Potential/SPB/tbasis_features（含 foreign_ligand/deranged_protein 对照臂）作 CIVS 负控与仪器；evaluate_qpsmp.py / build_double_cold_split.py 作契约；pilot20k 结构监督作方案 B 预训练源。
- **晋级门**：单 seed CIVS 1–7 全过 → ≥3 seeds 复现（bootstrap lo2.5>0 + LOFO 稳定）。不指标购物；UNRESOLVED 合法。

---

## 10. 成败判据与交付物

- **成功（核心任务一完成）**：晋级方案过 CIVS 1–7——冷靶点 **ordering** 增益（bootstrap+LOFO 稳定，**且 level 形式不单独成立**）、显著 load-bearing 消融 gap、反事实方向过阈值、nuisance/置换不变、指标优势不归因于重拟合 level。**逐柱报告 SUPPORTED/NOT SUPPORTED/UNRESOLVED。**
- **失败**：任一柱失败 = FALSIFIED-AS-TESTED；全部 CIVS-2 失败 = 当前约束下证伪/功率不足（合法终局）。
- **交付物/级**：PREREGISTRATION.md + SHA-256；阈值 addendum；data-audit JSON/MD（含 power table）；RESULT.json；REPORT.md；commands.jsonl；结构+数据契约测试；SHA256SUMS；FAILURES.md；append-only 同步 history.md/task.md/EVIDENCE_LEDGER.md。
- **生产约束**：研究期不改 model/、生产 scripts/；无 oracle 坐标入部署路径；CIIP potential 不入 BindingDB；% inhibition 不改标 Ki/Kd/pK；attention 热图永不作交互证据。

---

## 11. 参考文献（按主题；新增本轮检索）

### 11.1 冷靶点 / 零样本 DTA + 去偏（生物学）
[1] ZeroBind, Nat Commun 2023. https://www.nature.com/articles/s41467-023-43597-1
[2] DCI-SiteDTA, BMC Bioinformatics 2026. https://link.springer.com/article/10.1186/s12859-026-06446-8
[3] EBD-DTI, bioRxiv 2026. https://www.biorxiv.org/content/10.64898/2026.07.14.738384v1.full
[4] CompBind, J Chem Inf Model 2025. https://pubs.acs.org/doi/full/10.1021/acs.jcim.5c02451 — structure-free 推理（grounding B/C）
[5] Cross-Field Fusion DTI, arXiv 2024. https://arxiv.org/abs/2405.14545 — 字段交互（grounding F）
[6] Cross-scale MI + HGCL DTA, IEEE 2024. https://ieeexplore.ieee.org/document/11557120
[7] CoAff-DTI, J Biomed Inform 2026. https://www.sciencedirect.com/science/article/abs/pii/S1532046426001003
[8] AdaMBind (meta-learning task-adaptive DTA), Nat Commun 2026. https://www.nature.com/articles/s41467-026-70554-5 — R15 已审；训练机制借鉴（grounding D/E）
[9] Meta-learning inductive matrix completion (kinase). https://ouci.dntb.gov.ua/en/works/4bggWQNl/
[10] Task-conditioned DTI meta-learning, NeurIPS 2022. https://dev.neurips.cc/virtual/2022/57454
[11] TDC / DTI-DG benchmark, arXiv 2102.09548. https://arxiv.org/abs/2102.09548
[12] CS-DTA, Front Chem 2026. https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2026.1834317/full — 单向注意力变体消融实践
[13] GNN+Transformer few-shot nuclear receptor, J Cheminform 2024. https://link.springer.com/article/10.1186/s13321-024-00902-4
[14] Test-Time Adaptation w/o Source for OOD Bioactivity, ICLR 2026. https://proceedings.iclr.cc/paper_files/paper/2026/hash/0886e50806c3faf55e557bd63ba3e70c-Abstract-Conference.html
[15] KG-MACNF, PLOS ONE 2025. https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0331037
[16] LABind (ligand-aware binding site via protein-ligand interaction), Nat Commun 2025. https://www.nature.com/articles/s41467-025-62899-0 — pocket-grounding（grounding B）
[17] LaPro-DTA (latent dual-view drug + salient protein feature), arXiv 2026. https://ar5iv.labs.arxiv.org/html/2603.14792
[18] **TAPB: interventional debiasing for target prior bias in DTI, Nat Commun 2025.** https://www.nature.com/articles/s41467-025-66915-1 — **target-key 对症（grounding A）**
[19] DebiasedDTA: framework for generalizability of DTA, 2021. http://xxx.itp.ac.cn/abs/2107.05556v5 — 去偏（grounding A/F）

### 11.2 解释忠实性 / 捷径学习
[20] FiLM, AAAI 2018；Shortcut Learning in DNN, Geirhos et al., Nat Mach Intell 2020.

### 11.3 跨域机制（CV/NLP/RecSys/meta/因果/多任务）
[21] MAML/ANIL/CNP（grounding E）。
[22] Factorization Machines (Rendle ICDM 2010)/FFM (Juan RecSys 2016)（grounding F，本地已证伪朴素版，仅参照）。
[23] DML (Chernozhukov 2018)/R-learner (Nie & Wager 2021)/CFR (Shalit 2017)（剥离主效应）。
[24] Attention is not Explanation, Jain & Wallace, NAACL 2019. https://aclanthology.org/N19-1357/
[25] Concept Bottleneck Models, Koh et al., ICML 2020.
[26] Counterfactual augmentation, Kaushik et al., ACL 2020（grounding CIVS/A）。
[27] Antibody DomainBed, ICML 2023. https://icml.cc/virtual/2023/26300（不变表征，grounding A/F）。
[28] ActFound, Nat Mach Intell. DOI 10.1038/s42256-022-00581-6（成对差分，grounding B/E）。
[29] Marginal epistasis, Crawford et al., PLoS Genet 2017.
[30] PCGrad (Yu et al., NeurIPS 2020)/GradNorm (Chen et al., ICML 2018)/uncertainty weighting (Kendall, CVPR 2018)（**梯度/目标重平衡，grounding D**）。
[31] Interactive Graph IB / Iterative Substructure, ICLR 2025；Causal-spurious decoupling OOD, Expert Syst Appl。
[32] HyperNetworks (Ha, ICLR 2017)；Size-Guided Conditional MoE, IEEE. https://ieeexplore.ieee.org/document/11155130

---

## 12. 仓库内部冻结证据（本次实际审阅；grounding，非外部文献）
- report/CORE_TASK1_UNRESOLVED_TERMINAL_20260817.md — UNRESOLVED 终局 + 授权重启路径。
- report/EVIDENCE_LEDGER.md — level/ordering 不对称（L85-100,217,463）；配体侧成对 SAR（Δ-r+0.270）。
- report/CURRENT_MODEL_EVIDENCE.md / report/BOUNDARY_20260817_NIGHT.md — level wall。
- report/LITERATURE_R14_20260816.md / LITERATURE_R15_20260819.md — 既有文献审（AdaMBind 等，避免重复）。
- task.md — 当前契约；T2 基线；k=0 分解；双创新要求；10 篇 reopen 前必读。
- history.md — 全弧线（R0–R14、Stage A–X、CIIP-1A/2、collapse audit R_g=1081/rank 18/1700）。
- tools/research/stageS_sar_field / stageV_core_mmp / stageW0P_positive_control / stageX_csc_signal / stageCIIP2_olr_potential_20260820 / stageCIIP_context_propagation_20260820 / stageH_pocket — 各证伪报告。
- dataset/processed/meta_fewshot/*/manifest.json、dataset/processed/open_structures/pilot20k_* — 数据清单。
- model/interaction_grammar.py / qpsmp_meta.py / encoders.py — 可部署蛋白接口。

---

## 13. 下一步（具体）
1. 冻结 CT1-0 预注册 + SHA-256；用真实 component 数算 **power table / MDE**；审计 meta_test 靶结构可用性（方案 B 可行性）。
2. 跑 **CT1-1（方案 A：TP-ID 干预去偏）**——最小可部署、直接对症 R_g=1081 的 target-key；端到端验证 CIVS 管线。
3. 叠加 **方案 D（GORT 梯度重平衡）**，再推进 **方案 B（PG-PLG pocket 门控）** 作交互表征创新。

---

*本报告 v3 = 核心任务一的全局数据锚定完成程序。相对 v2：把 level/ordering 不对称升级为三重量化障碍（level 不对称 + R_g=1081 目标竞争 + rank 18/1700 表征退化）；用 TAPB/DebiasedDTA 把本地 target-key 对接文献"target prior bias"；方案围绕"证明交互的三必要条件"重构并显式映射任务契约的双创新。功率约束诚实标注，UNRESOLVED-by-power 为合法终局。最终核心任务（可信交互 × 训练机制创新 → 冷靶点零/少样本性能）为北极星。*