# CIIP-Potential Harness Execution Prompt

你是 `D:\MetaSieve` 项目的独立研究、审计和实现 agent。你的目标不是继续堆叠 Gate，也不是立即重构生产模型，而是验证一个能够直接进入 cold-target DTA 的统一函数：

```text
s_theta(P, L)                         # protein-conditioned SAR potential
g_theta(P, A -> B) = s_theta(P,B) - s_theta(P,A)
D_hat = [s(Pa,B)-s(Pa,A)] - [s(Pb,B)-s(Pb,A)]
y_hat_k = b_L(Lq) + b_P(P*) + s(P*,Lq) + b_S + delta_s(P*,Lq,S)
```

核心问题是：CIIP 识别的 differential signal 是否能成为最终 DTA 的 single-ligand score，并在真实 zero-shot/few-shot DTA 中产生可归因收益。不要把 functional inhibition 改名为 Ki/Kd/pK/DTA；不要把 scientific positive control 直接当作 BindingDB 性能结果。

## 必读资料与当前事实

先完整阅读：

- `task.md`
- `history.md`
- `report/EVIDENCE_LEDGER.md`
- `report/CONSOLIDATED_RESEARCH_RECORD.md`
- `report/research_ideas/ciip/CIIP_RESEARCH_PROPOSAL.md`
- `report/metasieve_research_programme.md`
- `tools/research/stageU_mmp_interaction/REPORT.md`
- `tools/research/stageV_core_mmp/REPORT.md`
- `tools/research/stageV_core_mmp/PHASE1_FINAL_DECISION.md`
- `tools/research/stageX_csc_signal/stageX0c_measurement_qualification_20260818/REPORT.md`
- 当前 Q2d 相关目录、报告和 `GOAL_ACTIVE.md`

当前必须保留的事实：

1. Stage U/V 已经基本耗尽 BindingDB exact-rectangle surface，不能换网络后重复同义实验。
2. Stage V 的结论是 exact interaction 在当前 surface 上 not estimable，不是 biological absence。
3. Q1 局部 ESM 表示具备突变敏感性；此前 Q2 one-hot harness 失败不能直接归因于真实生物学缺失。
4. ligand-side Tanimoto transport 是当前最强可复现 comparator，不能删除或用 protein 分支替代。
5. Q2d 当前链必须先完成终端整理；不得为了逃避 FAIL 继续无限创建 synthetic successor。
6. CIIP 不是最终 DTA 解。必须验证 `signal -> zero-shot utility -> few-shot utility`，才允许进入生产模型。

## 执行原则

- 先核对当前工作树、运行中的进程、Git 状态和已有结果；不得覆盖其他 agent 未提交文件。
- 不修改冻结预注册；若实现定义冲突，创建明确命名的 successor，并说明旧结果不受影响。
- 研究代码只能放在 `tools/research/` 对应新目录；`model/` 和生产 `scripts/` 在通过 utility 之前禁止修改。
- 不使用 ridge、闭式求解、test-time gradient、测试标签、query-panel transductive mean 或 free target ID。
- 所有训练必须是普通前向/反向联合训练；不得把研究 probe 伪装成最终模型性能。
- 每次实验只改变一个主要因素；先 CPU/synthetic smoke，再 GPU；失败立即停止该分支并写出原因。
- 所有命令写入 `commands.jsonl`，所有结果保存机器可读 JSON 和 Markdown 报告。
- 每个正结果必须同时运行 ligand-only、global protein、family-preserving shuffle、random protein、ligand-invariant shift 和 no-interaction 控制。
- 不得通过损坏 wrong arm 来制造 correct arm 的优势。

## 阶段 0：Q2d 终端归档

先确认当前 Q2d-1e 是否仍在运行。若运行，等待其自然结束并生成最终报告；若已结束，核对 gate、诊断、复现和 SHA。将其最终结论同步到 `history.md`、`task.md` 和 `report/EVIDENCE_LEDGER.md`，但不要修改旧 preregistration。

本阶段只做治理和归档，不重新训练新的 synthetic interaction family。

## 阶段 1：KiRHub 兼容性 census

审计 2026 年 KiRHub 数据是否适合作为 CIIP-1A/1B positive-control。主要来源：

`https://www.nature.com/articles/s41587-026-03090-8`

必须先得到：

- usable WT-variant pair 数；
- 每个 pair 的 identical ligand 数；
- parent、mutation、fusion 数量；
- construct、substrate、cofactor、ATP 条件完整率；
- duplicate 和饱和比例；
- functional inhibition 的端点方向、单位和可比较性；
- parent/pocket group connectivity；
- complete held-out parent 是否可构造；
- 每个 parent 的 ligand coverage 和有效 centered effect 方差。

任何字段缺失都要量化，不得用论文摘要中的总样本数替代 usable sample 数。若原始数据不能合法获取，停止模型训练，记录阻塞原因和可替代数据源。

## 阶段 2：统一 potential probe

先实现一个最小、可审计的 probe，不能直接上复杂 Transformer。提供两个诊断臂：

```text
P_pair: free pairwise predictor g(P,A,B)
P_potential: scalar potential s(P,L), with g=s(P,B)-s(P,A)
```

`P_potential` 必须通过以下结构测试：

- identity effect exactly zero；
- A->B 与 B->A 严格反号；
- A->B->C->A cycle sum 接近零；
- query 单点可以直接输出 `s(P,L)`；
- 不能读取 test labels 或 target ID；
- protein shuffle 真正改变 protein input，不改变 matched rows。

低容量 probe 可采用 `s=alpha(P)^T psi(L)`，但只作为诊断。不要把它自动升级成最终主干，因为 Stage S 已经显示 global protein compression 可能成为 target/family key。

只有在 positive control 显示低容量 potential 有真实信息后，才允许测试局部 interaction potential：

```text
residue-level protein states
atom/functional-group ligand states
local masked relevance
interaction tokens
scalar pooling -> s(P,L)
```

没有 common-frame complex 时，不得称为 contact probability、Cartesian interaction 或 atomic 3D recognition；只能称为 learned local residue relevance。

## 阶段 3：CIIP-1A 与 CIIP-1B

### CIIP-1A

在同一 parent 的 WT/variant panel 内，比较：

```text
ligand-only
global protein
correct local protein
family-preserving shuffle
ligand-invariant mutation shift
```

中心化目标为 mutation effect 减去该 mutation 的 ligand-wide mean。CIIP-1A 只证明 representation capacity，不得称 cross-parent transfer。

### CIIP-1B

完整留出 parent 或 pocket group。测试 correct local protein 是否在 unseen parent 上超过 chemical-only、global protein、family shuffle 和 ligand-invariant baseline。若只有 1A 通过而 1B 失败，结论必须是“有局部表示能力但没有迁移规则”，禁止进入 MetaSieve production。

## 阶段 4：instrument qualification

使用 observed graph 注入三类完全 synthetic labels：

1. protein main effect：`D = u_p + noise`；
2. family shortcut：`D = u_family(p,tau) + noise`；
3. true local interaction：`D = z_p^T M z_tau + noise`。

interaction probe 只能在第 3 类恢复，且必须在第 1/2 类拒绝。报告 Spearman、Pearson、dead-zone sign、scale recovery、false-positive rate 和 effect-size power curve。synthetic PASS 只证明 instrument 可用，不证明真实生物学。

## 阶段 5：Potential Bridge 到 native DTA

CIIP-1A/1B 未通过时禁止进入本阶段。通过后冻结 potential 参数初始化，但不要永久冻结，也不要让 absolute loss 单独重写它。

在 BindingDB-Ki native split 上比较：

```text
B0: ligand-only baseline
B1: ligand + global protein
B2: ligand + verified local potential
```

主要报告 within-target centered MSE/RMSE、CI、Spearman；absolute MSE 必须报告但不能单独归因 interaction。构造 oracle interaction ceiling：

- oracle 有效、learned 无效：representation/learner 不足；
- oracle 无效：interaction 不是当前 DTA 性能瓶颈；
- oracle 与 learned 都有效：才授权接入生产模型。

## 阶段 6：Few-shot utility

使用 nested support `k=1,2,3,5`，同一 query panel 比较：

```text
Z0: zero-shot
Z1: level-only support calibration
Z2: level + fixed Morgan/Tanimoto
Z3: level + potential-conditioned transport
```

核心不是 `MSE(k) < MSE(0)`，而是 Z3 是否超过 Z1 和 Z2，同时通过 wrong-target、label-shuffled、structure-only controls。

support 残差定义：

```text
r_i = y_i - y_hat_0(P*, L_i)
b_S = shrink_k(mean(r_i))
r_tilde_i = r_i - mean(r_i)
delta_s_q = eta_k * sum_i K(q_q, q_i) * r_tilde_i
```

必须满足：所有 `r_tilde_i=0` 时 `delta_s_q=0`。因此 `k=1` 初版只允许 level update；`k>=2` 才允许 centered shape update。固定 Tanimoto 是强基线，不得移除。

## 交付与终止

每个阶段产出：

- `PREREGISTRATION.md`（训练/评估前冻结）；
- `RESULT.json`；
- `REPORT.md`；
- `commands.jsonl`；
- 结构测试和数据契约测试；
- 更新 `history.md`、`task.md`、`report/EVIDENCE_LEDGER.md`。

只有以下全部成立，才允许修改 `model/` 和生产 `scripts/`：

```text
CIIP-1A interpretable PASS
CIIP-1B unseen-parent PASS
Potential Bridge zero-shot utility PASS
Few-shot Z3 > Z1 and Z2 for k>=2
no leakage/noise/censoring failure
component-level uncertainty supports the claim
```

若 KiRHub positive control 通过但 BindingDB bridge 失败，结论是“生物机制/表示存在，但 native DTA 数据或 estimand 不足”；若 KiRHub 1B 失败，停止 protein-conditioned interaction route；若 practical Z2 已经优于 Z3，则保留 Tanimoto 作为主性能路线，不强行加入 protein branch。

持续自主执行，但每个阶段必须在前一阶段完成并记录 PASS/FAIL/UNRESOLVED 后才进入下一阶段。目标是得到可重复的科学和性能结论，不是通过复杂实验制造正结果。

## Stage 1 失败后的强制故障归因流程

当前单 seed screening 已完成，`unified_local` 未通过冻结门槛，不能继续 CIIP-1B、BindingDB Bridge 或生产集成。你现在必须先执行只读 collapse audit，不得直接换 ESM、改 loss、增大 rank 或重新训练。

将问题拆成：

```text
输入信息 -> potential 坐标 -> 优化梯度
```

evaluation/统计问题必须单独审计。

### 1. 输出方差与可评价性

对于每个 WT/variant pair，计算：

```text
d_pred(L) = s(P_variant,L) - s(P_WT,L)
c_pred(L) = d_pred(L) - mean_L(d_pred(L))
```

报告 true variance、predicted variance、`sqrt(var_pred/var_true)` scale ratio、OLS slope、centered MSE、centered R2、sign accuracy、`N_target_informative`、`N_prediction_nonconstant`、`N_rank_evaluable`、`N_rank_evaluable/N_total` 和 parent-level 结果。

Spearman 对常数预测是数学上 undefined，不得将 NaN 当作真实零相关。zero-centered ligand-only baseline 使用 centered MSE、R2、sign 和 slope 比较；所有 rank 结果必须带有效 pair 分母。

### 2. Bilinear potential 的解析分解

当前 potential 近似：

```text
s(P,L) = alpha(P)^T psi(L)
Delta_P s_v(L) = [alpha(Pv)-alpha(Pw)]^T psi(L)
Var_L(Delta_P s_v) = Delta_alpha^T Cov(psi) Delta_alpha
```

从现有 checkpoint/feature cache 计算 protein input difference、`||alpha(Pv)-alpha(Pw)||`、ligand embedding covariance/effective rank、mutation difference matrix effective rank、potential variance，以及 ligand-invariant/ligand-dependent variance。

解释规则：input difference 很小表示 KLIFS representation bottleneck；input difference 大但 latent difference 很小表示 encoder/optimization collapse；protein/ligand latent 都有方差但 quadratic form 近零表示 cross-space misalignment 或 rank bottleneck；potential variance 健康但方向错误表示 interaction basis mismatch。

### 3. Gradient competition audit

只针对 interaction parameters `theta_s` 计算：

```text
g_abs = grad_theta_s(L_abs)
g_ctr = grad_theta_s(L_contrast)
R_g   = ||g_abs|| / (||g_ctr|| + eps)
C_g   = cosine(g_abs, g_ctr)
```

分别报告 `b_P`、`b_L`、`s_theta` 的梯度范数。若能读取训练轨迹，检查初始化、10%、25%、50%、75% 和最终阶段。

- `R_g >> 1` 且 `C_g < 0`：支持 absolute objective 干扰 centered learning；
- 两者梯度健康但 potential variance 近零：优先怀疑表示或参数化；
- `g_ctr` 从初始化就近零：优先怀疑输入或梯度路径。

### 4. 表示替换审计

当前 Stage 1 测的是 KLIFS pocket one-hot，不是此前 Q1 通过的 per-position local ESM。先做不训练比较 KLIFS 与 local ESM 的 WT/variant difference、mutation difference effective rank 和 mutation-centered information。只有审计显示 KLIFS 不足且 gradient conflict 不成立时，才允许创建 ESM-only successor；其他变量必须保持不变：ECFP4、rank、potential formula、optimizer、LR、steps、split、seed、metrics 和 gates。

### 5. Successor 决策树

审计完成后只能按以下规则选一个方向：

```text
representation bottleneck -> ESM-only successor
objective conflict         -> centered-only successor
两者都成立                 -> 预注册 2x2：KLIFS/ESM x joint/centered
两者都不成立               -> rank/capacity successor
```

不得同时替换 ESM、loss、rank、ligand features 和 optimizer。任何新实验都必须新建 preregistration，保留当前 Stage 1 原始结果，不得覆盖原 gate。

### 6. 对 free-pairwise 的解释

free-pairwise 只是 conditional expressivity ceiling，不是生产候选。只有当它在多个 parent 上稳定非恒定、centered R2/slope 良好，且 unified potential 在相同输入下仍失败时，才可以说当前 potential 表达不足。不得把 free-pairwise 优势解释为 biological non-integrability。

### 7. 交付和限制

在 read-only audit 完成并形成报告前，禁止 CIIP-1B、BindingDB Potential Bridge、生产代码修改、新大型 interaction architecture 和新 synthetic successor。生成：

- `STAGE1_COLLAPSE_AUDIT.json`；
- `STAGE1_COLLAPSE_AUDIT.md`；
- `commands.jsonl`；
- 新增结构/统计测试；
- `history.md`、`task.md`、`report/EVIDENCE_LEDGER.md` 同步。

最终报告必须明确：

```text
tested one-hot potential: FAIL
biological protein-conditioned signal: UNRESOLVED
primary diagnosed cause: representation / objective / capacity / unresolved
authorized successor: exactly one of ESM-only / centered-only / 2x2 / capacity
CIIP-1B: not authorized until successor is interpretable
```

## 统一修正版：所有阶段围绕同一个可部署函数对象

本节优先级高于本提示词中任何把 CIIP、Potential Bridge、few-shot transport 分成独立模型的旧表述。下一轮研究的核心对象不是三个模块，而是同一个标量蛋白–配体相互作用势：

```text
f_theta(P, L) = b_P(P) + b_L(L) + s_theta(P, L)
```

其中 `b_P` 和 `b_L` 是可解释的蛋白/配体主效应，`s_theta(P,L)` 是唯一允许承担 protein-conditioned interaction 的标量 potential。所有差分必须由这个同一函数导出：

```text
protein contrast:  s(P_variant,L) - s(P_WT,L)
ligand contrast:   s(P,L_B) - s(P,L_A)
CIIP double diff:   [s(Pv,LB)-s(Pv,LA)] - [s(Pw,LB)-s(Pw,LA)]
zero-shot:          b_P(P) + b_L(L) + s(P,L)
few-shot update:    s(P,L_q) + delta_s(P,L_q | support)
```

禁止分别训练一个 CIIP pair predictor、一个 DTA interaction head、一个 few-shot adapter，再用相关性拼接结果。自由 pairwise predictor 只能作为诊断对照，用于衡量“可积 potential 的表达约束是否过强”；它不能被包装成最终统一机制。

### 统一 potential 的最小实现要求

1. 第一版使用低容量可审计模型，优先采用 `s=alpha(P)^T psi(L)` 或局部状态的受限双线性形式；不得一开始使用大型 Transformer、复杂 hypernetwork 或闭式求解器。
2. 若使用局部表示，蛋白侧必须保留突变中心/口袋残基的局部 ESM 或经过核验的 residue state，配体侧使用 atom/functional-group state；无共同坐标时只能说 learned local relevance，不得称接触图或 3D recognition。
3. 明确区分 `b_P/b_L` 与 `s`。差分监督必须使用 cross-fitted nuisance、panel-centered target 或显式主效应头，不能让 potential 通过吸收 ligand main effect 或 target calibration 获得虚假增益。
4. 在同一模型中同时提供 `s(P,L)`、protein contrast、ligand contrast 和 double contrast 的接口，并用 identity、反对称和 cycle consistency 测试证明它们来自同一函数。
5. 低秩 global protein code 只能作为诊断或消融。若它在 unseen parent 上失败，必须转向 residue-local potential，而不是继续扩大 global code 容量。

### 真实数据执行顺序的重新定义

1. 完成并归档 Q2d 终局；不得再增加 synthetic successor。
2. 在 Duong-Ly 内先完成同一平台、同一配体的 CIIP-1A（同 parent）和 CIIP-1B（held-out parent）。使用 centered mutation contrast：

   ```text
   d_vl = y(Pv,L) - y(Pw,L)
   c_vl = d_vl - mean_L(d_vl)
   ```

   报告端点只能叫 `% inhibition`，不能转换为 pK/Ki/Kd。Davis 只能作为同平台 Kd 的独立复制或样本量不足说明；Anastassiadis 不能和 Duong-Ly 混合成同一 estimand。
3. CIIP-1A/1B 必须在同一 `s_theta` 上输出 contrast。若自由 pairwise predictor 胜过 potential，保留该结果作为“可积性约束造成的表达损失”，不要伪称 potential 已经成功。
4. 通过真实 CIIP 后才执行 BindingDB Potential Bridge。必须比较 ligand-only、固定 Morgan/Tanimoto、global protein、local potential；Tanimoto 是已验证的强基线，potential 只能以增量修正形式加入，不能替换它。
5. 通过 Bridge 后才执行 few-shot：

   ```text
   y_hat_q = b_P(P*) + b_L(L_q) + s(P*,L_q)
             + lambda_k * mean(r_i)
             + eta_k * sum_i [K_Tanimoto(q,i) + gamma*K_potential(q,i)] * r_tilde_i
   ```

   `r_i` 必须相对于完整 zero-shot base 计算；`r_tilde_i=r_i-mean(r)`。所有 centered residual 为零时，shape correction 必须严格为零。`k=1` 默认只做 level calibration；`k>=2` 才允许 shape transport。新 potential kernel 是 Tanimoto 的增量项，不得在没有 paired bootstrap 的情况下声称替代 Tanimoto。

### 训练创新的唯一主线

本项目的训练创新应命名为 **Centered Evidence Transport on an Integrable Protein–Ligand Potential**。训练不是闭式解、不是测试时求解器、不是多阶段部署优化，而是普通端到端梯度训练，包含：

- absolute endpoint loss，用于保持可用的 level/calibration；
- within-target ligand-pair ranking 或 centered difference loss，用于学习可迁移 SAR shape；
- exact rectangle 权重高于 contextual pseudo-pairs；
- protein/ligand main-effect nuisance 路由与 potential gradient 分离；
- label-permuted、family-preserving shuffle、ligand-invariant mutation shift 和 no-interaction 负控；
- potential cycle/identity regularization，但不得用无梯度的恒等正则冒充约束。

训练目标不能只优化 CIIP contrast，也不能只优化绝对 MSE。必须报告各损失项的梯度覆盖、`b_P/b_L/s` 的方差分解和 component-level bootstrap。若 shape 提升伴随 calibration 或 CI 系统性恶化，应判定为失败，而不是只报告较好的 Spearman。

### 允许进入生产代码的明确条件

只有当同一个 `s_theta` 同时满足以下条件，才允许把实现并入 `model/` 与生产 `scripts/`：

1. Duong-Ly CIIP-1A 通过，且不是由 ligand main effect、protein main effect、端点方向或删失伪影解释；
2. CIIP-1B 在 held-out parent 上仍保留，且优于 chemical-only 与 family-preserving shuffle；
3. BindingDB bridge 在 zero-shot 或 within-target centered 指标上不劣于 incumbent，并对固定 Tanimoto 有独立增益；
4. few-shot `k>=2` 的 potential-conditioned transport 在 paired target/component bootstrap 下优于 `level-only` 和 `Tanimoto-only`；
5. 所有收益在 ligand novelty、parent novelty、assay/platform 分层中方向一致，且没有 test-label 参与检索、标准化、checkpoint 选择或超参数选择；
6. 自由 pairwise 对照、local potential、global potential 的差异已解释清楚，不能把诊断模型的结果移植成生产模型的机制结论。

任一条件失败时，保留最强的 ligand-only/Tanimoto 实用路线，并明确关闭 protein-conditioned potential，不继续增加适配器、门控或 synthetic 变体。

### 研究效率与设备约束

研究阶段优先保证 estimand、split、数据契约和可复现性，而不是盲目提高 GPU 利用率。Q2d 诊断模型很小、arms 串行、数据在 CPU 生成会导致低利用率，这是实验性质决定的，不得为提高功耗而改变冻结预算。只有在 unified potential 进入真实面板训练后，才可采用 packed batch、预缓存 ligand/protein features、减少 CPU→CUDA 往返和并行 independent seeds；优化必须保持与基线相同的样本、split、步数语义。

### 最终交付

最终报告必须回答四个独立问题：

1. 是否存在可迁移的真实 protein-conditioned interaction signal；
2. 该 signal 是否由统一的 `s_theta(P,L)` 表示，而不是不可积的 pairwise shortcut；
3. 它是否改善 BindingDB 的 zero-shot 或 cold-target few-shot 性能；
4. 它是否在固定 Tanimoto、ligand-only、shuffle、删失和 component bootstrap 对照下仍成立。

报告必须分别标注 `SOLVED`、`PARTIAL`、`UNRESOLVED` 或 `FALSIFIED-AS-TESTED`，禁止把“CIIP 正控通过”直接写成“最终 DTA 性能提升”。
