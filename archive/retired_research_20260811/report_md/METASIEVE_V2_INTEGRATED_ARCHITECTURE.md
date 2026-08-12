# MetaSieve v2：保留已验证核心、仅替换失败模块的整体架构

日期：2026-08-11  
状态：`RESEARCH_ARCHITECTURE_APPROVED_IN_STAGES`  
生产迁移：`NOT_AUTHORIZED`

## 1. 结论

三份建议不能原样拼接，但可以收敛成一条可执行路线：

```text
保留：target-as-task + support-only adaptation + closed-form positive ridge
      + d<=5 linear arm + controls + frozen law operator

已测试并拒绝：uncentered section
              -> explicit support intercept + centered section
原因：centered 虽优于 intercept，但跨三个数据集均劣于 uncentered

条件替换：失败的 128-slot / 288D T-BASIS biological readout
          -> 通过结构和 measured-crossed Gates 的 exact-residue pair field

条件扩展：linear U kernel
          -> rich PSD pair kernel（只有 linear arm 已通过 biology Gate 后）
```

不把 PAREF、富核、partner anchor 和 CSMO 一次性加入，而是按
`K0 -> K1(已失败) -> R0 -> R1 -> K2 -> K3/K4 -> R2 -> z bridge` 顺序逐层解锁。
这样能够回答每个增益究竟来自校准、生物表示、核容量、跨靶点交互监督，
还是 support adaptation。

本轮已安全实现研究态的 `CenteredKernelSection`。它既可接原来的
`m=U^T phi`，也可在后续接高维 pair feature；未修改 `model/`、生产
CSMO 或训练脚本。

## 2. 当前模型真正面临的问题

### 2.1 Few-shot 增益大，但绝大部分是 target calibration

main-v0 的 Full MSE 为 `1.916`，显著优于 `d=0` 的 `8.711`，说明
support-only closed-form adaptation 在真实 Ki 数据上有效。但匹配的
support intercept 在 meta-validation target macro 上为 `1.441`，优于
Full 的 `1.578`；在已消费 meta-test 上，Full 相对 intercept 的
cluster-macro 特异增益仅约占 pair-support 总增益的 `2.1%`。

因此不能把“Full 大幅优于 d=0”全部解释成 ligand-specific meta-learning。
原模型没有先显式拿走 support 的均值校准，低维 section 很容易把其容量
用于 target offset。

### 2.2 绝对排序仍弱

main-v0 的 `R2=-1.244`、Pearson `0.097`、Spearman `0.086`。这说明
模型可以修正 target 的平均尺度，却没有可靠学会同一 target 内 ligand 的
亲和顺序。permuted support 与 correct support 很接近是同一问题的反映。

### 2.3 生物 partner specificity 未建立

33 个旧 test targets 只有 6 个 CD-HIT40 clusters，且一个 cluster 占
21/33。按 cluster sensitivity，correct protein 相对 wrong protein 的
MSE reduction 为 `-0.081`，LCB 为 `-0.227`。V1 的 2x2 干预进一步显示：
只换 support protein 或只换 query protein 会恶化，但 support/query 同时
换成同一个 wrong protein 后性能恢复。

这证明当前模型使用的是“support/query 坐标一致性”，不是必须对应真实
partner 的 affinity family。A0 又排除了“一组统一精确正交变换”作为完整
解释，所以不能把所有 wrong/wrong 恢复简化成单一 gauge 定理。

### 2.4 288D T-BASIS 不足以承担 partner-affinity 主张

在 physically redacted 的 11,278-row 2-core 上，T-BASIS 的观测方差
`98.07%` 可被 protein + ligand 加性项解释，fixed-ligand partner
dispersion 仅 `5.13%`。A1 的 dependency-component-held-out probe 中，
T-BASIS 还劣于 zero、ESM additive 和 rewired controls。

这只否定“当前 288D 表示中存在能被受限 probe 稳定恢复的、达到预注册
强度的 selectivity signal”；它不证明局部前端不存在任何非线性信息。

### 2.5 128 slots 是语义瓶颈，不只是 tensor 形状

现有 `MechanisticInteractionBridge` 的模块形状可以接 variable-length
residue states，但训练 supervision 把多个 residue 合并进 128 slots，
distance 取槽内最小值，contact 近似 OR。把旧 checkpoint 直接用于 exact
residue 输入会产生训练标签与部署语义错配，不能称 exact-residue model。

历史证据同时表明“恢复 exact residue”也不是充分条件：图感知 ligand
representation 确实增加过 residue-side signal，但 foreign ligand 仍几乎不
掉分；6 Å exact edge 的新增信息又几乎被 contact degree 吸收。因此新的
pair field 必须证明 ligand-specific atom×residue coupling，而不是只增加
分辨率或参数量。

### 2.6 V1 的 unrestricted pair prior 已失败

V1 pair-d0 MSE `3.806`，劣于 ligand-d0 `3.084`；V1-A/B 的 absolute MSE
也劣于 v0。所以下一版不能再次加入一个自由 `mu(P,L)` MLP，并指望它自动
成为 biological anchor。population 默认仍保留 `mu_L(L)`；任何 partner
anchor 必须先由真实 measured rectangles 和 absolute wrong/wrong Gate
独立准入。

### 2.7 数据供给能做 development，不能做 fresh confirmation

Metz、PDSP 和 BindingDB panel 有 crossed interaction development supply，
但 component 数、giant-component 占比以及与 main-v0 的 protein/ligand
重叠不满足 fresh confirmation。旧 main-v0 test、meta-validation 和 A1
都已消费，不能继续用于新架构的确认性选择。

### 2.8 冻结理论和上游 meta-learning 被混淆

权威冻结理论的唯一算子是：

\[
\mathsf A(F,z)=K_{\rm band}(B(z)F(z)).
\]

它明确不包含 support intersection、conditional fiber、confidence
coordinate、RKHS radius 或 architecture-specific training theorem；见
`FINAL_THEORY_COMPLETE.md` 第 90--96、470--491 行。因此：

- kernel ridge 的 `lambda` 不是冻结 simplex target 的 `mu`；
- kernel power 不是冻结理论 certificate；
- PAREF/kernel 与 law head 端到端联合训练不在当前 generalization proof 内；
- raw residue×atom tensor 不能进入 `z`。

要保留冻结 law guarantee，必须先在独立 development 数据上训练、冻结并
哈希上游；随后定义有界有限维 `z`，再用声明的 IID task sample 估计 `F`。

## 3. 保留、替换、后置

| 部分 | 决定 | 理由 |
|---|---|---|
| protein-as-task episodes | 保留 | 与研究问题和 AdaMBind 方法级协议一致 |
| frozen ESM residue states | 保留 | 尚未证明 PLM 本身失败 |
| ligand GINE atom states | 保留 | 图信息已证明优于简单 atom mean |
| P1B contact/distance bridge | 保留为 slot incumbent / geometry prior | 有已审计结构证据，但不能冒充 exact-residue teacher |
| 288D T-BASIS | 降为 legacy baseline | affinity/selectivity admission 失败 |
| learned `U`, `d<=5` | 首轮保留 | 已出现真实 support effect，尚未证明固定 d 是瓶颈 |
| positive closed-form ridge | 保留 | support-only、稳定、可微、线性核 primal/dual 严格等价 |
| uncentered section | 保留 predictor | K1 跨数据集均优于 centered；intercept另作控制 |
| free full-rank pair prior | 删除主线，仅保留失败 arm | V1 已实证有害 |
| exact-residue typed pair field | Gate 后替换前端 | 当前仅有候选资产，无可迁移 checkpoint |
| rich FC-MKS | Gate 后扩展 | 是严格 kernel 泛化，但不应和前端同时改变 |
| Q-PMA/MAML/support Transformer | 继续关闭 | 当前失败不是 support attention capacity |
| RFMS reserved block | 继续关闭 | `Xi(c-c_wrong)` 非零保证尚不存在 |
| CSMO/Band | 原封不动后置 | 当前没有 admitted biological `z` |

## 4. 整体模型

### 4.1 上游 pair representation

最终候选接口是全局 pair map，而不是只有同一蛋白内定义的相似度：

\[
\psi_\theta(P,L)\in\mathbb R^H,
\qquad
\kappa_\theta((P,L),(P',L'))
=\langle\psi_\theta(P,L),\psi_\theta(P',L')\rangle.
\]

这样 kernel 天然对称 PSD，而且 protein×support 2x2 干预中 support protein
与 query protein 不同仍有唯一、可复现的 cross-kernel 定义。

候选 exact-residue field 的数据流为：

```text
exact ESM residue states [R,Dp] ----\
                                      chunked atom-residue interaction field
ligand GINE atom states [A,Dl] -----/          |
P1B geometry prior / real typed labels --------|--> bounded pair feature psi(P,L)
```

实现必须 chunk 化，不能持久化大规模 `[B,A,R,H]` tensor。H-bond、salt、
pi、hydrophobic 等可共存，应使用 multi-label gates；没有真实 typed labels
或确定性物理规则的 latent channel 只能叫 interaction channel，不能命名为
具体物理作用，更不能叫 binding energy。

### 4.2 已测试候选：linear-U + explicit intercept + centered section

K1 曾测试以下候选：

\[
m(P,L)=U^T\psi(P,L),\qquad d\le5.
\]

对 support residual `r=y_S-mu_L(L_S)`：

\[
b_t=\frac1k\mathbf1^Tr,\quad
H=I-\frac1k\mathbf1\mathbf1^T,
\]

\[
M_c=HM,\qquad
\alpha=(M_cM_c^T+\lambda I)^{-1}Hr,
\]

\[
\hat y_q=\mu_L(L_q)+b_t+(m_q-\bar m_S)M_c^T\alpha.
\]

显式 intercept 使用一个 support-label 方向，centered section 的秩最多
`k-1`，合计 label-dependent 自由度仍不超过 `k`。常数 residual 只能进入
intercept，不再被误读为 ligand-specific section。该分解保留为诊断工具，
但因 K1 跨数据集预测 Gate 失败，不替换 uncentered predictor。

### 4.3 条件扩展：rich PSD kernel

只有上一步在 fresh biological Gate 上通过，且容量诊断显示 linear-U
不足时，才删除 `U` bottleneck：

\[
K_c=H\,\Psi_S\Psi_S^T\,H,
\quad
k_c(q,S)=(\psi_q-\bar\psi_S)\Psi_{S,c}^T,
\]

\[
\alpha=(K_c+\lambda I)^{-1}Hr,
\quad
\hat y_q=\mu_L(L_q)+b_t+k_c(q,S)\alpha.
\]

现有 ridge 是这个模型在线性 `m=U^T phi` kernel 下的精确特例：

\[
(M^TM+\lambda I)^{-1}M^T
=M^T(MM^T+\lambda I)^{-1}.
\]

因此升级的是共享函数空间，而不是另加一个 adaptation network。部署时
只有 `k` 个 support labels 进入闭式解；encoder、kernel、normalization、
`lambda` 都必须预先冻结且不得读 query labels。

### 4.4 coverage 只作外部 surrogate

对显式 intercept + fixed feature map，可报告：

\[
\rho^2(q;S)=\frac1k+
(\psi_q-\bar\psi_S)^T
(\Psi_{S,c}^T\Psi_{S,c}+\lambda I)^{-1}
(\psi_q-\bar\psi_S).
\]

它是在一个线性-Gaussian 参数模型下、以 observation variance 归一化的
mean-parameter variance；固定 feature map/`lambda` 时随 nested support
非增。可定义 bounded score `c=1/(1+rho^2)`。

它仍只是工程 coverage surrogate：不包含 model misspecification、测量外噪声
或总 predictive uncertainty，也不是冻结 CSMO theorem 的 radius。raw
kernel power 单独报告，因为 support recentering 后它本身不保证随 k 单调。

### 4.5 partner anchor 暂不接入

仅使用 pair kernel 不能数学上阻止 wrong/wrong。反例：

\[
\psi(P,L)=R_P\phi(L),\qquad R_P^TR_P=I
\]

会给所有蛋白同一个 Gram geometry。measured crossed loss可以降低这种风险，
但不是保证。由于 V1 pair prior 与 RFMS guarantee 均失败，partner anchor
不属于首轮模型。它只有在 R1 证明新 `psi` 含增量 measured interaction
signal 后，才能作为低容量独立 arm，并必须同时满足：

1. correct absolute loss 不劣于无 anchor；
2. `MSE_ww-MSE_cc` 的 component-level LCB 大于 0；
3. support solver 不能在已测 episode 上把 anchor 完全湮灭；
4. fresh donor map 未参与 anchor 训练。

### 4.6 下游 law operator

只有 R2 的 point/meta/biology Gates 全通过后，才定义有界有限维 statistic
map `z=g(frozen upstream episode)`, 冻结其版本、范围和缺失/abstention语义。
然后在独立 task sample 上训练 `F` 并使用原算子：

\[
\mathsf A(F,z)=K_{\rm band}(B(z)F(z)).
\]

point predictor 不等于 law operator；上游成功也不自动继承 law theorem。

## 5. 数据与损失的合法分工

### R0：结构 interaction admission

- 可用真实 holo/typed labels训练 exact-residue field；
- P1B slot checkpoint只能作 incumbent/geometry prior；
- InteractBind 的 Vina score 不是实验 Ki，只能用于结构定位开发；
- correct partner 必须优于 deranged protein/ligand、pocket prior 和
  composition baselines；
- inference unit 是 family/dependency component，不是 atom-residue rows。

### R1：measured crossed-affinity admission

只使用四格都实测、相同 endpoint/assay/panel 条件的 rectangle：

\[
\Delta_\times y=y_{11}-y_{12}-y_{21}+y_{22}.
\]

它严格消去 `a(P)+b(L)`，但不识别绝对 offset，必须和 absolute episodic
loss 联合。禁止把未测 wrong protein 当 non-binder。rectangles 高度依赖，
先 panel/group 聚合，再 dependency-component inference。

### R2：fresh unseen-target few-shot confirmation

- one protein = one task；CD-HIT40 cluster split；
- headline `k=5`，支持曲线 `k=1/2/3/5`；五 repeats；
- 同一 manifest/episodes 比较所有 arms；
- target-intercept 是永久 baseline；
- correct / zero / foreign / permuted support；
- `cc/cf/wc/wf` protein×support factorial；
- primary 仍含独立绝对 Gate `MSE_ww-MSE_cc>0`，不能只看交互项。

若 M 是 lower-is-better loss，正确的“correct-protein 下 support 特异增益”
方向定义为：

\[
\Delta_{\rm BioMeta}
=(M_{cf}-M_{cc})-(M_{wf}-M_{wc}).
\]

正值表示 correct support 的价值在 correct protein 下更大。但 wrong/wrong
恢复的例子也可能产生强 2x2 interaction，所以该量只能是机制指标，不能
替代 `M_ww-M_cc` 的绝对 Gate。

## 6. 分阶段消融与停止规则

| 阶段 | 唯一变化 | 解锁条件 |
|---|---|---|
| K0 | primal v0 ridge 改写成 dual linear kernel，零训练 | predictions/gradients float64 `<1e-8` |
| K1 | current `U,d<=5` + explicit intercept + centered section | **已执行并失败：不迁移** |
| R0 | exact-residue pair field 替换 T-BASIS | independent structural partner/localization Gate |
| R1 | frozen pair field 预测 measured `Delta_x` | component-held-out 增量优于 additive/T-BASIS/null |
| K2 | admitted pair field + original linear-U centered section | support、biology、absolute loss Gates |
| K3 | K2 仅放宽为 rich PSD kernel | K3 相对 K2 有配对增量 |
| K4 | K3 + measured crossed objective；可选低容量 anchor arm | K4 相对 K3 增量且 ww absolute Gate 通过 |
| R2 | 冻结模型在 fresh protein clusters 确认 | component LCB、coverage、noise Gates 全通过 |
| z bridge | bounded statistic admission | 独立数据上确认后才接 CSMO |

任一阶段失败，停止在该阶段；不得同时放宽 kernel、加入 anchor、改变前端
或打开 CSMO 来“救结果”。

## 7. 本轮代码落地与退役

K1 centered-section implementation and tests were research-only branches. They
were retired after the preregistered cross-dataset Gate failed; the
preregistration, JSON result, and written report remain as auditable evidence.

其已验证的实现内容包括：

- explicit support intercept；
- centered positive-ridge dual section；
- `k<=5` enforcement；
- rank、condition number、effective df；
- exact support-label sensitivity 与 bounded-noise radius；
- dimensionless information-variance surrogate；
- stable nonnegative raw kernel-power diagnostic；
- serialized ridge configuration、state/config consistency 和 dtype-aware
  ridge fail-closed；
- linear primal/dual equivalence witness；
- measured 2x2 `Delta_x` 与正确符号的 biology-support synergy helper。

验证覆盖：positive ridge、support budget、constant-residual separation、support
permutation invariance、primal/dual prediction及gradient equivalence、rank-
deficient/k<d/k>d、support-noise bound、orthogonal wrong/wrong counterexample、
kernel scale counterexample、float32 near-rank-deficient stability、nested
information variance、PSD/rank、centered end-to-end gradients、独立 support/query
protein 四臂、additive rectangle cancellation和 2x2 sign。

当前这些测试验证的是代数与接口，不是 biological PASS，也不授权真实新模型
训练或生产迁移。

## 8. 文献定位

- [ALPaCA](https://arxiv.org/abs/1807.08912) 已有 learned feature basis +
  closed-form Bayesian linear adaptation，说明这两个构件本身不是首创。
- [META-KEL](https://proceedings.mlr.press/v162/kassraie22a.html) 与
  [F-PACOH](https://arxiv.org/abs/2106.03195) 提供 meta-kernel/GP 近邻，
  但其 confidence theorem 不能转移给 MetaSieve。
- [PSICHIC](https://www.nature.com/articles/s42256-024-00847-1)、
  [Interformer](https://www.nature.com/articles/s41467-024-54440-6)、
  [PIGNet2](https://pubs.rsc.org/en/content/articlelanding/2024/dd/d3dd00149k)
  支持 interaction-aware/physicochemical pair representation 的研究价值；
  它们不是 MetaSieve 的 exact-residue few-shot protocol。
- [PBCNet](https://www.nature.com/articles/s43588-023-00529-9) 使用同一
  pocket 的 ligand-pair affinity difference，但不是 unseen-protein task。
- [InteractBind](https://arxiv.org/abs/2605.24045) 提供六类结构 interaction
  maps，但其 affinity 是 AutoDock Vina docking score，不是实验 Ki。
- [CNP](https://proceedings.mlr.press/v80/garnelo18a.html) 是从 context
  observations 推断函数的概念先例；没有 MetaSieve 的 support budget 或
  law-valued operator。

较稳妥的新颖性主张不是“发明 exact interaction field、双差或 kernel”，而是：

> 在 unseen-target few-shot DTA 中，以 support-label 信息预算约束闭式
> adaptation；用真实 measured protein×ligand rectangles 准入不可由主效应
> 解释的 biological interaction；再用 fresh protein×support 反事实证明该
> interaction 没有退化成 wrong/wrong 自洽坐标系，最后才桥接冻结的
> law-valued operator。

这是待验证的整合创新假设，不是已建立的首创结论。

## 9. 当前终局状态

```text
V0_META_ADAPTATION_RETAINED
UNCENTERED_CALIBRATION_MODULE_REPLACED_IN_RESEARCH
TBASIS_RETAINED_AS_LEGACY_BASELINE_ONLY
PAREF_REQUIRES_R0_AND_R1_ADMISSION
RICH_FCMKS_LOCKED_BEHIND_LINEAR_ARM
PARTNER_ANCHOR_NOT_YET_AUTHORIZED
CSMO_UNCHANGED_AND_DISCONNECTED
PRODUCTION_MIGRATION_NOT_AUTHORIZED
```

## 10. K1 跨数据集执行结果

本报告提出的首个 centered predictor replacement 已按预注册合同在公共
BindingDB、Davis、KIBA 上执行。centered arm 相对纯 intercept 的 MSE
reduction 在三个数据集均为正，但相对原 uncentered ridge 在三个数据集均
为负：`-0.0703/-0.0882/-0.0713`。equal-dataset pooled standardized LCB
为 `-0.1063`。

因此 centered section 不迁移。原 uncentered positive ridge 继续作为已保留
solver，support intercept 仅作为永久 calibration baseline。该结果把下一处
需要替换的失败模块进一步锁定为 biological pair representation，而不是
closed-form solver。完整结果见
`report/meta_fewshot/K1_CROSS_DATASET_CENTERED_SECTION_REPORT.md`。
