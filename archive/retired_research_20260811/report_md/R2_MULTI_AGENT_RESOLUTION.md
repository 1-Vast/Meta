# R2 多智能体终局复核：MetaSieve 当前问题与最小修复路线

日期：2026-08-11  
范围：Cowork R2 产物、冻结理论、E0–E3、RFMS 草图、本地 crossed-panel 数据与相邻文献。  
性质：已消费开发/测试集上的纠错性诊断；不打开新科学 Gate，不授权生产迁移。

```text
R2_H0_REGIME_FALSIFIED_BY_ITS_OWN_PREREGISTRATION
META_SECTION_PREDOMINANTLY_CALIBRATION_DESCRIPTIVELY
TBASIS_OBSERVED_DESIGN_NEAR_ADDITIVE_WITH_LOW_PARTNER_DISPERSION
LOCAL_CROSSED_TRAINING_SUPPLY_EXISTS_NO_FRESH_CONFIRMATION
RFMS_MATHEMATICALLY_UNDER_SPECIFIED_AND_NOT_AUTHORIZED
```

## 1. 结论

Cowork 找对了一个重要现象：当前 Meta-Section 的平均误差收益主要来自
support 对新 target 的整体 affinity 偏移校准，而不是稳定的 ligand-specific
排序。它随后给出的两个更强结论不成立：

1. 不能把 wrong/wrong 恢复归因于一个已被识别的 episode-wise `GL(2)`
   gauge。R2 自己预注册的 E1 规定 `gauge_ratio>1` 反证 H0；实际 meta-val
   为 `1.0306`，meta-test 为 `1.0987`，两次都触发反证。
2. `k=5,d=2` 不推出“表示和数据理论上都不可能有用”。该推论只在无噪、
   正确指定、support unisolvent、零正则等额外条件下成立；当前 BindingDB
   与 `ridge=1` 不满足这些条件，而且 Cowork 引用的 F20/CI-A3 来自旧开发
   handoff，并非当前权威冻结 CSMO operator。

所以当前真正的问题是三件事叠加，而不是一个单一 gauge bug：

- section 把 target-level calibration 与 ligand-specific correction 混在一起；
- 冻结 T-BASIS 在当前稀疏观测设计上几乎是加性的，固定 ligand 的 partner
  变化只占总特征离散的约 `5.13%`；
- 现有 wrong/wrong 对照同时改变 support/query 的 protein coordinate，无法
  单独证明完整预测器使用了绝对 partner identity。

最小可行修复不是 RFMS，而是先把校准显式分离：

```text
pair-specific population prior
+ explicit support-mean intercept
+ centered residual / centered-coordinate ridge section
```

该候选已经作为研究 helper 实现并通过数学测试，但 post-hoc 数字来自已消费
split，只能用于选择下一轮候选，不能作为新证据。

## 2. 当前完整预测器究竟是什么

现有 v0 不是 Cowork 报告中简化的纯 kernel。写成双形式：

```text
mu_q = f_L(L_q) + m_q c_pop
rho  = y_S - f_L(L_S) - M c_pop
H    = m_q M^T (M M^T + lambda I)^-1
yhat = mu_q + H rho
```

对固定 `rho`，`lambda=0` 且 `rank(M)=d` 时，`M->MG,m_q->m_qG` 对任意
可逆 `G` 保持 kernel；`lambda>0` 时只对正交 `G` 精确保持。完整预测器还含
固定 checkpoint 的 `c_pop`，wrong protein 替换也不会自动把它按 `G^-1`
同步变换。因此：

> kernel 的坐标不变性是必须审计的边界条件，但不是当前 wrong/wrong 恢复的
> 已证因果解释。

这与 representation identifiability 文献的边界一致：学习表示常只能识别到
某类线性变换，但“存在重参数化不变性”不等于“每次 partner 干预就是同一个
可转移变换”。参见 [Roeder et al., ICML 2021](https://proceedings.mlr.press/v139/roeder21a.html)
和带 task side-information 的非识别性分析
[Deng et al.](https://arxiv.org/abs/2201.07348)。

## 3. E0：校准确实占主导，但不是所有 family 都一致

E0 使用原五个 checkpoint 和原 episode draw。`pair_intercept` 保留 full
模型的正确 `(P,L)` population prior，只把 section 替换成 support residual
均值；因此它是隔离 section 形状信息的匹配 null，不是 ligand-only null。

| split / 统计单位 | full | pair intercept | pair zero | centered section |
|---|---:|---:|---:|---:|
| meta-val target macro | 1.5780 | **1.4408** | 2.7681 | 1.4721 |
| meta-val 9-cluster macro | 1.2550 | 1.1573 | 2.6978 | **1.1388** |
| meta-test target macro | 1.9162 | 1.8965 | 7.9260 | **1.8426** |
| meta-test 6-cluster macro | **1.9155** | 1.9608 | 4.0691 | 1.9145 |

meta-test cluster-macro 中，`pair_intercept-full=0.0453`，只占
`pair_zero-full` 总 support 增益的约 `2.10%`。这支持“平均增益主要是校准”。
但方向并不跨 cluster 一致：intercept 优于 full 的 cluster 数，meta-val 为
`4/9`，meta-test 为 `3/6`。因此禁止写成“所有 biology evidence 都只是
calibration”。另外，独立训练的 ligand-only intercept 仍明显差于 full：
meta-val cluster `1.5331`，meta-test cluster `2.5683`。

### 最小修复候选

令 `r=y_S-mu_S`，显式校准 `b=mean(r)`，再使用：

```text
M_c = M - mean_rows(M)
m_cq = m_q - mean_rows(M)
delta_q = m_cq M_c^T (M_c M_c^T + lambda I)^-1 (r-b)
yhat_q = mu_q + b + delta_q
```

它有三个优点：常数 residual 只能进入 `b`；section 只为 within-target 变化
负责；task-specific 连续状态仍不超过 `rank(M_c)<=k`。这不是声称性能已经
成立，而是把当前混杂变成可证伪的两个通道。代码位于
`research/meta_fewshot/r2_calibration_orthogonal_section.py`。

闭式 solver 进入 meta-learning 的构件有直接先例：
[R2-D2](https://openreview.net/forum?id=HyxnZh0ct7)。本项目新增的研究问题应是
“校准正交化后，生物坐标还能否产生跨 family 的 ligand-specific 增益”，而
不是把 ridge 本身称为创新。

## 4. E1：预注册 H0 被反证

| split | registered calibration share | gauge ratio | effective ridge ratio |
|---|---:|---:|---:|
| meta-val target macro | 0.6937 | **1.0306** | 0.3898 |
| meta-test target macro | 0.7793 | **1.0987** | 0.3428 |

预注册要求 `calibration_share>=0.9` 才支持近均匀权重，`gauge_ratio<=0.5`
才支持 near-gauge；并明确规定 `gauge_ratio>1` 反证 H0。结果既未过两个支持
阈值，又两次触发反证。`effective ridge ratio` 也不是可忽略的小量。

E1 的 support-weight 运算本身不使用 query label，但旧 `load_data()` 会把
整个 corpus 的 `y` 载入内存。因此准确措辞是：

```text
labels loaded by legacy loader; E1 structural statistic does not consume them
```

不是操作意义上的 sealed label-free experiment。

## 5. E2：当前 T-BASIS 近加性，但 Cowork 对前端的描述有误

修正版 E2 先生成物理删标的结构索引，核对 `cell_id` 顺序，并在 bipartite
2-core 上重算，避免大量 ligand singleton 自动制造零残差。

| 量 | 结果 |
|---|---:|
| cells | 17,717（不是旧报告的 21,473） |
| ligand singleton fraction | 0.6483 |
| crossed 2-core rows | 11,278 |
| crossed 2-core components | 47 |
| interaction degrees of freedom | 7,484 |
| additive explained fraction | 0.9807 |
| interaction residual fraction | 0.0193 |
| fixed-ligand partner dispersion fraction | 0.0513 |
| deranged/natural feature shift | 2.878 |

正确结论是：

```text
TBASIS_OBSERVED_CROSSED_DESIGN_NEAR_ADDITIVE
TBASIS_FIXED_LIGAND_PARTNER_DISPERSION_LOW
```

它不能升级为“T-BASIS 没有非线性交互容量”。backfitting 在 100 次迭代时
最终相对变化约 `5e-5`，且这是 observed sparse design 上的线性分解。

Cowork 还错误地写成“protein 只通过 6 类 composition 进入，没有 residue
geometry”。源码实际是 ESM residue states 先条件化 bridge distance logits，
随后才在 8 atom channel × 6 residue chemistry class × 6 radial shell 上聚合。
更合理的怀疑是：**最终聚合压弱了 partner variation**，而不是“前端从未看见
residue 信息”。`deranged/natural=2.878` 反而证明 statistic 会响应 partner
替换；缺的是已验证的 affinity direction，不是任何 partner sensitivity。

## 6. RFMS 为什么现在不能训练

设 adaptable support/query frame 为 `A,a`，reserved frame 为 `C,b`。共享
partner coefficient 改变量 `dc` 对预测的真实暴露为：

```text
Xi = b - a(A^T A + lambda I)^-1 A^T C
prediction gap = Xi dc
```

所以 `Var(c0)>0` 不是 partner specificity 证书。即使 `c0` 非常数，也可能
因 `Xi=0`、partner 差落入 `Xi` 零空间、query reserved energy 为零，或
reserved block 被 adaptable span 补偿而 wrong/wrong 恢复。该反例已经写入
单元测试。

RFMS 还存在四个理论/实验缺口：

1. ridge 解的最大 support residual 只是
   `inf_delta ||y-base-Cc-A delta||_inf` 的一个上界，不是 empty-fiber 证书；
2. hard `V_S/V_C` split 只对保持分块的正交变换自然不变，不是任意换基不变；
3. `d_support<=k<d` 与当前全局 `d<=5,k=5` 不能同时满足；
4. 用随机 wrong protein hinge 训练会把未测 partner 当 non-binder，并直接
   优化最终控制 Gate，不能作为生物证据。

因此 RFMS 目前只能保留为被阻止的研究草图。若以后重开，必须先证明
`Xi(c_correct-c_wrong)` 在实测 crossed panel 上非零，并用真正的
L-infinity feasibility solve，而不是 c0 方差或 ridge residual 替代。

部分适应、task-conditioned bias 和 side-information factorization 均有近邻：
[CAVIA](https://proceedings.mlr.press/v97/zintgraf19a)、
[ANIL](https://openreview.net/forum?id=rkgMkCEtPB)、
[MT-net](https://proceedings.mlr.press/v80/lee18a.html)、
[inductive matrix completion](https://www.microsoft.com/en-us/research/publication/provable-inductive-matrix-completion/)
和 [Macau](https://arxiv.org/abs/1509.04610)。因此 RFMS 若最终成立，创新也只能
放在可识别预算与可证伪 exposure certificate，不能说“reserved block”本身新。

## 7. 本地 crossed data：可开发，不可 fresh confirmation

E3 v2 按真实 schema 读取 design，完全不访问 `Y/pki`，并按 assay block
计算 interaction df；不同 panel 不相加伪造 confirmation PASS。

| panel | measured cells | interaction df | dependency components | largest share | 合法角色 |
|---|---:|---:|---:|---:|---|
| Metz XP2 | 32,849 | 31,775 | 1 | 1.000 | historical development |
| BindingDB panels | 6,363 | 3,188 | 8 | 0.912 | 当前不合格 |
| PDSP core | 10,701 | 8,313 | 1 | 1.000 | historical development |

这些包也不 fresh：BDB 的 129 条 protein 中 117 条与 main-v0 exact 重合，
2,845 个 ligand 中 825 个 SMILES exact 重合；PDSP 的 45 条 sequence 中 30
条 exact 重合；XP2 的 147 条中 25 条 exact 重合。必须先做完整
protein-cluster + scaffold + document closure，再重新数剩余独立 components。

Metz 是经典 dense kinase profiling panel
([Metz et al.](https://www.nature.com/articles/nchembio.530))，适合 source-only
interaction/PCM 开发；它只有一个 dependency closure component，不能充当
独立确认。现有 `esm2_t30_kinase_pocket85.npz` 也只覆盖 82 个 kinase；全 147
覆盖的是 full-sequence bank。虽然 XP2 内已有 147 条长度 85 的 KLIFS pocket
string，但在重新生成并 hash 前不能宣称 pocket bank 全覆盖。

## 8. 可执行的解决顺序

### R2-D0：冻结当前负结果

- 保留原 prereg，添加 post-execution invalidation，不回写假装预注册从未错；
- `r2_regime_audit.json` 标记为旧 v1，不再引用；
- meta-test/main-v0 已消费，不再用于架构选择；
- RFMS、Q-PMA、CSMO integration 和 production migration 全部停止。

### R2-D1：在新的 source-side development components 预注册最小消融

保持相同 encoder、optimizer、episode draws 与参数预算，只比较：

1. ligand population；
2. pair population；
3. pair population + support intercept；
4. 原 uncentered section；
5. explicit intercept + centered section。

primary 必须是 cluster/component-macro 的
`MSE(centered)-MSE(intercept)`，随后才检验 correct-permuted、correct-wrong
和 Pearson/Spearman；intercept 必须成为 standing baseline。不得从已消费
post-hoc 数字设阈值。

### R2-D2：先做 crossed-panel 的简单可证伪基线

在历史 development panel 上按 held-out protein family + held-out ligand
scaffold 比较：additive baseline、低秩 crossed interaction、PCM/IMC side-info。
只有低秩 cross-term 在 family/scaffold 双冷 split 稳定优于 additive，才值得
重建完整 147-pocket bank 或研究更复杂 partner anchor。

蛋白–化合物矩阵配 side information 的低秩建模已有成熟先例
([Macau](https://arxiv.org/abs/1509.04610))；同一化合物跨靶点的 selectivity
建模可参考 PCM 的问题定义
([van Westen et al.](https://pubs.rsc.org/en/content/articlelanding/2014/ib/c4ib00175c))。
这些是必须先过的简单基线，不应被更复杂 meta architecture 跳过。

### R2-D3：fresh confirmation 只能前瞻性建立

新确认集需与 meta-train/meta-val/main-v0 test 以及所有 R2 development panel
同时按 protein family、exact ligand/Murcko scaffold、document/assay closure
隔离。没有至少 30 个合格 targets 和足够独立 components，不声称 partner
specificity confirmation。

## 9. 已完成的代码纠错与验证

- E0/E1：split-specific 输出、真实 `protein_group_40` cluster macro、匹配
  pair/ligand intercept、checkpoint hashes；
- E2：物理删标结构索引、cell alignment、2-core、真实 cluster 字段、收敛信息；
- E3：显式 matrix/long schema、per-assay df、label-free design 读取；
- RFMS：删除 wrong-protein training hinge，加入 quotient exposure，修复单例
  batched protein shape，强制 CLI `d_support<=k`，保持训练硬阻止；
- 新增 8 个测试，覆盖 GL/orthogonal ridge 边界、校准正交化、E2/E3 synthetic
  算术与 RFMS 反例。

最终研究状态：MetaSieve 的少样本校准能力是真实且有用的，但当前没有足够证据
把主要收益解释成 ligand-specific、partner-specific biology。下一步应先用最小
centered section 与 crossed low-rank baselines 重新建立可识别的增量，再讨论
RFMS 或更强 biological frontend。
