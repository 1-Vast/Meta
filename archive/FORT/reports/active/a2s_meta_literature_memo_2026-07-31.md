# A2S-DTA 元学习与负迁移文献备忘录

日期：2026-07-31  
范围：abundant source target -> scarce recipient target；source/recipient target ID 严格不重叠；recipient query 标签只能用于最终评估。

## 结论先行

元学习适合 A2S-DTA，但不能把“target-as-task + MAML + task weighting + negative-transfer gate”整体宣称为新颖。以下已经有直接先例：

- MAML 已定义了跨任务学习一个可快速适配的初始化（Finn, Abbeel & Levine, ICML/PMLR 70, 1126-1135, 2017）。
- ANIL 显示 MAML 的大部分收益可来自共享特征复用，只在 task head 内循环更新（Raghu et al., ICLR 2020；arXiv:1909.09157）。
- FOMAML/Reptile 给出了只用一阶导的低成本替代（Nichol, Achiam & Schulman, arXiv:1803.02999, 2018）。
- CML 已把 MAML、drug-anchored/target-anchored task 与 contrastive block 用于 DTA（Li et al., IEEE BIBM 2022, DOI: 10.1109/BIBM55620.2022.9995372）。
- MetaDTA 已用 attention-based neural process 做 DTA few-shot（Lee et al., *MetaDTA: meta-learning-based drug-target binding affinity prediction*, ICLR 2022 Machine Learning for Drug Discovery workshop；该条目为 workshop，不应写成正式期刊证据）。
- ZeroBind 已将每个 protein 作为 task，并用 task-adaptive self-attention 对 task loss 加权（Wang et al., Nature Communications 14, 7861, 2023, DOI: 10.1038/s41467-023-43597-1）。
- AdaMBind 已在 DTA 上使用 MAML、query-loss/Support-query gradient similarity 调度和 label noise（Wan et al., Nature Communications 17, 2026, DOI: 10.1038/s41467-026-70554-5）。
- ATS 已把神经 scheduler 用于选择 meta-training tasks，并以 unseen-task 泛化优化 scheduler（Yao et al., NeurIPS 34, 7497-7509, 2021）。
- 更关键的是，Mera, Vogt & Bajorath 已在 *Scientific Reports* 15, 35236 (2025), DOI: 10.1038/s41598-025-22058-3 中提出面向药物设计的“meta-learning + transfer learning”负迁移控制：meta-model 根据 recipient 训练损失学习 source 实例权重，并优化迁移初始化。故“用元学习控制负迁移”本身不是空白。
- Tripuraneni, Jordan & Jin 的 NeurIPS 2020 理论（arXiv:2006.11650）表明共享表示的样本效率依赖训练任务多样性；这支持做 source-topology/diversity 审计，但不是 source routing 算法。
- PCGrad（Yu et al., NeurIPS 2020, arXiv:2001.06782）和 CAGrad（Liu et al., NeurIPS 2021, arXiv:2110.14048）处理的是共享多任务优化中的梯度冲突。它们可作为 pooled/multi-task 优化控制，但没有 recipient-conditioned donor 价值估计，也不能单独支撑 A2S 的 source selection 新颖性。

因此 A2S 的可辩护贡献应收窄为：在严格 target-disjoint、endpoint/provenance-closed 的长尾 DTA 评估中，学习 recipient-conditioned **source target** 选择，并以 cross-fitted recipient-level transfer gain 训练风险/弃权策略；不是单纯提出 MAML、ATS 或负迁移门控。

## 逐项可复用性与边界

### MAML / FOMAML / ANIL

可直接复用：target-as-task、support/query episodes、inner-loop adaptation、first-order 低成本实现；ANIL 适合作为 k=1/3/5 的默认结构，只更新 affinity residual head 或低秩 adapter，冻结 Transformer/Mamba body。

不能宣称：MAML 首次解决 scarce-target transfer，或“只更新 head”是 A2S 新颖点。必须和 pooled/no-adaptation、recipient-only calibration、FOMAML、ANIL、full MAML 对照。

### CML / MetaDTA / ZeroBind

CML 证明 DTA 可以构成 target-anchored 与 drug-anchored meta-task，并使用 task-inequality/contrastive block；它可作为 DTA-specific meta baseline。MetaDTA 是 few-shot attention-neural-process baseline，ZeroBind 的 task adaptive attention 是 task loss aggregation，不等于按 recipient 选择 donor target。

不能把“蛋白作为 task”“task weighting”“contrastive meta-learning”单独包装成 A2S 创新。A2S 需要显式输出 source target 的 recipient-conditioned score，并且 score 必须只依赖 recipient support、蛋白/化学表示和 source metadata；不能读取 query 标签。

### AdaMBind / ATS

AdaMBind 的 target-specific MAML、5/40-shot protocol、query-loss + support/query gradient-similarity scheduler 可作为强 baseline/骨架。ATS (Yao et al.) 的 scheduler 选择的是 meta-training task，不是推理时给 recipient 选择 donor source；二者在 A2S 中必须拆成 meta-task scheduler 与 source router。

AdaMBind 的公开结果使用 random task split 与 CD-HIT 40% novel-task split，主要报告 MSE/CI/R²/Spearman/Pearson，未定义 natural target-frequency tail、source/recipient disjoint transfer gain、provenance closure 或 negative-transfer rate。不能直接把其结果外推到 A2S。

### Mera et al. 2025：必须作为近邻先例

该方法在 kinase inhibitor 数据上，用 source 实例权重 meta-model：base model 在 source 上以预测权重训练，recipient training loss 作为外层 validation objective，再微调 recipient；作者在 source/target compound 不重叠时仍报告性能提升，并定义 negative-transfer index。它最适合作为 A2S 的“meta-weighted source transfer”强基线或消融。

但其任务是二分类 PKI、source 实例权重，不是 pKi/pKd 的 target-level donor routing；实验还包含 kinase-family 与 source/target 化学重叠场景。因此不能声称 A2S 首次将元学习用于药物设计负迁移，也不能直接复用其 binary/AUC 结论。

### Task diversity / gradient surgery 的定位

Task-diversity 理论可用于解释为何 source pool 不能只按数量扩张，应报告 source family、chemical-space 和 provenance 的覆盖。PCGrad/CAGrad 可在共享 encoder 的 pooled baseline 中做一个固定优化消融；若它们改善 pooled model，仍不能解释为学习了“哪个 source target 帮助哪个 recipient”。

## A2S 最小有力实验

### 1. 先做低成本闭合基线

固定当前 pKi primary roster（source `n_eff>=100`，recipient `n_eff<30`），k=`{1,3,5}`，只用 TRAIN。每个 recipient 运行：

1. B0 ligand-only / recipient calibration；
2. pooled source model + recipient fine-tune；
3. all-source average；
4. random-source；
5. protein-similarity-only；
6. chemistry-similarity-only；
7. support-compatible scalar router（现有 control）；
8. MAML/FOMAML/ANIL（source-only pseudo-recipient episodes）；
9. AdaMBind-style scheduler（仅训练 episode 采样）；
10. Mera-style meta-weighted source transfer（适配为 regression，作为近邻强 baseline）；
11. A2S router + risk/abstention gate。

### 2. 唯一允许学习的路由标签

对每个 abundant pseudo-recipient h，先留出 h，再从 `H\\{h}` 训练 source model。router feature 只能来自 protein relation、chemical compatibility、source depth/provenance、h 的 support residual/LOO utility；router label 才可用 h 的 query transfer gain。natural recipient 不得参与 router fitting 或阈值选择。

### 3. 必须报告的估计量

对每个 recipient r 和预算 k 定义：

`G_r(k) = loss(no-transfer recipient calibration) - loss(transfer)`。

报告 target-macro RMSE/MAE、within-target Spearman/pairwise accuracy、median G、target-bootstrap CI、negative-transfer rate `P(G_r<0)`、abstention rate 与 risk-coverage；分别给 pseudo-recipient 与 natural recipient，且 pKi/pKd 不合并。

### 4. 预先停止条件

- ANIL/FOMAML 与 full MAML 无差异而显著省算力：只保留 ANIL/FOMAML；
- Mera-style meta-weighted source transfer 不优于 pooled fine-tune：不再增加复杂 instance weighting；
- learned router 不优于 random/all-source 或只在 pseudo-recipient 有效：停止 A2S 神经路由；
- gate 只降低 coverage、不降低负迁移率：不能称为有效 abstention；
- natural pKi k=3/5 的 target-macro gain CI 不排除 0：不启动 Mamba/Transformer 扩容；
- 所有结果依赖 source/query 化学重叠、单一 homology bin 或 contamination：降级为诊断，不宣称迁移。

## 可引用链接

- Finn et al. MAML: https://proceedings.mlr.press/v70/finn17a.html
- Raghu et al. ANIL: https://arxiv.org/abs/1909.09157
- Nichol et al. FOMAML/Reptile: https://arxiv.org/abs/1803.02999
- Yao et al. ATS: https://papers.nips.cc/paper_files/paper/2021/hash/3dc4876f3f08201c7c76cb71fa1da439-Abstract.html
- CML: https://doi.org/10.1109/BIBM55620.2022.9995372
- ZeroBind: https://doi.org/10.1038/s41467-023-43597-1
- AdaMBind: https://doi.org/10.1038/s41467-026-70554-5
- Mera et al. negative-transfer meta-learning: https://doi.org/10.1038/s41598-025-22058-3
- Tripuraneni et al. task diversity theory: https://arxiv.org/abs/2006.11650
- PCGrad: https://arxiv.org/abs/2001.06782
- CAGrad: https://arxiv.org/abs/2110.14048
