# A2S-CMAL 实验方案

日期：2026-08-01  
状态：数据与模型框架已就绪；正式多种子训练未在本机运行。

## 1. 唯一科学问题

本项目只检验一个命题：能否从 abundant source targets 学到一个可迁移的
meta-adaptation operator，使它在严格未见的 recipient target 上，仅凭
`k={1,3,5}` 条 recipient support measurements，就对每个 query compound
产生 target-specific、query-dependent 的排序修正，并且优于不使用 support
的 DTA predictor 以及不学习适配机制的强基线。

最终模型不得退化为标量校准、插值、相似性检索、固定贝叶斯 posterior、
ridge/kernel 的闭式求解。这些方法只作为对照。

## 2. 证据标签：禁止把项目新增设置冒充论文设置

本文所有实验条目都带以下来源标签：

- **PAPER-EXACT**：按原论文 Methods/Figure protocol 复现；不修改数据集、
  split、support size、metrics 或论文 ablation。
- **PROJECT-FROZEN**：来自本项目研究问题、已给出的 Figure 2 意见或既有
  preregistration；不是“照搬论文”的内容。
- **ENGINEERING**：训练实现参数，只能在 source meta-validation 上选择，
  不得写成论文 protocol 或科学创新。
- **NOT RUN**：只完成设计或数据准备，尚无正式结果。

如果 PAPER-EXACT 与本项目严格泄漏规则冲突，分别报告两条 track；不得静默
改写论文 protocol 后仍称为复现。

## 3. 已准备的数据

### 3.1 主实验：ChEMBL 37 pKi（PROJECT-FROZEN）

正式主数据为
`dataset/formal_training/chembl37_pki_formal.v4`：157,613 条 exact pKi
measurements、1,098 个 human single-protein targets、82,646 个 parent
compounds。只纳入 `Ki`, `=`, `nM`, confidence 9, binding assay, human
wild-type single protein；pKi 与 pChEMBL 不一致、censored 和高噪声 assay
context 均隔离，不进入主标签。

冻结的 unseen-target roster 为
`dataset/formal_training/a2s_d0r_roster.v3`：

| 项目 | 数值 |
| --- | ---: |
| abundant source targets | 206 |
| recipient targets | 63 |
| recipient homology components | 55 |
| 每个 recipient 的 support draws | 5 |
| support budgets | 1, 3, 5 |
| source/recipient target、accession、document、parent、assay overlap | 0 |
| support/query parent、document overlap | 0 |
| homology-warm recipients | 9/63 |

因此主结论是严格 **target-ID unseen**，不是把 63 个 recipient 全部描述为
sequence-family unseen。除全体 63 targets 外，必须单独报告 54 个
homology-cold recipients；9 个 homology-warm targets 不得混入“严格
homology-cold”措辞。

roster v3 只从正式 `pki_measurements_context_main.parquet` 读取 metadata，
并为每个 support/query compound 冻结一个 `measurement_uid`；构建过程不读取
pKi。该修订防止 trainer 对同一 target-compound 的多 assay/time records
重新取 median，从而改变已冻结的 estimand。

模型就绪包为 `dataset/formal_training/a2s_cmal_episodes.v3`，content SHA-256
为 `2df5831bc8a51df93dc54531302327716fcca8900ec43f1aa37f16ed2fb9485a`。
它不含 affinity labels，共 30,123 episodes：23,127 meta-train、2,595
meta-validation、3,456 meta-test、945 recipient-test；每个 k 恰有 10,041
episodes。source target split 为 163/19/24，homology components 不跨 split。
support/query parent、document、measurement、ordered-time、nested-support 和三种
counterfactual same-target violations 均为 0。

### 3.2 AdaMBind 论文复现数据（PAPER-EXACT）

相关原始文件已经按 AdaMBind Git commit
`01a169a6d62fba0d6c003f47bfba539e55f5b344` 固定在
`tmp/adambind-data/data/raw`，并完成字节哈希审计：

| 文件 | SHA-256 | 本地审计 |
| --- | --- | --- |
| `bindingdb-full-data.csv` | `3ebd8dfabd2a20c0dbceba35cc59ba8e6dd44a90798667fef2c9059bab63fbba` | 42,203 rows；2 个 non-finite affinity，正式复现前按论文/官方代码的明确策略处理并报告 |
| `davis-full-data.csv` | `dc9331894d5eafa46787632cc0d9754406e5a96eb87980b27d4abe22308a6994e` | 30,056 rows |
| `kiba-full-data.csv` | `7b1e306a2344e38c5d5bbcda6f6112201440bbaa92d5081a4fc054ed83edca24` | 118,254 rows |

论文报告的 benchmark 规模为 Davis 68 drugs/442 targets、KIBA 2,111/229、
BindingDB 9,864/1,088；本地 CSV 的唯一序列/SMILES计数与论文表格不完全相等，
因此复现报告必须同时列出“paper-reported count”和“parsed-file count”，不能
用一个替代另一个。来源：[AdaMBind full text](https://www.nature.com/articles/s41467-026-70554-5)、
[固定代码](https://github.com/Moohyun-w/AdaMBind/tree/01a169a6d62fba0d6c003f47bfba539e55f5b344)、
[source data](https://doi.org/10.6084/m9.figshare.30963823.v1)。

### 3.3 LIT-PCBA（PAPER-EXACT，NOT RUN）

AdaMBind 使用 LIT-PCBA 的 15 targets、7,844 actives、407,381 inactives。
该 track 只有在数据按原论文版本与许可取得后才运行，不使用 ChEMBL 替代品。
原始 LIT-PCBA 基准来源：[JCIM paper](https://doi.org/10.1021/acs.jcim.0c00155)。

## 4. 模型：protein-conditioned learned adaptation operator

### 4.1 可识别的两阶段结构（PROJECT-FROZEN）

支持集为空时，基础预测器为

`f_theta(p, x) = pooled_source_prior(x) + neural_DTA_theta(p, x)`。

它只接收 recipient protein `p` 和 query ligand `x`，先在 source
meta-train 上训练。进入元适配阶段后冻结全部 base parameters，保证任何
support gain 都不是重新训练基础 DTA encoder 造成的。

对 recipient support
`S={(x_i,y_i)}_(i=1..k)`，模型先形成 measured residual token：

`t_i = pair_encoder(p,x_i) + residual_encoder(y_i-f_theta(p,x_i))`。

learned self-attention 将 `S` 编码成 adaptation state；每个 query
`x_q` 再对该 state 做 cross-attention，并由非线性 delta head 输出

`y_hat_q = f_theta(p,x_q) + Delta_phi(p,S,x_q)`。

`Delta_phi` 对 recipient protein、support measurements 和当前 query
同时有依赖。实现中没有 `torch.linalg.solve`、解析 ridge、GP posterior、
budget gate 或 deep-kernel branch。

### 4.2 反事实 support objective（PROJECT-FROZEN，来自所给 Figure 2 意见）

每个正 episode `(p,S+,Q)` 配三个 target-mismatched supports，recipient
protein 与 query 固定不变：

1. random：同 role/split/protocol/k 内随机选择其他 target；
2. protein-hard：选择 ESM-2 cosine 最接近的其他 target；
3. chemistry-matched：先最大化 support-only Murcko scaffold-set Jaccard，
   再以 ECFP4 support-centroid cosine 打破并列，选择其他 target；
4. label-swap：保持正确 support compounds/f0 不变，在线移植另一 target
   的 support residual，用于排除纯化学 arm 识别。

所有负样本 mining 均不读取 query labels；chemistry-matched 也不读取 query
chemistry。protein similarity 只用于挑选难负例，所有预测 arm 均输入原
recipient protein，避免模型从“wrong protein embedding”直接识别负例。

令 `L_rank(Q;S)` 为该 support 适配后在固定 query 上的 pairwise ranking
loss，则四分类 contrast 为：

`L_cf = -log exp(-L_rank(Q;S+)/tau) / sum_j exp(-L_rank(Q;S_j)/tau)`。

也就是说，contrast score 是 **post-adaptation query ranking performance**，
不是 support embedding 相似度。基础 query MSE/ranking loss 与 `L_cf` 联合
训练 adapter；正式 ablation 必须含 `lambda_cf=0`。

## 5. 正式实验矩阵

### E0. 数据与泄漏 gate（PROJECT-FROZEN）

正式训练前重新验证三个 package 的 size/hash、label firewall、nested k、
support/query parent/document/time、source/recipient target/accession/assay、
source meta-split component、三种 wrong-support target mismatch。任一硬 overlap
非零即停止，不运行模型。

### E1. 基础 DTA 性能（PROJECT-FROZEN）

目的：证明 support-free `f_theta(p,x)` 本身是合格的 DTA predictor，且
adapter 不以 support/query 化学泄漏弥补坏 base。

比较：pooled prior、protein-free base、完整 protein-ligand base；在 source
meta-test 与 recipient test 均报告 RMSE、MSE、MAE、Pearson、R²、CI、
Spearman、pairwise accuracy、NDCG@10、AUPR@7。recipient 标签只在正式外部
evaluation 阶段读取。base 完成后冻结；adapter 训练前后 support-free 输出
必须逐元素相同。

### E2. 核心 few-shot transfer（PROJECT-FROZEN）

在 63 个 unseen target、5 个冻结 draws、`k={1,3,5}` 上比较：

1. support-free base；
2. recipient-only intercept/slope calibration；
3. support kNN/similarity retrieval；
4. no-meta pooled ridge head；
5. fixed/closed-form deep-kernel posterior（现有 A2S-MDK）；
6. MAML；
7. ANIL；
8. A2S-CMAL（本文模型）。

1–5 只能作为非学习适配或闭式 baselines。6–7 是 learned meta baselines，
用于判断收益是否只是普通 gradient adaptation。所有方法必须使用相同
support/query rows、base representation、训练标签预算和五个 seeds：
`{1729,1731,1733,1741,1753}`。

primary endpoints 是 target-macro CI、within-target Spearman、NDCG@10 的
`adapted - support-free` gain；RMSE/MSE/MAE/Pearson/R²、pairwise accuracy、
AUPR@7 为 secondary。统计单位是 55 个 recipient homology components，
同时公开每 target/draw 结果、全体 63 target 汇总和 54-target homology-cold
汇总。

### E3. target-specificity / counterfactual support（PROJECT-FROZEN）

对同一个 `(recipient protein, query set)` 一次性计算：correct support、
random wrong target、protein-hard wrong target、chemistry-matched wrong target，
同化学组成的 label-swap，另做 support-label permutation。报告：

- `metric(correct) - metric(each wrong)`；
- `metric(correct) - metric(permuted labels)`；
- post-adaptation ranking-loss gap；
- 五 seed component-bootstrap interval。

只有 correct support 持续优于三种 wrong-target support，才能把增益解释为
“识别并使用 recipient measurements”。仅优于 random 而输给 hard negative，
只能解释为相似性捷径。

### E4. 机制 ablation（PROJECT-FROZEN）

固定其余条件，仅删除一个机制：

- `no-counterfactual`：`lambda_cf=0`；
- random-only、protein-hard-only、chemistry-match-only；
- embedding-similarity contrast（明确的 negative control）；
- zero/shuffled recipient protein；
- no support cross-attention；
- support labels permuted；
- no-meta ridge、A2S-MDK、MAML、ANIL。

已知 inert 的 budget gate 和 deep-kernel branch 不恢复为“创新模块”；它们只
保留历史失败基线。Figure 2 所述单 seed 结论只能当 pilot：counterfactual
contrast 当时是唯一 load-bearing component，MDK 在 k=5 未超过 ridge，且
ordered-vs-random 未被确认。正式结论必须来自上述五 seeds。

### E5. support protocol differential（PROJECT-FROZEN）

在 source meta-test 上比较 document-ordered 与 exchangeable-random episodes；
recipient 正式主结果只使用冻结 D0-R document-ordered episodes。由于 pilot
中 ordered 优于 random 的预测未成立，该项现在是机制诊断，不得作为既定
方向性证据。

## 6. 完整照搬的论文实验

### P1. AdaMBind 三数据集实验（PAPER-EXACT）

按 Wan et al. 2026 原样运行：

- target 是 task；每个 task 内分 support/query；
- random task split：meta-train/meta-validation/meta-test = 8:1:1；
- novel task split：CD-HIT 40% identity 聚类后，cluster 按 8:1:1 分配；
- majority support=40，few-shot support=5；
- support sensitivity：5、10、20、30、40，只用 R²画 sensitivity curve；
- datasets：BindingDB、KIBA、Davis；
- baselines：DeepDTA、HiSIF-DTA、ColdDTA、Co-VAE、PSICHIC、CML、
  MetaDTA、ZeroBind；non-meta baselines 按论文在 meta-test support 上 fine-tune
  5 steps；
- metrics：MSE、CI、R²、Spearman、Pearson；
- 每种方法 5 次 independent replication，报告 mean ± std；
- ablation：remove meta-learning module、remove adaptive task module、remove
  label-noise strategy；在 majority setting、5 seeds 下运行；
- cross-domain：novel-task split 下 KIBA → BindingDB。

这些结果是 published-protocol reproduction，不替代本项目的 target-specific
counterfactual test，也不把 random/CD-HIT split 描述成本项目的
document/parent/assay-closed D0-R。

### P2. cold-target representation control（PAPER-EXACT）

按 Nguyen et al., *Briefings in Bioinformatics* 2022 的 cold-target track
原样运行，不改 split：Davis train/valid/test = 15,708/3,877/4,964；PDBBind
v2019 = 9,134/2,282/2,595；validation/test targets 全部不出现在 train。
报告 RMSE、Pearson、Spearman、CI，超参数只由 validation 选择。该实验只检验
support-free cold-target representation，不证明 few-shot meta-adaptation。
来源：[paper](https://doi.org/10.1093/bib/bbac269)。

### P3. DeepDTA warm benchmark sanity check（PAPER-EXACT）

按 DeepDTA 使用 Davis/KIBA：随机分六份，一份 independent test，其余用于
5-fold hyperparameter selection；固定 test 上以五个训练组合报告平均 CI，
同时报告 MSE；baseline 为 KronRLS 与 SimBoost。该 track 只验证实现与经典
DTA benchmark 对齐，不支持 unseen-target claim。
来源：[DeepDTA](https://doi.org/10.1093/bioinformatics/bty593)。

### P4. LIT-PCBA early-recognition experiment（PAPER-EXACT，NOT RUN）

按 AdaMBind Figure 8 原样运行：15 protein tasks 中 ESR、TP53 为 test，其余
13 tasks 为 train；训练 task 丢弃 inactive，每 task 随机抽 15 actives 为
support，其余 actives 为 query；test task 同样以 15 actives 为 support，
其余 active+inactive 为 query。报告 EF@1%、EF@5%、Precision@10、
Precision@20、BEDROC(α=80.5)。不得把 ChEMBL pKi roster 替代此 screening
protocol。

## 7. 结论规则

不设新的任意效果量阈值。五 seed 正式结果必须满足以下逻辑，才能支持核心
论文命题：

1. A2S-CMAL 相对 support-free base 在 target-specific ranking primary
   endpoints 上给出正的 component-bootstrap interval；
2. A2S-CMAL 超过 no-meta ridge、fixed posterior/MDK，并与 MAML/ANIL 做
   equal-budget 比较；只在 k=1 胜出时，结论必须限定为 k=1；
3. correct support 同时优于 random、protein-hard、chemistry-matched wrong
   support、same-compound label-swap，并优于 label permutation；
4. shuffled/zero protein 显著削弱 correct-support advantage；
5. `no-counterfactual` 明显削弱 target-specific gain，才可称 counterfactual
   objective load-bearing；
6. support-free prediction 在冻结 base 前后完全不变；
7. 全体 recipient 与 homology-cold subgroup 结论分别陈述。

若只改善 RMSE 而不改善 target-specific ranking，结论是 calibration gain，
不是可迁移排序适配。若 correct/wrong support 无差异，结论是未识别 support
identity。若只在一个 seed 或一个 k 上成立，必须降级为 pilot，不写普遍机制。

既有 accuracy preregistration 的约 `0.12 RMSE` floor（当前 55 components
下 `sigma_delta=0.3` 对应 MDE80=0.1133）仍只
用于 accuracy claim；不新造替代阈值，也不把其用于 ranking claim。

## 8. 执行与算力边界

本机只允许：构建/校验 label-blind episode package、单元测试、source
meta-validation CUDA smoke 和 GPU profile。`--formal` 在未设置外部主机门禁
时会在读取 recipient labels 前退出。

当前 source validation/holdout mechanism gate 已失败，因而即使在外部主机
也不得开始 recipient formal run。最新失败诊断及下一轮全英文研究 handoff
见 `reports/active/CMAL_FAILURE_HANDOFF.md` 和
`reports/active/CMAL_EXTERNAL_AGENT_PROMPT.md`。在 source validation 与
source holdout 均以 paired component uncertainty 证实之前，recipient labels
继续封存。

```powershell
# 已存在 package 时只做 hash/audit；builder 默认拒绝覆盖
D:\anaconda\envs\drug\python.exe main.py prepare-cmal-data --device cuda

# 本机机制 smoke：不读 recipient labels
D:\anaconda\envs\drug\python.exe main.py a2s-cmal --smoke

# 指定外部 CUDA host；逐 seed 运行，不能在本机执行
$env:A2S_FORMAL_EXTERNAL='1'
D:\anaconda\envs\drug\python.exe main.py a2s-cmal --formal --seed 1729
```

正式外部运行必须保存：config、data content hash、seed、checkpoint、每
target/draw prediction、training log、GPU telemetry 和 bootstrap output。当前
`MECHANISM_SMOKE_ONLY` 文件不得进入论文性能表。

## 2026-08-01 Superseding Decision

上述 CMAL 正式训练方案已冻结，不得执行。平衡后的 v2 source-only
information gate 修复了旧设计 97% OOF 行集中于单折的问题，但 k=1/3/5 的
`Delta_label` 与 k=3/5 的 `Delta_assign` 均未得到正的 component-bootstrap
下界。决策为 `NO_GO_INFORMATION_NOT_ADMITTED`。完整证据与下一实验见
`reports/active/A2S_POST_REVIEW_V2_GATE_DECISION_2026-08-01.md`。在完成
label-free same-assay/MMP coverage/power census 并预注册局部门之前，不得训练
新 adapter、打开 locked-source/recipient labels 或执行 GitHub 发布。

## 9. 主要文献

- Wan et al., AdaMBind, *Nature Communications* 2026:
  <https://doi.org/10.1038/s41467-026-70554-5>
- Nguyen et al., cold-start interaction-knowledge transfer, *Briefings in
  Bioinformatics* 2022: <https://doi.org/10.1093/bib/bbac269>
- Öztürk et al., DeepDTA, *Bioinformatics* 2018:
  <https://doi.org/10.1093/bioinformatics/bty593>
- Nguyen et al., GraphDTA, *Bioinformatics* 2021:
  <https://doi.org/10.1093/bioinformatics/btaa921>
- Lee et al., MetaDTA, ICLR 2022 MLDD workshop:
  <https://openreview.net/forum?id=yzlif16IASM>
- Tran-Nguyen et al., LIT-PCBA, *Journal of Chemical Information and
  Modeling* 2020: <https://doi.org/10.1021/acs.jcim.0c00155>
