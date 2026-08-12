# Cold Target 少样本 DTA 元学习：数学理论审计与整合稿

> 状态：working theory，2026-08-12。本文不修改或覆盖
> `theory/FINAL_FROZEN_THEORY/`。冻结理论仍是仓库当前权威文本；本文记录其
> 可保留部分、已确认缺口，以及面向 Cold Target 的下一版统一理论。

## 0. 结论

现有理论不是整体错误，但尚未闭合为一个“可学习的 Cold Target 生物信息
理论”。目前可靠的部分是：

1. 旧 trace/section 理论给出了有限 support 下的精确 minimax 可辨识边界；
2. 新 band/simplex 理论给出了给定合法统计量 `z` 后的点态唯一正则目标；
3. 已消费的 main-v0 splits 显示 support-conditioned adaptation/calibration signal；
4. 在本次 governed BindingDB development panels 的 component bootstrap 下，X1 的
   observed-label rectangle double-difference magnitude 下界为正；它尚未从
   measurement noise 中识别 latent non-additivity。

当前没有被证明或被实验识别的是：

1. 从蛋白、配体和 support 到可迁移 interaction statistic 的充分性；
2. 对同源组件独立的 Cold Target 泛化；
3. 非加性交互的低秩性、descriptor 可迁移性或生物机制解释；
4. law-valued 输出到 RMSE、CI、Spearman、NDCG 的校准；
5. 当前 frozen theory 的完整 `continuity -> approximation -> consistency` 链。

因此，核心创新不应继续表述为“另一个 band/simplex 网络”，而应重新定位为：

$$
\boxed{
\text{component-level cold risk}
+\text{interaction quotient}
+\text{query-specific section radius}
+\text{cold-transport penalty}
+\text{selective-ranking guarantee}
}.
$$

其中，interaction representation 必须由部署可计算、在 source/deployment 有共享
支持的 descriptor 构成；descriptor 可以固定，也可以只用 source 学习。新靶标的
support-derived 自由状态必须受 `k` 维信息上限约束。共享深度表示的总维数不受
`k` 限制，只有从 `k` 个 support labels 识别的 task-specific 连续状态受限。

### 0.1 对 SAR-delta 核心困境的直接判定

最新理论能够解决“应当用什么数学对象归因”的问题，但尚未用数据解决
“该对象在 Cold Target 上是否存在且可学”的问题，更没有自动解决端到端 DTA。
严格分层如下：

| 层级 | 当前状态 | 理论作用 |
|---|---|---|
| 局部 SAR-delta 可预测性 | F-152/F-153 表观 PASS，后被 outcome-dependent construction 降级 | 尚无完整 G0-compliant 对称实验 |
| 排除 target-main、ligand-main、additive concat | 数学上可解决 | 使用 protein × ligand 的 rectangle quotient，而非原始 pair delta |
| Protein-conditioned transformation response | 尚未识别 | 必须通过新的 G2 cold quotient Gate |
| Support-conditioned SAR adaptation | 尚未识别 | 必须通过 G3a query-delta Gate；去掉 level 后可识别维数至多 $k-1$ |
| Scalar affinity / ranking 的端到端收益 | 尚未识别 | 必须由同一 scalar potential 通过 G3b，不能使用独立 edge head |
| Cold Target 泛化 | 只有条件性 bound | 还需 descriptor approximation、component generalization 和 transport 条件 |

因此最短合法研究链是

$$
\boxed{
\text{G2: cold protein-conditioned quotient}
\longrightarrow
\text{G3a: support-conditioned SAR delta}
\longrightarrow
\text{G3b: shared-potential scalar DTA}
\longrightarrow
\text{V1/CSMO preservation audit}.
}
$$

F-152/F-153 不能跳过 G2；G2 的 PASS 也不能跳过 G3a/G3b。故当前仍保持
`V1_INTEGRATION_AUTHORIZED=false` 和 `BIOLOGICAL_CLAIM_AUTHORIZED=false`。

## 1. 使用范围与基本假设

本文明确采用以下假设。任何一项不满足时，相应定理不得用于论文主张。

### A1. 观测类型

蛋白观测空间为 $\mathcal P$，配体观测空间为 $\mathcal L$，assay/context
空间为 $\mathcal A$。部署可见对象为蛋白序列或合法结构信息 $x_p$、配体分子图
$x_l$、assay 信息 $a$ 和 support labels；query label 在推理时不可见。

### A2. Affinity 可比性

潜在连续亲和力记为 $Y^*_{p,a,l}\in V=[v_{\min},v_{\max}]$。不同单位、方向、
端点和 censoring 必须先被类型化。如果不同 assay 之间没有已识别的数值变换，
绝对 RMSE 只能在同一 context 内定义；单调变换最多保留 within-task ranking。

### A3. 独立单位

在声明的同源阈值 $\rho$ 下，以完整依赖组件 $C$ 为潜在独立抽样单位。同一组件内
的多个蛋白、多个 support draw 和多个 query 允许相关，不能重复计作独立样本。
组件不相交只消除已登记的交叉依赖，不自动证明 superpopulation IID；第 6 节的
泛化定理会把 component-IID 作为额外假设。

### A4. Episode 合法性

一个任务为 $\tau=(p,a)$，episode 为

$$
E_\tau=(x_p,a,S_k,Q_m,Y_Q),\qquad
S_k=\{(x_{l_i},O_i)\}_{i=1}^k,
$$

其中 $O_i$ 是点值或合法 censoring interval。support/query 的 target-ligand
单元不交；重复测量先按预先声明的规则聚合；support 的选择不依赖隐藏 query
labels。任何 ligand pair 的 identity、inclusion、enumeration、subsampling/truncation
和方向都必须在读取 outcome 前冻结。方向由 canonical key 决定，或在训练和评估中
同时纳入 $(l_i,l_j,\Delta)$ 与 $(l_j,l_i,-\Delta)$。先按 affinity label 排序再定义
左右端点或截取前若干 pairs 都属于 label leakage，不能用于 SAR-delta admission。

### A5. 结构假设与证据分离

X1 的非零 raw rectangle RMS 只作为 observed-label double-difference magnitude 的
经验事实。即使 latent surface 完全 additive，四格测量噪声也会产生正 RMS，因此
latent non-additivity 还需要 replicate/noise-aware null 或去卷积。低秩、平滑、
蛋白 descriptor 充分性、跨组件机制不变性和 cold transport 都是待检验结构假设，
不是 X1 的推论。

### A6. Selection、确认集与证书目标

必须声明 affinity 行的 observation/missingness mechanism；否则 estimand 只是
observed-assay risk，不是所有潜在 DTA pairs 的总体风险。Censoring interval $O_i$
与额外误差界 $\varepsilon_i$ 必须有非重复语义。所有 architecture、descriptor、
dependency threshold、margin condition、Gate threshold 和 certificate cutoff 只在
source/calibration components 上选择；最终 confirmation components 不参与这些
选择。依赖图需明确区分 homology、document、assay、scaffold 和 repeated-measure
边。

所有 pair/rectangle 数据构造还必须登记 pair-selection/orientation function、可读字段
和 truncation 时点。Gate 运行前执行 reverse-pair audit；仅有架构形式上的
antisymmetry 不足以修复 outcome-dependent sample construction。

## 2. 现有理论审计

### 2.1 Blocking gap：`B(z)` 的变化在连续性证明中被遗漏

Frozen theory 定义

$$
J_\mu(z,p)=L_0(z,B(z)p)+\frac{\mu}{2}\lVert p\rVert^2,
\qquad
g_\mu^*(z)=\arg\min_{p\in\Delta_m}J_\mu(z,p).
$$

`S-CONT` 只控制固定 $\beta$ 时的

$$
|L_0(z,\beta)-L_0(z',\beta)|,
$$

但原证明把它直接用于

$$
L_0(z,B(z)p)-L_0(z',B(z')p).
$$

这里的 band argument 同时变化，因此证明缺一项。实现中的 context 又由离散
bucketization 产生，所以 $B(z)$ 的跳变是实际存在的，不是抽象边界情况。

一个严格反例是：令 band polytope 含标量线段 $[0,1]$，$m=1$，固定 anchor
$\beta_1=0$，并令

$$
\beta_0(z)=
\begin{cases}
0,&z<1/2,\\
1,&z\ge 1/2.
\end{cases}
$$

取 $L_0(z,\beta)=\beta$、$\mu=1$。此时 $L_0$ 完全不依赖 $z$，所以原
`S-CONT` 模数为零；但左侧唯一最小点是 $(1/2,1/2)$，右侧唯一最小点是
$(0,1)$。因此 $g_\mu^*$ 跳变。

结论：点态 existence/uniqueness 仍成立，measurability 也可保留；但目标连续性、
multilinear uniform approximation 以及依赖它的一致性结论当前不能由 frozen
assumptions 推出，现有证明链不成立。特定的 $L_0,B$ 仍可能恰好具有这些性质。

最小修复是直接假设 operative-risk continuity：

$$
\tag{OR-CONT}
\sup_{p\in\Delta_m}
|L_0(z,B(z)p)-L_0(z',B(z')p)|
\le \omega_J(d_Z(z,z')),
\qquad \omega_J(t)\to0\quad(t\downarrow0).
$$

此时原强凸性论证可正确给出

$$
\lVert g_\mu^*(z)-g_\mu^*(z')\rVert
\le \sqrt{\frac{2\omega_J(d_Z(z,z'))}{\mu}}.
$$

若分别假设

$$
\lVert B(z)-B(z')\rVert_{\rm op}\le\omega_B(d_Z(z,z')),
$$

则可取

$$
\omega_J\le\omega_\ell+L_{\rm Lip}\omega_B.
$$

另一种同样合法的修复是把
$Z=\bigsqcup_{c\in C_\kappa}Z_c$ 定义成 context 间有正距离的离散并，使 mesh
和插值永不跨 context 边界。

### 2.2 两代 operator 不是同一个对象

旧理论的 canonical object 是 support section 的 Chebyshev center 和 radius：

$$
A_{\rm sec}(\mathcal F,S,q,\varepsilon)
=\left(\operatorname{cen}I_q(S),
\frac12\operatorname{diam}I_q(S)\right).
$$

新 frozen theory 的唯一 operator 是

$$
\mathsf A(F,z)=K(B(z)F(z)),
$$

并明确排除 support intersection。两者可以分别成立，但不是同一个 theorem，也
不存在已证明的等价实现。

权威层级上没有形式冲突：`FINAL_FROZEN_THEORY` 高于 handoff。论文逻辑上仍缺
一个接口：旧理论应作为 information/admission layer，决定什么 support statistic
有资格进入下游；新理论只能作为给定该 statistic 后的 law decoder。

### 2.3 两个 positive ridge 必须分开

Frozen law layer 的 $\mu\lVert p\rVert^2/2$ 作用于 simplex coefficient $p$，只
保证强凸性和唯一 regularized target。该惩罚项在 simplex 上的唯一 minimizer，
也就是 regularization center，是 uniform mixture；当该项相对 base loss 占主导时
目标趋向它，而不是 population column $e_0$。若需要向 population prior 收缩，应另行定义
$\lVert p-e_0\rVert^2$，这会改变目标和定理。

Meta-Section 的 $\lambda\lVert c\rVert^2$ 作用于 task coefficient $c$，是
Tikhonov ridge。它是稳定估计器，但一般不等于 $L^\infty$ section 的 Chebyshev
center，也不会在 $\lambda>0$ 时保证 exact support interpolation。

### 2.4 `2h` 不是误差下界

Frozen transfer inequality

$$
d_H(K(\beta),K(\beta'))
\le D\lVert\beta-\beta'\rVert+2h
$$

中的 $2h$ 是该上界的 non-vanishing remainder。取 $F=g_\mu^*$ 时真实距离为零，
所以不能把 $2h$ 解释成所有模型不可突破的 error floor。

### 2.5 泛化单位与 Cold Target 部署不对齐

Frozen theory 把 episode $T_i=(S_i,Q_i,Y_i)$ 作为 IID 样本。实际同一蛋白或同一
同源组件下的多个 episodes 相关。用于 Cold Target 论文时，主样本量必须是独立
组件数 $N_C$，不是 episode 数、query 数或 target-ligand 行数。

Target-disjoint 不自动意味着 OOD：若完整组件从同一 superpopulation 随机抽取并
随机分配，source/deployment 可在组件层 IID。若测试条件显式要求离训练蛋白至少
某距离，或使用 novel-family stress，则是分布迁移，必须保留 transport term。

### 2.6 其他需修复的边界

1. 当前 $\Gamma_N$ 使用 $\log(\Lambda N)$；当 $\Lambda N\le1$ 时可为负或未定义，
   且文本未写 net 延拓项。若参数空间在所用 norm 下半径不超过 $R_\Omega$、
   ambient dimension 为 $D_N$、loss 有界于 $M_L$ 且对参数是 $\Lambda$-Lipschitz，
   一个合法的有限维 net bound 是

   $$
   \Gamma_N^{\rm net}
   =\inf_{\eta>0}\left\{
   2\Lambda\eta
   +M_L\sqrt{\frac{
   D_N\log(1+2R_\Omega/\eta)+\log(2/\delta)}{2N}}
   \right\}.
   $$

   这里 $2\Lambda\eta$ 是从 net 延拓到全参数空间的两侧误差。具体 norm 的
   covering number 可给出更紧常数，但不能省略该尺度项。

2. 若 $K(\beta)$ 指连续 $V$ 上仅受普通有限 CDF 点值约束的所有 laws，则 CDF
   upper constraints 对弱/$W_1$ 极限未必闭。可采用 archive 中的
   closed-lower/open-upper convention

   $$
   P([v_{\min},t_j])\ge l_j,
   \qquad P([v_{\min},t_j))\le u_j,
   $$

   并由 Portmanteau 证明闭性；更简单的计算路线是把 $K_h(\beta)$ 明确定义为固定
   grid-supported PMF polytope。二者必须选定一个，不能只写“satisfying the band”。
3. Full CSMO 只在结构上包含 frozen $\mathcal H_N$；其 learned gates 和多 view
   不能自动继承 frozen class 的 approximation/generalization bound。
4. Law-valued operator 没有唯一 scalar selector，因此不能自动推出 RMSE；逐 query
   marginal laws 也不能推出 CI、Spearman 或 NDCG。

## 3. 实验证据的严格含义

### 3.1 已识别

1. 在已消费的 main-v0 splits 上，存在描述性的 support-conditioned 低维
   adaptation/calibration signal：full-correct 明显优于 `d=0`、ligand-only、
   zero/foreign/permuted support。
2. 在这些 splits 上，该收益平均主要是 target-level calibration；cluster-robust
   correct-protein specificity 未被识别。
3. X1 在本次 governed development 的 21 panels、12 components 上得到正的
   component-bootstrap lower bound，说明该 split 的 observed labels 有
   double-difference magnitude。未做 replicate/noise correction，因此它不单独识别
   latent affinity surface 的 non-additivity，也不构成跨 split 复现。

### 3.2 未识别

1. F-159 只是一个可使用 target ID 的固定 oracle arm：
   `T-BASIS -> source PCA(d=5) -> ridge(lambda=1,k=5)`。它关闭这条预注册 recipe，
   不是该 feature family 的性能上界，也不是所有 nonlinear/descriptor model 的
   信息论上界。
2. X2 没有检验到低秩命题：intersection 前 train 有 17,433 个 exact
   transformations，development 有 163 个，但共享数为零，算法在拟合前即 fail
   closed。其结论只是 exact
   transformation ID 在当前 train/development split 上不可迁移。
3. 非零 rectangle 不推出低秩、双线性、因果耦合、接触机制或可迁移蛋白 encoder。
4. 历史 `p=1.0` rewiring diagnostic 已因 permutation 不合法和调参泄漏被降级，
   不能作为不可辨识证明。“一维状态必然不可辨”同样错误；可辨识性取决于函数族
   和 support design。

### 3.3 F-152 至 F-158 的重新归因

1. F-152/F-153 先按真实 `p_value/pK` 排序 ligand endpoints，再定义 label delta 和
   descriptor delta。其方向读取了 outcome，所以 target-main/intercept 可以预测被固定
   符号的差值幅度。两次 PASS 只能称为 outcome-oriented pair recipe 的表观预测性，
   不能作为合法 SAR transfer 或 target-conditioned interaction 证据。
2. F-153 在 BindingDB train 上重新拟合 ridge，并未把 ChEMBL 权重或模型迁移到
   BindingDB；只能说相同 recipe 在两个数据源分别获得表观结果，不能称 cross-corpus
   representation/model transfer。
3. F-156 正确发现 additive/target-main shortcut：`A/P` 很强但违反反对称性，bilinear
   `I` 虽结构反对称，却未胜 ligand-only、wrong-target 或 shuffled-target。
4. F-157 只对已选中的 pairs 做 forward/reverse 增广；上游仍先按 `pK` 排序并按组截取
   前 100 个 pairs。未截取/截取计数分别为 train `27,222/20,423`、development
   `1,122/1,033`，所以 inclusion/truncation 仍读取 outcome。增广后
   `P=Z=0.580072`、`A=L=0.601045`；这只说明固定选中样本上的 target-main 方向捷径被
   消除，不是完整 G0-compliant 对称实验，也不能把 `L` 未胜 zero 当作正式失败证据。
5. F-154 只关闭当前 1D neighbor-mean scalar recipe；F-155 只关闭当前八维
   unordered edge-summary 加 ridge recipe。它们没有否定所有 scalar potentials 或
   permutation-invariant set encoders。
6. F-158 只关闭固定 descriptors、linear ridge、当前监督与 split 的 U1
   score-difference arm。停止扩大 encoder 是项目 stop-tree 决策，不是 capacity
   impossibility theorem。

| Gate | 关键结果 | 本稿中的最终作用域 |
|---|---|---|
| F-152 | 5 assays/415 pairs；MSE `0.426781` vs `0.786496`；LCB `0.146371` | outcome-oriented protocol 的表观 PASS，降级 |
| F-153 | 1,033 dev pairs/8 components；`0.233714` vs `0.580072`；LCB `0.270172` | BindingDB 内重拟合且方向读 label；不是跨库模型迁移 |
| F-154 | 四个 CQ arms 均 `0.236361` | 只拒绝当前 1D neighbor-mean recipe |
| F-155 | correct = deranged `0.234537`；全部 LCB 失败 | 只拒绝当前八维 unordered summary recipe |
| F-156 | `I=1.241352`，`L=0.600399`；I 对 L/wrong/shuffle 全失败 | target-conditioned interaction 未识别 |
| F-157 | `P=Z=0.580072`，`A=L=0.601045`；全部 contrasts 失败 | 只在 outcome-selected/truncated pairs 上修复 direction；不是 G0-compliant Gate |
| F-158 | `R=1.428732` vs `L=0.601045`；LCB `-2.437710` | 关闭当前 U1 arm，不是 potential-class impossibility |

因此当前最强合法状态应改写为：

```text
OUTCOME_ORIENTED_SAR_DELTA_PREDICTABILITY_OBSERVED
G0_COMPLIANT_SYMMETRIC_SAR_DELTA_NOT_YET_RUN
TARGET_CONDITIONING_NOT_IDENTIFIED
CURRENT_NEIGHBOR_MEAN_AND_8D_EDGE_SUMMARY_RECIPES_REJECTED
CURRENT_PAIR_SCORE_DIFFERENCE_U1_REJECTED
END_TO_END_COLD_TARGET_UTILITY_OPEN
V1_INTEGRATION_AUTHORIZED=false
BIOLOGICAL_CLAIM_AUTHORIZED=false
```

## 4. 统一的 Cold Target 数学对象

### 4.1 Component-macro deployment risk

令 source、validation、deployment 组件集合两两不交：

$$
\mathcal C_{\rm src}\cap\mathcal C_{\rm val}
=\mathcal C_{\rm src}\cap\mathcal C_{\rm dep}
=\mathcal C_{\rm val}\cap\mathcal C_{\rm dep}=\varnothing.
$$

对 predictor $h$，定义回归风险

$$
R_{\rm reg}^{\rm dep}(h)
=\mathbb E_{C\sim\Pi_{\rm dep}}
  \mathbb E_{\tau\mid C}
  \mathbb E_{(S,Q,Y_Q)\mid\tau,C}
  \frac1{|Q|}\sum_{q\in Q}
  \ell_{\rm reg}(h(x_p,a,S,q),Y^*_{\tau q}).
$$

这里 $\Pi_{\rm dep}$ 是预先声明的组件分布；当前 benchmark 的主 estimand 取组件上
的 uniform empirical law。随后再抽 task/context 和 episode，防止巨型组件或高
测量密度 target 支配 estimand。若以 noisy observation 替代 $Y^*$，风险还包含
不可约测量噪声。

Within-target ranking risk 必须在同一 $(p,a)$ 的 query pair 上定义。令
$\Delta f=f_\tau(l)-f_\tau(l')$、
$\Delta h=h_\tau(l)-h_\tau(l')$。在 $\Delta f\ne0$ 的 pair 上采用

$$
\ell_{\rm rank}
=\mathbf1\{\Delta f\,\Delta h<0\}
+\frac12\mathbf1\{\Delta h=0\},
\qquad
R_{\rm rank}(h)=\mathbb E[\ell_{\rm rank}\mid\Delta f\ne0],
$$

即预测 tie 计半错；true ties 的纳入或排除必须另行声明。NDCG 还需要固定候选集
和 relevance map。

选择性预测输出 $(\widehat f,r,A)$，$A\in\{0,1\}$：

$$
R_\lambda^{\rm sel}
=\mathbb E[A\ell(\widehat f,Y^*)+(1-A)\lambda],
\qquad \operatorname{Cov}=\Pr(A=1).
$$

任何 selective 结果必须同时报告 coverage；否则 `always abstain` 是平凡解。

### 4.2 Additive 主效应与 interaction quotient

在每个 assay/context 内写成

$$
f_{p,a}(l)=\mu_a+u_{p,a}+v_a(l)+g_{p,a}(l).
$$

$u$ 只能改变绝对 location，在同靶标排序差分中严格抵消；$v$ 是 protein-blind
ligand shortcut；只有 $g$ 能产生 target-specific ligand reordering。可以用参考
测度下的零边际约束使分解唯一，但更稳健的对象是 gauge-invariant rectangle。
四个 cells 必须共享同一个可比 assay/context $a$：

$$
\mathcal R_a f(p,p';l,l')
=f_{p,a}(l')-f_{p,a}(l)-f_{p',a}(l')+f_{p',a}(l).
$$

它精确消去 $\mu_a+u_{p,a}+v_a(l)$。

**命题 1（interaction quotient 的完备性）。** 对固定 $a$，在完整 Cartesian
domain $\mathcal P_0\times\mathcal L_0$ 上，$\mathcal R_a f=0$ 对所有 rectangles
成立，当且仅当存在 $u_a,v_a$ 使
$f_{p,a}(l)=u_a(p)+v_a(l)$。

**证明。** “若”方向由代数消去得到。反向固定 $(p_0,l_0)$，令
$u_a(p)=f_{p,a}(l_0)$、
$v_a(l)=f_{p_0,a}(l)-f_{p_0,a}(l_0)$。对 rectangle
$(p,p_0;l_0,l)$ 使用 $\mathcal R_a f=0$ 即得
$f_{p,a}(l)=u_a(p)+v_a(l)$。$\square$

在稀疏观测图上，只有闭合 rectangles/cycles 所在 quotient 被识别；未观测区域的
interaction 仍需要 representation/transport 假设。

**命题 2（pair-delta、rectangle 与 shortcut 边界）。** 在命题 1 的完整 Cartesian
domain 上，定义 ligand difference 与
protein difference operator

$$
D_Lf(p;l,l')=f_{p,a}(l')-f_{p,a}(l),
$$

$$
D_P\delta(p,p';l,l')
=\delta(p;l,l')-\delta(p';l,l').
$$

则

$$
\mathcal R_af=D_PD_Lf.
$$

并且：

1. $D_Lf$ 对 ligand swap 反对称，$\mathcal R_af$ 对 protein swap 和 ligand swap
   分别反对称；
2. $\mathcal R_af\equiv0$ 当且仅当 $D_Lf$ 与 $p$ 无关，当且仅当
   $f_{p,a}(l)=u_a(p)+v_a(l)$；
3. 若 edge predictor 写成 $\delta(p;l,l')=a(p)+b(l,l')$，并在完整的正反 pair
   support 上满足反对称性，则 $a(p)$ 必为与 $p$ 无关的常数；该常数可通过 gauge
   吸收到 $b$，在 $a=0$ gauge 下 $b(l,l')=-b(l',l)$；
4. 在连通 ligand graph 上，antisymmetric edge field $\delta_p(l,l')$ 存在 scalar
   potential $s_p(l)$ 使
   $\delta_p(l,l')=s_p(l')-s_p(l)$，当且仅当它沿每个 ligand cycle 的 circulation
   为零；$s_p$ 只在加性常数意义下唯一。

**证明。** 第一式直接由两次差分展开。第 2 项由命题 1 得到。对第 3 项，正反 pair
相加给 $2a(p)=-b(l,l')-b(l',l)$，右侧不依赖 $p$，所以 $a$ 为常数；作 gauge
$a\mapsto0$ 后 $b$ 为 odd。第 4 项是连通图上 edge field 为 exact coboundary 的充要条件：
potential difference 的 cycle sum 必为零；反之固定根节点并沿任一路径累加，zero-cycle
保证路径无关。$\square$

这精确解释 F-156/F-157：对称化会消除 target-main，并使 additive concat 退化成
ligand-only odd edge，但反对称性本身不排除 ligand-only shortcut。还必须明确，
rectangle 与 scalar pair potential 不是互斥对象。一般
$s(p,l')-s(p,l)$ 的 rectangle 可以非零；当前第 5.2 节的双线性模型本身就来自
$s(p,l)=\rho(p)^\top M_ae(l)$。因此 F-158 只否定其冻结的低容量 score class/训练
recipe，绝不能推出所有 pair potentials 失败。

**命题 3（定量 quotient reconstruction）。** 固定 context $a$。对两个实值 surface
$f,g$，令 $e_a(p,l)=f_{p,a}(l)-g_{p,a}(l)$。假设存在参考 $(p_0,l_0)$，使对所有待
控制的 $(p,l)$，四个 cells $(p,l),(p,l_0),(p_0,l),(p_0,l_0)$ 都属于同一 $a$ 下的
reference-star domain。定义

$$
u_a(p)=e_a(p,l_0),
\qquad
v_a(l)=e_a(p_0,l)-e_a(p_0,l_0).
$$

则

$$
e_a(p,l)-u_a(p)-v_a(l)
=\mathcal R_a(f-g)(p,p_0;l_0,l).
$$

**证明。** 展开右侧并代入 $u,v$ 即得。$\square$

所以 reference rectangles 上的 uniform error $\varepsilon_R$ 可以控制模 additive
gauge 后的 interaction uniform error。但 X1 的 raw average rectangle MSE 或任意
平均风险不能自动推出逐 query error；这还需要 finite-design spectral gap、
norm-equivalence 或显式 coverage 条件。Additive gauge $u(p)+v(l)$ 也必须由共享
baseline 与新 target support 另行锚定。该命题是 rectangle head 到 scalar affinity
之间的第一个必要桥，而不是完整的端到端定理。

**推论 3.1（scalar bridge 的确定性误差分解）。** 在命题 3 的固定 context 和
reference-star domain 上，若

$$
\sup_{p,l}|u_a(p)+v_a(l)|\le\varepsilon_{\rm add},
\qquad
\sup_{p,l}|\mathcal R_a(f-g)(p,p_0;l_0,l)|\le\varepsilon_R,
$$

则

$$
\sup_{p,l}|f_{p,a}(l)-g_{p,a}(l)|
\le\varepsilon_{\rm add}+\varepsilon_R.
$$

**证明。** 对命题 3 的恒等式取绝对值并使用三角不等式。$\square$

在 MetaSieve 中，$\varepsilon_{\rm add}$ 必须由共享 scalar baseline、target-level
calibration 和 support section 共同控制；$\varepsilon_R$ 必须来自同一个 scalar
potential 的 rectangle error，而不是不相干的 edge head。Representation、finite-sample、
transport 和 measurement errors 还要分别进入第 5.4 节的 total certificate。

### 4.3 无结构时的 Cold Target 下界

**定理 1（no-free-lunch）。** 对新 target 的 query $q\notin S$，若存在两个数据
机制 $P_0,P_1$，它们对 source archive、所有合法 covariates 和 support 的联合分布
相同，但 $f_0(q)=v_{\min}$、$f_1(q)=v_{\max}$，则任意 estimator 满足

$$
\inf_{\widehat f}\sup_{j\in\{0,1\}}
\mathbb E_j|\widehat f(q)-f_j(q)|
\ge\frac{\operatorname{diam}(V)}2,
$$

$$
\inf_{\widehat f}\sup_j
\mathbb E_j(\widehat f(q)-f_j(q))^2
\ge\frac{\operatorname{diam}(V)^2}{4}.
$$

若两个 query 的真排序在两个不可区分机制中相反，则在上述 predicted-tie 计半错
的 loss 下，最坏 pairwise ranking error 至少为 $1/2$。

**证明。** 两个机制产生相同可见输入，故 estimator 在二者下的随机输出 $A$ 同分布。
对每个 $a$，$|a-v_{\min}|+|a-v_{\max}|\ge\operatorname{diam}(V)$，积分即得绝对
损失结论。令 $\bar v=(v_{\min}+v_{\max})/2$、$D=\operatorname{diam}(V)$，则

$$
\frac12\mathbb E[(A-v_{\min})^2+(A-v_{\max})^2]
=\mathbb E(A-\bar v)^2+\frac{D^2}{4}\ge\frac{D^2}{4}.
$$

相反排序的两机制中，严格预测至少在一个机制下错，预测 tie 则在两边各损失
$1/2$；随机化后同样由两机制平均风险不小于 $1/2$ 得证。$\square$

这一定理说明 meta-learning 不能从 task 数量本身创造 Cold Target 信息；必须声明
跨任务共享的 interaction structure。

### 4.4 Representation admission theorem

令完整合法输入为

$$
X=(x_p,x_l,a,S),
$$

在共同 probability space 上假设 $X$ 可测、$Z=z(X)$ 可测。定义

$$
R_Z^*=\inf_{g\;\sigma(Z)\text{-measurable}}\mathbb E(Y^*-g)^2,
\qquad
R_X^*=\inf_{g\;\sigma(X)\text{-measurable}}\mathbb E(Y^*-g)^2.
$$

$Y^*\in V$ 保证二阶可积。平方损失下有精确恒等式：

**定理 2（压缩的 Bayes 代价）。**

$$
R_Z^*-R_X^*
=\mathbb E\operatorname{Var}
\left(\mathbb E[Y^*\mid X]\mid Z\right)\ge0.
$$

**证明。** $Z$ 是 $X$ 的函数。分别写出平方损失 Bayes risk
$\mathbb E\operatorname{Var}(Y^*\mid Z)$ 和
$\mathbb E\operatorname{Var}(Y^*\mid X)$，再应用条件全方差公式。$\square$

因此，严格风险差成立的充要判据是
$\operatorname{Var}(\mathbb E[Y^*\mid X]\mid Z)>0$ 在正概率集合上成立。仅仅观察到
$Z$ 对 intervention “近似不变”并不足以推出严格下界：微小变化仍可能可解码，
intervention 也可能离开数据支持。Wrong-protein、shuffled-protein 和
permuted-support 应作为 intended biological mechanism 的必要 falsification/admission
diagnostics；定量下界还需不变性误差、decoder regularity 和 in-support intervention
假设。

## 5. 最小可学习结构

### 5.1 两种维度必须分开

设共享深度 encoder 的维数为 $r_{\rm shared}$，task coefficient 的设计维数为
$r_{\rm task}$。给定 support design，真正从 labels 可识别的连续维数是

$$
d_{\rm id}(S)=\operatorname{rank}(\Phi_S)
\le\min\{r_{\rm task},k\}.
$$

完整恢复整个 task coefficient 需要
$\operatorname{rank}(\Phi_S)=r_{\rm task}$，因而要求 $r_{\rm task}\le k$；但只对
特定 query 作预测时，即使 $r_{\rm task}>k$ 也可能可识别。理论从不要求
$r_{\rm shared}\le k$，也不要求整个统计量 $z$ 的维数不超过 $k$。蛋白、query
ligand 和 assay covariates 不由 support labels “创造”，可以高维；只有实际从
support 提取的 member-specific 连续信息维数受 `k` 限制。

以下 $k-1$ 结论有明确前提：$k\ge1$ 是原始 support cells 的数量（不是已经形成的
delta observations），fixed support design 下 residual observation model 为
$r_S=\mu\mathbf1+\Phi_Sc$，且未知自由 nuisance intercept $\mu$ 与 $c$ 都由同一组
$k$ labels 估计。若把该 target-level calibration 与 SAR-specific adaptation 分开，令

$$
H_k=I_k-\frac1k\mathbf1\mathbf1^\top,
\qquad
\Phi_S^\Delta=H_k\Phi_S.
$$

则

$$
d_{\rm SAR,id}(S)
=\operatorname{rank}(\Phi_S^\Delta)
\le\min\{r_{\rm task},k-1\}.
$$

这是因为 $\operatorname{rank}(H_k)=k-1$。因此在这些前提下，`k=5` 时，在与 support level
正交的部分最多识别四个连续 SAR coordinates；把五维非截距 task state 宣称为由五个
labels 完整识别，会重新混入 calibration。具体 query 是否可识别仍取决于 centered
row-space，而不是只看维数。若 $\mu$ 已知，或输入本身是独立 delta observations，
不能机械使用该 $k-1$ 上界。

### 5.2 Descriptor-based crossed interaction

Exact transformation ID 应替换为部署时可计算、在 source/deployment 有共享支持的
antisymmetric descriptor：

$$
\psi(l,l')=-\psi(l',l).
$$

$\psi$ 的 orientation 必须由 canonical ligand key、固定 graph edit direction 或其他
outcome-independent 规则决定；另一合法方案是完整 forward/reverse augmentation。
F-156 中 `L/I` 的零 antisymmetry error 只是无截距线性结构在输入翻号下的恒等式，
F-157 中 `P` 的反对称性来自预测塌缩为零、`A` 来自退化为 `L`。这些结构审计都不
等于 protein-conditioned mechanism admission。

令 $\rho(p)$ 为蛋白 descriptor。一个低容量 inductive rectangle 候选 family 是

$$
\widehat{\mathcal Rf}(p,p';l,l')
=\bigl(\rho(p)-\rho(p')\bigr)^\top
M_a\psi(l,l'),
\qquad \operatorname{rank}(M_a)\le r_{\rm int}.
$$

若 $\psi(l,l')=e(l')-e(l)$，它来自 pair potential
$g(p,l)=\rho(p)^\top M_a e(l)$，并自动满足两种交换的反对称性。

这是一个可证伪的结构假设，不是 X1 的结论。其可识别至少要求：

1. source rectangles 在 protein-contrast 与 transformation descriptor 方向上有足够
   Gram/covariance rank；
2. development descriptor 位于 source descriptor 的覆盖区，而不是只比较 exact ID；
3. 拆分按依赖组件进行；
4. correct model 同时优于 additive、ligand-delta、wrong-target 和 shuffled-target；
5. 对反向 ligand pair 与反向 protein pair 的反对称性由结构保证；
6. 对稀疏/不精确 side information 明确保留 approximation residual。
7. 若最终需要 scalar affinity，edge predictor 必须来自共享 potential，或在声明的
   ligand domain 上结构性满足 zero-circulation；observed-cycle audit 只能提供部分经验
   检查。仅有 unordered edge distribution 不足以给出 query score。

这与 inductive matrix completion 的数学范式一致：低秩只有在 side information
覆盖真实 row/column space 且观测设计充分时才可恢复。当前 X2 的 one-hot exact ID
不满足 deployment coverage。

### 5.3 可迁移 baseline 加 support section

用 source components 学共享对象 $b_\theta$ 和 $\phi_\theta$：

$$
f_{\tau,c}(l)
=b_\theta(x_p,x_l,a)
+\phi_\theta(x_p,x_l,a)^\top c_\tau,
\qquad \lVert c_\tau\rVert_2\le R,
$$

其中：

1. $b_\theta$ 是可包含 crossed-interaction pretraining 的 zero-shot pair predictor；
2. $\phi_\theta$ 的非截距坐标必须随 ligand 变化，否则只能学习 calibration；
3. $\phi_\theta$ 必须依赖正确蛋白且通过 derangement Gate；
4. $c_\tau$ 是唯一的 task-specific state channel；support 通常只识别其 affine
   fiber 或 row-space projection，只有 full-rank 时整个 $c_\tau$ 才唯一；
5. 截距 coordinate 与 centered interaction coordinates 分开报告。

对 SAR 路线进一步要求一个 single shared scalar-potential family。把共享 baseline 分成
不参与 SAR arm 差异的 support-level/ligand additive 部分
$b_{\theta,S}^{\rm add}$ 和 protein-conditioned potential：

$$
b_{\theta,S}^{\rm add}(p,l,a)
=\mu_\theta(a)+u_\theta(p,a)+v_\theta(l,a)
+\kappa_\theta(T_{\rm level}(S),p,a),
$$

其中 $T_{\rm level}(S)$ 只能是对 ligand identity 与 label binding 的置换不变 level
统计量，且 $\kappa_\theta$ 不依赖 query ligand。因而对固定 $(p,S,a)$，该 baseline
不含 protein $\times$ ligand interaction，并在 protein-contrast rectangle 中消去；
ligand difference 中的 $v_\theta(l,a)$ 不消去，但它在所有 SAR arms 中固定相同。

$$
s_{\theta,c}(p,l,a)
=s_\theta(p,l,a)+\phi_\theta(p,l,a)^\top c,
$$

$$
f_{\theta,c}(p,l,a\mid S)
=b_{\theta,S}^{\rm add}(p,l,a)+s_{\theta,c}(p,l,a).
$$

G2 的 frozen zero-shot arm 取 $c=0$；G3a 只从 support 得到
$\widehat c_p(S)$；G3b 使用 $f_{\theta,\widehat c_p(S)}(\cdot\mid S)$。所有 edge、delta 和 rectangle
必须由这一 family 的 scalar outputs 相减，不允许另设不共享该 potential 的 head。
$b_{\theta,S}^{\rm add}$ 可读取合法 support 以完成 level calibration，但不得使用 query
labels；它的算法、参数和输入在所有 SAR arms 中完全相同，并先于 SAR arm 比较冻结。
因此 arm 间差异只能来自 $s_{\theta,c}$，不能由 support-level calibration 或
ligand-only baseline 冒充。这里的 “shared” 是冻结 gauge 和 channel parameterization
后的架构约束，不是 potential 或生物分解的唯一性定理。

### 5.4 精确 section theorem

本节先给 model-conditional 结果。条件于固定 support $S$，定义 c-free baseline

$$
\widetilde b_{\theta,S}(p,l,a)
=b_{\theta,S}^{\rm add}(p,l,a)+s_\theta(p,l,a).
$$

固定 $\widetilde b_{\theta,S},\phi_\theta$，并假设真实 residual
function 属于声明的有界线性族，即存在 $\lVert c_\tau\rVert\le R$ 使所有相关 query
满足
$f_\tau(l)=\widetilde b_{\theta,S}(p,l,a)+\phi_\theta(p,l,a)^\top c_\tau$。若该 well-specification
不成立，下面的 radius 只量化 affine fiber ambiguity，不是对真实 affinity 的完整
证书。

对 exact support，定义

$$
\Phi_S=
\begin{bmatrix}
\phi_\theta(x_p,x_{l_1},a)^\top\\
\vdots\\
\phi_\theta(x_p,x_{l_k},a)^\top
\end{bmatrix},
\qquad
r_S=y_S-\widetilde b_{\theta,S},
\qquad
c^\dagger=\Phi_S^+r_S.
$$

**定理 3（有界线性残差族的精确 query section）。** 若
$r_S\in\operatorname{range}(\Phi_S)$ 且 $\lVert c^\dagger\rVert\le R$，则 query
$q$ 的实数值可行集合是区间，其绝对损失下的精确 minimax center/radius 为

$$
\widehat f_q
=\widetilde b_{\theta,S,q}+\phi_q^\top c^\dagger,
$$

$$
r_q
=\sqrt{R^2-\lVert c^\dagger\rVert^2}
\left\lVert P_{\ker\Phi_S}\phi_q\right\rVert.
$$

特别地，

$$
r_q=0
\iff
\phi_q\in\operatorname{row}(\Phi_S)
$$

（除去可行球已退化为单点的边界情况）。在非退化情形
$\lVert c^\dagger\rVert<R$ 下，完整恢复 $c_\tau$ 需要
$\operatorname{rank}(\Phi_S)=r_{\rm task}$，所以完整 task-state recovery 要求
$r_{\rm task}\le k$；即使 $r_{\rm task}>k$，个别 query 仍可能因 row-space 条件
而可识别。一般的 coefficient 唯一性充要条件是

$$
\ker\Phi_S\cap
\{n:\lVert c^\dagger+n\rVert\le R\}=\{0\}.
$$

**证明。** 所有 exact solutions 写成
$c=c^\dagger+n$，$n\in\ker\Phi_S$。Moore-Penrose 解与 kernel 正交，所以可行
$n$ 构成半径 $\sqrt{R^2-\lVert c^\dagger\rVert^2}$ 的 Euclidean ball。线性泛函
$\phi_q^\top n$ 在该球上的极值为正负该半径乘
$\lVert P_{\ker\Phi_S}\phi_q\rVert$。$\square$

该公式首先针对 $\mathbb R$-值函数族。若只对当前 query $q$ 施加 pointwise
$V$ 约束，可把该 interval 与 $V$ 相交后再取 midpoint/radius。若要求整条函数在
所有相关 queries 上均为 $V$-值，则必须把所有这些 range constraints 加入 coefficient
feasible set 后重新投影；一般不能只截断当前 $q$。另一种足够条件是预先假设所有
$\lVert c\rVert\le R$ 的相关预测均已落在 $V$，此时无需截断。

对 bounded noise 或 censoring，要求 $O_i$ 是非空闭 interval、
$\varepsilon_i\ge0$，且二者语义不重复。定义 compact convex 可行集

$$
\mathcal C_S=
\left\{c:\lVert c\rVert\le R,\;
O_i\cap[\widetilde b_{\theta,S,i}+\phi_i^\top c-\varepsilon_i,
\widetilde b_{\theta,S,i}+\phi_i^\top c+\varepsilon_i]\ne\varnothing,
\;i\le k\right\}.
$$

分别最小化和最大化 $\widetilde b_{\theta,S,q}+\phi_q^\top c$ 得到实值端点；pointwise 或 global $V$ 约束
必须按上一段分别处理。所得 midpoint/radius 是标量 absolute loss 下的精确
model-conditional minimax 决策和 certificate。空集必须 fail closed。

实际 learned representation 的总证书还需

$$
r_q^{\rm total}
\ge r_q^{\rm section}
+\eta_{\rm repr}(q)+\eta_{\rm meta}(q)
+\eta_{\rm transport}(q)+\eta_{\rm obs}(q).
$$

这些项分别控制模型族失配、source 表示估计、条件机制迁移和观测误差。若没有
可证明的上界，只能将 $r_q^{\rm section}$ 称为 conditional radius，并经验检验其
校准，不能声称真实 affinity coverage。

Positive ridge

$$
\widehat c_\lambda
=\arg\min_c\lVert\Phi_Sc-r_S\rVert^2+\lambda\lVert c\rVert^2
$$

可作为稳定、可微的工程近似，但必须与上述 certified section 分开命名和评估。

### 5.5 最小 meta-learning objective

对 SAR route 的 point-valued support，内层必须实际使用第 5.1 节的 centered design，
而不是只在解释中声称 $k-1$。令
$H_k=I_k-k^{-1}\mathbf1\mathbf1^\top$、
$r_{\theta,S}=y_S-\widetilde b_{\theta,S,S}$，并取 $\lambda>0$：

$$
\widehat c_\theta(S)
=\arg\min_{\lVert c\rVert\le R}
\frac1k\left\lVert
H_k(r_{\theta,S}-\Phi_Sc)
\right\rVert_2^2
+\lambda\lVert c\rVert^2,
$$

该目标连续、强凸，所以在闭球上有唯一 minimizer；作为数据和冻结参数的连续
piecewise-smooth 映射，它给出可测 adapter；这里假设 baseline 与 feature maps 均为
Borel measurable。它只拟合去 level 后的 contrasts，故实际
SAR state 的可识别 design 是 $H_k\Phi_S$。Censoring support 若要进入该 adapter，必须
另行预注册 convex profiled interval loss，并重新证明存在、唯一和 centered rank；在此
之前只进入第 5.4 节的可行集证书，不能机械继承 $k-1$ 声明。

对合法 source zero-shot rectangles，明确定义 component-macro auxiliary loss

$$
\mathcal L_{\rm rect}(\theta)
=\frac1{N_C}\sum_{C\in\mathcal C_{\rm src}^{\rm rect}}
\mathbb E_{r\mid C}
\left[
D_r-\mathcal R_a f_{\theta,0}
(p,p';l^-,l^+\mid\varnothing)
\right]^2.
$$

两个 protein 在这里各自取 $c=0$、空 support level statistic；rectangle labels 只作
outer source supervision，绝不进入任一 support。由于
$b_{\theta,\varnothing}^{\rm add}$ 是 additive，
$\mathcal R_a f_{\theta,0}=\mathcal R_a s_{\theta,0}$，但前式明确表明它仍来自同一
full scalar surface。

定义按组件等权的训练准则

$$
\mathcal J(\theta)
=
\frac1{N_C}\sum_{C\in\mathcal C_{\rm src}}
\frac1{|\mathcal T_C|}\sum_{\tau\in\mathcal T_C}
\mathbb E_{S,Q\mid\tau}
\frac1{|Q|}\sum_{q\in Q}
\ell_{\rm reg}(f_{\theta,\widehat c_\theta(S)}(q\mid S),Y_q)
+\lambda_{\rm rect}\mathcal L_{\rm rect}(\theta).
$$

外层深度目标不保证存在全局 argmin；因此把 $\widehat\theta$ 定义为预注册、可测的
训练算法 `Train` 对 $\mathcal J$ 的输出，而不是写成未经证明存在的 argmin。第 6 节
使用该输出实际达到的 empirical regression regret $\gamma_{\rm reg}$。

$\mathcal L_{\rm rect}$ 只在合法 source rectangles 上监督同一个 scalar predictor。
若目标是端到端 DTA，不允许另设一个与 scalar affinity
head 无共享 potential 的独立 edge/rectangle head；否则可能再次出现 F-153 式的
delta head 成功而 DTA 无用。所有 query deltas 必须由

$$
\Delta\widehat f_{p,S}(l^-,l^+)
=\widehat f_{p,S}(l^+)-\widehat f_{p,S}(l^-)
$$

计算。它不使用 deployment labels。
Correct/wrong/shuffled controls属于 admission/evaluation，不应通过把未测 wrong pair
伪造为 non-binder 来训练。

第 6 节的 regression bound 不会因为上述 penalized objective 自动成立。对最终
predictor $\widehat h$ 必须另行计算或上界其 empirical regression regret

$$
\gamma_{\rm reg}
=\widehat R_{\rm src}(\widehat h)
-\inf_{h\in\mathcal H}\widehat R_{\rm src}(h)\ge0.
$$

若采用第二阶段 regression ERM，它可退化为 optimization tolerance；若直接使用
rectangle-penalized 解，$\gamma_{\rm reg}$ 还包含 auxiliary regularizer 引入的
suboptimality，不能省略。

## 6. Cluster generalization 与 cold transport

令固定、pointwise measurable/separable 的 hypothesis class $\mathcal H$ 中，每个
$\bar\ell_h(C)\in[0,M_L]$ 是先在组件内平均 task/episode 后的 loss。定义

$$
R_{\rm src}(h)=\mathbb E_{C\sim\Pi_{\rm src}}\bar\ell_h(C),
\qquad
R_{\rm dep}(h)=\mathbb E_{C\sim\Pi_{\rm dep}}\bar\ell_h(C),
$$

$$
\Delta_\mathcal H
=\sup_{h\in\mathcal H}|R_{\rm dep}(h)-R_{\rm src}(h)|.
$$

对 observed IID source blocks $C_1,\ldots,C_{N_C}$，定义

$$
\widehat R_{\rm src}(h)
=\frac1{N_C}\sum_{i=1}^{N_C}\bar\ell_h(C_i),
$$

以及 expected Rademacher complexity，其中 $\sigma_i$ 是与 $C_i$ 独立的 IID
Rademacher signs：

$$
\mathfrak R_{N_C}(\bar\ell\circ\mathcal H)
=\mathbb E_{C,\sigma}
\sup_{h\in\mathcal H}
\frac1{N_C}\sum_{i=1}^{N_C}\sigma_i\bar\ell_h(C_i).
$$

以下定理把每个 observed component block 及其内部抽样看成一个原子随机对象；
组件内增加 episodes 可降低该 block loss 的测量噪声，但不会把独立样本数从
$N_C$ 增加为 episode 数。若还要把有限 episode 平均逼近组件条件期望，需另加
within-component sampling term。

**定理 4（组件级 ERM 加 transport）。** 对 $N_C$ 个 IID source components，若
$\widehat h$ 的 empirical regression regret 不超过 $\gamma_{\rm reg}$，则以至少
$1-\delta$ 的概率，

$$
R_{\rm dep}(\widehat h)-\inf_{h\in\mathcal H}R_{\rm dep}(h)
\le
4\mathfrak R_{N_C}(\bar\ell\circ\mathcal H)
+2M_L\sqrt{\frac{\log(2/\delta)}{2N_C}}
+\gamma_{\rm reg}
+2\Delta_\mathcal H.
$$

**证明。** 标准 Rademacher uniform deviation 给出 source ERM excess 的前三项；
在 $\widehat h$ 和任意 deployment $\epsilon$-optimal comparator 两端各使用一次
$|R_{\rm dep}-R_{\rm src}|\le\Delta_\mathcal H$，最后令 $\epsilon\downarrow0$。
$\square$

更细地，令 $D$ 是部署可计算、具有 metric $d_{\rm mech}$ 的组件 descriptor，
$\Pi_d^D$ 是 domain $d\in\{\mathrm{src},\mathrm{dep}\}$ 中 $D$ 的分布。对每个 $h$，
固定 conditional risk 的指定版本
$m_h^d(D)=\mathbb E_d[\bar\ell_h(C)\mid D]$。假设两分布有有限一阶矩，并假设
$m_h^{\rm src}$ 存在一个在
$\operatorname{supp}(\Pi_{\rm src}^D)\cup\operatorname{supp}(\Pi_{\rm dep}^D)$
上 everywhere-defined、统一 $L_{\rm mech}$-Lipschitz 的版本；下面始终使用这个指定
版本（或其固定 McShane extension）。定义条件机制变化项

$$
\Delta_{\rm mech}
=\sup_{h\in\mathcal H}
\mathbb E_{D\sim\Pi_{\rm dep}^D}
|m_h^{\rm dep}(D)-m_h^{\rm src}(D)|.
$$

Kantorovich-Rubinstein duality 给出

$$
\Delta_\mathcal H
\le L_{\rm mech}
W_1(\Pi_{\rm src}^D,\Pi_{\rm dep}^D;d_{\rm mech})
+\Delta_{\rm mech}.
$$

Source/deployment 共享同一 conditional risk mechanism 时
$\Delta_{\rm mech}=0$。40% sequence identity split 本身既不证明 $W_1$ 小，也不证明
该项为零。若组件被同分布随机分配，可声明相应 component law 相同；若是
novel-family stress，两个 shift terms 及其估计误差必须显式报告为经验量或
sensitivity parameters。

令 $R_{\rm dep}^*$ 是允许所有合法部署信息的 Bayes risk，并定义 descriptor/model
class 的部署 approximation gap

$$
\varepsilon_{\rm app}^{\rm dep}
=\inf_{h\in\mathcal H}R_{\rm dep}(h)-R_{\rm dep}^*\ge0.
$$

把定理 4 与 transport bound 合并，立即得到：

**推论 4.1（Cold Target descriptor bridge）。** 在定理 4 及其指定 Lipschitz
conditional-risk version 的条件下，以至少 $1-\delta$ 的概率，

$$
R_{\rm dep}(\widehat h)-R_{\rm dep}^*
\le U_{\rm gen}
+2\left[
L_{\rm mech}W_1(\Pi_{\rm src}^D,\Pi_{\rm dep}^D;d_{\rm mech})
+\Delta_{\rm mech}
\right]
+\varepsilon_{\rm app}^{\rm dep},
$$

其中

$$
U_{\rm gen}
=4\mathfrak R_{N_C}(
\bar\ell\circ\mathcal H)
+2M_L\sqrt{\frac{\log(2/\delta)}{2N_C}}
+\gamma_{\rm reg}.
$$

**证明。** 定理 4 先控制 $\widehat h$ 相对 class-optimal predictor 的 excess risk，
再用 transport 上界代替 $\Delta_\mathcal H$，最后加上 class optimum 相对 Bayes risk
的差。$\square$

$\varepsilon_{\rm app}^{\rm dep}$ 正是“蛋白/化学 transformation descriptor 和模型类
是否足够”的不可省略项。固定 $\mathcal H,M_L,\delta$ 并控制 $\gamma_{\rm reg}$ 时，
增加独立 components 通常按标准 rate 减小统计项；它不保证 $U_{\rm gen}$ 逐点单调，
更不能自动消除 $\varepsilon_{\rm app}^{\rm dep}$。
若训练的是 observed rectangles 而结论指向 latent affinity，还要把 measurement
mechanism discrepancy 加入 $\Delta_{\rm mech}$ 或 approximation term。若要进一步从
rectangle predictor 得到 scalar Cold Target query，还必须同时满足：edge zero-cycle
可积性、命题 3 的 quotient coverage、additive gauge 的 baseline/support 锚定，以及
第 5.4 节的 total certificate 条件。

## 7. Regression、ranking 与 abstention 的统一连接

假设在一个 joint event $\mathcal E_{q,q'}$ 上，两条 accepted queries 的 total
certificates 同时给出

$$
|\widehat f(q)-f(q)|\le r(q).
$$

若一对 query 产生非零 pairwise loss，也就是被严格错序或被预测为 tie，则必有

$$
|f(q)-f(q')|\le r(q)+r(q').
$$

因此得到：

**定理 5（selective ranking）。** 令
$A_{\rm pair}$ 表示两条 query 都满足 $r(q)\le r_0$，并假设
$\Pr(A_{\rm pair},\Delta f\ne0)>0$。若真实 margin 满足

$$
\Pr(0<|f(q)-f(q')|\le t\mid A_{\rm pair},\Delta f\ne0)
\le C_{\rm mar}t^\alpha,
$$

并且在同一个条件事件下，joint certificate failure probability
$\Pr(\mathcal E_{q,q'}^c\mid A_{\rm pair},\Delta f\ne0)$ 不超过
$\eta_{\rm pair}$，则 accepted non-tie pair 的 ranking error 至多

$$
\min\{1,C_{\rm mar}(2r_0)^\alpha+\eta_{\rm pair}\}.
$$

**证明。** 在 $A_{\rm pair}\cap\{\Delta f\ne0\}$ 中，严格错序或 predicted tie 都
蕴含 true margin 不超过两个证书半径之和；在 $\mathcal E_{q,q'}$ 上应用 margin
条件，在其补集上用 loss 不超过 1。
$\square$

如果可行值集合恰好是对称 interval
$[\widehat f(q)-r(q),\widehat f(q)+r(q)]$，且 $r$ 是它的 exact Chebyshev radius，
那么对 absolute loss 和 abstention cost $\lambda$，点态 minimax 决策为

$$
A(q)=\mathbf1\{r(q)\le\lambda\},
\qquad
\inf\sup L_\lambda=\min\{r(q),\lambda\}.
$$

平方损失下阈值改为 $r(q)^2\le\lambda$。对一般 conservative certificate，这些
阈值只是安全的充分接受规则，$\min\{r,\lambda\}$ 只是 worst-case upper bound，
不保证等号。这条路线直接连接 regression、ranking、coverage 和 abstention，
不需要把 law band 的 midband readout 错当成 ranking theorem。

## 8. 可选 law-valued decoder 的合法位置

Law layer 只能位于 certified/learned section statistic 之后：

$$
(x_p,a,S,q)
\xrightarrow{\text{interaction + section}}
z_{\rm adm}
\xrightarrow{F}
p
\xrightarrow{B}
\beta
\xrightarrow{K_h}
\text{law class}.
$$

为使本整合稿中的 $K_h$ 成为闭合对象，先固定有限 affinity grid
$a_1<\cdots<a_G$，定义 mesh
$h=\max_{1\le s<G}(a_{s+1}-a_s)$，并固定 $J$ 个 CDF thresholds
$t_1<\cdots<t_J$。令
$\beta=(\ell_1,u_1,\ldots,\ell_J,u_J)$，以及

$$
\Delta^{G-1}
=\left\{\pi\in\mathbb R_+^G:\sum_{s=1}^G\pi_s=1\right\},
$$

$$
K_h(\beta)
=\left\{\pi\in\Delta^{G-1}:
\ell_j\le\sum_{s:a_s\le t_j}\pi_s\le u_j,
\ j=1,\ldots,J\right\}.
$$

定义有效 band 域

$$
\mathcal B_h=\{\beta:K_h(\beta)\ne\varnothing\},
$$

并要求 decoder 的全部输出满足 $B(z)p\in\mathcal B_h$。于是 $K_h(\beta)$ 对每个
有效 $\beta$ 都是非空 compact convex polytope；将 $pi$ 识别为离散 law
$\sum_{s=1}^G\pi_s\delta_{a_s}$。若改用 continuum laws，则必须明确采用
closed-lower/open-upper CDF 约定并单独证明所选 topology 下的闭性。

要把该层整合回 frozen archive，仍需完成两项证明修复，并采用本 working 稿已给出的
两项接口修复：

1. 使用 `OR-CONT` 或 context-disjoint topology 修复 $g_\mu^*$ 连续性；
2. 对具体 scalar selector 另证 regression calibration；ranking 使用第 7 节的 joint
   query/certificate 路线。

接口上应采用上文已定义的 nonempty grid $K_h$（或另证 continuum 版本）以及定理 4
的 component-level bound；二者尚未写回 frozen archive。

在 archive 完成上述合并与证明前，law layer 可以作为结构有效的 uncertainty decoder，
但不是 Cold Target 生物创新的主定理。

## 9. 另一 agent 问题清单的最终判定

| 问题 | 判定 | 精确边界 |
|---|---|---|
| `B(z)` 连续性漏项 | 成立，blocking | 中断 uniform approximation/consistency，不影响点态唯一性 |
| 仓库有两套 operator | 数学对象非等价；权威层级不冲突 | 缺少 admission-to-decoder composition theorem |
| Meta-Section 等于 Chebyshev section | 不成立 | ridge 只是稳定近似 |
| positive ridge 收缩到 population | 不成立 | 当前收缩中心是 simplex uniform mixture |
| `2h` 是不可突破 floor | 不成立 | 只是 proved upper bound 的余项 |
| S-IID 推出 homology-cold | 当前不能 | 需组件层 IID；novel-family 需 transport term |
| F-159 oracle 否定全部 SAR headroom | 不成立 | 只关闭固定 `T-BASIS -> PCA(d=5) -> ridge(lambda=1,k=5)` recipe |
| X2 证明低秩失败 | 不成立 | exact transformation 共享数为零，命题未被检验 |
| X1 证明二阶/低秩生物机制 | 不成立 | 只识别 observed-label double-difference magnitude；latent non-additivity 需 noise-aware null |
| F-152/F-153 证明 SAR-delta transfer | 不成立 | pair construction 读取 outcome；F-153 还在 BindingDB 重拟合，并非跨库模型迁移 |
| F-157 symmetry PASS 证明 interaction | 不成立 | 只修复已选 pairs 的 direction；inclusion/truncation 仍读 outcome，未通过 G0 |
| F-154/F-155 否定 scalar/unordered 表示 | 不成立 | 只关闭当前 1D neighbor-mean 与八维 summary recipes |
| F-158 否定 pair-score potentials | 不成立 | 只关闭固定 descriptors、linear ridge 和当前 supervision 的 U1 arm |
| 标量或 `p=1` 必然不可辨 | 不成立 | 一维也可辨；历史 `p=1` diagnostic 已降级 |
| Meta-learning 已贡献 Cold biology | 尚未识别 | 已消费 splits 上有 adaptation/calibration signal，未识别 cluster-robust specificity |
| Full CSMO 继承 frozen bound | 不成立 | 只保证 simplex/band 结构，full class 需另证 |
| Law calibration 推出 RMSE/CI/Spearman | 不成立 | 需 selector-risk 与 joint-ranking theorem |

## 10. 可证伪的研究 Gate

“transferable descriptor + rectangle quotient + certified section”是当前资源优先检验的
候选路线，不是排除其他 interaction/meta-learning 路线的数学结论。
每个 Gate 使用独立组件 bootstrap 或组件级 concentration，不使用行级伪重复。

### 共同预注册契约

除 G0 外，每个 Gate 在查看 confirmation outcomes 前必须冻结七项：模型类
$\mathcal M_G$、primary estimand（或预声明的 co-primary 向量）$\Theta_G$、实际独立单位 $C$、最小有意义效应
$\tau_G$、错误率 $\alpha_G$、单侧 confidence/lower-bound 构造、以及 `PASS/FAIL` 后
允许采取的动作。若同一 Gate 有多个必要 contrasts，需预先声明 intersection rule 和
multiplicity correction。默认成功条件是

$$
\operatorname{LCB}_{1-\alpha_G}(\Theta_G)>\tau_G.
$$

数值阈值必须由 assay noise、positive control 或临床/筛选容差在运行前确定，不能从
当前 development result 反推。secondary metrics 不能挽救 primary failure。除预注册的
确定性有效性或 coverage Gate 外，效果 Gate 失败只表示在冻结模型、数据域、程序、
样本量和阈值下未识别预注册效应；可按 stop-tree 关闭该路线，但不构成函数类或真实
效应不存在的不可能性结论。

### G0. 数据与 estimand

1. 冻结 target × assay task 定义、affinity 转换和 censoring 规则；
2. 冻结 homology、document/time/publication、assay/panel、shared-ligand 和 scaffold
   依赖边，并逐类声明它用于 hard split、block bootstrap 还是 overlap stratification；
3. 报告 source/validation/deployment 的组件数和最大组件占比；
4. 注册为 hard-separation 的边跨 split 违规数必须为零；
5. preprocessing、descriptor/metric learning 和 threshold selection 只读 source/calibration；
6. confirmation/deployment labels 在正式开封前不得被 loader materialize；
7. pair orientation function 不得读取 affinity outcome；若无 canonical direction，必须
   同时物化 forward/reverse pair，并验证 label、descriptor 和 prediction 均翻号；
8. pair identity、inclusion、enumeration 和 subsampling/truncation 必须在读取 label 前
   冻结，并用 label-redacted construction test 验证；
9. support/query unit 严格不交，query labels 不进入 representation 或调参；
10. 在读取 confirmation outcomes 前冻结最小独立 protein components、每组件最小有效
    rectangles、有效 transformations、最大 component/panel 权重和 multiway cluster 数；
    不足即 data-supply fail closed。小 cluster 数只能使用预注册的 exact/randomization
    procedure，不能事后套用普通渐近区间。

G0 是确定性有效性 Gate；任一项失败即
`DATA_ESTIMAND_GATE_FAIL_CLOSED`，后续效果检验无权覆盖它。

### G1. Descriptor coverage

1. transformation descriptor 必须在 train/development 有共享数值支持；
2. 报告最近邻距离、convex-hull/Gram leverage 和有效 rank；
3. planted descriptor positive control 必须可恢复；
4. exact ID overlap 只作诊断，不再作为迁移表示；
5. 若声称 chemical transformation descriptor transfer，confirmation 必须 hold out
   exact transformation ID，同时保持 descriptor coverage；若 exact edge 可复用，
   最多主张 protein-cold transfer。

运行前冻结 coverage、leverage、effective-rank 和 positive-control recovery 的可接受
阈值；五项全部通过才 admission。失败只关闭该 descriptor route。

### G2. Cold protein-conditioned quotient admission

只使用同一可比 assay/context $a$ 的合法有向 transformation
$e=(l^-,l^+)$ 与两个 protein，目标为

$$
D_r=
[Y_{p,a}(l^+)-Y_{p,a}(l^-)]
-[Y_{p',a}(l^+)-Y_{p',a}(l^-)].
$$

若最终主张 scalar DTA，模型必须来自第 5.3 节的 shared scalar-potential family；G2
固定 $c=0$，故

$$
\widehat D_r
=\mathcal R_a f_{\theta,0}(p,p';l^-,l^+\mid\varnothing)
=\mathcal R_a s_{\theta,0}(p,p';l^-,l^+),
$$

否则只能主张 edge/quotient prediction。Primary risk 固定为

$$
R_{\rm rect}(h)
=\mathbb E_C\mathbb E_{r\mid C}(D_r-\widehat D_r)^2.
$$

rectangle 是同时接触两个 protein nodes、两个 ligand edges 和四个 measurement cells 的
dyadic/graph observation。依赖处理必须预注册以下之一：按共享 cell、protein、document、
transformation 和 panel 的传递闭包形成完整 dependency components；对 protein nodes
使用 dyadic/graph bootstrap 并嵌套 transformation/panel blocks；或在有效独立单位很少时
使用 exact/randomization inference。普通的行级 bootstrap 或“每条 rectangle 每一维只
给一个 cluster ID”的 multiway cluster 不能自动处理这种 overlap。每个 confirmation
rectangle 的 $p,p'$ 都必须来自 source-disjoint
homology components，训练不得读取二者的任何 affinity labels；二者彼此是否属于同一
依赖闭包必须进入 $C$。若 $p'$ 是 source protein，只能称 source-to-cold contrast，
不能通过本 cold-by-cold Gate。

稀疏 observed design 中，protein-pair、transformation、panel、missingness 和 noise
variance 仍可能预测 rectangle。简单加性且双反对称的 nuisance class 只含零函数，
不能作为有效 comparator。故在 source/calibration 上冻结 design-nuisance class

$$
h_{\rm nuis}(p,p',e,a)
=q_P(p,p',a)^\top W_{\rm nuis}q_E(e,a),
$$

其中 $q_P(p,p',a)=-q_P(p',p,a)$、$q_E(e^{\rm rev},a)=-q_E(e,a)$ 只包含在
outcome-redacted construction 中冻结的 degree、availability、panel、missingness、
selection-probability 和 noise-level 等设计特征；它们不得包含主模型的 protein 或
chemical transformation descriptors $\rho,\psi$。该 tensor class 可以捕获双反对称的
观测设计 interaction，但不能冒充待检验的 biological descriptor interaction。仅用
calibration 选定 confirmation 上唯一 comparator，其候选集合至少包含 zero quotient。
所有 comparator 使用与主模型相同的 canonical orientation、完整 forward/reverse
augmentation 和双重反对称约束，不能自行恢复 outcome-order shortcut。
冻结三个
co-primary contrasts：

$$
\Theta_{\rm rect}=R_{\rm best\ design\ nuisance}-R_{\rm correct},
$$

$$
\Theta_{\rm partner}=R_{\rm matched\ wrong}-R_{\rm correct},
\qquad
\Theta_{\rm shuffle}=R_{\rm shuffled}-R_{\rm correct}.
$$

三者的 multiplicity-adjusted LCB 均须超过阈值。Matched wrong protein 必须
in-support、来自另一同源组件，并匹配 assay、descriptor leverage 和 transformation
coverage。另有三个必要但不计作效果证据的结构审计：完整 forward/reverse
augmentation、双重反对称、edge zero-cycle 可积性。Replicate/noise-aware null 或重复
测量去卷积仍为 latent-interaction 主张的必要条件。

G2 PASS 最多授权：“在冻结模型族和数据域内，source-disjoint protein components
上的 observed-affinity transformation response 含有可预测的 protein-conditioned
quotient。”它不授权因果耦合、接触机制或普适生物 SAR。任一 co-primary 或结构
审计失败，都关闭该冻结 descriptor/potential 路线。

### G3a. Few-shot SAR-delta admission

只有 G2 的 shared-potential 分支通过并冻结
$\Delta\widehat f_{p,0}(e)
=f_{\theta,0}(p,l^+,a\mid S)-f_{\theta,0}(p,l^-,a\mid S)$，才可进入 G3a；quotient-only
G2 PASS 没有单蛋白 delta baseline，不授权继续。G3a 专门检验 support 是否改变
held-out query transformation，而不是一般 absolute affinity。对 $e_q=(q^-,q^+)$，定义

$$
R^\Delta_{u,v}
=\mathbb E_C\mathbb E_{S,e_q}
\left[
\Delta Y_p(e_q)-\Delta\widehat f_{u,v}(e_q\mid S)
\right]^2,
$$

其中 $u\in\{\mathrm{correct},\mathrm{matched\ wrong}\}$ 只选择 adaptation state 的
support protein，$v\in\{\mathrm{real},\mathrm{permuted}\}$ 表示 support
transformation-label binding。所有 arm 固定正确 query protein 的 zero-shot potential：

$$
\Delta\widehat f_{u,v}(e_q\mid S)
=\Delta\widehat f_{p,0}(e_q)+A_{u,v}(S,e_q),
$$

$$
A_{u,v}(S,e_q)
=\bigl[\phi_\theta(p,q^+,a)-\phi_\theta(p,q^-,a)\bigr]^\top
\widehat c_u(S_v).
$$

因此 wrong-protein control 只腐坏 support-derived state，不同时替换 frozen G2 baseline
或正确 query basis。$\Delta\widehat f_{p,0}$ 与
$b_{\theta,S}^{\rm add}$ 始终用同一个 correct real support 的 level statistic 计算；
该统计量对 binding permutation 不变。冻结

$$
\Theta_{\rm adapt}
=R^\Delta_{\rm frozen\ G2}-R^\Delta_{\rm correct,real},
$$

$$
\Theta_{\rm binding}
=R^\Delta_{\rm correct,permuted}-R^\Delta_{\rm correct,real},
$$

$$
\Theta_{\rm spec}^{\Delta}
=R^\Delta_{\rm wrong,real}-R^\Delta_{\rm correct,real},
$$

$$
\Theta_{\rm SAR-DID}
=\bigl(R^\Delta_{\rm wrong,real}-R^\Delta_{\rm correct,real}\bigr)
-\bigl(R^\Delta_{\rm wrong,permuted}-R^\Delta_{\rm correct,permuted}\bigr).
$$

另报告 prediction-update interaction

$$
\Xi_A=
(A_{\rm correct,real}-A_{\rm correct,permuted})
-(A_{\rm wrong,real}-A_{\rm wrong,permuted}),
$$

及其 component-level magnitude、与 query residual 的方向一致性和对上述 risk DID 的
贡献；仅有 $\Xi_A\ne0$ 而没有 risk improvement 不构成 PASS。

若论文主张 transformation-conditioned 而非一般 chemical-neighborhood adaptation，
再加入必要 co-primary

$$
\Theta_{\rm trans}
=R^\Delta_{\rm molecule\ kernel}
-R^\Delta_{\rm transformation\ adapter}.
$$

所声明的全部 co-primary 的 multiplicity-adjusted LCB 均须超过阈值。
Transformation-binding permutation
在 support compounds 间置换 labels，保持 compound 集合、label multiset 和 protein
不变，再按冻结的 canonical directions 重算 deltas，从而只破坏
$\psi(e_s)\leftrightarrow\Delta Y(e_s)$ 配对。Support/query 原始 cells 严格不交，
query exact edges 不得与 support 重合；报告 query descriptor leverage 和 centered
row-space coverage。Support selection 按 transformation Gram rank 预注册或分层，并与
容量匹配的 molecule-kernel adapter 比较。SAR-specific state 使用
$\Phi_S^\Delta$，故可识别维数至多 $k-1$。

Foreign-support 只作复合 stress test。若不检验 $\Theta_{\rm trans}$，G3a PASS 只能授权
source-disjoint confirmation components 上 observed-assay risk 的 few-shot
protein-conditioned chemical-neighborhood adaptation；检验并通过后才可称
SAR-transformation adaptation。两者都不授权 scalar DTA 或 homology-cold
superpopulation 泛化。

### G3b. Shared-potential scalar DTA bridge

G3b 要求 G2、G3a 和 scalar affinity 共用完全相同的 potential：

$$
\widehat f_{p,S}(l)
=b_{\theta,S}^{\rm add}(p,l,a)
+s_{\theta,\widehat c_p(S)}(p,l,a).
$$

至少设置三个 arms：`no SAR`、`frozen G2 potential without support adaptation`、以及
`full G2 + G3a adaptation`。采用冻结的分阶段训练：先在 source/calibration 上得到同一个
$b_{\theta,S}^{\rm add}$ checkpoint，并在全部 arms 中禁止更新；再训练一次 G2
$s_{\theta,0}$ checkpoint，供 `frozen G2` 与 `full` 共同载入并禁止更新；最后只在
`full` 中训练 $\phi_\theta$ 与 centered adapter。每个新增 channel 必须在其训练阶段
真实优化，禁止只在 inference 时置零。初始化、优化预算和 channel capacity 的规则
预先冻结。所有 arms 使用同一 additive checkpoint 并在相同 support 输入上计算，所以
support-level calibration 或 ligand-only baseline 不能贡献 arm 间差异。Query deltas
直接由相应 full scalar surface 相减。冻结两个
scalar co-primary：

$$
\Theta_{\rm E2E}
=R_{\rm no\ SAR}-R_{\rm full},
\qquad
\Theta_{\rm scalar\ adapt}
=R_{\rm frozen\ G2}-R_{\rm full}.
$$

二者的 multiplicity-adjusted LCB 都必须超过阈值；第二项专门识别 support/meta-learning
超越 zero-shot potential 的 scalar 收益。另冻结 within-target ranking noninferiority
margin，以及 SAR-channel-specific 的 correct-vs-matched-wrong protein 和
real-vs-transformation-permuted support contrasts：这些 controls 全程固定正确 query
protein、正确 real support 所产生的 $b_{\theta,S}^{\rm add}$，只替换
$s_{\theta,c}$ 中的 protein descriptor 或 support state。还需冻结其阈值、scalar
selector、intersection rule 和最小 confirmation components。命题 3 的 gauge 锚定与
quotient coverage 也必须可计算。只有全部通过，
query-specific potential correction 及其 candidate/conditional radius 才可进入 G4a；
只有 G4a PASS 后，calibrated total radius 才有资格进入实际 acceptance/readout 的 `z`。

进入 V1/CSMO 后还必须用最终 scalar/law readout 重放同一组 contrasts；handoff 后任一
contrast 消失即 fail closed。没有 coverage 的 ridge adapter 只能称 adaptation，不能称
certified section；G2/G3a 单独 PASS 只能说明中间对象已识别。G3b 的授权仍限于
source-disjoint confirmation components 的 observed-assay risk；G5 通过前不得外推
homology-cold superpopulation 或 latent affinity。

### G4a. Certificate calibration

在独立 calibration components 上冻结 $r_0$，以及最小 query acceptance coverage
$c_q$、目标 certificate coverage $1-\alpha_{\mathrm{cert}}$ 和最小有效 confirmation
component 数。令 $A_q=\{r^{\rm total}(q)\le r_0\}$。Confirmation PASS 要求
simultaneous component-level confidence bounds 同时满足

$$
\operatorname{LCB}\Pr(A_q)\ge c_q,
\qquad
\operatorname{LCB}\Pr(|\widehat f(q)-f(q)|\le r^{\rm total}(q)\mid A_q)
\ge1-\alpha_{\rm cert}.
$$

Synthetic/latent truth 可检验 conditional section radius；若使用 noisy query
observation，必须把 representation、transport、measurement-noise/censoring 项加入
total certificate。若只有 $r_q^{\rm section}$ 而没有其余误差上界，PASS 只能称为
finite-cohort empirical calibration，不能称现实 affinity certificate validity。

### G4b. Selective ranking

在 G4a 通过后，用同一个已冻结 $r_0$。预注册最小 accepted-pair coverage
$c_{\mathrm{pair}}$、最大 pairwise ranking error $e_{\mathrm{pair}}$ 和最小有效
confirmation component 数。令 $A_{\rm pair}$ 表示两条 query 都在 $A_q$ 中。PASS
要求 $\Pr(A_{\rm pair}\mid\Delta f\ne0)$ 的 LCB 不低于
$c_{\mathrm{pair}}$，同时在 $A_{\rm pair}\cap\{\Delta f\ne0\}$ 上 ranking error 的
UCB 不高于 $e_{\mathrm{pair}}$；二者之一
失败即关闭当前 selective-ranking claim。若报告整条 curve，
需 simultaneous confidence band；否则只检验冻结的单点。CI、Spearman、NDCG 各自是
单独 estimand，不能从 marginal law 自动继承。

### G5. Cold transport

1. 区分组件随机分配与 novel-family stress；
2. 在声明的 mechanistic descriptor metric 下报告 source/deployment 距离；
3. 做 transport-distance 分层风险曲线；
4. 任何“cold generalization theorem”都必须显示 $N_C$ 与
   $\Delta_\mathcal H$，不能用 episode 数替代。

运行前冻结 $\Delta_{\rm mech}$ sensitivity set $\mathcal S$、Lipschitz 上界、各估计量的
置信规则和相对 Bayes risk 的 deployment excess-risk budget $B_{\rm dep}$。记
$U_{\rm gen}$ 为定理 4 中 source generalization 加 empirical regret 的上界，
$U_{W_1}$ 为 descriptor transport 距离的上置信界，$U_{\rm app}$ 为
$\varepsilon_{\rm app}^{\rm dep}$ 的预注册上界或 sensitivity value。Sensitivity-conditional
PASS 定义为

$$
U_{\rm gen}
+2\left(L_{\rm mech}U_{W_1}
+\sup_{\delta\in\mathcal S}\delta\right)
+U_{\rm app}
\le B_{\rm dep}.
$$

若无法由数据识别 $\Delta_{\rm mech}$，必须报告使不等式刚好失效的 tipping point。
Descriptor distance 不能验证 conditional mechanism invariance。预算超限时，结论是该
deployment domain 不在当前 sensitivity-conditional scope 内；它不是机制不变性的检验。

## 11. 论文可主张与不可主张

在 G0-G5 完成前，可主张：

1. 有限 support 下 query-specific identifiability radius 的一般理论；
2. component-level Cold Target risk 与 transport decomposition；
3. 本次 governed BindingDB development split 的 component bootstrap 显示
   observed-label rectangle double-difference magnitude 的下界为正；其中 latent
   interaction 与 measurement/assay noise 尚未分离，也尚未跨 split 复现；
4. 在已消费 main-v0 splits 上，当前 v0 adapter 的收益平均主要表现为 support
   calibration；cluster-robust protein-specific SAR 未识别；
5. descriptor-based crossed interaction 加 bounded adaptation 是由失败证据约束出的
   可预注册候选模型类。
6. F-152/F-153 在 outcome-oriented pair protocol 上出现过表观预测性，但该证据已被
   F-156/F-157 降级；完整 G0-compliant symmetric SAR-delta Gate 尚未运行。

不可主张：

1. 已证明在 homology-cold targets 上一致；
2. X1 已识别低秩或生物机制；
3. X2/F-159 构成一般不可能性定理；
4. ridge section 是旧 canonical minimax operator；
5. band law theorem 自动保证 scalar RMSE 或 ranking；
6. 更大的 encoder 会自然解决信息不足；
7. F-154/F-155 已否定所有 scalar potential 或 permutation-invariant set model；
8. F-158 已否定所有 UniPert/pair-score-difference model。

上述主张是 scope statement，不是充分的论文新颖性声明。这里的 quotient、线性
section、Rademacher bound 和 transport inequality 都是已知数学工具的 DTA 特化。
真正尚待完成的理论桥是：在 component-level source sampling 下，把 learned
$\phi_\theta$ 的表示估计/失配误差传入 $r_q^{\rm total}$，并由同一有限样本结果同时
控制 deployment regression、selective ranking 和 abstention。没有这条桥，模型只能
称为受理论约束的候选实现，不能称为完整 Cold Target 一致性理论。

## 12. 与相关数学工作的关系

共享表示文献在共享真实表示、任务多样性和样本条件下给出新任务收益或表示恢复
结果；这些条件不能直接外推到当前非线性 DTA。
参见 [Tripuraneni, Jin, Jordan (ICML 2021)](https://proceedings.mlr.press/v139/tripuraneni21a.html)
和 [Maurer, Pontil, Romera-Paredes (JMLR 2016)](https://jmlr.org/beta/papers/v17/15-242.html)。

Descriptor-based rectangle model 属于带 row/column side information 的 inductive
low-rank estimation 范式；side information 覆盖真实子空间和观测设计充分是定理
前提，不是训练自动产生的事实。参见
[Jain and Dhillon, Provable Inductive Matrix Completion](https://www.microsoft.com/en-us/research/publication/provable-inductive-matrix-completion/)。

生物信息上，sequence-derived pair interaction fingerprints 已有经验先例，例如
[PSICHIC](https://www.nature.com/articles/s42256-024-00847-1)；这只说明表示路线值得检验，
不证明其现成 embedding 对 MetaSieve 的 support-conditioned interaction 充分。
AdaMBind 也明确把每个 target 构造成 support/query task 并在 unseen targets 上适配，
但其经验结果不能替代本文的组件独立性、可辨识半径和 transport 条件：
[AdaMBind](https://www.nature.com/articles/s41467-026-70554-5)。

## 13. 最终停止判据

当前数学研究可以在以下结论上停止继续抽象扩张：

1. 旧 identifiability layer 保留；
2. frozen law layer 降为可选 decoder，并登记第 2.1、2.6 节修复；
3. 项目当前优先候选对象是 descriptor-based interaction quotient 加 bounded section；
4. 所有统计结论改用 component-macro 风险；
5. Cold Target 泛化统一写成 source generalization 加 transport penalty；
6. regression、ranking 和 abstention 通过 query certificate 连接；
7. 在 G1 或 G2 通过前，不再扩大模型容量；只有 G3a 与 G3b 均通过后，才可进入
   V1/CSMO preservation audit；
8. 若 G2 在合法 descriptor coverage 和 positive control 下仍失败，则结论限定为：
   在冻结的数据域、descriptor、模型、训练程序和效应阈值下，未识别预注册效应，
   关闭该具体路线，而不是继续重命名同类特征。

这构成当前证据支持的紧凑候选框架：命题 1-3、推论 3.1、定理 1-5 和推论 4.1 在
所列假设下可证，模型类可训练且 Gates 可预注册；具体 interaction representation
的充分性、可学性、总证书和 Cold deployment scope 仍待 G1-G5 验证。

## 14. 仓库证据索引与复核方法

本稿由三个并行只读审计角色交叉复核：axiom/theorem、Cold Target/meta-learning
formulation、implementation/experiment/external claims。任何角色提出的反例只有在主稿
重新推导并通过末轮交叉检查后才纳入。关键仓库证据如下：

| 结论 | 最小证据位置 |
|---|---|
| Frozen theory 是唯一权威，handoff 仅为支持材料 | `theory/README.md:5-14`；`theory/PROJECT_HANDOFFS/README.md:3-6` |
| `S-CONT` 只控制固定 $\beta$，且没有 shift theorem | `theory/FINAL_FROZEN_THEORY/chapters/01_CORE_FOUNDATIONS/01_FOUNDATIONS.md:63-69` |
| 连续性证明直接用 `S-CONT` 控制变化后的 composite objective | `theory/FINAL_FROZEN_THEORY/chapters/01_CORE_FOUNDATIONS/03_STRONG_CONVEXITY_AND_REGULARITY.md:23-36` |
| Approximation/consistency 依赖该连续性链 | `theory/FINAL_FROZEN_THEORY/chapters/02_LEARNING_THEORY/05_APPROXIMATION_THEORY.md:27-43`；`06_CALIBRATION_AND_GENERALIZATION.md:75-95` |
| 实现中的 context 和 $B(z)$ 是 bucketized table lookup | `model/meta_operator.py:53-68,630-635` |
| Full CSMO 不继承 frozen bounds | `model/meta_operator.py:19-26` |
| V1 adapter 是 positive dual ridge | `model/metasieve_v1.py:93-116` |
| Scalar law readout 是 frozen operator 外的工程选择 | `model/mathematical.py:569-577` |
| 旧 handoff canonical object 是 section center/radius | `theory/PROJECT_HANDOFFS/FINAL_THEORY_TO_MODEL_HANDOFF.md:38-42,186-193` |
| F-121 support adaptation 与 cluster-specificity failure | `history.md:5174-5205` |
| F-123 calibration-dominant 的 consumed-split 诊断 | `history.md:5317-5363` |
| F-152/F-153 的原始 pair protocol 与表观 PASS | `history.md:7942-8122` |
| Outcome-dependent endpoint ordering 的实现证据 | `research/source_affinity/train_chembl_assay_sardelta.py:118-127`；`research/crossed_interaction/train_bindingdb_sardelta_attribution.py:118-146` |
| F-153 在 BindingDB 内重新拟合而非跨库权重迁移 | `research/crossed_interaction/train_bindingdb_sardelta_cq_bridge.py:154-167` |
| F-154/F-155 recipe-level lift failure | `history.md:8124-8294` |
| F-156 attribution、F-157 symmetry 与 F-158 U1 failure | `history.md:8399-8644` |
| F-159 固定 target-ID oracle arm 及其失败边界 | `history.md:8647-8728` |
| X1 rectangle census/result 与原始作用域 | `history.md:8731-8805` |
| X2 当前 split 的 exact-ID intersection 为零 | `history.md:8808-8847` |

末轮机器检查包括：Markdown/LaTeX display delimiter 配对、禁止控制字符扫描、
`git diff --check`、过强/陈旧措辞扫描，以及 frozen archive 工作区零改动检查。
