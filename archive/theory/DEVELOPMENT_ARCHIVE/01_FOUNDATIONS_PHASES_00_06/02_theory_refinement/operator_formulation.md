# Operator Formulation

> **Status:** theory-refinement layer, 2026-08-02. Source of record: `../00_raw_outputs/identifiability_treatise_raw.md`; formal base: `theorem_formalization.md` (F-numbers); conceptual layer: `core_principle.md` (CP-numbers). All propositions OP-1 … OP-11 were adversarially refereed; the record below is the corrected version (notably: the nonlinearity witness required $k\ge3$, and several scope clauses were added).

**Question.** Can the theory's result be expressed as an abstract operator $A(S)$, or another mathematical mapping? The answer is *yes, in a stratified and provably forced sense* — but the operator cannot have the shape $A(S)$ with $S$ the sample alone. This file determines the exact shape, the canonical choice, its uniqueness, and its structural properties.

Setting as in F-Part 0. $T=T_{D,x}\subseteq\mathbb R^{k+1}$; section $S_\varepsilon(\tilde y)$; modulus $\omega=\omega_{x,D}$; $R^*=\tfrac12\omega(2\varepsilon)$.

---

## 1. Three candidate readings of "$A(S)$"

1. **Unary in the sample:** a fixed map $\mathbb R^k\to\mathbb R$ (or of the labeled sample $\{(x_i,\tilde y_i)\}$ and query) independent of the family. — *Impossible* (§2).
2. **Binary:** $A(T,\tilde y)$ — family-level argument and member-level argument. — *Exists, canonical, unique in a precise sense* (§3).
3. **Set-valued:** the correspondence $(T,\tilde y)\mapsto S_\varepsilon(\tilde y)$. — *Forced and assumption-free; the binary operator is its canonical selection* (§3–§4).

---

## 2. No unary operator exists

**OP-10 (Proposition — necessity of the family argument). [refereed: confirmed]**
There is no map $\Phi:\mathbb R^k\to\mathbb R$, chosen independently of the family, with finite worst-case error at any $x\notin D$ uniformly over families: already against the single family $\mathcal F=\mathbb R^{\mathcal X}$, $\omega_{x,D}(0)=\infty$ at every $x\notin D$, so any fixed $\Phi$ has infinite worst-case error (F4). Hence any operator formulation is **at least binary**: the family dependence is not a nuisance parameter but a mathematical necessity.

By CP-1/CP-2, the family argument can be taken to be exactly $T$ — no more (sufficiency) and no less (minimality, with the $\varepsilon=0$ clause; at any fixed $\varepsilon>0$ alone, $T$ is determined only up to dense-in-section modifications, e.g. $\mathbb Q\times\{0\}$ vs $\mathbb R\times\{0\}$).

---

## 3. The canonical binary operator

**OP-1 (Proposition — domain of the family argument). [refereed: confirmed]**
Every nonempty $T\subseteq\mathbb R^{k+1}$ is the trace set of some family: for $(u,t)\in T$ let $f_{(u,t)}(x_i)=u_i$, $f_{(u,t)}(x)=t$, $0$ elsewhere (uses only that $x_1,\dots,x_k,x$ are $k+1$ distinct points, guaranteed by $x\notin D$). Hence the family argument ranges over **all** of $2^{\mathbb R^{k+1}}\setminus\{\emptyset\}$; no regularity of $T$ may be presumed by the operator.

**Definition (the central operator).** For nonempty bounded $S\subseteq\mathbb R$ let $\operatorname{cen}(S)=\tfrac12(\inf S+\sup S)$. Define, for realizable $\tilde y$ (i.e. $S_\varepsilon(\tilde y)\ne\emptyset$) with bounded section,
$$A(T,\tilde y)\;=\;\operatorname{cen}\big(S_\varepsilon(\tilde y)\big),$$
extended arbitrarily off the realizable set (immaterial: no risk term lives there; under the observation model realizability is automatic, and for misspecified data the convention of record is projection onto the realizable set with the doubling bound $\tfrac12\omega(2\varepsilon+2\eta)$, treatise §10.4).

**OP-Thm A (optimality).** $A(T,\cdot)$ attains the minimax value $R^*$ for every family realizing $T$ (F1). This holds with no assumptions beyond SA1–SA4.

**OP-2 (Proposition — the window of optimal rules). [refereed: confirmed with caveats]**
Assume $R^*<\infty$ (then optimal rules exist — the central rule is one, by F1). Every minimax-optimal $\Phi$ satisfies, at **every** realizable $\tilde y$:
$$\Phi(\tilde y)\ \in\ \big[\sup S_\varepsilon(\tilde y)-R^*,\ \ \inf S_\varepsilon(\tilde y)+R^*\big],$$
a window of width $\omega(2\varepsilon)-\operatorname{diam}S_\varepsilon(\tilde y)\ge0$; conversely any rule selecting inside the window everywhere is optimal. (The per-data constraint is genuinely forced: each realizable $\tilde y$ contributes its own conditional supremum as a term of the worst-case risk.) At data attaining the supremal diameter $\operatorname{diam}S=\omega(2\varepsilon)$ — *when such data exist; F1 guarantees only a supremum* — the window degenerates and all optimal rules coincide with the midpoint.

**OP-3 (Theorem — universal property: unique conditional minimaxity). [refereed: confirmed with scope]**
Assume $\omega(2\varepsilon)<\infty$. For each realizable $\tilde y$, $a\mapsto\sup_{t\in S}|a-t|=\max(\sup S-a,\ a-\inf S)$ is strictly decreasing left of $\operatorname{cen}(S)$ and strictly increasing right of it, so the midpoint is the **unique** minimizer, with value $\tfrac12\operatorname{diam}S$. Hence $A(T,\cdot)$ is the unique rule that is minimax **conditionally at every data point** (on the realizable set; values elsewhere immaterial), and it pointwise-dominates the guarantee of every other minimax-optimal rule, strictly wherever they differ. This is the precise sense in which the operator is *canonical rather than chosen*.

---

## 4. Structure theorems

**OP-decomposition.** $A=\operatorname{cen}\circ\,\mathrm{sec}$, where $\mathrm{sec}:(T,\tilde y)\mapsto S_\varepsilon(\tilde y)$ is the (assumption-free) correspondence and $\operatorname{cen}$ the selection. The regularity of the two factors separates cleanly:

**OP-5 (Proposition — the selection is 1-Lipschitz). [refereed: confirmed, constant sharp]**
On nonempty bounded subsets of $\mathbb R$ with the Hausdorff distance, $S\mapsto\sup S$ and $S\mapsto\inf S$ are $1$-Lipschitz, hence $\operatorname{cen}$ is $1$-Lipschitz; the constant $1$ is attained (translations). ($d_H$ is only a pseudometric on non-closed sets; harmless — $\sup,\inf,\operatorname{cen}$ are closure-invariant.)

**OP-6 (Proposition — cross-data spread of the operator). [refereed: confirmed; constant $2\varepsilon+\delta$, not $2\varepsilon+2\delta$]**
For realizable $\tilde y,\tilde y'$ with $\delta=\|\tilde y-\tilde y'\|_\infty$: every $a\in S_\varepsilon(\tilde y)$, $a'\in S_\varepsilon(\tilde y')$ satisfy $|a-a'|\le\omega(2\varepsilon+\delta)$ (one triangle through the two data vectors), hence $|A(T,\tilde y)-A(T,\tilde y')|\le\omega(2\varepsilon+\delta)$. This is a spread bound, not a continuity claim (it tends to $\omega(2\varepsilon^+)$, not $0$, as $\delta\to0$). At $\varepsilon=0$ under identifiability it is exact and optimal: $A(T,\cdot)=\Psi$, and $\omega$ *is* the modulus of continuity of $\Psi$ on $\mathcal F|_D$ — no rule has a better modulus (F1/F3).

**OP-7 (Proposition — equivariance). [refereed: confirmed, including reflections]**
(i) *Permutations:* $A$ commutes with permutations of the design coordinates (the $\ell_\infty$ ball is permutation-invariant); consequently $A(T,\cdot)$ is well-defined on the sample as a finite **set** $\{(x_i,\tilde y_i)\}$ — the natural formulation is a map from finite labeled samples to values at queries.
(ii) *Affine action on values:* for $\alpha\ne0,\beta$, transforming members $f\mapsto\alpha f+\beta$ maps $T\mapsto\alpha T+\beta$ (all value coordinates), $\tilde y\mapsto\alpha\tilde y+\beta\mathbf1$, $\varepsilon\mapsto|\alpha|\varepsilon$, and then $A(\alpha T+\beta,\alpha\tilde y+\beta\mathbf1)=\alpha A(T,\tilde y)+\beta$ — valid also for $\alpha<0$, where $\inf$ and $\sup$ swap but the midpoint is reflection-equivariant.
(iii) *Gauge invariance:* $A$ depends on the family only through $T$ (CP-1); reparametrization acts trivially.

**OP-8 (Proposition — monotonicity of the guarantee, not of the estimate). [refereed: confirmed]**
$T'\subseteq T$ (subfamily) or design refinement (extra observations) shrinks every section, so the conditional guarantee $\tfrac12\operatorname{diam}S$ is non-increasing. The estimates themselves are not ordered: nested sections $[0,10]\supset[9,10]\supset[9,9.2]$ have midpoints $5,\ 9.5,\ 9.1$.

**OP-9 (Theorem — factorization / minimal sufficiency). [refereed: confirmed with the $\varepsilon=0$ clause]**
The inference problem at $(D,x)$ — section correspondence, risk of every estimator, optimal-rule set — is a function of $T$ alone; and $T$ is recovered from the exact-data section map as its graph. Therefore the assignment (family) $\mapsto$ (problem, all $\varepsilon\ge0$, **including $\varepsilon=0$**) has exactly the fibers of (family) $\mapsto T$: the operator's family argument is minimally sufficient. With $\varepsilon$ ranging over $(0,\infty)$ only, minimality fails (dense-in-section modifications are invisible; §2).

---

## 5. Linearity analysis: what is linear, and what is not

For linear families, F7 provides a minimax-optimal **linear** rule $\tilde y\mapsto\hat w^\top\tilde y$. The relation between that rule and the canonical operator is delicate and was refereed to the following corrected statements:

**OP-4a (worst-case optimality does not imply pointwise optimality). [refereed: confirmed, all numbers check]**
$V=$ constants, $D=\{x_1,x_2\}$, $\ell_\infty$ noise: $\Lambda_*=\min\{\|w\|_1:w_1+w_2=1\}=1$; minimax $=\varepsilon$ (consistent with F1: $\omega(t)=t$). The rule $\tilde y\mapsto\tilde y_1$ (weights $(1,0)$) is minimax-optimal; yet at $\varepsilon=1$, $\tilde y=(0,2)$ the section is the singleton $\{1\}$ — the data determine the answer — and this rule outputs $0$, erring by $\varepsilon$ where zero error is achievable. The central rule outputs $1$ with error $0$ and, by OP-3, pointwise-dominates. *Optimal selections can waste information the data have already paid for; the canonical operator cannot.*

**OP-4b (the canonical operator is nonlinear — first at $k=3$). [corrected per referee]**
At $k=2$ the central rule for the constants family is the two-point midrange, and for two numbers $\max+\min=y_1+y_2$: the central rule **is** the linear rule $w=(\tfrac12,\tfrac12)$, everywhere. The nonlinearity claim at $k=2$ was false. At $k=3$ it is real: $V=$ constants, $D=\{x_1,x_2,x_3\}$, $S(\tilde y)=[\max_i\tilde y_i-\varepsilon,\ \min_i\tilde y_i+\varepsilon]$, central rule $=$ the three-point midrange $\tfrac12(\max_i\tilde y_i+\min_i\tilde y_i)$, which violates additivity: $\operatorname{mid}(1,0,0)=\operatorname{mid}(0,1,0)=\tfrac12$ but $\operatorname{mid}(1,1,0)=\tfrac12\ne1$ (all three data vectors realizable for $\varepsilon\ge\tfrac12$). Hence:

**Conclusion (linearity is a property of selections, not of the canonical operator).** For linear families, *some* optimal selection is linear (F7 — a theorem, not a design choice), and for $k\ge3$ the *canonical* (unique conditionally-minimax) operator is genuinely nonlinear even for a one-dimensional linear family. The two agree on worst-case guarantee always, and disagree pointwise wherever the window of OP-2 is nondegenerate.

---

## 6. The total operator and its partiality

The full data-to-value map is a **composition**:
$$\mathbb A(\text{archive},\ \tilde y,\ x)\;=\;A\big(\widehat T(\text{archive})_{D,x},\ \tilde y\big),$$
where $\widehat T$ is the archive-to-trace-set map. This factorization is exact (CP-1/CP-2: the archive can contribute only through $T$), but $\widehat T$ is **partial**:

- it is defined, within the linear exactly-$d$-dimensional class, exactly under the rank conditions of F17 (necessary and sufficient there);
- it is undefined at queries off the covered set, where F18's dichotomy bounds every possible extension;
- **[remark, not a theorem]** with an inexactly identified constraint (archive noise), the composed error splits into a noise term through $\varepsilon\Lambda_*(x)$ and a model term through the perturbed trace set; a quantitative additive bound requires fixing a metric on subspaces (e.g. largest principal angle), the norm on coefficients, and a smallness regime for the perturbation (the perturbation moves $\Lambda_*$ itself), and is recorded here only as a first-order expectation, not an established result. The refereed status of this item is: factorization confirmed, quantitative constant **open**.

---

## 7. Verdict

The result of the theory **is** expressible as a mathematical mapping, in exactly the following stratified sense, each stratum carrying its certificate:

1. **Forced, assumption-free:** the set-valued correspondence $(T,\tilde y)\mapsto S_\varepsilon(\tilde y)$ — this much exists for arbitrary $\mathcal F$ and is the invariant content of the problem (CP-1, OP-1).
2. **Canonical, uniquely determined:** the single-valued selection $A=\operatorname{cen}\circ\mathrm{sec}$ — minimax-optimal always (F1), and the *unique* conditionally-minimax rule on realizable data when $\omega(2\varepsilon)<\infty$ (OP-3). Its structure: $1$-Lipschitz selection over a correspondence whose regularity is exactly the trace modulus (OP-5/OP-6); permutation- and affine-equivariant, gauge-invariant (OP-7); guarantee-monotone (OP-8); nonlinear in general, including for linear families once $k\ge3$ (OP-4b), while linear optimal selections exist in the linear case (F7).
3. **In the identifiable exact-data regime:** a genuine function $\Psi=A(T,\cdot)|_{\varepsilon=0}$, unique on realizable traces, with optimal modulus $\omega$ (F2, OP-6).
4. **Never unary:** no formulation $A(S)$ with $S$ the sample alone exists (OP-10); the second argument is necessary and is minimally-sufficiently the trace set (OP-9).
5. **At the archive level, partial:** the total map is the composition $A\circ\widehat T$ with $\widehat T$ defined exactly on F17's conditions and bounded by F18's dichotomy off the covered set.

**Answer to the posed question:** the correct abstract form is not $A(S)$ but
$$A:\ \big(2^{\mathbb R^{k+1}}\setminus\{\emptyset\}\big)\times\mathbb R^k\ \dashrightarrow\ \mathbb R,\qquad A(T,\tilde y)=\operatorname{cen}\big(S_\varepsilon(\tilde y)\big),$$
a partial binary operator — correspondence followed by canonical selection — whose two arguments carry, respectively and irreducibly, the constraint and the selection identified in `core_principle.md`.
