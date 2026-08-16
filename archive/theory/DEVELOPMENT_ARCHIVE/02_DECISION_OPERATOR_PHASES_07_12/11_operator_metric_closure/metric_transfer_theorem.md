# Metric Transfer Theorem (§1)

> **Status:** Phase-11 (operator metric closure), 2026-08-03. Phases 0–7 frozen; Phase-10 objects cited, with its operator-consistency *sentence* retracted (the audit's counterexample is valid) and replaced by the theorems below. Audit of record: `../11_final_meta_operator_audit/FINAL_VERDICT.md` (`META_OPERATOR_LEARNABILITY_INCOMPLETE`). New results carry **OM-** numbers, tagged **[proved] / [declared] / [retracted]**.

---

## 0. The defect, exactly

Phase-10's LC-15 controls $d_{\mathrm{desc}}$ (uniform endpoint deviation) and then *asserted* consistency in the operator metric $d_{\mathbb M}$ (uniform Hausdorff-TV on feasible law sets). The audit's counterexample — dyadic queries $Q_k$ with $|\Omega_{Q_k}|=2^k$, records revealing $U\sim\mathrm{Unif}[0,1]$ exactly, singleton-cell atlases — has finite-VC indicators (dyadic intervals), hence $d_{\mathrm{desc}}\to0$, while for every finite $N$ the empirical law sits in the estimated polytope with $\mathrm{TV}(\hat p_n,p_n)=1-N/n\to1$: $\sup_Q$ Hausdorff-TV stays near $1$. **The assertion is retracted.** The transfer needs a geometric hypothesis; this file supplies it.

## 1. Route A (adopted): restricted operator class, transfer proved

**Assumptions (Route A), each declared:**
- **(RA1) Bounded query outcome complexity.** $\sup_{Q\in\mathcal Q_0}|\Omega_Q|\le\bar n$ and $\sup_{\iota\in\mathcal I}|\mathcal E_\iota|\le\bar e$: outcome spaces and per-index event atlases are uniformly finite. *(This is what the counterexample violates: $|\Omega_{Q_k}|=2^k$ is unbounded.)*
- **(RA2) Compact/effectively finite constraint geometry.** Follows from (RA1) [proved, OM-1(i)]: the constraint matrices range over a finite set of patterns.
- **(RA3) Finite effective dimension.** The Phase-10 VC declaration $d^*$ for the induced indicator class, retained (with the §3 context-allocation fix).

**Lemma OM-1 (uniform Hoffman constant). [proved]**
(i) Under (RA1), each index's polytope is $P(b)=\{p\in\mathbb R^{\Omega}:A p\le b\}$ where $A$ stacks the simplex rows and the $0/1$ event-incidence rows: at most $\bar n$ outcomes and $\bar e$ events, so $A$ ranges over a **finite** family $\mathcal A(\bar n,\bar e)$ of matrices.
(ii) For each fixed $A$, Hoffman's error bound gives $H(A)<\infty$ with $\mathrm{dist}_{\ell_1}\!\big(x,P(b)\big)\le H(A)\,\big\|(Ax-b)_+\big\|_\infty$ for every $b$ with $P(b)\ne\emptyset$ and every $x$. Set $\bar H=\max_{A\in\mathcal A(\bar n,\bar e)}H(A)<\infty$ — finite as a maximum over a finite set, and **independent of the index and of the endpoint vectors**. $\square$

**Theorem OM-2 (metric transfer). [proved]**
Let $M_1,M_2\in\mathbb M$ have, at every index, nonempty values defined by the same constraint pattern with endpoint vectors differing by at most $\varepsilon=d_{\mathrm{desc}}(M_1,M_2)$. Under (RA1):
$$d_{\mathbb M}(M_1,M_2)\ \le\ \tfrac12\,\bar H\,\varepsilon\qquad\text{(Hausdorff distance in TV }=\tfrac12\ell_1\text{)}.$$
*Proof.* Fix an index with matrix $A$ and vectors $b_1,b_2$, $\|b_1-b_2\|_\infty\le\varepsilon$. For $x\in P(b_1)$: $(Ax-b_2)_+\le(b_1-b_2)_+\le\varepsilon$ componentwise, so $\mathrm{dist}_{\ell_1}(x,P(b_2))\le\bar H\varepsilon$ by OM-1; symmetrically for $x\in P(b_2)$; hence the Hausdorff-$\ell_1$ distance is $\le\bar H\varepsilon$, i.e. Hausdorff-TV $\le\tfrac12\bar H\varepsilon$; take the sup over the index class — the constant is uniform. Nonemptiness: the estimator's polytope contains the empirical law by construction (its coordinates are the empirical constraint centers), and the target's contains the population law; both sides of the comparison are nonempty whenever the values are defined (else the failure flag fires upstream). Corollary: $d_{\mathrm{desc}}\to0\ \Rightarrow\ d_{\mathbb M}\to0$, with an explicit linear modulus. $\square$

**Remark (elementary special case).** For singleton-cell atlases, TV between laws on $\le\bar n$ outcomes obeys $\mathrm{TV}\le\tfrac{\bar n}{2}\max_\omega|p(\omega)-p'(\omega)|$ — a direct $\bar H\le\bar n$-type bound showing concretely how (RA1) kills the counterexample: the audit's construction needs $n\to\infty$ to defeat the transfer, and (RA1) forbids it.

## 2. Route B (declared alternative): unrestricted atlas with stability assumptions

Keep unbounded outcome spaces and the Hausdorff-TV metric; **declare** the stability inequality
$$d_{\mathbb M}(M_1,M_2)\ \le\ C\ d_{\mathrm{desc}}(M_1,M_2)$$
directly, via geometric sufficient conditions, each of which must then be verified per instantiation: **(RB1)** uniform Hoffman/conditioning bound $\sup_\iota H(A_\iota)\le C'<\infty$ (no longer automatic — the matrix family is infinite); or **(RB2)** a TV-determining atlas with uniformly bounded dual constants: for every index, every law in the value set is determined within TV-$\varepsilon$ by its atlas evaluations within $\varepsilon/C''$ (equivalently, the atlas contains a basis whose inverse has uniformly bounded $\ell_1$-operator norm). Route B is mathematically legitimate but shifts the burden to an unbounded verification; **Route A is adopted** as the contract default because its hypotheses are finitely checkable declarations and its constant is proved, not assumed. The choice is echoed in the ledger; the audit's counterexample stands as the permanent witness that *some* such hypothesis is irreducible — no assumption-free transfer exists. **[proved by the counterexample]**
