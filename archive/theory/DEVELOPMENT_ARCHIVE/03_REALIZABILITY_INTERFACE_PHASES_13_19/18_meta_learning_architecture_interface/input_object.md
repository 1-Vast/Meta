# Input Object (§1)

> **Status:** Phase-18 (architecture interface derivation), 2026-08-03. Base: the Phase-17 reconstruction (`../17_meta_operator_reconstruction/`, MR-numbers), unmodified. New results carry **MI-** numbers, tagged **[proved] / [declared] / [forbidden]**. Role: mathematics only — the interface a future model must approximate, never the model.

---

## 1. The exact input

**Definition MI-1 (support object).** The current-task input is
$$S_T\ =\ \Big(\ \{(x_i,\tilde y_i)\}_{i\le k}\ \text{(finite, }k\le5\text{)},\ \ \varepsilon,\ \ c\ \text{(optional declared auxiliary label)},\ \ (\kappa\text{-context is a function of the above)}\ \Big),$$
together with the declared index arguments (query object $Q$, specification $\gamma$) — which are *arguments of the operator*, not part of $S_T$. At meta-training time the input is additionally the history $H_N$ (sequence/multiset-typed, MC-2/3).

- **Uncertainty information:** $\varepsilon$ is a **mandatory input** (frozen F1 Rem. 1.4: no rule is optimal across noise levels); the partiality flags are outputs, not inputs.
- **Task identity information:** **[forbidden]** — a task ID is a non-transferable lookup key reducing to case-(i) of MC-7; the only admissible task-level side information is the *declared observable* auxiliary label $c$ (the Phase-6 fiber datum), which enters through $\kappa$.

## 2. The mathematical type of $S_T$ — decided with proofs

**Theorem MI-2 (the support is a finite set; the history is a multiset; nothing else is admissible). [proved]**
(i) **Not a sequence.** The identified set $I(O)=\bigcap_i\{f:|\tilde y_i-f(x_i)|\le\varepsilon\}$ is a conjunction — invariant under every reordering of the observations; so is the canonical band $b_{\mathrm{can}}$ and hence every family member (MR-3). An order-sensitive input type would carry structure the operator provably cannot use; typing it away is exact, not simplifying.
(ii) **Set, not multiset, within a task.** The within-task constraints are **idempotent**: an exact duplicate $(x,\tilde y)$ repeats a constraint and changes nothing ($\cap$ is idempotent); two *distinct* values at the same $x$ are distinct elements of the set and both bind (their intersection tightens the interval at $x$) — so set semantics is information-preserving and information-complete for the identification channel. Crucially, there is **no within-task frequency semantics to preserve**: the frozen noise model is adversarial-bounded, not stochastic (the standing firewall) — which is also why the input is **not a distribution or a measure**: typing $S_T$ as an empirical measure would silently inject within-task sampling semantics that the frozen theory rejects (DE-H4), and would make duplicate observations carry mass they provably do not carry.
(iii) **Multiset across tasks.** $H_N$ retains multiplicities — frequencies are the population learning target (MC-3); the set/multiset split *is* the two-channel separation, one level down.
(iv) **Not a graph.** No relational structure among observations is declared by the theory; any edge structure would be an undeclared assumption injecting information the ledger cannot echo. If a deployment declares one (e.g. structured $\mathcal X$), it enters through the declared family/closure class — i.e. through the *identification* side — never through the input type of the operator. $\square$

## 3. Invariance requirements (the contract on any consumer of $S_T$)

**MI-3 (invariance ledger). [proved / declared]**
1. **Permutation invariance** in the observation list — required, proved in `permutation_invariance.md`.
2. **Duplicate idempotence** within a task — required by MI-2(ii): $A_\theta(S\cup\{s\})=A_\theta(S)$ when $s\in S$.
3. **Gauge invariance** — no meaning may attach to internal coordinates of any representation of $S_T$ (MC-7; frozen CP §4.5); only the emitted $\mathbb B$-point and side channels are semantic.
4. **Channel typing** — the consumer may read $S_T$ only through the two declared channels: identification (set-image; produces $b_{\mathrm{can}}$ and the current identified objects) and population (through $\kappa(S_T)$ and, at training time, multiset statistics of $H_N$). Any other read is leakage (audited in `failure_audit.md`).
