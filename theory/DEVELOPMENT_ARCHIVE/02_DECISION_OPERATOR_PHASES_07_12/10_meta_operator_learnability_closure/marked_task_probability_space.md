# Marked Task Probability Space (§1)

> **Status:** Phase-10 (meta-operator learnability closure), 2026-08-03. Phases 0–7 frozen; repaired decision statements cited, not modified. Audit of record: `../11_meta_learning_final_audit/` (`META_LEARNING_OPERATOR_INVALID`, finding 2: the ideal target invoked $P(g(f_T)\mid T)$ from a law on the *observable* task space, which does not contain $f_T$ — MC-16/ML-L1 false as written, **retracted here and replaced**). New results carry **LC-** numbers, tagged **[proved] / [conditional] / [declared]**.

---

## 1. The defect, and the repair strategy

The observable task $T=(O,S,Q,\gamma)$ does not contain the latent member $f_T$; a law $\Pi$ on the observable space $\mathbb T$ cannot define any conditional involving $f_T$ — two different joint lifts share the same observable law with different latent conditionals (the audit's point, adopted as LC-3(iii) below). The repair: **mark tasks with their latent decision object** — not with the full member, which would break measurability for uncountable $\mathcal X$ — and let every latent conditional live on the marked space by declaration.

## 2. The marked task space

**Definition LC-1 (latent decision object and marked space).**
Fix a declared **countable** query-specification atlas $\mathcal Q_0$ (query tuples $Q$ with their pushforward maps $g_Q$; countability is the separability declaration consumed in §3/§5). The **latent decision object** of a task with member $f$ is
$$Y\ =\ \big(g_Q(f)\big)_{Q\in\mathcal Q_0}\ \in\ \mathbb Y\ :=\ \prod_{Q\in\mathcal Q_0}\mathbb V_Q,$$
where $\mathbb V_Q$ is the (standard Borel) value space of $g_Q$ ($\mathbb R^{m_Q}$, $\mathbb R$, or the finite order space $\Omega_{m_Q}$). $\mathbb Y$ is a countable product of standard Borel spaces, hence standard Borel — the measurability that $\mathbb R^{\mathcal X}$ lacks. The **marked task space** is
$$\mathbb T^\bullet\ =\ \mathbb T\times\mathbb Y,\qquad T^\bullet=(T,\,Y_T),$$
standard Borel as a product. (Marking with $(T,f_T)$ itself is the declared alternative when $\mathcal X$ is countable; the $Y$-marking is the general form — "equivalent latent decision object", and it is exactly the part of $f_T$ the decision layer ever consumes.)

**Definition LC-2 (the two laws — explicitly distinguished).**
- **Marked latent task law $\Pi^\bullet$:** a probability law on $\mathbb T^\bullet$, **declared** (as part of the population model), required to satisfy the frozen support coupling: $Y_T$ is a.s. consistent with the observable record — $g_Q(f)\in g_Q(I(O_T))$ for all $Q\in\mathcal Q_0$, a.s. (the mark of a task always lies in that task's identified image; this is a *constraint linking mark and record*, not an assumption — it holds by definition of $I$).
- **Observable task law $\Pi_{\mathrm{obs}}=\mathrm{proj}_{\mathbb T\,*}\,\Pi^\bullet$:** the pushforward to $\mathbb T$ — the only law the history ever samples.

**Typing rule LC-3 (no unanchored latent conditionals). [declared, enforced; (iii) proved]**
(i) Every expression of the form $P(Y\in\cdot\mid\cdot)$, $P(g_Q(f)\mid\cdot)$, or any conditional touching the mark, is well-formed **only** relative to a declared $\Pi^\bullet$ (regular conditionals exist on the standard Borel $\mathbb T^\bullet$ — this, and only this, is the corrected existence route; the retracted MC-16/ML-L1 asserted it from $\Pi_{\mathrm{obs}}$).
(ii) Statements estimable from data are typed over $\Pi_{\mathrm{obs}}$ (plus the per-task identified objects, which are $\mathbb T$-measurable).
(iii) *The gap between the two is real and exactly the identification tier:* two lifts $\Pi^\bullet_1\ne\Pi^\bullet_2$ with $\mathrm{proj}_*\Pi^\bullet_1=\mathrm{proj}_*\Pi^\bullet_2$ and different latent conditionals exist whenever some task's record leaves its decision object undecided — witness: records whose identified image contains two points $\{y,y'\}$; put the conditional mass on $y$ under lift 1 and $y'$ under lift 2; both lifts satisfy the support coupling and project identically. **[proved]** Hence: latent conditionals are declared-model objects; observable data identify them only up to the lift class (§4). $\square$

## 3. What the marked space repairs downstream

**LC-4 (relocations). [proved]**
(i) The ideal meta-target is redefined on $\mathbb T^\bullet$: $M^\star_{\Pi^\bullet}(c,Q,\gamma)=\Pi^\bullet\big(g_Q\text{-component of }Y\in\cdot\ \big|\ \kappa(O)=c\big)$ — existing by regular conditionals **given a declared $\Pi^\bullet$** (`existence_identification_learning.md`).
(ii) The Phase-9 conditioning theorem (MC-11) already lived on a joint (member, record) law; it is re-anchored verbatim on $\mathbb T^\bullet$ with the kernel-indexed sufficiency unchanged — no content change, the ambient space is now declared rather than informal.
(iii) The support coupling of LC-2 is what makes the frozen ceiling automatic on the marked space: any $\Pi^\bullet$-conditional given the record is supported in the record's identified image, by construction — DE-H2 becomes a property of the space rather than a separate argument. $\square$
