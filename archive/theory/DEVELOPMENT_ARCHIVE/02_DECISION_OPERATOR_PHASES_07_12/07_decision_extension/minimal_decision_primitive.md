# The Minimal Decision Primitive (Part II)

> **Status:** Phase-7, 2026-08-03. Frozen corpus cited, not modified. New results carry **DE-P** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. Mandate discipline: probability measures, priors, expected losses, preference orderings, minimax-regret rules, and ambiguity sets are treated as *candidates*, not premises; the primitive is **derived** from what a well-defined selection requires.

---

## 1. What a decision rule must do

Fix a query set $Q$ and a loss $L:\mathcal A\times\mathbb R^m\to\mathbb R$. Each action $a$ induces its **loss profile** on the identified set:
$$P_a:\ J_Q(O)\to\mathbb R,\qquad P_a(v)=L(a,v).$$
Identification hands the decision layer exactly the indexed family $\{P_a\}_{a\in\mathcal A}$ of profiles — nothing else (DE-S1). A *selector* $\mathsf D$ assigns to every decision problem (a nonempty identified set $J$, an available action menu $\mathcal B\subseteq\mathcal A$) a nonempty choice set $\mathsf D(J,\mathcal B)\subseteq\mathcal B$. Requirements, from weakest to strongest:

- **(R0) Well-definedness.** $\mathsf D$ is a function: same $(J,\mathcal B)$, same choice.
- **(R1) Dominance-respect.** If $P_a\le P_{a'}$ pointwise on $J$ with strict inequality somewhere, then $a'\notin\mathsf D(J,\mathcal B)$ whenever $a\in\mathcal B$. (Violating R1 contradicts the identified set itself — the only unconditional knowledge available.)
- **(R2) Menu-coherence** (Sen $\alpha$+$\beta$ / WARP): if $\mathcal B'\subseteq\mathcal B$ and $\mathsf D(J,\mathcal B)\cap\mathcal B'\ne\emptyset$ then $\mathsf D(J,\mathcal B')=\mathsf D(J,\mathcal B)\cap\mathcal B'$; and choices depend on actions only through their profiles.

---

## 2. What identification already provides: the dominance preorder

**Proposition DE-P1. [proved]**
The pointwise order $P_a\le P_{a'}$ on profiles — the **dominance preorder** — is determined by $J_Q(O)$ and $L$ alone, and it is the *finest* action comparison so determined: any comparison of two dominance-incomparable actions varies across decision contexts consistent with the same $J_Q(O)$ (witness: DE-S3's two weightings reverse the comparison of $a=\tfrac13$ vs $a=\tfrac23$ on $J=\{0,1\}$, squared loss). $\square$

So the question "what is absent from identification but necessary for selecting an action" has a first, exact answer: **the ability to compare dominance-incomparable loss profiles.**

---

## 3. The derivation: completions of dominance

**Theorem DE-P2 (dichotomy). [proved]**
(i) If only R0 is demanded, the minimal primitive is a bare **choice function** on undominated sets — always sufficient, but carrying no structure: it is itself an unexplained selection, is unfalsifiable, and (Part VIII) not learnable.
(ii) If R0–R2 are demanded on a domain rich in finite menus, then $\mathsf D$ is **rationalized by a total preorder $\preceq$ on loss profiles** ($\mathsf D(J,\mathcal B)=\{a\in\mathcal B:P_a\preceq P_b\ \forall b\in\mathcal B\}$), and R1 forces $\preceq$ to be **monotone**: $P\le P'$ pointwise $\Rightarrow P\preceq P'$, strictly when strict somewhere on $J$.
*Proof.* (ii) is the classical revealed-preference rationalization (Arrow–Sen: $\alpha+\beta$ over finite menus $\iff$ rationalizability by a total preorder), applied to profiles because R2 makes choice profile-determined; monotonicity is R1 read into $\preceq$. (i) is immediate. $\square$

**Theorem DE-P3 (existence, non-uniqueness, and the exact residual). [proved]**
(i) *Existence:* the dominance preorder always admits monotone total-preorder completions (Szpilrajn's extension theorem, order-extension form).
(ii) *Non-uniqueness:* the completion is unique (as a preorder on achievable profiles) **iff** dominance is already total on them — which for the scalar prediction problem of DE-S2 happens iff $J_Q(O)$ is a singleton, i.e. iff identification already decided everything.
(iii) Hence the **residual decision ambiguity is exactly the multiplicity of monotone completions of the dominance preorder**; picking one completion *is* the decision primitive. $\square$

$$\boxed{\ \textbf{Minimal decision primitive}\ =\ \text{a monotone total-preorder completion of the dominance order on loss profiles.}\ }$$

This answers the mandated question. Identification supplies the partial order; decision requires — and requires only — one of its completions. Everything on the candidate list is a *generator* of such completions, not the primitive itself.

---

## 4. The candidates, located as completions

| Candidate object | Induced completion of dominance | Extra declared content |
|---|---|---|
| single law (prior + likelihood) $\mu$ | $P\preceq P'\iff\int P\,d\mu\le\int P'\,d\mu$ | a measure on $J$ (or $I(O)$) — DE-S5: never canonical |
| ambiguity set of laws $\mathcal Q$, worst-case expected loss | compare $\sup_{\mu\in\mathcal Q}\int P\,d\mu$ | a set of measures |
| minimax (sup-order) | compare $\sup_J P$ | **none** — set-definable; the frozen theory's own rule |
| minimin, Hurwicz$_\alpha$ | compare $\alpha\sup_J P+(1-\alpha)\inf_J P$ | one scalar $\alpha$ (already a $\Delta$) |
| lexicographic (e.g. $\sup$, then $\inf$) | lexicographic on $(\sup_J P,\inf_J P)$ | none, but see DE-P5 |
| minimax regret | see DE-P4 | a **reference function** + menu dependence |

**Proposition DE-P4 (minimax regret is not a profile preorder). [proved]**
Minimax regret compares $\sup_{v\in J}[P_a(v)-h_{\mathcal B}(v)]$ with $h_{\mathcal B}(v)=\inf_{b\in\mathcal B}P_b(v)$. Since $h_{\mathcal B}$ depends on the menu, the induced comparison of two fixed profiles can reverse as $\mathcal B$ varies — regret violates contraction consistency (Sen $\alpha$), a classical fact. Consequently regret is **not** rationalized by any single preorder on profiles; its primitive is the *pair* (sup-order, declared reference function $h$), with $h=h_{\mathcal B}$ one conventional choice. It satisfies R0–R2 only at a fixed action universe. This does not disqualify regret; it prices it: one more declared object than minimax. $\square$

**Proposition DE-P5 (the primitive is strictly weaker than any real-valued risk functional). [proved]**
There exist monotone total preorders on profiles representable by **no** real-valued functional: on two-point $J$, profiles live in $\mathbb R^2$, and the lexicographic order is total, monotone, and non-representable (Debreu). Hence demanding "an expected loss" or any numerical risk index is *strictly stronger* than the minimal primitive; the preorder itself is the floor. $\square$

**Proposition DE-P6 (canonical completions exist but are exhausted by set-symmetric statistics). [proved, scoped]**
Completions definable from $J$ and $L$ alone must be invariant under all automorphisms of the identified structure; on profiles this confines them to functions of member-symmetric statistics — in the bounded scalar case, of $(\inf_J P,\sup_J P)$ and the value set of $P$ — yielding the minimax / minimin / Hurwicz-with-declared-$\alpha$ / lexicographic band. Among these, minimax is singled out by the frozen optimality theory (Theorem 1: it is the exact worst-case-optimal selection, with the certificate) and by DE-S4 (the only intrinsic equivariant behavior on symmetric sets). Every completion *outside* this band — in particular every non-uniform weighting — requires breaking member symmetry, i.e. a declared $\Delta$. $\square$

---

## 5. The answer to the mandated question

**What information is absent from identification but necessary for selecting an action?**

A **comparability structure across the admissible members** — formally, a monotone completion $\preceq$ of the dominance order on loss profiles. Identification is invariant under every relabeling of $I(O)$ that preserves the set; a selection among dominance-incomparable actions must break exactly that invariance. The candidates on the mandate's list are sufficient generators of completions, ordered by declared content:

$$\text{choice function (no structure)}\ \prec\ \text{monotone total preorder (the minimum)}\ \prec\ \text{real risk functional}\ \prec\ \text{ambiguity set of laws}\ \prec\ \text{single law}.$$

The frozen minimax rule is the unique inhabitant of the "no added content" tier with an optimality certificate; the single law is the maximal tier. Part III–IV determine which tier the data (historical members included) can actually justify.
