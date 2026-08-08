# Irregular-Archive Theory (Part II)

> **Status:** Phase-6, 2026-08-03. Removes the common-core assumption of F17. Sources: frozen corpus (F17, DM-3, DM-7). New results **CI-B**; tags **[proved] / [conditional] / [impossible] / [open]**. Adversarially refereed (PASS-WITH-FIXES); the load-bearing corrections — the real constant-rank replacement of the projective intersection theorem in CI-B1, the explicit holonomy witness in CI-B4, the nonzero-subgraph fix in CI-B7, and the honest B8 literature status — are incorporated.

**Setting.** Exactly-$d$ linear class ($d$ known); candidate objects are subspaces $V\subseteq\mathbb R^U$ of dimension exactly $d$, $|U|=N>d$. Members are elements of the true $V$, observed with **exact** values (noise deferred to Part III) on arbitrary designs $D_a\subseteq U$, $|D_a|=k_a$. Incidence pattern $P=\{(a,x):x\in D_a\}$. "Generic position" of observed values unless stated.

**The question.** What structure of $P$ is necessary/sufficient to identify the finite-window information $V|_U$ that the frozen operator needs (windows of size $\le k+1$, DM-3)?

---

## 1. The correct object is a sheaf-like consistency system, derived — not a graph

Reading identifiability as "recover the global section $V|_U$ from locally observed restrictions" makes the object a **presheaf of restrictions with a gluing rule**: to each $S\subseteq U$ the (partially known) window $V|_S$, restriction maps $V|_{S'}\to V|_S$ ($S\subseteq S'$), and a compatibility/gluing law across overlaps. Simple graph connectivity is provably **not** the invariant (CI-B3); the invariant is whether local sections glue to a unique global section, which is a sheaf-theoretic — not graph-theoretic — condition. The derivation below produces this object rather than naming it in advance, and isolates its **holonomy** obstruction (CI-B4).

---

## 2. Necessity

**Theorem CI-B1 (general counting necessity — extends DM-7 beyond the core pattern). [proved, with the real-analytic scoping]**
For any pattern, if $\sum_a(k_a-d)_+ < d(N-d)$ then generic identifiability of $V|_U$ fails.
*Proof (refereed — uses the real constant-rank theorem, NOT the complex projective intersection theorem, which does not transfer to real point sets).* Each member's consistency locus $\Sigma_a=\{V\in\operatorname{Gr}(d,N):y_a\in\operatorname{proj}_{D_a}V\}$ has, at the true $V$, codimension equal to the **corank of the restriction map** $\operatorname{proj}_{D_a}:V\to\mathbb R^{D_a}$, i.e. $(k_a-d)_+$ on the unisolvent cell (smaller off it — which only strengthens necessity). The consistent set is cut in the smooth real manifold $\operatorname{Gr}_{\mathbb R}(d,N)$ (dimension $d(N-d)$) by $c=\sum_a(k_a-d)_+$ real-analytic equations. At the true $V$ — a *known smooth real solution* — the constant-rank / real implicit-function theorem gives the solution set local real dimension $\ge d(N-d)-c$. When $c<d(N-d)$ this is $>0$: a real-analytic arc of genuinely distinct consistent subspaces through $V$, so identification fails. $\square$
*Scope discipline (refereed):* the bound is **necessity only** (a dimension lower bound); it must not be read as sufficiency (the equations can be dependent and the true component larger, never provably smaller). The pathology that a real variety's point set can undershoot its complex scheme dimension is avoided precisely by evaluating local dimension *at the known smooth real point*.

**Theorem CI-B2 (per-point necessity, no core). [conditional: non-pivot genericity]**
If some $x\in U$ is observed by $|J_x|<d$ members, and $V|_{U\setminus\{x\}}$ has dimension $d$ (i.e. $x$ is a **non-pivot** coordinate: the restriction to $U\setminus\{x\}$ is injective on $V$), then $V|_U$ is not identifiable: the value at $x$ is a linear functional $a_x\in\mathbb R^d$ subject to $|J_x|<d$ observed linear equations, leaving a realized $(d-|J_x|)$-dimensional affine freedom (this is exactly the F17 ambiguity computation with submatrix rank $\le|J_x|<d$). **Hence every point needs $\ge d$ observing members** — *under the standing genericity that every coordinate is non-pivot*.
*Pivot caveat (refereed):* if $\dim V|_{U\setminus\{x\}}=d-1$, $x$ carries an intrinsic degree of freedom (a member supported only at $x$) and the count changes — a genuinely different regime where a pivot coordinate can sometimes be pinned by fewer observers. State the non-pivot hypothesis explicitly; it is generic for $N>d$ but is a hypothesis, not a theorem.

**Theorem CI-B3 (connectivity is NOT sufficient — the decisive counterexample). [proved]**
Any pattern in which every member observes exactly $d$ points — *however connected the incidence structure* — leaves a co-null set of consistent candidates. Each such member has $(k_a-d)_+=0$, so $\sum(k_a-d)_+=0<d(N-d)$ (for $N>d$), and CI-B1 gives a positive-dimensional consistent set; independently, each $k_a=d$ member excludes only a null set (established), and a finite union of null sets has continuum complement in the positive-dimensional $\operatorname{Gr}(d,N)$. **Connectivity is orthogonal to identifiability** — the intuition that reaches for it is wrong; the invariant is the excess count plus its distribution, not connectivity.

---

## 3. Transport through overlaps: the gluing lemma and its holonomy

**Lemma CI-B4a (gluing). [proved]**
If $V|_{S_1}$ and $V|_{S_2}$ are known (both dimension $d$) and the overlap is **unisolvent** ($\dim V|_{S_1\cap S_2}=d$, so restriction is injective on $V$, requiring $|S_1\cap S_2|\ge d$ with unisolvent columns), then
$$V|_{S_1\cup S_2}=\{(u,w):u\in V|_{S_1},\,w\in V|_{S_2},\ u,w\text{ agree on }S_1\cap S_2\}.$$
*Proof.* A compatible pair is $(g|_{S_1},h|_{S_2})$ for members $g,h\in V$ agreeing on the unisolvent overlap; a member of $V$ is determined by its restriction to any unisolvent set, so $g=h$ and they glue. The unisolvent overlap makes each restriction an isomorphism onto its image, upgrading agreement-on-overlap to agreement-of-members. $\square$

**Theorem CI-B4b (holonomy — the sheaf obstruction, with explicit witness). [proved]**
When an overlap has **rank $<d$**, compatibility is weaker than membership: each rank-1 overlap collapses pairwise compatibility to a *single* functional, and around a cycle these functionals can fail to span the member space, leaving a **realized residual gauge** (holonomy) — the glued object strictly contains $V|_{\bigcup S_i}$.
*Explicit $d=2$ witness (refereed).* $U=\{1,\dots,6\}$, $V=\operatorname{span}(e,f)$ with $e=(1,1,1,1,1,1)$, $f=(0,1,2,0,1,2)$; patches $S_1,S_2,S_3$ with pairwise overlaps $\{1,4\},\{2,5\},\{3,6\}$ chosen so $f$ takes **equal values on each overlap pair** ($f_1=f_4=0$, $f_2=f_5=1$, $f_3=f_6=2$). Every member $\alpha e+\beta f$ then takes equal values on both overlap coordinates, so each overlap has rank $1$ and imposes one scalar equation; the three overlap-functionals are affinely dependent around the 3-cycle, leaving a residual 1-parameter $\beta$-gauge — **all three pairwise gluings individually consistent, global identification fails.**
*Sharp statement (refereed):* rank-1 overlaps are **necessary but not sufficient** for holonomy; holonomy occurs **exactly when the transported partial identifications fail to jointly determine the member space around some cycle.** This is the sheaf-cohomological content: local consistency $\ne$ global sectionability.

---

## 4. Sufficiency

**Theorem CI-B5 (chaining sufficiency). [conditional on the stated local conditions]**
If (i) the patches **cover $U$**, (ii) on each patch the local archive satisfies F17-type conditions and identifies $V|_{S_i}$ to **full local dimension $d$**, and (iii) the patch-overlap graph is **connected with every overlap unisolvent** ($\ge d$ points, rank $d$), then $V|_U$ is identified, by induction along a spanning tree via CI-B4a.
*Cycles need no separate check (refereed proof).* With every overlap unisolvent, each transport is a genuine isomorphism recovering the *literal* restriction of the one true $V$; the result is tree-independent, so cycle consistency is automatic — the contrapositive of CI-B4b (unisolvence everywhere $\Rightarrow$ trivial holonomy).

**Theorem CI-B6 (truncation suffices — coherence with DM-3). [proved]**
For the exactly-$d$ class, knowing all windows of size exactly $d+1$ determines $V|_U$: pick a unisolvent $d$-set $B\subseteq U$ (exists — a basis matrix of $V|_U$ has an invertible $d\times d$ minor); let $v\in V$ match $w$ on $B$; for each $x$, the hypothesis on $B\cup\{x\}$ plus unisolvence of $B$ forces $w(x)=v(x)$; so $w=v|_U\in V|_U$. Hence the operator's needs (windows of size $\le k+1$, $k\ge d$) are met exactly by the identified $V|_U$ and nothing more — the finite-window truncation of Part I/DM-3 is the right target here too.

**Theorem CI-B7 ($d=1$ characterization). [proved, with the nonzero-subgraph fix]**
For $d=1$, identifiability of the line $V|_U$ holds **iff the bipartite member–point incidence subgraph restricted to pairs with nonzero observed value is connected across all coordinates that are nonzero in $V$.** Ratios $v_x/v_{x'}$ transport along member paths; cycles are automatically consistent (exact data telescopes the ratio to $1$); disconnection leaves a free relative scale between components (realized $(\#\text{components}-1)$-dimensional ambiguity).
*Zero handling (refereed — the raw graph is the wrong object):* a coordinate that is **zero in $V$** is forced to $0$ and needs only one observer; the danger is a **nonzero** coordinate whose only observers pass through zero coordinates — transport breaks there even when the raw incidence graph is connected. Hence connectivity must be of the **nonzero-value** subgraph.

---

## 5. Literature anchor and the honest open item

**CI-B8 [proved status of the literature; one sub-item open].**
The relevant characterization is **Pimentel-Alarcón, Boston, and Nowak**, "A Characterization of Deterministic Sampling Patterns for Low-Rank Matrix Completion," *IEEE JSTSP* 10(4):623–636, 2016 (arXiv:1503.02596; plus the published Corrections) — note **Boston**, not "Boumal." Their model is the transpose of ours (fixed column space, many drawn columns ↔ our one subspace $V$, members observed on designs). Findings, mapped:
- **Finite completability is exactly characterized** (necessary and sufficient): each column $\ge r+1$ entries, plus a Laman/Hall-type matching condition ($|N(S)|\ge|S|+r$ for every column subset of an $r(N-r)$-column submatrix). **[proved in the literature]**
- **Unique completability has only *sufficient* combinatorial conditions**; **no exact necessary-and-sufficient combinatorial characterization for general $d$ is known.** **[open in the literature]**
- **Consistency with CI-B1.** A column with $r+1$ entries contributes exactly $1$ excess ($=(k_a-d)_+$ at $k_a=d+1$); needing $r(N-r)$ such columns $\iff$ needing $\sum(k_a-d)_+\ge d(N-d)$ — **CI-B1 is exactly the dimension-count necessity**, matching the literature's entry count. The literature adds that the count must be **supplemented by a Hall-type matching condition** to become sufficient (excess concentrated on too few coordinates does not identify others — the CI-B2/CI-B3 phenomenon), and even then yields *finite*, not *unique*, completability.

---

## 6. Ledger of the irregular-archive theory

| Result | Statement | Tag |
|---|---|---|
| CI-B1 | $\sum(k_a-d)_+\ge d(N-d)$ necessary (real constant-rank proof) | **proved** |
| CI-B2 | every point needs $\ge d$ observers (non-pivot genericity) | **conditional** |
| CI-B3 | connectivity not sufficient (zero-excess continuum) | **proved** |
| CI-B4a | gluing across unisolvent overlaps is exact | **proved** |
| CI-B4b | rank-$<d$ overlaps → holonomy; explicit $d=2$ witness | **proved** |
| CI-B5 | connected patches + unisolvent overlaps + local F17 → identified; cycles automatic | **conditional** |
| CI-B6 | size-$(d{+}1)$ windows determine $V|_U$ | **proved** |
| CI-B7 | $d=1$ iff nonzero-value subgraph connected | **proved** |
| CI-B8 | count-necessity matches literature; **exact unique-completability characterization for general $d$** | **open (literature)** |

**Net:** the necessary structure is the counting bound CI-B1 plus the per-point bound CI-B2; a clean sufficient structure is the unisolvent-overlap chaining CI-B5; the obstruction between them is holonomy CI-B4b. The exact necessary-and-sufficient combinatorial line for general $d$ is the one genuinely open item — and (see `final_theory_stopping_criterion.md`) it is **non-blocking**: CI-B1 + CI-B5 sandwich it, and Part III's validity never requires it.
