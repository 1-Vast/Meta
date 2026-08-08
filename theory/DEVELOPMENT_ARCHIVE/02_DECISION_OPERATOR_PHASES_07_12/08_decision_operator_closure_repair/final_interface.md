# Final Interface (Closure Target 5)

> **Status:** Phase-8.2, 2026-08-03. The closed decision-operator interface, superseding the interface sections of Phases 8 and 8.1. Domains and codomains complete; every validity predicate cites a proved (or explicitly declared) result; the three scoping corrections of the audit are absorbed. New results **DC-I1–I2**. No architectures.

---

## 0. Spaces

$\mathcal O=(\mathcal X\times\mathbb R)^{\le k}\times[0,\infty)$ (current record); $\mathcal H$ = finite **sets** of historical task records (see $V_I$(iii)); $\mathcal K_m$ = nonempty closed subsets of $\mathbb R^m$; $\mathcal K^{\mathrm{fin}}_m$ = finite subsets (witness lists); $\Omega_m=S_m$ (declared tie convention); $\mathfrak Q(\Omega)$ = constraint polytopes of laws on the relevant outcome space (DC-R1: full space, events as constraints); $\mathrm{Rung}=\{1{:}\varnothing,\,2{:}\text{marginal},\,3{:}\text{conditional},\,4{:}\text{posterior}\}$ (DC-C3); $\mathrm{Ctx}$ = declared decision contexts $\gamma=(\mathcal A\ni a_{\mathrm{abs}}\text{ with declared cost},\,L,\,\text{criterion},\,g,\,T\text{ tolerance},\,\eta,\,\delta)$; $\mathrm{Flags}$ = frozen partiality flags $\cup$ {invalid-certificate, empty-class, untypable-demand, tolerance-infeasible, branch-switch-proximity}; $\mathrm{Ledger}$ as below.

## 1. Identification: $I_\theta(O)$

$$I_\theta:\ \mathcal H\times\mathcal O\times\mathcal X^m\ \longrightarrow\ \mathcal K_m\times\mathcal K^{\mathrm{fin}}_m\times\mathrm{Flags},\qquad (H,O,Q)\mapsto(\widehat J,\ \widetilde J,\ \mathrm{flags}).$$
*(The mandate's $I_\theta(O)$ abbreviates fixed archive + declared closure class.)* **$V_I$:** (i) **outer:** $\widehat J\supseteq J_Q(O)$ (order projection $\widehat\Sigma\supseteq\Sigma(J_Q(O))$) under the declared closure class; (ii) **witnesses:** every element of $\widetilde J$ is $e_Q$ of an explicitly verified admissible member ($\widetilde J\subseteq J_Q(O)$; floor device only, never a feasibility report); (iii) **set semantics (audit scoping 3):** $I_\theta(H,\cdot)=I_\theta(H',\cdot)$ whenever $H,H'$ contain the same set of distinct records — permutation and duplicate-multiplicity invariance, enforced as an axiom so the feasibility channel cannot carry frequency information.

## 2. Population conditioning: $M_\phi(H,O,Q)$

$$M_\phi:\ \mathcal H\times\mathcal O\times\mathcal X^m\times\mathrm{Ctx}\ \longrightarrow\ \mathfrak Q(\Omega)\times(0,1]\times\mathrm{Rung}\times\mathrm{Tags},\qquad (H,O,Q,\gamma)\mapsto(\widehat{\mathcal Q},\,1-\delta,\,r,\,\mathrm{tags}).$$
**$V_M$:** (i) **typed dependence on $O$:** $M_\phi(H,O,Q,\gamma)=M_\phi(H,\kappa(O),Q,\gamma)$ — the current record enters *only* through the declared context map $\kappa$ (DC-C2); any finer dependence would smuggle current-task information into the population channel; (ii) **rung-typed coverage:** with probability $\ge1-\delta$ under the tagged declarations ((C-IID-$\kappa$), (SUFF-$\kappa$), transport $\rho$, union-bound family $\mathcal E$), $\widehat{\mathcal Q}$ contains the true object *of rung $r$* (DC-C3/C4); ranking classes are DC-R1 constraint polytopes on the full $\Omega_m$; (iii) no $M_\phi\to I_\theta$ edge (DE-H2/H3).

## 3. Decision: $D_\psi(J_Q,\Delta,L,\tau)$

$$D_\psi:\ \big(\mathcal K_m\times\mathcal K^{\mathrm{fin}}_m\times\mathrm{Flags}\big)\times\big(\mathfrak Q(\Omega)\times(0,1]\times\mathrm{Rung}\big)\times\mathrm{Ctx}\times\mathcal T\ \longrightarrow\ \Big(\big(2^{\mathcal A}\!\setminus\!\{\emptyset\}\big)\ \cup\ \mathrm{FailureReports}\Big)\times\mathrm{Ledger},$$
where $\mathcal T$ = declared tie-breaks (possibly absent) and a singleton output is the special case $|\cdot|=1$ (via DR-S2 uniqueness or $\tau$). **$V_D$:**
1. **rung check:** $\widehat{\mathcal Q}$ consumed only at its tagged rung (DC-C3); rung $<3$ ⇒ conditional rows suppressed or degraded to the DE-T4 endpoint;
2. **guarantee validity (audit scoping 2):** the unconditionally valid guarantee row is the **selected action's own outer risk** $\sup_{v\in\widehat J}L(\hat a,v)$ — always well-defined; $G_{\mathrm{cert}}=R_{\mathrm{set}}(\widehat J,\mathcal A,L)$ is quoted as *achieved* only with an attainment certificate (compact/l.s.c., DR-S1), else as an infimum with $+\eta$ slack (DR-S5);
3. **floor honesty (audit scoping 1):** witness floor statements are policy-typed — deterministic rules are bounded below by $R_{\mathrm{set}}(\widetilde J,\mathcal A,L)$, randomized rules by $R_{\mathrm{rand}}(\widetilde J,\mathcal A,L)\le R_{\mathrm{set}}$; "no rule" claims cite $R_{\mathrm{rand}}$; outer values are never floors (DR-F4-R);
4. **selection:** $\mathcal A^*_\eta$ or $\tau(\mathcal A^*_\eta)$; no hidden measures (DR-S3); branch-switch proximity flagged per DC-S1/S5;
5. **abstention/failure split:** abstention selected iff criterion-optimal (DC-A1); the five failure cases of DC-A4 emit $\mathrm{FailureReports}$, never actions; tolerance $T$ is a post-selection feasibility check (DC-A4(5));
6. Phase-7 honesty axioms H1–H6 carried.

**Ledger rows:** identification (declaration-invariant); bracket ($R_{\mathrm{rand}}/R_{\mathrm{set}}$ on $\widetilde J$ $\big|$ selected action's outer risk, policy-typed per $V_D$3); conditional rows with rung + axiom tags (LP brackets of DC-R5 for ranking); tolerance/tightness ($\eta$, bracket width); echo (all declarations: closure class, $\kappa$, rungs, $\mathcal E$, $\tau$, $a_{\mathrm{abs}}$-cost, $T$).

---

## 4. Closure theorems

**Theorem DC-I1 (closed composition). [proved, given the cited results]**
For $\mathbb D=D_\psi\circ(I_\theta\times M_\phi)$ with $V_I,V_M,V_D$: every emission is valid under its stated type — guarantee rows by outerness ($\widehat J\supseteq J$, definitional sup-domination); floor rows by witness membership + monotonicity, policy-typed; conditional rows by DC-C2 at rung 3 (the conditioning now *proved from declared assumptions*, not asserted) and by DC-R5's LPs on the correctly-typed polytope; selection rows by DR-S2/S3/S5 with DC-S1's corrected discontinuity scope; abstention/failure by DC-A1/A4. Degradation is monotone, loss-typed, and rung-typed: each invalidated input collapses exactly its own rows to the conservative endpoint (flags; DE-T4 minimax under the declared $(\mathcal A,L)$). No step invokes a refuted statement: the DR-S4-R claim is replaced by DC-S1; the 8.1 abstention clause by DC-A1/A4; the subset-law object by DC-R1; the marginal-as-conditional step by DC-C2–C4. $\square$

**Theorem DC-I2 (nothing left to invent). [proved as an audit of obligations]**
Every judgment the operator must make is either (a) a proved theorem cited above, or (b) an explicitly declared input slot (closure class, $\kappa$, rungs, transport, $\mathcal E$, criterion, $\tau$, $a_{\mathrm{abs}}$-cost, $T$, $\eta$, tie convention) surfaced in the echo row. The audit's list of theorem-level choices a builder would have had to invent — conditioning route, joint ranking structure, abstention-failure behavior, randomized-vs-deterministic floors, approximate attainment — is now covered by (a) DC-C2/C3, DC-R1–R5, DC-A1/A4, $V_D$3, $V_D$2 respectively. $\square$
