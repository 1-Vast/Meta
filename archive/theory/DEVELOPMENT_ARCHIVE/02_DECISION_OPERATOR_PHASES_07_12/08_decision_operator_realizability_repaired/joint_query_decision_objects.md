# Joint Query Decision Objects — REPAIRED (Part III)

> **Status:** Phase-8.1, 2026-08-03. Supersedes `../08_decision_operator_realizability/joint_query_decision_objects.md`. The audit **passed** this file; the sole repair is the certificate re-typing of one sentence in DR-J3(iii) (T1 propagation). DR-J1, DR-J2, DR-J3(i)–(ii), DR-J4 and the representation constraints are carried verbatim by citation; the corrected clause and the audit's own added witness are recorded below.

---

## 1. Carried results (audit: pass)

- **DR-J1 [carried].** Pushforward sufficiency per loss class: if $L(a,v)=\ell(a,g(v))$, everything factors through $g(J)$.
- **DR-J2 [carried].** Strict quotient lattice: marginals $\nrightarrow$ differences $\nrightarrow$ order set $\nrightarrow$ $J$; none recoverable from a coarser one. The audit's compact witness is adopted into the record: $J_{\mathrm{diag}}=\{(0,0),(1,1)\}$ vs $J_{\mathrm{anti}}=\{(0,1),(1,0)\}$ — identical marginal intervals; the first always ties the two queries, the second admits either strict order and never ties.
- **DR-J3(i)–(ii) [carried].** Listwise sufficiency of $\Sigma(J)\subseteq S_m$; strict over-admission of the pairwise-compatible set $\Sigma^{\mathrm{pair}}\supseteq\Sigma$, witness $J=\{(0,1,2),(2,1,0)\}$: all pairwise signs bilateral, $\Sigma^{\mathrm{pair}}=S_3$, $\Sigma(J)=$ the two reversals.
- **DR-J4 [carried].** Minimal decision-sufficient object for the ranking application: $\Sigma(J_Q(O))$ per queried $Q$ (+ the declared population pushforward on $S_m$ when preferences beyond minimax are demanded; + magnitudes only for graded losses).

---

## 2. DR-J3(iii)-R — the over-admission cost, correctly typed

On the reversal witness, with Kendall discordance loss:
- computed on the **exact** order object $\Sigma(J)$: $R_{\mathrm{set}}=2$ (triangle bound $\lceil 3/2\rceil$, attained by the middle orders) — this is the **true information floor** of the listwise problem (DR-F2 applies: both reversal members are admissible);
- computed on the **outer proxy** $\Sigma^{\mathrm{pair}}=S_3$: the value $3$ is $G_{\mathrm{cert}}$ — a **conservative surrogate guarantee** (DR-F4-R(b)), valid as an upper bound on what the proxy-using system can promise, **not** "a floor of 3": no admissible member forces discordance $3$ against the middle orders. *The Phase-8 sentence calling the proxy value a conservatively-reported floor is retracted; the inclusion $\Sigma\subseteq\Sigma^{\mathrm{pair}}$ and the action-set discrepancy stand.*

The witness therefore now exhibits, in one example, the entire repaired type discipline: exact object → floor $2$; outer proxy → guarantee $3$; an inner witness pair (the two reversal members, explicitly constructible) → certified floor lower bound $2$, closing the bracket exactly: $[2,3]$ reported honestly, with the gap = the declared price of carrying only pairwise structure.

---

## 3. Representation constraints (carried, re-typed)

Outer semantics for orders is unchanged as a **feasibility/validity** rule: $\widehat\Sigma\supseteq\Sigma(J_Q(O))$ under the declared class; a spuriously excluded order is a false Tier-1 certificate. What outer $\widehat\Sigma$ buys is re-typed per DR-F4-R: valid **guarantees** and valid **Tier-1 checks** (outer singleton $\Rightarrow$ true singleton), never floor values. Floor values for ranking cite witness pairs/tuples (DR-F5-R): two admissible members realizing opposite signs certify the pairwise floor $\min(1,c)$ under the declared $(\mathcal A,L)$ — loss-typed per DR-F1-R(iii). Size discipline (pairwise data + joint realizability predicate; $\Sigma^{\mathrm{pair}}$ as the always-available outer proxy, flagged non-tight) carried unchanged.
