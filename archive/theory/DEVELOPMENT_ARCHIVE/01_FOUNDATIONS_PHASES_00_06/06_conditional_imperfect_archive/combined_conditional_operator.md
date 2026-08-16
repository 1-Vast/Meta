# The Combined Conditional Operator (Parts IV & V)

> **Status:** Phase-6, 2026-08-03. Composes Parts I–III into a single conditional operator and derives the mandated failure/substitution tests. Sources: `auxiliary_information_theory.md` (CI-A), `irregular_archive_theory.md` (CI-B), `imperfect_archive_theory.md` (CI-C), frozen corpus. New results **CI-D**; tags **[proved] / [conditional] / [impossible] / [open]**. The composition-commutation and no-coupling hypotheses are the refereed ones.

---

## IV. The combined conditional operator

**Definition CI-D1.** For a new member with known auxiliary $c_b$ and support $S_b$ (design $D$, data $\tilde y$, noise $\varepsilon$), query $x$:
$$A_c(\text{archive},\,c_b,\,S_b,\,x,\,\varepsilon)\;:=\;\text{the canonical base operator applied to the }c_b\text{-fiber of the union family},$$
i.e. center $+$ radius of $S_\varepsilon(\tilde y\mid c_b)$ taken in $\big(\mathcal F_{\mathrm{union}}\big)_{c_b}$, where:
- **archive noise** enters via the union family $\mathcal F_{\mathrm{union}}=\bigcup_{W\in\mathcal W_{\mathrm{arch}}}\mathcal F_W$ (Part III, CI-C1);
- **irregular pattern** enters via CI-B's identification conditions determining $\mathcal W_{\mathrm{arch}}$ (which windows are consistent, and to what residual ambiguity);
- **auxiliary information** enters via the fiber $(\cdot)_{c_b}$ (Part I, CI-A1); and $\mathcal W_{\mathrm{arch}}$ is computed in the **augmented (labeled) class** so archive labels $c_a$ fold into the $\delta_a$-consistency test.

**Theorem CI-D2 (well-posed composition and optimality). [proved under well-specification + no-coupling + augmented-class $\mathcal W_{\mathrm{arch}}$]**
(i) $A_c$ is minimax-optimal, conditional on the archive, by composition of the three reductions — no new optimality proof: Part I makes conditioning a fiber restriction, Part III makes archive noise a union, and Theorem 1 applies to any nonempty family, so it applies to the union-fiber family $(\mathcal F_{\mathrm{union}})_{c_b}$.
(ii) **Fiber and union commute** (refereed set identity): fibering is intersection with the label slice $\{c=c_b\}$, and $\big(\bigcup_W \mathcal F_W\big)\cap\{c=c_b\}=\bigcup_W\big(\mathcal F_W\cap\{c=c_b\}\big)$; $c_b$ does not enter archive consistency (archive members carry their own labels inside $\mathcal W_{\mathrm{arch}}$'s definition), so union-then-fiber $=$ fiber-then-union. Candidates whose $c_b$-fiber is empty drop out automatically — **this is the misspecification-detection channel**. Requires labels intrinsic to candidate members (guaranteed by the augmentation reduction).
(iii) The irreducible remainder is exactly $\tfrac12\,\omega^{\mathrm{union}}_{x,D}(2\varepsilon\mid c_b)$ (Theorem 1 on the union-fiber family), conditional on that family being nonempty; the empty branch reports detected inconsistency, no certificate.

**Exact information ledger. [proved]** For the combined problem:

| Source | Supplies | Certificate |
|---|---|---|
| **Archive** | the consistent set $\mathcal W_{\mathrm{arch}}$ of augmented, truncated (size $\le k+1$), covered-region window systems — identified under CI-B's pattern conditions, to residual ambiguity $h$ under CI-C | CI-B1/B5, CI-C1 |
| **$c_b$** | fiber selection within each candidate ($\mathcal F_W\to(\mathcal F_W)_{c_b}$) **plus** misspecification detection (empty fiber) | CI-A1, CI-A3 |
| **$S_b$** | the section cut — at most $k\le5$ continuous dimensions of member identity ($2k$ inequalities cutting the section) | F20/CP-3, DM-3 |
| **$x$** | selection of the window and the validity flag (covered? in the identifiable region of each $W$?) | F18, O4 |
| **Irreducible remainder** | exactly $\tfrac12\,\omega^{\mathrm{union}}_{x,D}(2\varepsilon\mid c_b)$ | Theorem 1 on the union-fiber family |

**No hidden interaction (refereed).** The three reductions compose without cross-terms under the stated hypotheses: the new member's label does not feed back into archive consistency; union commutes with fiber; and conditioning commutes with the section operation. The only couplings that would break exact optimality (member choices constrained jointly with archive realizations) are absent in the exactly-$d$ linear class and are flagged as the no-coupling hypothesis for the general-class statement — where the union operator remains **valid** (outer) but exact optimality may weaken.

---

## V. Failure and substitution tests

Each test names a decisive quantity that **must change** if the corresponding source is genuinely informative — a falsification protocol at the level of the operator, not of any implementation.

**Test T1 — replace $c_b$ by an unrelated $c'$. [proved]**
Decisive quantity: the **realizability flag** = emptiness of $S_\varepsilon(\tilde y\mid c')$ (CI-A5, refereed). If $c'$ is genuinely wrong and informative, the flag fires with positive frequency over realizable data; in the silent regime the reported interval still lies within $\omega_{x,D}(2\varepsilon\mid c',c_b)$ of the truth. **What should change:** the fired-flag rate, and (when silent) the interval location by up to the cross-fiber modulus. **What should *not* change if $c$ is useless:** nothing — by CI-A3, uselessness $\iff$ $T_{D\cup\{x\}}(c')=T_{D\cup\{x\}}(c_b)$, in which case substitution is provably harmless.

**Test T2 — remove $c_b$. [proved]**
Decisive quantity: the **radius**. It must increase from $\tfrac12\omega^{\mathrm{union}}(2\varepsilon\mid c_b)$ to $\tfrac12\omega^{\mathrm{union}}(2\varepsilon)$; the increase is exactly the conditional information $\Gamma_c\ge0$ (CI-A, determination 4). If removing $c_b$ leaves the radius unchanged, $c_b$ was worthless at $(D,x)$ in the worst case (useful-iff, CI-A3).

**Test T3 — replace $S_b$ by observations from another member. [proved]**
Decisive quantity: the **realizability flag and the section**. Foreign data $\tilde y'$ from a member $f'$ in a *different* fiber or a different consistent window either fires the emptiness flag (if $\tilde y'\notin$ the $\varepsilon$-tube of $(\mathcal F_{\mathrm{union}})_{c_b}$) or shifts the center to $f'$'s value — detectable against the genuine member's value whenever the two members differ at $x$ beyond the noise floor. If the section is unchanged by the swap, the two members are indistinguishable at $(D,x,\varepsilon)$ (sub-resolution, base theory) — a certified equivalence, not a bug.

**Test T4 — break the archive overlap structure. [proved]**
Decisive quantity: the **residual ambiguity dimension** of $\mathcal W_{\mathrm{arch}}$ (hence the radius, on the affected covered points). Removing a member so some point falls below $d$ observers (CI-B2), or downgrading a unisolvent overlap to rank $<d$ (CI-B4b holonomy), enlarges $\mathcal W_{\mathrm{arch}}$ by a realized affine/gauge freedom — the union sections widen at exactly the points that lose identifiability. If breaking overlap leaves all radii unchanged, that overlap carried no identifying information (its excess was redundant — the CI-B counting slack).

**Test T5 — increase archive uncertainty $\delta$. [proved]**
Decisive quantity: the **radius**, which is nondecreasing in $\delta$ (CI-C4). Under $\sigma_0$-conditioning the increase is $O(\delta)$ (CI-C2); near degeneracy it can be $\Theta(\delta/\sigma_{\min})$ or $\kappa^L$ (CI-C3) — a diverging radius under small $\delta$ is the operational signature of ill-conditioned transport, not of a wrong operator. If the radius is insensitive to $\delta$, the archive constraints were not binding (saturation, base theory).

**Common structure of the tests.** Every genuine information source has a **certificate quantity** whose change is forced when the source is perturbed: fiber $\to$ flag/radius (T1,T2), support $\to$ section (T3), pattern $\to$ ambiguity dimension (T4), archive noise $\to$ radius (T5). A source whose perturbation moves *no* certificate quantity is, by the corresponding useful-iff / saturation theorem, provably non-informative at that configuration — the tests are decisive in both directions.

---

## Ledger

| Result | Statement | Tag |
|---|---|---|
| CI-D1/D2 | combined operator well-posed; optimal by composition; fiber/union commute; exact ledger | **proved** (well-spec + no-coupling + augmented class) |
| T1–T5 | decisive certificate quantity per source, forced-to-change iff informative | **proved** |
| — | exact optimality under member/archive coupling (general non-linear classes) | **open (non-blocking; validity holds)** |

**Net.** There is a mathematically justified conditional operator $A_c(\text{archive},c_b,S_b,x,\varepsilon)$: it is the base canonical operator on the union-fiber family, optimal by composition, with an exact five-way information ledger and five decisive substitution tests. Nothing new had to be invented — the three additions are three reductions to Theorem 1.
