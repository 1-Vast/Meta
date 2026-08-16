# Input Typing Repair (Item 1)

> **Status:** Phase-22.1 (target alignment minimal repair), 2026-08-03. Phases 0–21 unmodified; the operator is not redesigned, no ranking theory is added, no generality is expanded. Audit of record: `../22_target_alignment_audit/FINAL_AUDIT.md` (`TARGET_ALIGNMENT_INVALID`). This file resolves the `z_H` ambiguity the audit flagged (finding 2). New/repaired results carry **-R** suffixes. Tagged **[declared] / [proved]**.

---

## The choice: **A — $z_H$ is a declared deployment constant**

$$\boxed{\ g^\star(z),\qquad z=z(S_T,Q_T,\gamma)\in Z,\qquad z_H\ \text{fixed by the deployment (SKEL).}\ }$$

**Declaration DT-A (deployment constant). [declared]**
The deployment state $z_H$ — trained parameters, fiber counts, per-context population bands $b^{\mathrm{pop}}$, anchors $b_0,\dots,b_m$ — is **frozen once at deployment** and is **not** an argument of the statistic $z$, of the target $g^\star$, or of any theorem in this phase. Every result below is a statement about **one fixed $z_H$**: the theory defines and learns *one* risk-optimal operator *per declared deployment*.

**What is therefore claimed, and what is not (audit finding 2, resolved without overclaim). [proved trivially / declared]**
- *Claimed:* for each fixed $z_H$, a single risk-optimal target $g^\star=g^\star_{z_H}$ on the support/query-conditioned statistic $z(S_T,Q_T,\gamma)$ (so this is genuinely support-and-query conditioned, not support-ignoring regression — the audit's own concession).
- *Not claimed:* one operator conditioned on a **varying** observable $(z_H,S,Q,\gamma)$. Option B (folding $z_H$ into the statistic) would be a *different, larger* object; the audit correctly noted no equivalence between fixing $z_H$ and including it was ever proved, and this phase **does not assert one**. A varying-$z_H$ operator is explicitly out of scope here.

**Consistency note.** Because $z_H$ (hence $b^{\mathrm{pop}}$ and the anchors) is a fixed constant, the assembly's basis points $\beta_0=b^{\mathrm{pop}},\beta_1=b_1,\dots,\beta_m=b_m$ are fixed vectors — a fact used essentially in `convex_parameterization.md` to obtain a genuinely linear (not merely bilinear) assembly. Fixing $z_H$ is thus not only the honest typing choice but the one that makes the strong-convexity repair available.
