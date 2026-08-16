# Stopping Criterion (Target Alignment Minimal Repair)

> **Status:** Phase-22.1 terminal decision, 2026-08-03. Sources: the five repair files (DT-A, RP-1–RP-10) and the audit `../22_target_alignment_audit/FINAL_AUDIT.md` (`TARGET_ALIGNMENT_INVALID`). Phases 0–21 unmodified; operator not redesigned; no ranking theory; no expanded generality. The repair touches exactly the six items the mandate named, and nothing else.

---

## The six repairs, item by item

| Item | Audit finding | Repair | Status |
|---|---|---|---|
| **1. Input typing** | $z_H$ neither an argument nor declared constant | **Route A chosen:** $z_H$ is a **declared deployment constant**; every theorem is about one fixed deployment; the varying-$(z_H,S,Q,\gamma)$ operator is explicitly out of scope, with no false equivalence asserted (DT-A) | resolved |
| **2. Strong convexity** | assembly bilinear in $c=(\lambda,w)$; ridge did not guarantee $\mu$-strong convexity ($(\lambda w-a)^2$ witness) | **Barycentric reparameterization** $p\in\Delta_m$, $\mathsf{asm}(p)=Bp$ **linear** (RP-1); convex band loss ∘ linear map is convex, ridge gives modulus **exactly $\mu$** (RP-2); measurability fixed by including A-CONT in the hypotheses (RP-3) | proved |
| **3. Continuity** | linear argmin modulus $\varpi_\ell/\mu$ false | **Square-root modulus** $\|g^\star(z)-g^\star(z')\|\le\sqrt{2\varpi_\ell/\mu}$ adopted as the theorem (RP-4); linear rate available **only** under an explicitly declared gradient assumption (A-GRAD), not default (RP-4.1) | proved / corrected |
| **4. Approximation** | interpolation rate used the false modulus | Restated at the square-root modulus (RP-5); $\varepsilon_{\mathrm{approx}}\to0$ survives (continuity suffices for convergence; only the rate changes), and nothing downstream needed the linear rate (RP-6) | proved |
| **5. Calibration** | mesh floor omitted from $\Phi$; $\Phi(t)\to0$ overclaimed | Floor written **additively and separately**: $\|d_{\mathbb M}(F,g^\star)\|_{L^2}\le\Phi(\mathcal E(F))+\varepsilon_{\mathrm{design}}$, $\Phi(t)=D_V\sqrt{2t/\mu}\to0$, total $\to0$ **only if** $\varepsilon_{\mathrm{design}}=2h\to0$ (RP-8); end-to-end consistency on an **explicit joint sieve schedule** $\dim\Omega_N\ln N/N\to0,\ h_N\to0$ (RP-9) | proved / conditional on the declared schedule |
| **6. Scope** | end-to-end scope overclaimed; no ranking theory present | **Explicit restriction:** continuous point-valued affinity regression only; **no ranking guarantees**; ridge target is the regularized regressor; ranking stays a separate Route-A objective, not derived (RP-10) | declared |

## The minimal obstruction, removed

The audit's *minimal* obstruction was PT-6's false Lipschitz modulus, used by PT-9 and the chain. It is removed: RP-4 proves the square-root modulus (the audit's own bound) and marks the linear rate as requiring a declared gradient hypothesis; RP-5 restates approximation at that modulus; RP-8/RP-9 carry the honest additive floor and the explicit sieve schedule so the chain's vanishing is genuine. The independent secondary issues (z_H typing, non-affine assembly, Route-B floor, sieve growth, ranking scope) are each addressed above.

## Residual, named (not defects)

Square-root (not linear) approximation rate unless (A-GRAD) is declared; a positive operator-error floor at any fixed mesh $h$, vanishing only along $h\to0$; consistency conditional on the declared sieve schedule; single-deployment ($z_H$-constant) scope; regularized target under a ridge; regression-only, ranking excluded. All declared; none is a hidden target or an unproved regularity claim.

## Verdict

$$\boxed{\textbf{TARGET\_ALIGNMENT\_REPAIRED}}$$

The single risk-optimal target survives on a corrected footing: genuinely strongly convex in barycentric coordinates (so uniqueness and calibration are sound), continuous at the honest square-root modulus (so approximation converges, at the correct rate), calibrated with the mesh floor written explicitly and a joint sieve schedule making the end-to-end error vanish, and scoped in writing to continuous affinity regression with no ranking claim. Every one of the audit's six items is repaired without redesigning the operator, adding ranking theory, or expanding generality.

**Verdict: `TARGET_ALIGNMENT_REPAIRED`.**
