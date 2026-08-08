# Stopping Criterion (Final Self-Contained Theory Closure)

> **Status:** Phase-25.1 terminal decision, 2026-08-03. Sources: the four files of this folder (FC-1–FC-19, symbol index) and the audit `../25_final_fixed_theory_audit/FINAL_AUDIT.md` (`FIXED_DEPLOYMENT_THEORY_INVALID`; sole obstruction: the freeze package was not self-contained — calibration and consistency contained undefined symbols). Phases 0–24 unmodified. No theory redesigned, no claim strengthened, no continuum/ranking/varying-$z_H$ content restored. Only definitional closure and the $B\to B(z)$ notation fix performed.

---

## The obstruction, removed

The audit **passed** scope, target identity, strong convexity, and approximation, and rejected only self-containedness: calibration (finding 5) and consistency (finding 6) invoked undefined risk, excess-risk, metric, transfer-constant, hypothesis-class, empirical-objective, and estimator symbols, unverifiable from the freeze folder alone. This phase defines all of them **in-folder**:

| Mandate item | Delivered | Where |
|---|---|---|
| **1. Complete calibration definitions** | $R_\mu(F)$ (FC-7), $\mathcal E_\mu(F)=R_\mu(F)-R_\mu(g^\star_\mu)$ (FC-8), $d_{\mathbb M}$ (FC-9), transfer constant $D_V=D_V^{\mathrm{val}}\kappa_B$ **with $\|B(z)\|$ explicitly inside via $\kappa_B$** (FC-9, the audit's precise open question), $\Phi$ (FC-10); calibration restated with every symbol defined: $\|d_{\mathbb M}(F,g^\star_\mu)\|_{L^2}\le\Phi(\mathcal E_\mu(F))+2h$ (FC-11) | `calibration_closed.md` |
| **2. Complete meta-learning contract** | $P_T$, $T_i=(S_i,Q_i)$ (FC-12), $\mathcal H_N$/$\Omega_N$ (FC-13), $\widehat R_{\mu,N}$ (FC-14), $\hat\omega_N$ (FC-15), $R_\mu$ (FC-7), $\Gamma_N$ (FC-17), $\gamma^{\mathrm{opt}}_N$ (FC-18), $\varepsilon_{\mathrm{approx}}$ (FC-16), $L_p$ (FC-18); consistency restated with every symbol defined (FC-19) | `meta_learning_contract_closed.md` |
| **3. Notation fix $B\to B(z)$** | $B(\cdot)$ a fixed deployment-determined matrix-valued rule, constant on each finite context cell, pointwise linear $p\mapsto B(z)p$; $\kappa_B=\sup_z\|B(z)\|_{\mathrm{op}}$ (FC-4) | `preliminaries.md` |
| **4. Scope unchanged** | fixed deployment $\mathcal D$; continuous point-valued affinity regression only; no ranking; no continuum refinement; no varying $z_H$ (restated throughout; carried, not modified) | all files |

## Self-containment check

The symbol index (`symbol_index.md`) lists **every** symbol in the two retained theorems (FC-11 calibration, FC-19 consistency) with an in-folder definition location; **no entry is external**. The retained mathematical content is identical to the Phase-24.1 fixed-mesh results the earlier audit passed — nothing is added or strengthened; the sole change is that the definitions those results relied on are now transcribed into this package, and the $B(z)$ notation is corrected. The design floor $2h$ remains a fixed positive constant; consistency is stated as $\limsup\le2h$ (the FD-2 correction, carried).

## Retained scope limits (declared, not gaps)

No continuum-refinement guarantee; no ranking guarantee; no varying-$z_H$ guarantee. Each is a declared boundary of a complete, self-contained theory on $\mathcal D$ — unchanged from Phase 24.1.

## Verdict

$$\boxed{\textbf{FINAL\_FIXED\_THEORY\_CLOSED}}$$

Every symbol appearing in every retained theorem is defined within `25_1_final_theory_closure/`: the calibration inequality and the consistency theorem are now auditable from this folder alone, with the transfer constant carrying the assembly norm $\|B(z)\|$ explicitly, the excess risk and empirical objective written out, the estimator and hypothesis class defined, and the $B(z)$ notation made consistent. The theory content — one regularized target $g^\star_\mu$, strong convexity, approximation, calibration, fixed-deployment consistency, regression-only scope — is unchanged and was already passed; only its self-containedness, the audit's single obstruction, is now established.

**Verdict: `FINAL_FIXED_THEORY_CLOSED`.**
