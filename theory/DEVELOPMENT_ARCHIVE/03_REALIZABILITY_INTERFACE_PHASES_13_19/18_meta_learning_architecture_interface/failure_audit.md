# Failure Audit (§6)

> **Status:** Phase-18, 2026-08-03. Mandated falsification pass over the five named risk areas. Results **MI-14–MI-18**, tagged **[attack refuted] / [attack found — repaired] / [attack found — scoped]**.

---

## 1. Representation sufficiency

**Attack.** $r(S)=(b_{\mathrm{can}},z)$ discards within-task information that some operator outside the declared family could exploit (e.g. fine geometry of the observation configuration beyond the canonical bands and context).
**Outcome: [attack found — scoped, provably tight].** True — and exactly quantified: MI-6 proves $r$ is minimal-sufficient *for the declared family*; sufficiency for a richer family is neither claimed nor possible without enlarging the family first (which is a new declared skeleton/deployment). The scoping is not a loophole: the family was derived in Phase 17 as the minimal valid object, so information beyond $r$ is information the *valid* operator class provably cannot consume without new theory. A second probe — does $r$ include everything the family needs? — verified: $\varepsilon$ and the auxiliary label enter $b_{\mathrm{can}}$ and $\kappa$ respectively; the index arguments $(Q,\gamma)$ are operator arguments; nothing consumed is outside $r$ and the declared arguments. $\square$

## 2. Finite dimension

**Attack.** The fixed-$p$ claim hides growth: refine accuracy, grow the skeleton, $p$ grows.
**Outcome: [attack found — scoped, with the impossibility already proved].** Carried verbatim from MR-4/MR-5: no fixed $p$ across resolutions (proved impossible); fixed $p$ within a declared skeleton (proved); the relaxation is the weakest possible. The interface adds nothing new to attack: the coefficient set $\mathcal C$, the latent dimensions of MI-8's poolings (bounded via $k\le5$), and $p$ are all skeleton-constants. No silent leakage found beyond the declared, verdict-qualifying relaxation. $\square$

## 3. Hidden leakage

**Attack battery.** (a) Task identity smuggled through the representation; (b) query answered from population instead of identification (marginal-to-conditional leaps); (c) training loss reading latent marks; (d) $\kappa$ reading unidentified quantities; (e) history multiplicities leaking into the identification channel.
**Outcome: [attack refuted, all five].** (a) $r$ contains no task index; MC-7(ii)'s exclusion is structural (no slot exists). (b) Population values enter only rung-typed and support-restricted; the rung ladder makes the leap a type violation (MC-14/ML-K). (c) The loss consumes identified query information only (MI-10; censored tasks score against compatible regions — observable). (d) $\kappa$ is a declared function of the observable record (MC-11's $\kappa$-DESIGN). (e) The set/multiset channel typing (MC-4, MI-2) is enforced at the input type. Each check is a re-verification of a standing theorem against the *new* interface surface; none failed. $\square$

## 4. Non-identifiability

**Attack.** (a) $\theta$-gauge: many parameters, one operator; (b) representation gauge: many poolings, one $r$; (c) elicitation masquerade: calibrated bands presented as identified bands.
**Outcome: (a),(b) [attack found — harmless by design]:** only extensional operator values are semantic (MC-7, MR gauge discipline); the convex program's minimizing face makes gauge explicit rather than hidden. **(c) [attack found — repaired by standing firewall]:** MI-12(ii)/MR-13 — certificates are $\theta$-invariant, calibration statements rung-tagged; the masquerade is unrepresentable in the output schema (no slot exists in which a learned band can be emitted as a certificate). $\square$

## 5. Task distribution assumptions

**Attack.** The generalization theorem (MI-11) silently assumes what DE-T3 forbids assuming: population stability.
**Outcome: [attack refuted — the assumptions are loud, not silent].** MI-11 is tagged (IID)/(C-IID-$\kappa$) with fiber counts and inherits the missing-fiber term; the transport ladder (DE-T2) prices declared shift; undeclared shift triggers the proved adversarial reversal (DE-T3) and the graceful degradation to the frozen minimax endpoint — a failure mode with a theorem, a flag, and a fallback, not an unexamined premise. The echo row carries every tag. $\square$

## Verdict input

Two findings survive as declared scope (family-relative sufficiency; skeleton-relative finiteness — both proved tight), one as a standing typed firewall (elicitation). No attack invalidates the interface; none requires a new object.
