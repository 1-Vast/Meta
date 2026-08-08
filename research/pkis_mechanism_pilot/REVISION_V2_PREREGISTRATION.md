# PKIS mechanism pilot: one permitted representation revision (v2)

Status: **FROZEN BEFORE ANY V2 SCORE WAS COMPUTED**  
Date: 2026-08-08  
Formal status: external exploratory falsification only; this is **not** E-AFF-X1,
does not authorize a biological `z`, and cannot change a historical Gate.

## 1. Why exactly one revision is allowed

The preregistered v1 statistic failed on the PKIS1 -> PKIS2 dual-cold transfer:
its correct-protein arm did not beat the additive, zero-interaction, or deranged-
protein arms.  Its ligand representation compressed each molecule to global
pharmacophore counts in four Murcko-distance shells.  The registered diagnosis
is therefore a loss of *local chemical environment identity*, not evidence that
a larger generic encoder, attention block, orientation module, or frozen
mathematical change is warranted.

V2 makes one and only one change: it replaces the v1 linear outer-product tensor
with a five-channel, local-environment product kernel.  V1 and all of its
negative artifacts remain immutable.

## 2. Read firewall and data roles

| Data | Role in v2 | Outcome access before this freeze |
|---|---|---|
| PKIS1 | sole training/source panel | already consumed |
| PKIS2 | consumed development transfer; same fixed v2 is reported on it | consumed by v1 |
| Anastassiadis 2011 | external cross-assay transfer | schema/identifiers and a small printed header preview only; no outcome statistic, fit, contrast, or selection |
| ChEMBL X1 labels | forbidden | zero |
| DAVIS | forbidden | zero |
| recipient/support query labels | forbidden | zero |

The Anastassiadis workbook is the official PMC Supplementary Table 3.  Its
frozen SHA-256 is
`cd756bf2b6ad541a1781508c563caf0da6da876dfb71f2546fbff02e13d98684` and its
publisher-supplied MD5 is `1c8e90029491e8081c826675fefec23a`.

The schema inspection established only that the workbook contains 178 compounds
(CAS identifiers), 300 assay constructs, and percent remaining catalytic
activity.  It also established identifier-only feasibility: 184 constructs map
unambiguously to one human 85-position KLIFS pocket; 63 of those HGNC targets
are absent from the PKIS1 source and 20 are also KLIFS-family-cold.  Chemical
structure resolution and scaffold-cold counts are still preconditions, not
guaranteed outcomes.

## 3. Biological representation

For each ligand, RDKit assigns atoms to five mechanistic channels and emits a
radius-2 Morgan bit vector whose centers are restricted to those atoms:

1. hydrogen-bond channel: donor or acceptor atoms;
2. ionic channel: positive- or negative-ionizable atoms;
3. aromatic channel: aromatic atoms;
4. hydrophobic channel: hydrophobe/lumped-hydrophobe atoms;
5. steric channel: all heavy atoms.

For each protein, the 85 aligned KLIFS pocket positions are encoded with the
public SiteAlign physicochemical table.  Channel-specific vectors are:

1. residue donor and acceptor capacity;
2. positive and negative charge capacity;
3. aromatic capacity;
4. aliphatic capacity;
5. normalized side-chain size and free-space proxy.

No target, family, assay, compound, CAS, or dataset identifier enters a feature
vector.  KLIFS family/group is used only for splitting, controls, and reporting.

## 4. Five product kernels

For channel \(c\), define a ligand Tanimoto kernel and an aligned-pocket RBF
kernel

\[
 k^L_c(L,L')=\operatorname{Tanimoto}(\ell_c(L),\ell_c(L')),\qquad
 k^P_c(P,P')=\exp\{-d_c(P,P')/\tau_c\},
\]

where \(d_c\) is mean squared distance across the declared 85-position channel
vector and \(\tau_c\) is the strictly positive median nonzero source distance.
The pair kernel is the tensor product

\[
 k_c((L,P),(L',P'))=k^L_c(L,L')\,k^P_c(P,P').
\]

This construction says that a source interaction may transfer only when *both*
the relevant local ligand environment and the aligned protein chemistry are
similar.  It is permutation invariant, pair local, bounded, identifier free,
and has the same five biological coordinates on every panel.

## 5. Source-only estimator

Let \(Y\) be PKIS1 activity on \([0,1]\), and let

\[
 R=Y-\bar Y_{L\cdot}-\bar Y_{\cdot P}+\bar Y
\]

be the source double-centered interaction residual.  Each channel fits
separable kernel ridge regression.  With eigendecompositions
\(K^L_c=U_c\operatorname{diag}(s_c)U_c^\top\) and
\(K^P_c=V_c\operatorname{diag}(t_c)V_c^\top\), its coefficient matrix is

\[
 C_c=U_c\left[\frac{U_c^\top R V_c}{s_ct_c^\top+\lambda_c}\right]V_c^\top.
\]

Prediction on a new rectangle is

\[
 G_c(L_*,P_*)=K^L_c(L_*,L)\,C_c\,K^P_c(P,P_*).
\]

The regularization grid is frozen as
\(\{10^{-2},10^{-1},1,10,10^2,10^3,10^4\}\).  Selection uses only PKIS1 and
three deterministic dual-cold folds: validation cells have both a held-out
generic Murcko scaffold fold and a held-out KLIFS group fold.  Ties choose the
largest \(\lambda\).

Out-of-fold predictions of the five channels are combined by weights
\(w_c\ge0,\ \sum_c w_c\le1\), minimizing squared error on the union of the
dual-cold validation cells.  The zero vector is feasible.  This six-vertex
convex constraint prevents a weak channel from manufacturing scale and keeps
the biological mechanism dimension at five.

The five signed contributions are

\[
 q_c(L,P)=w_cG_c(L,P),\qquad
 z_c=\tfrac12\{1+\tanh[2(q_c-a_c)/(b_c-a_c)-1]\},
\]

where \(a_c,b_c\) are the source 1st and 99th percentiles, with a fixed
unit-span fallback.  Thus \(z_{\rm bio}\in[0,1]^5\).  These coordinates are
candidate section coordinates for the frozen law-valued operator, not a
declared energy decomposition and not yet admitted to `model/config.py`.

## 6. Nuisance and control arms

Population, ligand-only, protein-only, and additive nuisance predictors use the
same source-only ridge protocol as v1.  The interaction estimator is added only
in the diagnostic correct and deranged arms:

\[
 \widehat Y_{\rm correct}=\widehat Y_{\rm additive}+\sum_c q_c(L,P),\qquad
 \widehat Y_{\rm wrong}=\widehat Y_{\rm additive}+\sum_c q_c(L,\pi(P)).
\]

The derangement \(\pi\) is deterministic and within KLIFS group wherever at
least two eligible targets exist; it has no fixed point.  There is no query-label,
target-ID, assay-ID, or scalar bypass.  Raw scalar values are scored only as a
falsification diagnostic before any law-interface admission.

## 7. Transfer construction

PKIS2 uses the original continuous activity scale.  Anastassiadis reports
percent remaining catalytic activity, so its prespecified target variable is

\[
Y=\operatorname{clip}(1-\text{remaining activity}/100,0,1).
\]

For each transfer panel:

- retain only unambiguous wild-type human KLIFS mappings;
- retain targets whose HGNC gene is absent from PKIS1;
- resolve structures from the workbook CAS identifiers through PubChem and
  cache exact CID/SMILES responses;
- retain ligands with neither exact canonical-SMILES overlap nor generic Murcko-
  scaffold overlap with PKIS1;
- require at least 100 ligands and 20 targets before any score is interpreted;
- score finite cells only and cluster uncertainty by target.

The primary stratum is exact-target-cold plus scaffold-cold.  KLIFS-family-cold
is secondary and requires at least 20 targets.  No group-cold claim is made when
that stratum is unavailable.

## 8. Frozen estimands and pass rule

Interaction admission requires all of the following in the primary stratum:

1. ligand-only positive-control MSE reduction over population has 95% target-
   bootstrap lower confidence bound (LCB) above zero;
2. correct interaction residual beats zero interaction, LCB \(>0\);
3. correct interaction residual beats deranged protein, LCB \(>0\);
4. target-macro interaction Pearson correlation, LCB \(>0\).

Location admission separately requires correct raw prediction to beat each of
population, ligand-only, protein-only, additive, and deranged arms in target-
macro MSE, with every 95% target-bootstrap LCB above zero.

There are 10,000 deterministic bootstrap draws, seed 20260808.  No p-value,
fold, arm, family, or metric may substitute for a failed registered comparison.
The external cross-assay result is considered a candidate solution only if the
Anastassiadis primary stratum passes the corresponding rule.  PKIS2 success
alone is insufficient.  Any failure produces `REVISION_V2_NOT_VALIDATED`, keeps
`admission_to_biological_z_authorized=false`, and ends the permitted revision.

## 9. Mathematical interface and abstention

If and only if a biological statistic later passes a formal Gate, its five
bounded coordinates may be registered as a declared view of the existing CSMO
input.  The frozen map remains

\[
F(z)\in\Delta_m,\qquad B(z)F(z)\in\text{Band},\qquad
A(F,z)=K(B(z)F(z)),
\]

so simplex membership, positive ridge, Band feasibility, and the law-valued
output are unchanged.  Coverage is reported from nearest source ligand
Tanimoto and nearest source KLIFS-pocket similarity; unsupported pairs route to
the existing broad/abstention law rather than receive an invented point value.

This pilot does not claim a valid interval, minimax radius, end-to-end DTA model,
non-kinase generalization, or \(k\le5\) adaptation.
