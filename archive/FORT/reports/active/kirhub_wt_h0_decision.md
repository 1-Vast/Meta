# KirHub WT H0 decision

Verdict: `H0_PASS_AUTHORIZE_FROZEN_PROTEIN_NECESSITY_PROBE`.

The 409×92 WT matrix is adequately powered for a low-capacity protein-necessity probe under strict
family/scaffold isolation. The label-blind registry maps 408/409 preparations to 95 KLIFS families
and all 92 inhibitors to 87 Bemis–Murcko components. After 5–95% saturation filtering, 358 genes
remain across 92 represented family components; 46.22% of cells are usable. Every one of the 25
balanced family-fold × scaffold-fold combinations retains at least 579 cells and 59 target profiles
with five query ligands. Family-macro MDE80 is +0.0292 at paired SD 0.10.

The leave-family-out global-ligand rank has macro Spearman 0.3924
[0.3713, 0.4136], but median rho-squared is only 0.1858 and the residual-energy fraction is 0.8336.
Residual profiles are protein-family structured: within-family minus same-group cross-family
Spearman is +0.2860 [0.2370, 0.3364] over 74 family units; the exposure- and family-size-preserving
family-label permutation gives p=0.0005.

The attachment-requested null based on permuting target labels against total residual/cycle energy
was replaced because that energy is invariant to row-label permutation. The valid inferential null
permutes KLIFS family labels within kinase groups while preserving family counts. This is a
measurement correction, not a model revision.

H0 authorizes only a frozen-protein necessity probe. It does not authorize a neural encoder,
Transformer, posterior, support mechanism or multi-seed training.

