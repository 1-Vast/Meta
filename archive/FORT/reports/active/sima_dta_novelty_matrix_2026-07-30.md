# SIMA-DTA Bounded Novelty Matrix

**Audit scope:** primary sources named in the transition instruction, plus
directly relevant foundational sources. This audit bounds claims; it does not
claim exhaustive literature coverage or establish performance.

## Primary-Source Audit Record

| Method/family | Primary source checked | Relevant prior art | SIMA-DTA claim boundary |
| --- | --- | --- | --- |
| Mamba | Gu and Dao, arXiv:2312.00752 | Selective state-space long-context sequence model | Hybridization alone is not novel |
| Mamba-2 | Dao and Gu, arXiv:2405.21060 | Structured state-space duality and efficient sequence modeling | Do not claim new SSM principle |
| Mamba-3 | Bounded search did not verify a canonical primary source under this exact name | Name cannot support a novelty claim until source/version is pinned | No claim |
| Jamba | Lieber et al., arXiv:2403.19887 | Transformer-Mamba hybrid with MoE | Hybrid protein adapter is an implementation choice |
| Samba/hybrid SSM studies | Bounded audit found hybrid SSM-Transformer work; exact selected source must be version-pinned before FSA-B0 | Attention plus SSM control families exist | Require matched pure Mamba and Transformer controls |
| ProtMamba | Exact primary source/version was not verified in the bounded audit | Protein SSM/PLM applications are relevant prior art | No protein-SSM novelty claim |
| LC-PLM | Exact primary source/version was not verified in the bounded audit | Long-context protein language modeling is relevant prior art | No long-context PLM novelty claim |
| MetaDTA | Few-shot DTA meta-learning family; source pin required before baseline execution | Target-wise support adaptation exists | Include equal-budget baseline |
| CML | Acronym is ambiguous; no claim is based on it | Contrastive/cross-modal learning is established | Do not claim CML novelty |
| AdaMBind | AdaMBind release, doi:10.6084/m9.figshare.30963823.v1 | DTA target adaptation with 5 or 40 supports | Compatible k=5 baseline, not strict-cold proof |
| Neural processes | Garnelo et al., arXiv:1807.01622 | Set-conditioned latent task representations | Set encoding is not novel |
| In-context few-shot molecular prediction | Bounded family audit | Support-conditioned molecular prediction is established | No generic few-shot novelty claim |
| Active meta-learning/design | Classical optimal design and active meta-learning families | Information-based support selection exists | Do not claim D/V-optimality novelty |
| Low-rank/hypernetwork adaptation | Ha et al., arXiv:1609.09106; Hu et al., arXiv:2106.09685 | Hypernetworks and low-rank adaptation are established | FiLM/low-rank adaptation is not novel |

## Novelty Matrix

| Method | Protein backbone | Adapted parameters | Support encoding | Support selection | Calibration separation | Wrong-support objective | Scaffold-cold query | Exact nested null | Efficiency claim |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Standard MAML | Varies | Broad inner-loop parameters | Support labels/features | Usually random | Not required | Usually absent | Varies | Varies | Not primary |
| MetaDTA family | DTA encoder | Task/head or broader parameters | Target support | Usually random | Often not explicit | Often absent | Must be audited | Varies | Not primary |
| AdaMBind | DTA encoder | Target adaptation | 5 or 40 supports | Reported protocol | Not sufficient for this task | Not primary | Not closed in its protocol | Not this null | Not primary |
| Neural process family | Varies | Latent task code | Permutation-invariant set | Usually random | Optional | Optional | Varies | Optional | Not primary |
| Hybrid Mamba-Transformer | Sequence backbone | Backbone or adapters | Not inherently support-conditioned | Not inherent | Not inherent | Not inherent | Not inherent | Not inherent | Common efficiency motivation |
| SIMA-DTA | Frozen cached protein tokens plus hybrid adapter | a_t, b_t, c_t and lightweight train-time adapters | Set of ligand, B0, and residual triples | Query-span under strict closures | Explicit [1,B0] residualization | Correct/wrong/permuted/calibration | Required by task definition | c_t=0 | Must beat matched memory/throughput control |

## Allowed Paper-Level Statement

SIMA-DTA is an identifiability-constrained meta-learning framework in which a
scaffold-diverse support set must identify a low-dimensional target-specific
ligand-reordering operator, implemented with an efficient support-conditioned
hybrid protein adapter.

## Disallowed Novelty Claims

Do not claim novelty for Mamba-Transformer hybridization, MAML, orthogonal
projection, D/V-optimal design, wrong-target controls, set encoders, FiLM,
low-rank adaptation, or few-shot DTA individually. Novelty requires the
combined identifiability contract to survive the registered controls.

