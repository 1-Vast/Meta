# CITA, SSNM, and CRISP multi-agent feasibility review

Date: 2026-07-31

## Decision

`STOP_ALL_THREE_CURRENT_DATA_ADMISSION_FAILED`

All three proposals are implementable in principle, but none is scientifically
admissible for a new training run on the current FORT information state. This is
not a CUDA or software failure. Each proposal's own cheapest-first STOP rule is
already triggered by the existing TRAIN-only data and mechanism audits.

No candidate model was implemented or trained in this review. Starting a larger
run would repeat a failed protein-specificity question with more capacity.

## Scope and frozen inputs

| proposal | SHA-256 |
| --- | --- |
| CITA-DTA attachment | `D262C8BD0A025EDABA24E5B00F454FAF355F551FB33552EE914C1FBC8AAC1951` |
| SSNM-DTA attachment | `07FC2934DF86F61AE0BC0FB219362C28DFAB5EF2DE53CFC4348E121B60200C24` |
| CRISP-Mamba-DTA attachment | `3C8FB809BCE5CC725AF5275F070400AB949390B6B606F9C6D37E544217ABAA59` |

The review preserved target/homology-cold roles, support-query scaffold and
chemical-component closure, endpoint separation, provenance fields, episode
hashes, and the confirmation/sealed firewall. Existing user worktree changes
were not modified.

## Environment verification

- Conda environment: `D:\anaconda\envs\drug`.
- PyTorch `2.6.0+cu124`; CUDA runtime `12.4`.
- GPU: NVIDIA GeForce RTX 4060 Laptop GPU, 8188 MiB.
- `mamba-ssm 2.2.4` imports successfully.
- `python -m pytest -q tests`: `76 passed` in 10.76 s.

## Route decisions

| route | direct literature support | decisive current evidence | decision |
| --- | --- | --- | --- |
| CITA-DTA | PCM, matched molecular pairs, activity cliffs, mixed double-mutant cycles, PSICHIC, TAPB, and AdaMBind support individual ingredients. No paper establishes cross-target affinity double differences as a sequence-DTA pretraining objective. | ChEMBL has zero verified same-assay crossed rectangles; Papyrus is one aggregated row per parent-target and has zero document-resolved replicated cells; BindingDB lacks a source-native assay identifier. The smallest correct-protein probe is worse than protein-free. | Stop at D0/P0. Retain only as a new-data-triggered protocol. |
| SSNM-DTA | ms-Mamba validates multiple internal SSM sampling rates for time-series forecasting. Nested Learning and delta-rule fast weights support general memory concepts. Protein Mamba papers support sequence feasibility, not multi-scale few-shot DTA. | Existing memory and protein-specificity gates fail. Runtime `Mamba.forward` does not accept Delta; `dt_scale` is initialization-only. The proposed memory formula does not itself enforce the claimed rank bound. | `STOP_FOR_DATA_IDENTIFIABILITY`; do not implement T1/T2/T3. |
| CRISP-Mamba-DTA | MMR-Mamba, MambaPro, and Shuffle Mamba are peer-reviewed cross-modal image methods. The Topological Trouble paper is a recurrence perspective, not a static DTA validation. | The current model already performs ligand-dependent residue pooling, so the diagnosis is not simple concatenation-only late fusion. FORT has no ligand token sequence, and the protein-specificity gate fails before recurrent expansion. | No-go; do not implement T1-T5. |

## High-value literature findings

### CITA

- Papyrus: Bequignon et al., Journal of Cheminformatics 2023,
  [10.1186/s13321-022-00672-x](https://doi.org/10.1186/s13321-022-00672-x).
  It supports curated aggregated bioactivity, not reconstruction of independent
  raw document observations.
- Rigorous PCM evaluation: Avdiunina et al., JCIM 2025,
  [10.1021/acs.jcim.5c00395](https://doi.org/10.1021/acs.jcim.5c00395).
  Protein permutation tests directly motivate protein-free and wrong-protein
  admission gates.
- MMP cliffs: Hu et al., JCIM 2012,
  [10.1021/ci3001138](https://doi.org/10.1021/ci3001138), and van Tilborg et al.,
  JCIM 2022, [10.1021/acs.jcim.2c01073](https://doi.org/10.1021/acs.jcim.2c01073).
- Mixed double-mutant cycles already instantiate a protein-variant by ligand-edit
  double difference: Hudspith et al., OBC 2019,
  [10.1039/C9OB01558B](https://doi.org/10.1039/C9OB01558B). CITA's novelty is
  therefore the large-scale sequence-DTA transfer, not the four-cell algebra.
- PSICHIC: Koh et al., Nature Machine Intelligence 2024,
  [10.1038/s42256-024-00847-1](https://doi.org/10.1038/s42256-024-00847-1).
  It supports a sequence-based physicochemical interaction operator, but the
  public XL training membership and row-level source lineage are incomplete.
- TAPB and AdaMBind are real peer-reviewed sources, but TAPB is binary DTI and
  AdaMBind does not close scaffold, assay, document, or source overlap:
  [TAPB](https://doi.org/10.1038/s41467-025-66915-1),
  [AdaMBind](https://doi.org/10.1038/s41467-026-70554-5).

### SSNM

- ms-Mamba: Karadag et al., Neurocomputing 680 (2026),
  [10.1016/j.neucom.2026.133226](https://doi.org/10.1016/j.neucom.2026.133226).
  The evidence is time-series-only; no residue-index or DTA transfer is shown.
- Nested Learning: Behrouz et al., NeurIPS 2025,
  [arXiv:2512.24695](https://arxiv.org/abs/2512.24695). Its Hope system does not
  establish SSNM's permutation-invariant one-write support memory.
- Delta-rule fast weights are a closer mathematical precedent: Schlag et al.,
  ICML 2021, [PMLR 139](https://proceedings.mlr.press/v139/schlag21a.html).
- Caduceus, PTM-Mamba, and ProtMamba support bidirectional or long-context
  biological sequence processing, but none tests multi-Delta few-shot DTA:
  [Caduceus](https://proceedings.mlr.press/v235/schiff24a.html),
  [PTM-Mamba](https://doi.org/10.1038/s41592-025-02656-9),
  [ProtMamba](https://doi.org/10.1093/bioinformatics/btaf348).
- Mastropietro et al., Nature Machine Intelligence 2023,
  [10.1038/s42256-023-00756-9](https://doi.org/10.1038/s42256-023-00756-9),
  provides strong evidence that affinity models can rely on ligand memorization.

### CRISP

- MMR-Mamba: Zou et al., Medical Image Analysis 2025,
  [10.1016/j.media.2025.103549](https://doi.org/10.1016/j.media.2025.103549).
  Its aligned MRI modalities and unshared stacked TCM blocks do not validate
  shared-weight DTA recurrence.
- MambaPro: Wang et al., AAAI 2025,
  [10.1609/aaai.v39i8.32879](https://doi.org/10.1609/aaai.v39i8.32879).
  Its prompt is globally learned and is not written from k labeled supports.
- The Topological Trouble With Transformers, [arXiv:2604.17121](https://arxiv.org/abs/2604.17121),
  argues for state propagation across sequential input chunks. Repeating the
  same static protein-ligand pair R times is a separate hypothesis.
- Shuffle Mamba: Cao et al., IEEE TCSVT 2026,
  [10.1109/TCSVT.2026.3668923](https://doi.org/10.1109/TCSVT.2026.3668923).
  A fusion-order stability diagnostic is reasonable, but is not evidence for a
  better protein-ligand interaction mechanism.
- Merrill et al., ICML 2024,
  [PMLR 235](https://proceedings.mlr.press/v235/merrill24a.html), warns that an
  SSM state should not be assumed to imply successful state tracking.

## Decisive local measurements

1. ChEMBL-37 registry-closed audit: 274,400 TRAIN rows and 16,126,069 algebraic
   rectangles collapse to 35,301 target-pair/document units. Verified
   `same_assay_rectangles = 0`. pKi and pKd DD/noise(q90) are about 0.80 and
   0.30; the largest document-source family contributes about 69.2% of units.
2. Papyrus 05.7++: 707,461 aggregated parent-target rows; 147,434 strict rows;
   `document_replicated_parent_target_cells = 0`.
3. IDG-RBP admission: 58 strict pKi episodes. Protein-free
   RMSE/Spearman/pairwise are `1.3550/0.0865/0.5334`; correct protein is
   `1.4266/0.0520/0.5182`. Correct-minus-protein-free RMSE gain is
   `-0.0709 [-0.1221, -0.0200]`; `promote=false`.
4. The earlier trainable AnchorDelta arm has correct and wrong-protein outputs
   that are effectively identical. The memory gate also stops: correct does not
   beat wrong protein or the support-only ranking controls.

The low-capacity IDG-RBP run already used CUDA and cost 126.8 s with 1861 MiB
peak Torch memory. Repeating the question with three Mamba scales, recurrent
Cross-Mamba, or support memory is not the next cheapest falsification.

## Reopening rule

Reopen architecture work only after either:

1. a source/assay/protocol-comparable crossed affinity panel with independent
   provenance blocks is registered; or
2. a new authorized interaction label passes the same target, homology,
   scaffold, document, assay, and source closure.

Then rerun a frozen low-rank correct/protein-free/homology-matched-wrong probe.
Only simultaneous positive component-bootstrap lower bounds for one error
metric and one ranking metric, above a preregistered MDE, authorize one
architecture kill test. CITA should be evaluated first because it changes the
supervision signal; SSNM memory and CRISP recurrence remain secondary capacity
hypotheses.
