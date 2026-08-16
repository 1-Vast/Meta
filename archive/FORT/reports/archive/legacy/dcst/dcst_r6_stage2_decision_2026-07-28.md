# DCST-R6 Stage-2 decision

Date: 2026-07-28  
Decision: `STOP_R6_ADVANCE_TO_CONTENT_ADDRESSED_ROLE_TRANSPORT`

## Admissible result

`dcst_r6_two_stage_seed1729.json` reproduced the source result inside the
complete run: privileged SMB certified `2/4` bands and SMB-NoPriv certified
`0/4`; all source gates passed.

On ChEMBL-37 strict dual-cold development:

| arm | mean target Spearman | RMSE |
| --- | ---: | ---: |
| B0 | 0.0982 | 1.4560 |
| Scratch | 0.0971 | 1.464 |
| NaiveFT | 0.1031 | 1.457 |
| FrozenEncoderFT | 0.0976 | 1.455 |
| FullTransferResidual | 0.0901 | 1.460 |
| DCST | 0.0956 | 1.462 |
| DCST-CertShuffle | 0.0934 | 1.453 |
| DCST-NoPriv | 0.0945 | 1.460 |

The primary paired `DCST - B0` effect was `-0.0025`, 95% bootstrap interval
`[-0.0110, 0.0038]`, versus the frozen MDE `0.0586`. Naive fine-tuning's
small positive point effect (`0.0034`) also crossed zero
(`[-0.0030, 0.0101]`). DCST did not beat scratch, naive fine-tuning, or
no-privileged transfer. Target- and ligand-destroyed effects were
`-0.0041` and `-0.0003` versus B0 and therefore did not establish downstream
mechanism removal. Only the no-material-RMSE-loss gate passed among the
downstream gates.

Wall time was `661.574 s`; peak allocated CUDA memory was `1721.4 MiB`.
Confirmation and sealed test were not scored.

## Diagnosis

R6 solved information entry but not information transport. Its source
interaction-energy matrix indexes 32 normalized sequence positions. Those
positions are exact within a protein but are not homologous coordinates across
the diverse PLINDER and ChEMBL target families. The fitted source matrix was
materially position-specific: its 32 row norms had coefficient of variation
`0.371`, with the largest rows at zero-based segments 2, 19, 27, 15, 6, 21,
26, and 29.

The active downstream gates stayed close to their source certificate priors,
so the failure is not explained by an optimizer erasing them. The more direct
limitation is that a source energy attached to "segment 19" has no invariant
biochemical meaning on an unrelated cold target.

R6 establishes a useful intermediate result: privileged structure creates
held-source affinity directions, but absolute relative-position coordinates
do not transport them across target families. The next route must canonicalize
the protein side by content before spectral transfer.

