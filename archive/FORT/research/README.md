# Research Packages

Research code is grouped by the experimental substrate rather than by run date.
Each runner remains executable with `python -m <module>` from the repository
root. Reports retain historical source paths; use the package map below for the
current source location.

| Package | Contents |
|---|---|
| `research.a2s` | A2S-DTA data, baseline, gate, CMAL, and transfer-object runners. |
| `research.psep` | PSEP substrate, measurement, operator, transfer, and replication runners. |
| `research.adambind` | AdaMBind data audit and smoke runners. |
| `research.shared` | Shared data processing and earlier meta-learning diagnostics. |
| `research.fable` | Existing separately packaged FABLE work. |

New isolated experiments enter the matching package only after a written
hypothesis and test plan. Failed work is preserved in `reports`, recorded in
`history.md`, and is not promoted to `model` or `scripts`.
