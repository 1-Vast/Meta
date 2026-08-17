# Tools workspace

`tools/` groups development infrastructure outside the admitted model runtime:

- `research/`: unadmitted experiments. Passing work is promoted into `model/`
  and `scripts/`; failed implementations are removed only after their compact
  verdict and required reproducibility artifacts are recorded.
- `tests/`: maintained pytest contracts and test-only fixtures.
- `runtime/`: ignored local executables/downloads such as CD-HIT and MMseqs2.

`main.py` orchestrates only admitted code from `model/` through entry points in
`scripts/`. It never dispatches code directly from `tools/research/`.
