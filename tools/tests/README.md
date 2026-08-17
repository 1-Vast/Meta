# Tests

The suite covers retained production and comparator contracts: data governance,
split seals, QPSMP/BPSF, similarity, level-shape, relative transport, Cartesian
algebra, R14 diagnostics, checkpoint loading and stage execution.

Tests dedicated solely to deleted research modules were removed with those
modules. Their historical pass counts are evidence records, not active regression
contracts. Use `pytest tools/tests -q`; set `RUN_SLOW=1` only for full-corpus
subprocess smokes and `RUN_RESEARCH_GATES=1` for retained opt-in training gates.
