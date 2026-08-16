# UBSE-P0A Attempt-2 Detached-Execution Amendment

**Frozen:** 2026-07-29 20:04 Asia/Shanghai  
**Scope:** execution transport only  
**Scientific protocol change:** none

## Attempt-1 disposition

Attempt 1 started at 2026-07-29 18:56:29 with Windows PID 7976.
At 19:53:49, after 57.3 minutes, it was still active with:

- 100% reported GPU utilization;
- 7,817/8,188 MiB reported framebuffer use;
- 4,667.9 process CPU seconds;
- no stderr or scientific result.

By approximately 20:02, the process no longer existed and GPU utilization
was 0%. No seed checkpoint, validation ledger, result JSON, temporary output,
Python traceback, Windows Application Error, or Windows Error Reporting event
was present.

The foreground tool execution session also became unavailable across the
continuation boundary. This is consistent with, but does not prove, an
execution-session lifetime termination. It is not evidence about P0A model
quality.

Binding disposition:

`ABORT_UBSE_P0A_ATTEMPT1_EXECUTION_SESSION_LOSS_NO_SCIENTIFIC_DECISION`

No partial state exists and nothing from attempt 1 may enter the accepted
result.

## Frozen restart

Attempt 2 must use the identical:

- `research/ubse_p0a.py`;
- source and panel manifests;
- seeds `(1729, 1730, 1731)`;
- four epochs and evaluation epochs `(1, 2, 4)`;
- optimizer, batching, controls, gates, and thresholds;
- CUDA-only execution and outcome firewall.

The P0A code SHA-256 remains:

`c1732e6d4cfacd8e9ee42ab03bde20d14d98a5ea9d9dff0c6f247cf998211ad7`

The already frozen shuffle-correction protocol SHA-256 remains:

`bdc9636a711189cd52aab090b509e55175abda417951b0824fe108c3e5711d94`

The only execution change is to launch the same module as a hidden detached
Windows process and redirect unbuffered stdout/stderr to:

- `tmp/ubse_p0a_attempt2.stdout.log`;
- `tmp/ubse_p0a_attempt2.stderr.log`.

The detached process is monitored by PID, GPU telemetry, log growth, and
atomic scientific artifacts. It must not be interpreted from logs alone.
After a successful raw run, the fixed-model-seed shuffle correction remains
mandatory and only the corrected ledger may drive the P0A decision.

If attempt 2 exits without a complete raw result, the logs must be preserved
and diagnosed before any further run. Gates and thresholds may not be changed
as a response to an execution failure.
