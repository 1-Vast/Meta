# LEXOR L1E Environment Overall Report

Date: 2026-07-27

## Scope

This report concerns isolated provider engineering instruments only. No LEXOR
scientific stage ran; no paper, supplement, raw measurement, structure, FORT
identifier, development label, confirmation label, or sealed label was sent to
the provider. `L1_LOCKED` remains the governing scientific state.

## Frozen Results

| Instrument | Request scope | Result | What it establishes |
| --- | --- | --- | --- |
| L1E | One synthetic abstention chat request | `fail`, no response, one call consumed | The original adapter path could not complete a request; this result is sealed and was not retried. |
| L1E-A | One authenticated `GET /v1/models` request | `pass`, HTTP 200, configured model listed | The configured HTTPS endpoint, credential, and configured model catalog entry are usable. |
| L1E-D | One synthetic chat request using `max_tokens` | `fail`, HTTP 400 | Chat-completions reached the server, but the request was rejected. |
| L1E-E | One synthetic chat request using documented `max_completion_tokens` | `fail`, HTTP 400 | Replacing the deprecated output-token parameter alone does not resolve the rejection. |
| L1E-F | One synthetic Kimi K2.6 profile with omitted `temperature`, disabled thinking, and `max_completion_tokens=16` | `pass`, HTTP 200, 15 tokens | The configured API accepts a bounded no-data chat request under a documented K2.6-compatible profile. |

The two chat failures have the same secret-free response digest and byte count:
`6ec6fe8351887771fd202acc9d4ab46f06cd54ff49ba4ac970558ac1ea7bfccf`,
108 bytes. That supports a stable server-side rejection condition rather than
a stochastic output-contract failure.

## Evidence-Based Interpretation

The Kimi K2.6 model parameter reference narrows the cause: it permits explicit
`temperature` only at `1.0` with thinking enabled or `0.6` with thinking
disabled, and recommends omitting it otherwise. Both failed chat probes sent
`temperature: 0.0`; their identical response digest is consistent with that
documented rejection. `max_tokens` is separately deprecated in favor of
`max_completion_tokens`.

L1E-F froze a documented compatible profile before dispatch: it omitted
`temperature`, set `thinking: {"type":"disabled"}`, retained an empty
synthetic firewall, and sent only `Return exactly: OK`. It passed with HTTP 200
and the exact `OK` response contract; its result contains only digests, status,
and token accounting.

This proves that the configured credential, endpoint, model, and an
appropriately constrained K2.6 chat request are usable. It does not isolate
the causal contribution of temperature alone because L1E-F also selected the
documented non-thinking mode and a smaller fixture. That isolation is not
needed for the engineering objective, so no further live call is authorized or
required.

## Decision

```text
LEXOR_L1E_ENGINEERING_CHAT_COMPATIBLE
L1_LOCKED
```

This is an engineering compatibility pass, not an L1 extraction-reliability
result and not evidence for the LEXOR scientific hypothesis. It authorizes
neither scientific L1 nor any L1-L6/F1-F4 transition.

## Next Permitted Exploration

No additional engineering call is needed. A future scientific L1 run remains
blocked by L0/L0C and must separately preregister its model profile. For Kimi
K2.6, it must not send `temperature: 0.0`; it must either omit `temperature`
or use the provider's documented temperature/thinking combination.
