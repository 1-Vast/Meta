# L1E-D Protocol Correction

`L1E-D` reached the configured OpenAI-compatible chat-completions endpoint but
received HTTP 400. The preceding `L1E-A` result had already established that
the HTTPS endpoint, credential, and configured model catalog entry were valid.

The official Kimi API chat-completions documentation was checked on
2026-07-27 at `https://platform.moonshot.cn/docs/api/chat`. Its parameter
reference states that `max_tokens` is deprecated and `max_completion_tokens`
must be used. The same reference permits `temperature` from 0 to 1, so the
frozen value 0.0 is retained.

The next instrument may change only the output-token parameter name from
`max_tokens` to `max_completion_tokens`. It remains a one-call, no-data chat
compatibility probe and cannot authorize L1 or any scientific stage.
