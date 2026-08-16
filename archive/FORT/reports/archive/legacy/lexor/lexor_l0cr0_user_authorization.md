# LEXOR L0C-R0 Direct User Authorization Receipt

The user instruction to continue and confirm whether the operation is effective
authorizes one local `L0C-R0` Reinecke calibration run. The authorization is
limited to the already acquired local S1/S2/S3 bytes bound by the frozen
calibration manifest. It does not authorize downloading a new source, calling
an API, reading `.env`, L1 extraction, label access, model training, retries,
or a scientific-stage transition.

<!-- LEXOR_L0CR0_SOURCE_MANIFEST_SHA256: c406208350e83c05b6ec372e627e7546f6385b5ae39f925f03f6a14e1acc972d -->
<!-- LEXOR_L0CR0_DIRECT_USER_AUTHORIZATION: true -->
<!-- LEXOR_L0CR0_REQUIRES_CLI_EXECUTE: true -->
