# Remote Operations

The active endpoint is declared in `hosts.json`; `ssh_config` supplies the
stable alias `metasieve-remote`. The project stores only the public key and its
fingerprint. The private key remains at
`C:/Users/59964/.ssh/metasieve_remote_ed25519_v3`, outside this project. Passwords
are never persisted.

Connect and verify key-only authentication:

```powershell
ssh -F config/remote/ssh_config -o BatchMode=yes metasieve-remote
```

The canonical remote workspace is `/root/autodl-tmp/MetaSieve`. External raw
datasets and model weights must be downloaded there first:

```bash
/root/miniconda3/bin/python3.12 -m pip install -r config/remote/requirements.txt
/root/miniconda3/bin/python3.12 config/remote/bootstrap_remote.py \
  --hf-endpoint https://hf-mirror.com
```

`bootstrap_remote.py` downloads DAVIS and KIBA from Harvard Dataverse and
checks their SHA-256 values. It downloads all declared ESM-2 snapshots from
Hugging Face into `weights/hub`. Datasets, generated caches, and weights are not
uploaded from the local machine. A failed direct download remains an explicit
remote-setup failure until the authoritative source succeeds.

On the currently declared host, `huggingface.co:443` is unreachable while the
Hugging Face-compatible `https://hf-mirror.com` endpoint is reachable. The
explicit `--hf-endpoint` argument records that provenance instead of relying on
an undeclared environment override.

`sync_manifest.txt` is the allow-list for project synchronization. Raw data,
weights, generated caches, Python bytecode, and smoke outputs are intentionally
not part of general code synchronization.

## Hardware or Key Rotation

1. Add the new endpoint to `hosts.json` and change `active` only after it has
   been verified.
2. Add a separate `Host` block to `ssh_config`; do not reuse a key across
   providers unless explicitly required.
3. Generate the private key under the user SSH directory. Add only its `.pub`
   file and fingerprint here.
4. Install the public key on the new host and test with `BatchMode=yes` before
   removing the old endpoint.
5. Update observed hardware metadata after running `nvidia-smi`; treat it as an
   observation, not a permanent capability guarantee.

Validated research implementations remain under `research/` until their
declared experiments pass. Only then may the corresponding core changes be
promoted into `model/` or `scripts/`.
