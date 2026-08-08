# Manual Validation

This directory contains manually invoked validation and audit commands. The
automated pytest suite remains in `tests/`.

Run the end-to-end source-view invariant check with
`python -m test.smoke --sealed-dir <sealed-dir> --protein-cache <cache>`.
Run the runtime-bound audit with `python -m test.audit <dataset> --sealed-dir
<sealed-dir> --protein-cache <cache>`. `audit_sealed_dataset.py` verifies the
sealed source/metaval artifact boundary before either label view is mounted.

`compiled_dataset_tools.py` contains offline fixtures used by regression tests
to verify that query labels remain outside model inputs. It is not part of the
production training runtime.
