# Batch 03 validation — 2026-10-11

Base main: `06f85f36a694315faf8483ff62c0b930f948660b`.

Environment: Python 3.11.15, PyTorch 2.13.0, PyYAML 6.0.3, pytest 9.1.1, CPU.
All checks used already installed dependencies and ran offline.

- `python -m pytest -p no:cacheprovider` in this module: **11 passed**.
  Existing tests cover flight drag/gravity, spin direction, table bounce,
  direct/swept racket contact, inverse contact recovery, bounded contact fitting,
  deterministic nominal sampling, randomization ranges and spin magnitude.
- `python -m unittest discover -s mocap/tests -v` at repository root:
  **22 passed** (existing adapter regression suite).
- `PIP_NO_INDEX=1 python -m pip wheel --no-deps --no-build-isolation .`:
  wheel successfully built with the YAML included. Build outputs are excluded
  from the commit.
- Extracted-wheel check from an independent directory, without source package
  import: imports and default YAML lookup passed. All four profiles ran with
  seed 37, 128 environments and 200 steps of 0.005 s, repeated twice.
  Position/velocity/spin were finite and exactly reproducible, inactive-state
  position/velocity/spin were preserved, and each profile produced 127 bounces.
- Source-to-release manifest hashes verified; Python syntax and import closure
  checked. New runtime imports are standard library, PyTorch, PyYAML and local
  package modules.
- Scope scan checked own-code brand names, private paths, addresses, credentials,
  key material, generated files, datasets, videos, weights and file sizes.
  Remaining names are necessary original API identifiers in the transformation
  manifest and the preserved research attribution. No excluded artifacts added.
- `git diff --check`: passed.

At the pre-publication check, main was unprotected and the repository had zero
GitHub Actions runs and zero commit check runs. There are no workflow files in
this repository. Local test success does not imply remote CI ran.

No complete training, checkpoint, restore, environment evaluation, ONNX export,
GPU backend or hardware test was performed for this standalone physics batch.
The simulation coefficients are not empirically validated by these checks.
