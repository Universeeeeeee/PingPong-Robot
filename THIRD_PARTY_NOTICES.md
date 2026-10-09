# Third-party notices for training foundations, batch 01

These notices apply to the files added by this batch. They do not relicense the
repository's pre-existing code or imply that the full upstream project is included.

## Upstream project and contributed evaluation/configuration code

- Source: https://github.com/hitchopen/HOPE
- Pinned revision: `d14886823fc72304c611d99ad6cdc105c51d3aca`
- Upstream project license: Apache-2.0, retained verbatim in
  `LICENSES/HOPE-Apache-2.0.txt`.
- The copied evaluator retains its original header:
  `Copyright (c) 2025, Intelligent Racing Inc. (dba Hitch Interactive).`
  Its SPDX identifier is `Apache-2.0`.
- The evaluator adds `BALL_PHYSICS_CONFIG` while retaining the original
  environment variable as a compatibility alias. Physics YAML edits change
  explanatory/provenance text only. The original metric test is unchanged.

## whole_body_tracking framework

- Upstream's notices describe `hope_training/whole_body_tracking/` as a
  fork/derivative of BeyondMimic (`HybridRobotics/whole_body_tracking`), an Isaac Lab
  motion-tracking reinforcement-learning extension.
- Upstream framework license: MIT.
- Copyright: `Copyright (c) 2024, The Isaac Lab Project Developers.`
- The original MIT text is retained verbatim in
  `training/whole_body_tracking/LICENCE`.
- Only the observation registry, standalone evaluator, and selected tests are
  included here; no complete training framework or third-party robot assets are
  redistributed in this batch. Specific Apache-2.0 evaluator headers are preserved.

The source/target mapping and hashes are recorded in
`training/whole_body_tracking/SOURCE_MANIFEST.json`. Necessary source names,
copyright attribution, research references, and compatibility identifiers are
preserved. NumPy, PyYAML, and pytest are runtime/test dependencies installed from
their own distributions; their source and binaries are not vendored here.
