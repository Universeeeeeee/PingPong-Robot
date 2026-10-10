# Third-party notices for training source releases

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

## Batch 02: A3 PPO source dependency closure

- Source and pinned revision remain https://github.com/hitchopen/HOPE at
  `d14886823fc72304c611d99ad6cdc105c51d3aca`; actual source/release SHA-256 hashes
  are in `training/whole_body_tracking/SOURCE_MANIFEST_BATCH_02.json`.
- The copied `whole_body_tracking` framework is a derivative of the
  MIT-licensed BeyondMimic / Isaac Lab tracking project. Its license remains
  verbatim in `training/whole_body_tracking/LICENCE`, including
  `Copyright (c) 2024, The Isaac Lab Project Developers.`
- Upstream HOPE-specific contributions are described upstream under Apache-2.0;
  the full text remains in `LICENSES/HOPE-Apache-2.0.txt`. Existing source
  copyright/SPDX headers are retained. This notice does not relicense vendor
  materials or third-party dependencies.
- Robot configuration transcribes upstream Agibot A3 interface names, gains and
  limits. The vendor URDF, meshes, CAD, deployment SDKs and weights are excluded;
  operators must obtain authorized copies externally. The source asset
  preparation script's historical claim of a bundled vendor asset does not
  apply to this release.
- Shared table geometry constants and ball-model helper source are included as
  transitive dependencies of PPO return shaping. No table USD/mesh asset or
  standalone match-simulation environment is redistributed.
- `configs/ball_physics_venue.yaml` retains upstream numerical parameters while
  stripping comments containing personal research paths and unshipped dataset
  references. Upstream experimental fit claims have not been reproduced here.
  Motion data, tuple-bank data/receipts, raw recordings and private paths are
  excluded from this source release.
- PyTorch, Isaac Sim/Lab, rsl_rl, ONNX and other runtime dependencies are installed
  separately under their own terms. No third-party runtime binaries are vendored.

Runtime training, checkpoint creation/resume, environment evaluation and ONNX
export are unverified for batch 02. The user authorized source publication
without those checks; this is not a claim of training or deployment readiness.
