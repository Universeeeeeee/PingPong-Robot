# Third-party notices for source releases

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

## Batch 03: standalone ball physics

- The local project-owned simulation source is released under the user's
  explicit publication authorization. This batch does not relicense existing
  code or introduce a new blanket license. Source and transformed release
  hashes are recorded in `simulation/ball_physics/SOURCE_MANIFEST.json`.
- The original source cites Duerr et al. (Nature, 2026), equations 1–8, for
  aerodynamic and contact structure. Necessary research attribution is retained;
  no paper text, external implementation, vendor asset, fitted residual network
  or weight is redistributed. The provisional contact coefficients are documented.
- PyTorch, PyYAML, pytest and setuptools remain separately installed dependencies;
  their source and binaries are not vendored. Existing third-party license and
  copyright notices elsewhere in the repository are unchanged.

## Batch 04: standalone C++ trajectory and racket target planning

- Source: `hope_ws/src/hope_planner_cpp/` in
  https://github.com/hitchopen/HOPE at
  `d14886823fc72304c611d99ad6cdc105c51d3aca`. The upstream package declares
  Apache-2.0, and its project notices identify planner contributions as
  Apache-2.0. The original full license remains verbatim in
  `LICENSES/HOPE-Apache-2.0.txt`.
- Only two solver implementations, their headers, the required shared types
  and the associated regression test file are derived from that package.
  Source/release hashes and modifications are recorded in
  `planner/cpp_core/SOURCE_MANIFEST.json`. Modified files carry a notice.
- Namespace/include paths and person-derived model API names are neutralized
  consistently across code and tests. Necessary Nakashima/Ace scientific
  attribution remains in comments. No formulas, numerical defaults or numerical
  test assertions are retuned. No external paper text is redistributed.
- Eigen and GoogleTest are separately installed build/test dependencies with
  their own licenses. No dependency source, header tree or binary is vendored.
  ROS nodes/messages, hardware controls, raw data, checkpoints and weights are
  excluded. This notice does not relicense other existing repository code.
