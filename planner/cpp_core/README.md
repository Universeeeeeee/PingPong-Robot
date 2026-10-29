# Standalone C++ planner core — batch 04

This source batch contains trajectory prediction and racket target solving in
the `pingpong_planner` namespace. It predicts a ball's crossing of a strike
plane, including drag, optional spin and table bounce, and solves an outgoing
velocity and racket target for a requested landing point and flight time.
It includes six upstream regression tests with their numerical assertions.

The library takes `BallPhysics`, `TableParams` and `PlannerConfig` values in
code. It does not load YAML automatically. Upstream numerical defaults and
formulas are unchanged; their empirical accuracy has not been established here.
Coordinates and time must match the supplied model: the default table spans
`x=[0, length]`, `y=[-width, 0]`, with its surface at `z=0`; positions use metres
and time uses seconds. Spin inputs use radians per second.

## Validation status

Validated on macOS ARM64 on 2026-10-29 with AppleClang 21, CMake 4.4.5,
Eigen 3.4.1 and GoogleTest 1.11.0: the Release library build and all six C++
regression tests passed. Repository regressions also passed: 22 adapter tests
and 13 foundation tests. These are local results; no remote CI runs were
available. See [VALIDATION.md](VALIDATION.md) for commands and limits.

## Build

Use CMake 3.16+, a C++17 compiler and a separately installed Eigen 3.3+.
GoogleTest is also required when tests are enabled. From this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build --target pingpong_planner_core
```

To configure the included tests on an equipped machine:

```sh
cmake -S . -B build-tests -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build-tests
ctest --test-dir build-tests --output-on-failure
```

Use your dependency installation's `Eigen3_DIR` and `GTest_DIR` if CMake cannot
locate them. No dependency is downloaded by this project's CMake configuration.

## Scope

This is a deterministic mathematical planning component. No ROS node, estimator,
robot controller, hardware configuration, replay data, trained policy, model
weight or complete reinforcement-learning runtime is included. Training,
checkpoint/resume, evaluation/export and hardware behavior are not validated by
this batch. The earlier PPO source release retains its own unverified status.

## Source and compatibility

Derived from the Apache-2.0 planner in
[the pinned upstream revision](https://github.com/hitchopen/HOPE/tree/d14886823fc72304c611d99ad6cdc105c51d3aca/hope_ws/src/hope_planner_cpp).
The [source manifest](SOURCE_MANIFEST.json) records each source/release hash and
the modifications. The full [license](../../LICENSES/HOPE-Apache-2.0.txt) and
[third-party notices](../../THIRD_PARTY_NOTICES.md) apply.

Public includes and the namespace use `pingpong_planner`. Person-derived bounce
identifiers use `Coulomb` / `coulomb`; the associated scientific attribution is
retained in comments. Shared declarations are limited to the seven types needed
by these solvers. This new standalone API does not promise binary compatibility
with upstream ROS packages. Implementation, declarations and test references
were updated together; formulas, parameter values and numerical assertions were
preserved. Modified source files carry explicit modification notices.

Eigen and GoogleTest remain external build/test dependencies under their own
licenses. Neither their source nor their binaries are uploaded in this batch.
