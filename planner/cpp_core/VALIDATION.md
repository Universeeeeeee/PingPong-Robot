# Standalone C++ core validation — 2026-10-29

Environment: macOS ARM64, AppleClang 21.0.0.21000101, CMake/CTest 4.4.5,
Eigen 3.4.1, GoogleTest 1.11.0. Python regressions used Python 3.13.13,
NumPy 2.5.4, PyYAML 6.0.3 and pytest 9.1.1.

From the repository root, set `Eigen3_DIR` and `GTest_DIR` to the installed
packages' CMake configuration directories, and `GTEST_LIBRARY_DIR` to the
GoogleTest shared-library directory when that installation needs an RPATH:

```sh
cmake -S planner/cpp_core -B build/cpp-core \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
  -DEigen3_DIR="$Eigen3_DIR" -DGTest_DIR="$GTest_DIR" \
  -DCMAKE_BUILD_RPATH="$GTEST_LIBRARY_DIR"
cmake --build build/cpp-core --parallel 2
ctest --test-dir build/cpp-core --output-on-failure
python -m unittest discover -s mocap/tests -v
python -m pytest -q -p no:cacheprovider \
  training/whole_body_tracking/tests/test_actor_observation_contract.py \
  training/whole_body_tracking/tests/test_ball_physics_config.py \
  training/whole_body_tracking/tests/test_success_metric.py
```

Results:

- Release library and test executable built successfully. Six C++ tests passed:
  scalar-reference prediction/solve, invalid no-crossing, legacy spin dispatch,
  Coulomb topspin, venue-grip scalar-reference crossing and effective bounce gain.
- Existing adapter regressions: 22 passed. Foundation regressions: 13 passed.
- Six imported source hashes and six transformed release hashes matched the
  source manifest. The Apache-2.0 license remains unchanged. Attribution and
  modification notices are preserved.
- Scoped naming, credential/private-path and artifact scans found no excluded
  content. Local quoted include paths resolve. `git diff --check` passed.
- No GitHub Actions runs, check runs or status contexts existed at the
  pre-merge review. Local results do not represent remote CI success.

Only Eigen was newly installed; CMake and GoogleTest were reused. Build outputs,
dependency source/binaries and local paths are excluded from this release.
The initial GoogleTest discovery could not find its shared library; setting the
local build RPATH resolved that environment issue before all tests ran.

These results validate the standalone mathematical core and its regression
cases. They do not establish empirical coefficient accuracy, ROS integration,
GPU/other-platform behavior, hardware readiness, or complete RL training,
checkpoint/resume, environment evaluation or ONNX export.
