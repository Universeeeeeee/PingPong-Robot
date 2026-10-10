# A3 PPO source release

Batch 02 adds the A3 PPO source dependency closure. **Training, checkpoint
creation/resume, environment evaluation and ONNX export have not been run for
this release.** The user explicitly waived those runtime checks for this batch.
This is a source release, not a ready-to-run asset or pretrained-policy bundle.
The earlier host checks were performed on an upstream review snapshot; they do
not establish that the adapted release code passes tests or trains successfully.

The default task is `PingPong` / `PingPong-Hitter-AgibotA3-v0`: the formal
`hitter_pure` actor consumes 110 values and emits 31 joint action values at 50 Hz.
The 29 active channels inside that 31-value interface are not a 107-to-29
legacy network. Historical source classes, checkpoint keys, the
`HOPE-HitterPingPong-AgibotA3-v0` Gym alias and the `HOPEPingPong` task recipe
remain for compatibility. Imports register only the A3 training task and its
motion-tracking baseline. Independent match-simulation tasks are not registered.

## External runtime and assets

Use Linux with an NVIDIA CUDA GPU and mutually compatible Isaac Sim, Isaac Lab
and `rsl_rl` installations. Use that installation's Python 3.10+ interpreter.
The package additionally imports PyTorch, Gymnasium, Hydra/OmegaConf, NumPy,
PyYAML, TensorDict, ONNX and TensorBoard-related runner dependencies. Its
`setup.py` does not install the simulator or GPU stack. W&B-related compatibility
code remains, but the provided training entry uses local logs/checkpoints.
No paid storage or external tracking service is enabled by this release.

These required inputs are **not included**:

- An authorized racket-equipped Agibot A3 URDF and meshes. Obtain these from
  the vendor or your authorized asset provider. This release does not establish
  permission to redistribute the vendor CAD files.
- Forehand and backhand motion NPZ files, in that order, matching the 31-joint
  canonical order in `../config/joint_order_agibot_a3.yaml`. Supply both paths
  explicitly; the recipe's historical default filenames are placeholders here.
- The recipe's default-enabled correlated tuple bank and JSON receipt. Supply
  an authorized matching pair via the task overrides below. The recipe requires
  schema 2 and at least 50,000 rows per side; source and hashes must agree with
  its loader. This release does not include a generator or that data bank.
- A real checkpoint if resuming, playing, evaluating or exporting. The upstream
  `model_21800.pt` found locally was only a 132-byte Git LFS pointer. Neither that
  pointer nor any weights are uploaded. Training from scratch needs no old weight.

The default recipe also loads `../../configs/ball_physics_venue.yaml` by upward
search. Its parameters are retained from upstream; their fit claims were not
reproduced here. Canonical action adapter constants are in
`../../configs/action_adapter.yaml`.

## Commands for an equipped machine (not executed for this release)

From `training/whole_body_tracking`, select your existing Isaac Python and,
if needed, `ISAACLAB_ROOT`, then source the launcher:

```bash
export ISAAC_PYTHON=/path/to/authorized/isaacsim/python.sh
export ISAACLAB_ROOT=/path/to/IsaacLab
source setup_train_env.sh
pingpong_isaac_py -m pip install -e source/whole_body_tracking
```

Prepare an externally supplied URDF package (`urdf/*.urdf`, `meshes/*`):

```bash
python3 scripts/prepare_a3_isaac_asset.py --source-root /path/to/authorized/a3-package
```

Alternatively set `PINGPONG_A3_URDF` to an already prepared, authorized URDF
whose mesh references resolve locally. This variable is read before the robot
configuration is imported. Vendor assets are never downloaded automatically.

```bash
pingpong_isaac_py scripts/train.py task=PingPong algo=ppo headless=true \
  motion_file=/path/to/forehand.npz motion_file_2=/path/to/backhand.npz \
  task.racket.venue_tuple_bank_path=/path/to/tuple-bank.npz \
  task.racket.venue_tuple_bank_receipt_path=/path/to/tuple-bank.receipt.json
```

For a bounded future smoke run, add `num_envs=12 max_iterations=2
algo.runner.save_interval=1`. The formal recipe uses six PPO minibatches, so
select a small environment count compatible with that setting. For ordinary
weight/optimizer/iteration resume, add `checkpoint_path=/path/to/model.pt`.
The entry calls `runner.load`; it does not claim exact environment/RNG replay.
Keep all motion and task overrides identical when resuming. Strict model loading
rejects incompatible tensor shapes; old 107-to-29 checkpoint migration is not
promised or verified.

```bash
pingpong_isaac_py scripts/evaluate.py --checkpoint /path/to/model.pt \
  --motion-file /path/to/forehand.npz --motion-file-2 /path/to/backhand.npz
pingpong_isaac_py scripts/export_onnx.py --checkpoint /path/to/model.pt \
  --motion-file /path/to/forehand.npz --motion-file-2 /path/to/backhand.npz
```

Evaluation/export rebuild registered defaults, rather than replaying every
training Hydra override or loading the saved `env.yaml`. Check their effective
environment recipe against the training recipe before interpreting results.
The ONNX joint-order gate and manifest code are supplied, but no exported graph
or loading result is delivered. Optional rigid-ball truth audits refer to
standalone match-simulation configuration omitted from this batch; keep the
formal `physical_ball: false` setting. Those audits are not supported here.
FFS code/weights are not bundled or established as ready by this release.

See [SOURCE_MANIFEST_BATCH_02.json](SOURCE_MANIFEST_BATCH_02.json) for exact
source/release hashes and exclusions. Numerical PPO and task parameters were
not retuned. Changes adapt repository paths, expose neutral public aliases,
limit automatic registration, and allow an external URDF location. Necessary
upstream names, copyrights and checkpoint metadata identifiers remain.

## Batch 01 foundations (historical scope)

This batch publishes a standalone observation-contract registry and a no-spin
return-success evaluator. It does not include PPO training, Isaac Lab environments,
robot assets, motion data, checkpoints, hardware control, or ONNX export.

The `hitter_pure` observation contract has 110 entries, including 31 previous
action entries. Other historical/experimental layouts remain in the registry for
explicit compatibility checks; they must not be mixed with `hitter_pure`.
Contract validation rejects reordered terms even when total dimensions match.
No policy action decoder is delivered in this batch.

`success_metric.py` evaluates contact, net clearance, and the first bounce in
the opponent half using quadratic drag and gravity. It is a no-spin approximation,
not an RL training algorithm or a hardware accuracy certification. The original
synthetic return tests are retained.

## Run the scoped tests

Use Python 3.10+ and NumPy, PyYAML, and pytest. From the repository root:

```bash
python -m pip install -r training/whole_body_tracking/requirements-batch-01.txt
python -m pytest -q training/whole_body_tracking/tests
```

The legacy direct test entry remains available:

```bash
python training/whole_body_tracking/tests/test_success_metric.py
```

The YAML configuration is found by walking upward to `configs/ball_physics.yaml`.
Set `BALL_PHYSICS_CONFIG` to override it. The legacy `HOPE_BALL_PHYSICS_CONFIG`
variable remains supported, with the neutral name taking precedence. The upstream
fallback behavior is unchanged: a missing explicit file or unavailable PyYAML
returns defaults. Install PyYAML and run the configuration tests to verify that
your intended file is actually loaded.

## Source, changes, and compatibility

Imported from [the upstream repository at the pinned revision](https://github.com/hitchopen/HOPE/tree/d14886823fc72304c611d99ad6cdc105c51d3aca).
See [SOURCE_MANIFEST.json](SOURCE_MANIFEST.json) for source/target mappings and
SHA-256 hashes, and [third-party notices](../../THIRD_PARTY_NOTICES.md) for licenses.

The registry, original metric tests, and license texts are copied without changes.
The metric adds a neutral environment-variable name and keeps the legacy alias.
The physics YAML changes only explanatory/provenance text and references to
upstream-only material; physical parameters and units remain unchanged. The
capture/fitting assertions are upstream provenance, not newly reproduced results.
New focused tests cover observation-order rejection and actual YAML loading.

Academic contract names, source links, copyright holders, and license wording
remain intact for attribution and compatibility. This release makes no changes
to existing repository modules, public message schemas, or checkpoint formats.
