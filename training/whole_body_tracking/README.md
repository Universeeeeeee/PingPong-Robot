# Training foundations: batch 01

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
