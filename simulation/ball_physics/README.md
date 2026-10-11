# Standalone ball physics — batch 03

This package publishes the existing local tensor-only table-tennis physics
implementation and its dependency closure. It adds drag, gravity, Magnus force,
spin-coupled table bounce, batched racket contact with swept contact detection,
and inverse / velocity-box-constrained contact fitting. `DomainRandomizer`
provides four existing YAML profiles; `BallState` and `BallEvents` are the two
unchanged state contracts needed by this module.

## Run

Python 3.10–3.11, PyTorch 2.x and PyYAML 6.x are required. From this directory:

```bash
python -m pip install -e '.[dev]'
python -m pytest
```

Or use existing installed dependencies without installation:

```bash
PYTHONPATH=src python -m pytest
```

```python
import torch
from pingpong_physics import BallState, BallPhysicsWorld, DomainRandomizer

parameters = DomainRandomizer().sample("nominal", 1)
state = BallState(
    position=torch.tensor([[0.0, 0.0, 1.0]]),
    velocity=torch.tensor([[5.0, 0.0, 0.0]]),
    angular_velocity=torch.zeros(1, 3),
    active=torch.ones(1, dtype=torch.bool),
)
next_state, events = BallPhysicsWorld().step(state, parameters, 0.01)
```

## Source and compatibility

`SOURCE_MANIFEST.json` records the original working-tree file hashes, release
hashes, source repository revision and exact transformations. This is the local
project's own simulation implementation, published under the user's explicit
authorization; this batch does not assign a new license to existing project code.
No vendor implementation, runtime binary, robot asset, dataset or weight is
included. The original source attributes the world equations to Duerr et al.
(Nature, 2026); that research attribution remains in the code and manifest.

The own-code public names are `BallPhysicsConfig` and `BallPhysicsWorld`, under
`pingpong_physics`. All shipped callers/tests use these neutral names. This is a
new package; previous repository APIs are unchanged. The source development
package is untouched. The randomization YAML is included as package data, so
default configuration lookup works after wheel installation.

## Validation and limits

See `VALIDATION.md` for this batch's offline checks. This is a kinematic ball
world, using the source's center-height convention and fixed-step integration.
Racket restitution/tangential coefficients are provisional. The research-fitted
restitution and learned residual model are unavailable and are not reproduced.
Tests demonstrate software behavior, not measured real-ball accuracy.

This module contains no training loop or robot integration. Batch 02's training,
checkpoint creation/resume, environment evaluation and ONNX export remain
unverified. There is no training-success or hardware-readiness claim here.
