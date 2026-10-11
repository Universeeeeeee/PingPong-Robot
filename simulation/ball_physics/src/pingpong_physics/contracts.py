from __future__ import annotations

from dataclasses import dataclass

import torch


@dataclass
class BallState:
    """Privileged world state. It must never be passed directly to the planner."""

    position: torch.Tensor
    velocity: torch.Tensor
    angular_velocity: torch.Tensor
    active: torch.Tensor

    def clone(self) -> "BallState":
        return BallState(
            self.position.clone(),
            self.velocity.clone(),
            self.angular_velocity.clone(),
            self.active.clone(),
        )


@dataclass
class BallEvents:
    table_bounce: torch.Tensor
    table_contact_position: torch.Tensor
