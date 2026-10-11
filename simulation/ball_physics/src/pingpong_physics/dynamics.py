from __future__ import annotations

from dataclasses import dataclass
import math

import torch

from .contracts import BallEvents, BallState
from .randomization import DomainParameters


@dataclass(frozen=True)
class BallPhysicsConfig:
    mass_kg: float = 2.7e-3
    radius_m: float = 0.02
    air_density_kg_m3: float = 1.204
    drag_coefficient: float = 0.55
    gravity_m_s2: float = 9.81
    table_height_m: float = 0.02
    table_half_length_m: float = 1.37
    table_half_width_m: float = 0.7625
    restitution_min: float = 0.75
    restitution_max: float = 0.99
    post_bounce_linear_sigma_m_s: tuple[float, float, float] = (0.03, 0.03, 0.06)
    post_bounce_angular_sigma_rad_s: tuple[float, float, float] = (3.0, 3.0, 6.0)
    racket_restitution: float = 0.8
    racket_tangential_transfer: float = 0.2
    racket_contact_distance_m: float = 0.08
    racket_ball_separation_m: float = 0.04


def racket_contact_frame(normal: torch.Tensor) -> torch.Tensor:
    """Return the stable local frame used by the simulated racket contact."""
    normal = normal / torch.linalg.vector_norm(
        normal, dim=-1, keepdim=True
    ).clamp_min(1e-8)
    reference = torch.tensor(
        (0.0, 0.0, 1.0), device=normal.device, dtype=normal.dtype
    ).expand_as(normal)
    parallel = torch.abs(torch.sum(reference * normal, dim=-1)) > 0.95
    fallback = torch.tensor(
        (0.0, 1.0, 0.0), device=normal.device, dtype=normal.dtype
    ).expand_as(normal)
    reference = torch.where(parallel[:, None], fallback, reference)
    tangent_x = torch.cross(reference, normal, dim=-1)
    tangent_x = tangent_x / torch.linalg.vector_norm(
        tangent_x, dim=-1, keepdim=True
    ).clamp_min(1e-8)
    tangent_y = torch.cross(normal, tangent_x, dim=-1)
    return torch.stack((tangent_x, tangent_y, normal), dim=-1)


def racket_contact_outgoing_velocity(
    incoming_velocity: torch.Tensor,
    angular_velocity: torch.Tensor,
    racket_velocity: torch.Tensor,
    racket_normal: torch.Tensor,
    *,
    config: BallPhysicsConfig | None = None,
) -> torch.Tensor:
    """Evaluate the simulated linear racket-contact velocity map."""
    cfg = config or BallPhysicsConfig()
    rotation = racket_contact_frame(racket_normal)
    relative = incoming_velocity - racket_velocity
    relative_local = torch.bmm(
        rotation.transpose(1, 2), relative.unsqueeze(-1)
    ).squeeze(-1)
    omega_local = torch.bmm(
        rotation.transpose(1, 2), angular_velocity.unsqueeze(-1)
    ).squeeze(-1)
    k = cfg.racket_tangential_transfer
    outgoing_relative_local = relative_local.clone()
    outgoing_relative_local[:, 0] = (
        (1.0 - k) * relative_local[:, 0]
        + k * cfg.radius_m * omega_local[:, 1]
    )
    outgoing_relative_local[:, 1] = (
        (1.0 - k) * relative_local[:, 1]
        - k * cfg.radius_m * omega_local[:, 0]
    )
    outgoing_relative_local[:, 2] = (
        -cfg.racket_restitution * relative_local[:, 2]
    )
    return (
        torch.bmm(rotation, outgoing_relative_local.unsqueeze(-1)).squeeze(-1)
        + racket_velocity
    )


def invert_racket_contact_velocity(
    outgoing_velocity: torch.Tensor,
    incoming_velocity: torch.Tensor,
    angular_velocity: torch.Tensor,
    racket_normal: torch.Tensor,
    *,
    config: BallPhysicsConfig | None = None,
) -> torch.Tensor:
    """Solve the simulated linear contact map for racket velocity."""
    cfg = config or BallPhysicsConfig()
    k = cfg.racket_tangential_transfer
    if not 0.0 < k <= 1.0:
        raise ValueError("racket_tangential_transfer must be in (0, 1]")
    rotation = racket_contact_frame(racket_normal)
    world_to_local = rotation.transpose(1, 2)
    incoming_local = torch.bmm(
        world_to_local, incoming_velocity.unsqueeze(-1)
    ).squeeze(-1)
    outgoing_local = torch.bmm(
        world_to_local, outgoing_velocity.unsqueeze(-1)
    ).squeeze(-1)
    omega_local = torch.bmm(
        world_to_local, angular_velocity.unsqueeze(-1)
    ).squeeze(-1)

    racket_local = outgoing_local.clone()
    racket_local[:, 0] = (
        outgoing_local[:, 0]
        - (1.0 - k) * incoming_local[:, 0]
        - k * cfg.radius_m * omega_local[:, 1]
    ) / k
    racket_local[:, 1] = (
        outgoing_local[:, 1]
        - (1.0 - k) * incoming_local[:, 1]
        + k * cfg.radius_m * omega_local[:, 0]
    ) / k
    racket_local[:, 2] = (
        outgoing_local[:, 2]
        + cfg.racket_restitution * incoming_local[:, 2]
    ) / (1.0 + cfg.racket_restitution)
    return torch.bmm(rotation, racket_local.unsqueeze(-1)).squeeze(-1)


def fit_racket_contact_velocity(
    outgoing_velocity: torch.Tensor,
    incoming_velocity: torch.Tensor,
    angular_velocity: torch.Tensor,
    racket_normal: torch.Tensor,
    minimum_velocity: torch.Tensor,
    maximum_velocity: torch.Tensor,
    *,
    config: BallPhysicsConfig | None = None,
    iterations: int = 24,
) -> torch.Tensor:
    """Fit the simulated contact map inside the planner velocity box."""
    cfg = config or BallPhysicsConfig()
    if iterations < 0:
        raise ValueError("iterations must be non-negative")
    if bool((minimum_velocity > maximum_velocity).any()):
        raise ValueError("minimum_velocity must not exceed maximum_velocity")

    command = invert_racket_contact_velocity(
        outgoing_velocity,
        incoming_velocity,
        angular_velocity,
        racket_normal,
        config=cfg,
    )
    command = torch.maximum(
        torch.minimum(command, maximum_velocity), minimum_velocity
    )

    rotation = racket_contact_frame(racket_normal)
    local_response = torch.tensor(
        (
            cfg.racket_tangential_transfer,
            cfg.racket_tangential_transfer,
            1.0 + cfg.racket_restitution,
        ),
        device=command.device,
        dtype=command.dtype,
    ).expand(command.shape[0], -1)
    response = torch.bmm(
        rotation * local_response[:, None, :], rotation.transpose(1, 2)
    )
    step_size = 1.0 / float(local_response.max().square())
    for _ in range(iterations):
        achieved = racket_contact_outgoing_velocity(
            incoming_velocity,
            angular_velocity,
            command,
            racket_normal,
            config=cfg,
        )
        error = achieved - outgoing_velocity
        gradient = torch.bmm(
            response.transpose(1, 2), error.unsqueeze(-1)
        ).squeeze(-1)
        command = command - step_size * gradient
        command = torch.maximum(
            torch.minimum(command, maximum_velocity), minimum_velocity
        )
    return command


class BallPhysicsWorld:
    """Batched kinematic world using Duerr et al. (2026), equations 1-6."""

    def __init__(self, config: BallPhysicsConfig | None = None):
        self.cfg = config or BallPhysicsConfig()

    def acceleration(self, state: BallState, parameters: DomainParameters) -> torch.Tensor:
        velocity = state.velocity
        omega = state.angular_velocity
        speed = torch.linalg.vector_norm(velocity, dim=-1, keepdim=True)
        omega_norm = torch.linalg.vector_norm(omega, dim=-1, keepdim=True)
        safe_speed = torch.clamp(speed, min=1e-6)
        c_m = 0.1 * self.cfg.radius_m * omega_norm / safe_speed - 0.001

        area = math.pi * self.cfg.radius_m**2
        drag_force = (
            -0.5
            * self.cfg.drag_coefficient
            * parameters.drag_scale[:, None]
            * self.cfg.air_density_kg_m3
            * area
            * speed
            * velocity
        )
        magnus_force = (
            -(4.0 / 3.0)
            * c_m
            * parameters.magnus_scale[:, None]
            * self.cfg.air_density_kg_m3
            * math.pi
            * self.cfg.radius_m**3
            * torch.cross(velocity, omega, dim=-1)
        )
        acceleration = (drag_force + magnus_force) / self.cfg.mass_kg
        acceleration[:, 2] -= self.cfg.gravity_m_s2
        return torch.where(state.active[:, None], acceleration, torch.zeros_like(acceleration))

    def _table_bounce(
        self,
        velocity: torch.Tensor,
        omega: torch.Tensor,
        parameters: DomainParameters,
        mask: torch.Tensor,
        generator: torch.Generator | None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        if not mask.any():
            return velocity, omega
        ids = torch.nonzero(mask, as_tuple=False).flatten()
        v_minus = velocity[ids]
        w_minus = omega[ids]
        mu = parameters.table_mu[ids]
        eps = parameters.restitution_intercept[ids] - parameters.restitution_slope[ids] * v_minus[:, 2]
        eps = eps.clamp(self.cfg.restitution_min, self.cfg.restitution_max)

        v_t = torch.stack(
            (v_minus[:, 0] - self.cfg.radius_m * w_minus[:, 1],
             v_minus[:, 1] + self.cfg.radius_m * w_minus[:, 0]),
            dim=-1,
        )
        v_t_norm = torch.linalg.vector_norm(v_t, dim=-1).clamp_min(1e-8)
        sliding_alpha = mu * (1.0 + eps) * torch.abs(v_minus[:, 2]) / v_t_norm
        alpha = torch.minimum(sliding_alpha, torch.full_like(sliding_alpha, 0.4))

        v_plus = v_minus.clone()
        v_plus[:, 0] = (1.0 - alpha) * v_minus[:, 0] + alpha * self.cfg.radius_m * w_minus[:, 1]
        v_plus[:, 1] = (1.0 - alpha) * v_minus[:, 1] - alpha * self.cfg.radius_m * w_minus[:, 0]
        v_plus[:, 2] = -eps * v_minus[:, 2]

        w_plus = w_minus.clone()
        coupling = 3.0 * alpha / (2.0 * self.cfg.radius_m)
        w_plus[:, 0] = -coupling * v_minus[:, 1] + (1.0 - 1.5 * alpha) * w_minus[:, 0]
        w_plus[:, 1] = coupling * v_minus[:, 0] + (1.0 - 1.5 * alpha) * w_minus[:, 1]

        linear_sigma = torch.tensor(
            self.cfg.post_bounce_linear_sigma_m_s, device=velocity.device, dtype=velocity.dtype
        )
        angular_sigma = torch.tensor(
            self.cfg.post_bounce_angular_sigma_rad_s, device=velocity.device, dtype=velocity.dtype
        )
        stochastic = parameters.position_noise_sigma[ids] > 0.0
        if stochastic.any():
            v_noise = torch.randn(v_plus.shape, device=velocity.device, dtype=velocity.dtype, generator=generator)
            w_noise = torch.randn(w_plus.shape, device=velocity.device, dtype=velocity.dtype, generator=generator)
            v_plus = v_plus + v_noise * linear_sigma * stochastic[:, None]
            w_plus = w_plus + w_noise * angular_sigma * stochastic[:, None]

        new_velocity = velocity.clone()
        new_omega = omega.clone()
        new_velocity[ids] = v_plus
        new_omega[ids] = w_plus
        return new_velocity, new_omega

    def step(
        self,
        state: BallState,
        parameters: DomainParameters,
        dt: float,
        *,
        generator: torch.Generator | None = None,
    ) -> tuple[BallState, BallEvents]:
        previous_position = state.position
        acceleration = self.acceleration(state, parameters)
        next_velocity = state.velocity + acceleration * dt
        next_position = state.position + 0.5 * (state.velocity + next_velocity) * dt
        next_omega = state.angular_velocity.clone()

        crossed_surface = (
            state.active
            & (previous_position[:, 2] > self.cfg.table_height_m)
            & (next_position[:, 2] <= self.cfg.table_height_m)
            & (next_velocity[:, 2] < 0.0)
        )
        denominator = next_position[:, 2] - previous_position[:, 2]
        fraction = (self.cfg.table_height_m - previous_position[:, 2]) / torch.where(
            denominator.abs() > 1e-9, denominator, torch.ones_like(denominator)
        )
        contact = previous_position + fraction[:, None].clamp(0.0, 1.0) * (next_position - previous_position)
        inside_table = (
            (contact[:, 0].abs() <= self.cfg.table_half_length_m)
            & (contact[:, 1].abs() <= self.cfg.table_half_width_m)
        )
        bounce = crossed_surface & inside_table
        next_velocity, next_omega = self._table_bounce(
            next_velocity, next_omega, parameters, bounce, generator
        )
        if bounce.any():
            next_position = next_position.clone()
            next_position[bounce, 2] = self.cfg.table_height_m + 1e-6

        next_position = torch.where(state.active[:, None], next_position, state.position)
        next_velocity = torch.where(state.active[:, None], next_velocity, state.velocity)
        next_omega = torch.where(state.active[:, None], next_omega, state.angular_velocity)
        events = BallEvents(table_bounce=bounce, table_contact_position=contact)
        return BallState(next_position, next_velocity, next_omega, state.active.clone()), events

    def apply_racket_contact(
        self,
        state: BallState,
        racket_position: torch.Tensor,
        racket_velocity: torch.Tensor,
        racket_normal: torch.Tensor,
        already_hit: torch.Tensor,
        *,
        previous_ball_position: torch.Tensor | None = None,
        previous_racket_position: torch.Tensor | None = None,
    ) -> tuple[BallState, torch.Tensor, torch.Tensor]:
        """Apply the base linear racket-contact structure from equations 7-8.

        The paper-fitted restitution and residual MLP are unavailable. The
        configured restitution and tangential coefficient are provisional and
        intentionally exposed as world parameters for later calibration.
        """

        normal = racket_normal / torch.linalg.vector_norm(racket_normal, dim=-1, keepdim=True).clamp_min(1e-8)
        relative = state.velocity - racket_velocity
        normal_speed = torch.sum(relative * normal, dim=-1)
        wrong_side = normal_speed > 0.0
        normal = torch.where(wrong_side[:, None], -normal, normal)
        normal_speed = torch.sum(relative * normal, dim=-1)
        contact_racket_position = racket_position
        if previous_ball_position is not None and previous_racket_position is not None:
            relative_start = previous_ball_position - previous_racket_position
            relative_end = state.position - racket_position
            relative_delta = relative_end - relative_start
            closest_fraction = -torch.sum(relative_start * relative_delta, dim=-1) / torch.sum(
                relative_delta * relative_delta, dim=-1
            ).clamp_min(1e-12)
            closest_fraction = closest_fraction.clamp(0.0, 1.0)
            closest_relative = relative_start + closest_fraction[:, None] * relative_delta
            distance = torch.linalg.vector_norm(closest_relative, dim=-1)
            contact_racket_position = previous_racket_position + closest_fraction[:, None] * (
                racket_position - previous_racket_position
            )
        else:
            distance = torch.linalg.vector_norm(state.position - racket_position, dim=-1)
        contact = (
            state.active
            & (~already_hit)
            & (distance <= self.cfg.racket_contact_distance_m)
            & (normal_speed <= -0.05)
        )
        if not contact.any():
            return state, contact, distance

        ids = torch.nonzero(contact, as_tuple=False).flatten()
        n = normal[ids]
        rel = relative[ids]
        # Build a stable local contact frame whose z axis is the racket normal.
        rotation = racket_contact_frame(n)
        rel_local = torch.bmm(rotation.transpose(1, 2), rel.unsqueeze(-1)).squeeze(-1)
        omega_local = torch.bmm(
            rotation.transpose(1, 2), state.angular_velocity[ids].unsqueeze(-1)
        ).squeeze(-1)

        k = self.cfg.racket_tangential_transfer
        radius = self.cfg.radius_m
        alpha = torch.full((ids.shape[0],), k, device=n.device, dtype=n.dtype)
        v_plus_local = rel_local.clone()
        v_plus_local[:, 0] = (1.0 - alpha) * rel_local[:, 0] + alpha * radius * omega_local[:, 1]
        v_plus_local[:, 1] = (1.0 - alpha) * rel_local[:, 1] - alpha * radius * omega_local[:, 0]
        v_plus_local[:, 2] = -self.cfg.racket_restitution * rel_local[:, 2]
        coupling = 3.0 * alpha / (2.0 * radius)
        omega_plus_local = omega_local.clone()
        omega_plus_local[:, 0] = -coupling * rel_local[:, 1] + (1.0 - 1.5 * alpha) * omega_local[:, 0]
        omega_plus_local[:, 1] = coupling * rel_local[:, 0] + (1.0 - 1.5 * alpha) * omega_local[:, 1]

        next_state = state.clone()
        next_state.velocity[ids] = (
            torch.bmm(rotation, v_plus_local.unsqueeze(-1)).squeeze(-1) + racket_velocity[ids]
        )
        next_state.angular_velocity[ids] = torch.bmm(
            rotation, omega_plus_local.unsqueeze(-1)
        ).squeeze(-1)
        next_state.position[ids] = contact_racket_position[ids] + n * self.cfg.racket_ball_separation_m
        return next_state, contact, distance
