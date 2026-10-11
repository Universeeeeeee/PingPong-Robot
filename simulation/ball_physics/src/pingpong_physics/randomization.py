from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import torch
import yaml


@dataclass
class DomainParameters:
    drag_scale: torch.Tensor
    magnus_scale: torch.Tensor
    table_mu: torch.Tensor
    restitution_intercept: torch.Tensor
    restitution_slope: torch.Tensor
    initial_spin_magnitude: torch.Tensor
    position_noise_sigma: torch.Tensor
    velocity_noise_sigma: torch.Tensor
    latency_s: torch.Tensor
    timestamp_jitter_s: torch.Tensor
    dropout_probability: torch.Tensor
    position_bias: torch.Tensor

    def to_log_dict(self, env_id: int = 0) -> dict[str, float | list[float]]:
        result: dict[str, float | list[float]] = {}
        for key, value in vars(self).items():
            sample = value[env_id].detach().cpu()
            result[key] = sample.tolist() if sample.ndim else float(sample)
        return result


class DomainRandomizer:
    """Samples episode-static parameters. Physics values never change per tick."""

    def __init__(self, config_path: str | Path | None = None):
        if config_path is None:
            config_path = Path(__file__).resolve().parent / "domain_randomization.yaml"
        with Path(config_path).open("r", encoding="utf-8") as stream:
            self._profiles = yaml.safe_load(stream)["profiles"]

    @property
    def profiles(self) -> tuple[str, ...]:
        return tuple(self._profiles)

    @staticmethod
    def _uniform(
        bounds: list[float],
        shape: tuple[int, ...],
        *,
        device: torch.device,
        dtype: torch.dtype,
        generator: torch.Generator | None,
    ) -> torch.Tensor:
        low, high = float(bounds[0]), float(bounds[1])
        if low == high:
            return torch.full(shape, low, device=device, dtype=dtype)
        return low + (high - low) * torch.rand(shape, device=device, dtype=dtype, generator=generator)

    def sample(
        self,
        profile: str,
        num_envs: int,
        *,
        device: str | torch.device = "cpu",
        dtype: torch.dtype = torch.float32,
        generator: torch.Generator | None = None,
    ) -> DomainParameters:
        if profile not in self._profiles:
            raise ValueError(f"unknown profile {profile!r}; expected one of {self.profiles}")
        cfg = self._profiles[profile]
        device = torch.device(device)

        def scalar(name: str) -> torch.Tensor:
            return self._uniform(cfg[name], (num_envs,), device=device, dtype=dtype, generator=generator)

        bias_sigma = torch.full((num_envs, 1), 0.003 if profile != "nominal" else 0.0, device=device, dtype=dtype)
        bias = torch.randn((num_envs, 3), device=device, dtype=dtype, generator=generator) * bias_sigma
        return DomainParameters(
            drag_scale=scalar("drag_scale"),
            magnus_scale=scalar("magnus_scale"),
            table_mu=scalar("table_mu"),
            restitution_intercept=scalar("restitution_intercept"),
            restitution_slope=scalar("restitution_slope"),
            initial_spin_magnitude=scalar("spin_magnitude_rad_s"),
            position_noise_sigma=scalar("position_noise_m"),
            velocity_noise_sigma=scalar("velocity_noise_m_s"),
            latency_s=scalar("latency_s"),
            timestamp_jitter_s=scalar("timestamp_jitter_s"),
            dropout_probability=scalar("dropout_probability"),
            position_bias=bias,
        )

    @staticmethod
    def sample_spin(
        parameters: DomainParameters,
        *,
        generator: torch.Generator | None = None,
    ) -> torch.Tensor:
        """Mixture: 70% top/back spin, 20% side spin, 10% mixed."""

        magnitude = parameters.initial_spin_magnitude
        num_envs = magnitude.shape[0]
        device, dtype = magnitude.device, magnitude.dtype
        family = torch.rand((num_envs,), device=device, dtype=dtype, generator=generator)
        signs = torch.where(
            torch.rand((num_envs,), device=device, dtype=dtype, generator=generator) < 0.5,
            -torch.ones((num_envs,), device=device, dtype=dtype),
            torch.ones((num_envs,), device=device, dtype=dtype),
        )
        spin = torch.zeros((num_envs, 3), device=device, dtype=dtype)
        top = family < 0.7
        side = (family >= 0.7) & (family < 0.9)
        mixed = family >= 0.9
        spin[top, 1] = signs[top] * magnitude[top]
        spin[side, 2] = signs[side] * magnitude[side]
        if mixed.any():
            direction = torch.randn((int(mixed.sum()), 3), device=device, dtype=dtype, generator=generator)
            direction = direction / torch.clamp(torch.linalg.vector_norm(direction, dim=-1, keepdim=True), min=1e-8)
            spin[mixed] = direction * magnitude[mixed, None]
        return spin
