import torch

from pingpong_physics.randomization import DomainRandomizer


def test_nominal_is_deterministic() -> None:
    randomizer = DomainRandomizer()
    params = randomizer.sample("nominal", 8)
    assert torch.all(params.drag_scale == 1.0)
    assert torch.all(params.magnus_scale == 1.0)
    assert torch.all(params.table_mu == 0.25)
    assert torch.all(params.position_bias == 0.0)


def test_dr_v0_stays_in_declared_ranges() -> None:
    generator = torch.Generator().manual_seed(7)
    params = DomainRandomizer().sample("dr_v0", 4096, generator=generator)
    assert torch.all((params.drag_scale >= 0.85) & (params.drag_scale <= 1.15))
    assert torch.all((params.magnus_scale >= 0.70) & (params.magnus_scale <= 1.30))
    assert torch.all((params.table_mu >= 0.20) & (params.table_mu <= 0.30))
    assert torch.all((params.dropout_probability >= 0.0) & (params.dropout_probability <= 0.02))


def test_sample_spin_is_bounded() -> None:
    generator = torch.Generator().manual_seed(9)
    randomizer = DomainRandomizer()
    params = randomizer.sample("dr_v0", 1024, generator=generator)
    spin = randomizer.sample_spin(params, generator=generator)
    assert spin.shape == (1024, 3)
    assert torch.allclose(torch.linalg.vector_norm(spin, dim=-1), params.initial_spin_magnitude, atol=1e-4)


def test_low_spin_profile_only_reduces_initial_spin_range() -> None:
    randomizer = DomainRandomizer()
    low_generator = torch.Generator().manual_seed(11)
    full_generator = torch.Generator().manual_seed(11)
    low = randomizer.sample("dr_low_spin", 1024, generator=low_generator)
    full = randomizer.sample("dr_v0", 1024, generator=full_generator)

    assert torch.all(low.initial_spin_magnitude <= 60.0)
    torch.testing.assert_close(low.drag_scale, full.drag_scale)
    torch.testing.assert_close(low.magnus_scale, full.magnus_scale)
    torch.testing.assert_close(low.table_mu, full.table_mu)
