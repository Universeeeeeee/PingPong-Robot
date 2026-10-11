import torch

from pingpong_physics.contracts import BallState
from pingpong_physics.dynamics import (
    BallPhysicsWorld,
    fit_racket_contact_velocity,
    invert_racket_contact_velocity,
    racket_contact_outgoing_velocity,
)
from pingpong_physics.randomization import DomainRandomizer


def test_flight_has_drag_and_gravity() -> None:
    params = DomainRandomizer().sample("nominal", 1)
    state = BallState(
        position=torch.tensor([[0.0, 0.0, 1.0]]),
        velocity=torch.tensor([[5.0, 0.0, 0.0]]),
        angular_velocity=torch.zeros((1, 3)),
        active=torch.ones(1, dtype=torch.bool),
    )
    next_state, events = BallPhysicsWorld().step(state, params, 0.01)
    assert next_state.velocity[0, 0] < 5.0
    assert next_state.velocity[0, 2] < 0.0
    assert not events.table_bounce[0]


def test_table_contact_reverses_vertical_velocity() -> None:
    params = DomainRandomizer().sample("nominal", 1)
    state = BallState(
        position=torch.tensor([[0.0, 0.0, 0.025]]),
        velocity=torch.tensor([[1.0, 0.0, -1.0]]),
        angular_velocity=torch.zeros((1, 3)),
        active=torch.ones(1, dtype=torch.bool),
    )
    next_state, events = BallPhysicsWorld().step(state, params, 0.01)
    assert events.table_bounce[0]
    assert next_state.velocity[0, 2] > 0.0
    assert next_state.position[0, 2] > 0.02


def test_spin_changes_flight_direction() -> None:
    params = DomainRandomizer().sample("nominal", 2)
    state = BallState(
        position=torch.tensor([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
        velocity=torch.tensor([[5.0, 0.0, 0.0], [5.0, 0.0, 0.0]]),
        angular_velocity=torch.tensor([[0.0, 0.0, 200.0], [0.0, 0.0, -200.0]]),
        active=torch.ones(2, dtype=torch.bool),
    )
    acceleration = BallPhysicsWorld().acceleration(state, params)
    assert acceleration[0, 1] * acceleration[1, 1] < 0.0


def test_racket_contact_marks_hit_and_reverses_normal_component() -> None:
    state = BallState(
        position=torch.tensor([[0.0, 0.0, 0.0]]),
        velocity=torch.tensor([[-2.0, 0.0, 0.0]]),
        angular_velocity=torch.zeros((1, 3)),
        active=torch.ones(1, dtype=torch.bool),
    )
    next_state, hit, distance = BallPhysicsWorld().apply_racket_contact(
        state,
        racket_position=torch.zeros((1, 3)),
        racket_velocity=torch.zeros((1, 3)),
        racket_normal=torch.tensor([[1.0, 0.0, 0.0]]),
        already_hit=torch.zeros(1, dtype=torch.bool),
    )
    assert hit[0]
    assert distance[0] == 0.0
    assert next_state.velocity[0, 0] > 0.0


def test_swept_contact_catches_between_tick_crossing() -> None:
    state = BallState(
        position=torch.tensor([[0.1, 0.0, 0.0]]),
        velocity=torch.tensor([[-2.0, 0.0, 0.0]]),
        angular_velocity=torch.zeros((1, 3)),
        active=torch.ones(1, dtype=torch.bool),
    )
    _, hit, distance = BallPhysicsWorld().apply_racket_contact(
        state,
        racket_position=torch.zeros((1, 3)),
        racket_velocity=torch.zeros((1, 3)),
        racket_normal=torch.tensor([[1.0, 0.0, 0.0]]),
        already_hit=torch.zeros(1, dtype=torch.bool),
        previous_ball_position=torch.tensor([[-0.1, 0.0, 0.0]]),
        previous_racket_position=torch.zeros((1, 3)),
    )
    assert hit[0]
    assert distance[0] < 1e-6


def test_sim_contact_inverse_recovers_desired_outgoing_velocity() -> None:
    incoming = torch.tensor([[-2.2, 0.3, -0.6], [-1.4, -0.2, 0.7]])
    outgoing = torch.tensor([[3.1, -0.1, 2.0], [2.2, 0.4, 1.5]])
    omega = torch.tensor([[10.0, -30.0, 80.0], [-20.0, 15.0, -40.0]])
    normal = outgoing - incoming
    racket = invert_racket_contact_velocity(
        outgoing, incoming, omega, normal
    )
    recovered = racket_contact_outgoing_velocity(
        incoming, omega, racket, normal
    )
    torch.testing.assert_close(recovered, outgoing, atol=1e-5, rtol=1e-5)


def test_bounded_sim_contact_fit_improves_over_clipped_inverse() -> None:
    incoming = torch.tensor([[-2.2, 0.3, -0.6]])
    outgoing = torch.tensor([[3.1, -0.1, 2.0]])
    omega = torch.zeros((1, 3))
    normal = outgoing - incoming
    minimum = torch.tensor([[0.9, -0.3, 0.01]])
    maximum = torch.tensor([[2.0, 0.3, 0.3]])
    inverse = invert_racket_contact_velocity(
        outgoing, incoming, omega, normal
    )
    clipped = torch.maximum(torch.minimum(inverse, maximum), minimum)
    fitted = fit_racket_contact_velocity(
        outgoing, incoming, omega, normal, minimum, maximum
    )
    clipped_outgoing = racket_contact_outgoing_velocity(
        incoming, omega, clipped, normal
    )
    fitted_outgoing = racket_contact_outgoing_velocity(
        incoming, omega, fitted, normal
    )
    clipped_error = torch.linalg.vector_norm(clipped_outgoing - outgoing)
    fitted_error = torch.linalg.vector_norm(fitted_outgoing - outgoing)

    assert fitted_error < clipped_error
    assert bool((fitted >= minimum).all())
    assert bool((fitted <= maximum).all())
