"""Validate observation layouts without importing Isaac Lab or the training package."""

import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest


path = (Path(__file__).resolve().parents[1] / "source" / "whole_body_tracking"
        / "whole_body_tracking" / "tasks" / "tracking" / "actor_observation_contract.py")
spec = importlib.util.spec_from_file_location("batch_actor_contract", path)
contracts = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = contracts
spec.loader.exec_module(contracts)


def environment(layout, total):
    return SimpleNamespace(observation_manager=SimpleNamespace(
        active_terms={"policy": [name for name, _ in layout]},
        group_obs_term_dim={"policy": [(dim,) for _, dim in layout]},
        group_obs_dim={"policy": (total,)},
    ))


def test_hitter_pure_contract_accepts_canonical_layout():
    contract = contracts.resolve_actor_observation_contract("hitter_pure")
    assert contract.total_dim == 110
    assert sum(dim for _, dim in contract.layout) == 110
    env = environment(contract.layout, 110)
    assert contracts.validate_actor_observation_contract(env, "hitter_pure") is contract
    assert contracts.infer_actor_observation_contract(env) is contract


def test_equal_size_reordered_observations_are_rejected():
    contract = contracts.HITTER_PURE
    layout = list(contract.layout)
    layout[1], layout[2] = layout[2], layout[1]
    env = environment(layout, 110)
    with pytest.raises(ValueError, match="contract mismatch"):
        contracts.validate_actor_observation_contract(env, "hitter_pure")
    assert contracts.infer_actor_observation_contract(env) is None


def test_wrong_total_dimension_is_rejected():
    with pytest.raises(ValueError, match="contract mismatch"):
        contracts.validate_actor_observation_contract(
            environment(contracts.HITTER_PURE.layout, 109), "hitter_pure")


def test_legacy_alias_and_unknown_contract():
    assert contracts.resolve_actor_observation_contract("real_sensor_only") is contracts.DEPLOY_PARITY
    with pytest.raises(ValueError, match="Unknown actor observation contract"):
        contracts.resolve_actor_observation_contract("unknown_contract")
