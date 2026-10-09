"""Exercise shipped YAML and the neutral/legacy environment override boundary."""

import importlib.util
from pathlib import Path
import sys

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[3]
path = (Path(__file__).resolve().parents[1] / "source" / "whole_body_tracking"
        / "whole_body_tracking" / "utils" / "success_metric.py")
spec = importlib.util.spec_from_file_location("batch_success_metric_config", path)
metric = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = metric
spec.loader.exec_module(metric)


@pytest.fixture(autouse=True)
def clean_overrides(monkeypatch):
    monkeypatch.delenv("BALL_PHYSICS_CONFIG", raising=False)
    monkeypatch.delenv("HOPE_BALL_PHYSICS_CONFIG", raising=False)


def test_shipped_yaml_is_discovered_and_loaded():
    assert Path(metric._find_config_path()) == ROOT / "configs" / "ball_physics.yaml"
    cfg = metric.load_ball_physics_config()
    expected = yaml.safe_load((ROOT / "configs" / "ball_physics.yaml").read_text())
    assert cfg == expected
    physics = metric.BallPhysics.from_config()
    table = metric.TableGeometry.from_config()
    assert physics.drag_k == expected["drag"]["k"]
    assert physics.ball_radius == expected["ball"]["radius"]
    assert table.net_height == expected["net"]["height"]


def test_neutral_override_precedes_legacy_and_merges_defaults(tmp_path, monkeypatch):
    new = tmp_path / "new.yaml"
    old = tmp_path / "old.yaml"
    new.write_text("gravity: 8.0\ndrag:\n  k: 0.25\n")
    old.write_text("gravity: 7.0\n")
    monkeypatch.setenv("BALL_PHYSICS_CONFIG", str(new))
    monkeypatch.setenv("HOPE_BALL_PHYSICS_CONFIG", str(old))
    cfg = metric.load_ball_physics_config()
    assert cfg["gravity"] == 8.0
    assert cfg["drag"]["k"] == 0.25
    assert cfg["drag"]["velocity_clip"] == 50.0


def test_legacy_override_remains_compatible(tmp_path, monkeypatch):
    path = tmp_path / "legacy.yaml"
    path.write_text("gravity: 7.0\n")
    monkeypatch.setenv("HOPE_BALL_PHYSICS_CONFIG", str(path))
    assert metric.BallPhysics.from_config().gravity == 7.0


def test_invalid_explicit_override_keeps_existing_default_fallback(tmp_path, monkeypatch):
    monkeypatch.setenv("BALL_PHYSICS_CONFIG", str(tmp_path / "missing.yaml"))
    assert metric._find_config_path() is None
    assert metric.load_ball_physics_config()["gravity"] == 9.81
