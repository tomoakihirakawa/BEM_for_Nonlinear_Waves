"""Tests for pycable.params.CableParams."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pycable.params import CableParams


def test_defaults_match_cable_solver_defaults():
    """The dataclass defaults must match cable/cable_solver.cpp getX(...) defaults."""
    p = CableParams()
    assert p.point_a == (500.0, 0.0, -58.0)
    assert p.point_b == (0.0, 0.0, 0.0)
    assert p.cable_length == 522.0
    assert p.n_segments == 40
    assert p.line_density == 348.5
    assert p.EA == 1.4e9
    assert p.damping == 0.5
    assert p.diameter == 0.132
    assert p.gravity == 9.81
    assert p.mode == "equilibrium"
    assert p.max_equilibrium_steps == 500_000
    assert p.equilibrium_tol == 0.01
    assert p.snapshot_interval == 10_000


def test_to_dict_produces_natural_json_types():
    """Numbers should be unquoted floats/ints, arrays should be python lists."""
    p = CableParams()
    d = p.to_dict()
    assert isinstance(d["point_a"], list)
    assert all(isinstance(x, float) for x in d["point_a"])
    assert isinstance(d["cable_length"], float)
    assert isinstance(d["n_segments"], int)
    assert isinstance(d["mode"], str)


def test_dict_round_trip():
    p1 = CableParams(cable_length=100.0, n_segments=20, EA=1e9)
    d = p1.to_dict()
    p2 = CableParams.from_dict(d)
    assert p1 == p2


def test_from_dict_fills_missing_keys_with_defaults():
    p = CableParams.from_dict({"cable_length": 100.0})
    defaults = CableParams()
    assert p.cable_length == 100.0
    # Other fields fall back to defaults.
    assert p.n_segments == defaults.n_segments
    assert p.EA == defaults.EA


def test_write_and_read_json_round_trip(tmp_path: Path):
    p1 = CableParams(
        point_a=(1.0, 2.0, 3.0),
        point_b=(4.0, 5.0, 6.0),
        cable_length=42.0,
        n_segments=10,
    )
    path = tmp_path / "input.json"
    p1.write_json(path)
    assert path.is_file()

    # The on-disk form must be human-readable JSON with unquoted numbers.
    raw = json.loads(path.read_text())
    assert raw["cable_length"] == 42.0
    assert raw["point_a"] == [1.0, 2.0, 3.0]

    p2 = CableParams.read_json(path)
    assert p1 == p2


def test_write_json_creates_parent_directory(tmp_path: Path):
    p = CableParams()
    path = tmp_path / "sub" / "dir" / "input.json"
    p.write_json(path)
    assert path.is_file()


def test_vec3_field_is_a_tuple_after_from_dict():
    """from_dict must coerce list inputs back to tuples."""
    p = CableParams.from_dict({"point_a": [1.0, 2.0, 3.0]})
    assert p.point_a == (1.0, 2.0, 3.0)
    assert isinstance(p.point_a, tuple)
