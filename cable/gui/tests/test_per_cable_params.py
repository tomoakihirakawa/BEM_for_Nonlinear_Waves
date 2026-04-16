"""Tests for PerCableParams (per-cable JSON format)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pycable.params import CableParams, CableSpec, PerCableParams


# ---------------------------------------------------------------------------
# Basic round-trip
# ---------------------------------------------------------------------------

def test_per_cable_defaults():
    p = PerCableParams()
    assert p.name == "cable"
    assert p.end_a_position == (0.0, 0.0, 0.0)
    assert p.end_a_body == ""
    assert p.end_a_motion == "fixed"
    assert p.n_points == 41
    assert p.mode == "equilibrium"


def test_per_cable_to_dict_minimal():
    p = PerCableParams(name="C01",
                       end_a_position=(-104.1, 0, 0),
                       end_b_position=(0, 0, 45),
                       cable_length=113.41,
                       n_points=31)
    d = p.to_dict()
    assert d["name"] == "C01"
    assert d["end_a_position"] == [-104.1, 0, 0]
    assert d["end_b_position"] == [0, 0, 45]
    assert d["cable_length"] == 113.41
    assert d["n_points"] == 31
    # body keys omitted when empty
    assert "end_a_body" not in d
    assert "end_b_body" not in d
    # motion keys omitted when fixed
    assert "end_a_motion" not in d
    assert "end_b_motion" not in d


def test_per_cable_to_dict_with_body():
    p = PerCableParams(name="C01",
                       end_a_position=(-104.1, 0, 0),
                       end_a_body="deck",
                       end_b_position=(0, 0, 45),
                       end_b_body="tower",
                       cable_length=113.41,
                       n_points=31)
    d = p.to_dict()
    assert d["end_a_body"] == "deck"
    assert d["end_b_body"] == "tower"


def test_per_cable_to_dict_dynamic_motion():
    p = PerCableParams(name="D01",
                       end_a_position=(0, 0, 95),
                       end_b_position=(30, 0, 20),
                       end_b_motion="sinusoidal",
                       end_b_motion_dof="heave",
                       end_b_motion_amplitude=2.0,
                       end_b_motion_frequency=0.1,
                       end_b_motion_phase=0.0,
                       cable_length=80,
                       n_points=31,
                       mode="dynamic",
                       dt=0.01,
                       t_end=10.0)
    d = p.to_dict()
    assert d["end_b_motion"] == "sinusoidal"
    assert d["end_b_motion_dof"] == "heave"
    assert d["end_b_motion_amplitude"] == 2.0
    assert d["dt"] == 0.01
    assert d["t_end"] == 10.0
    # end_a motion is fixed, so no motion keys emitted
    assert "end_a_motion" not in d


def test_per_cable_from_dict_round_trip():
    p_in = PerCableParams(
        name="C03",
        end_a_position=(-75.708, 0, 0),
        end_a_body="deck",
        end_b_position=(0, 0, 36.2),
        end_b_body="tower",
        cable_length=83.917,
        n_points=31,
        line_density=63.2514,
        EA=1.530925e9,
        damping=0.5,
        diameter=0.10129,
    )
    d = p_in.to_dict()
    p_out = PerCableParams.from_dict(d)
    assert p_out.name == p_in.name
    assert p_out.end_a_position == p_in.end_a_position
    assert p_out.end_b_position == p_in.end_b_position
    assert p_out.end_a_body == p_in.end_a_body
    assert p_out.end_b_body == p_in.end_b_body
    assert p_out.cable_length == p_in.cable_length
    assert p_out.n_points == p_in.n_points
    assert p_out.EA == p_in.EA


def test_per_cable_json_round_trip(tmp_path: Path):
    p_in = PerCableParams(
        name="M01",
        end_a_position=(500, 0, -58),
        end_b_position=(0, 0, 0),
        end_b_body="floater",
        cable_length=522,
        n_points=41,
    )
    fp = tmp_path / "M01.json"
    p_in.write_json(fp)

    raw = json.loads(fp.read_text())
    assert "end_a_position" in raw
    assert raw["end_b_body"] == "floater"

    p_out = PerCableParams.read_json(fp)
    assert p_out.name == "M01"
    assert p_out.end_b_body == "floater"
    assert p_out.cable_length == 522


# ---------------------------------------------------------------------------
# Conversion from legacy formats
# ---------------------------------------------------------------------------

def test_from_cable_params():
    legacy = CableParams(
        point_a=(500, 0, -58),
        point_b=(0, 0, 0),
        cable_length=522,
        n_segments=40,
        line_density=348.5,
        EA=1.4e9,
    )
    p = PerCableParams.from_cable_params(legacy, name="M01")
    assert p.name == "M01"
    assert p.end_a_position == (500, 0, -58)
    assert p.end_b_position == (0, 0, 0)
    assert p.cable_length == 522
    assert p.n_points == 41  # n_segments + 1
    assert p.EA == 1.4e9
    assert p.end_a_body == ""
    assert p.end_a_motion == "fixed"


def test_to_cable_params():
    p = PerCableParams(
        name="C01",
        end_a_position=(-104.1, 0, 0),
        end_a_body="deck",
        end_b_position=(0, 0, 45),
        end_b_body="tower",
        cable_length=113.41,
        n_points=31,
        line_density=63.2514,
        EA=1.530925e9,
    )
    legacy = p.to_cable_params()
    assert legacy.point_a == (-104.1, 0, 0)
    assert legacy.point_b == (0, 0, 45)
    assert legacy.cable_length == 113.41
    assert legacy.n_segments == 30
    assert legacy.EA == 1.530925e9


def test_from_cable_spec():
    spec = CableSpec(
        name="C01",
        point_a=(-104.1, 0, 0),
        point_b=(0, 0, 45),
        cable_length=113.41,
        n_points=31,
        line_density=63.2514,
        EA=1.530925e9,
        damping=0.5,
        diameter=0.10129,
    )
    p = PerCableParams.from_cable_spec(spec, gravity=9.81)
    assert p.name == "C01"
    assert p.end_a_position == (-104.1, 0, 0)
    assert p.n_points == 31
    assert p.EA == 1.530925e9


# ---------------------------------------------------------------------------
# Read actual per-cable example file
# ---------------------------------------------------------------------------

def test_read_bridge_c01_example():
    p = Path(__file__).resolve().parents[1] / "examples" / "bridge_C01.json"
    if not p.exists():
        pytest.skip(f"bridge_C01.json not found at {p}")
    pc = PerCableParams.read_json(p)
    assert pc.name == "C01"
    assert pc.end_a_body == "deck"
    assert pc.end_b_body == "tower"
    assert pc.cable_length == 113.41
    assert pc.n_points == 31
