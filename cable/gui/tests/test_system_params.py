"""Tests for LumpedCableSystemParams (Phase 4 multi-cable support)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pycable.params import (
    CableParams,
    CableSpec,
    LumpedCableSystemParams,
    DEFAULT_SINGLE_CABLE_NAME,
)


# ---------------------------------------------------------------------------
# CableSpec basic
# ---------------------------------------------------------------------------

def test_cable_spec_n_segments_property():
    spec = CableSpec(n_points=41)
    assert spec.n_segments == 40


def test_cable_spec_to_flat_array_round_trip():
    spec = CableSpec(
        name="C03",
        point_a=(-75.708, 0.0, 0.0),
        point_b=(0.0, 0.0, 36.2),
        cable_length=83.917,
        n_points=31,
        line_density=63.2514,
        EA=1.530925e9,
        damping=0.5,
        diameter=0.10129,
    )
    flat = spec.to_flat_array()
    assert len(flat) == 13
    assert flat[0] == "C03"
    assert flat[7] == 83.917
    assert flat[8] == 31

    decoded = CableSpec.from_flat_array(flat)
    assert decoded == spec


def test_cable_spec_from_flat_array_rejects_wrong_length():
    with pytest.raises(ValueError, match="13 elements"):
        CableSpec.from_flat_array([1, 2, 3])


def test_cable_spec_from_cable_params():
    legacy = CableParams(cable_length=100.0, n_segments=20, EA=1e9)
    spec = CableSpec.from_cable_params(legacy, name="C99")
    assert spec.name == "C99"
    assert spec.cable_length == 100.0
    assert spec.n_points == 21
    assert spec.EA == 1e9


def test_cable_spec_to_cable_params_round_trip():
    legacy_in = CableParams(
        point_a=(1.0, 2.0, 3.0),
        point_b=(4.0, 5.0, 6.0),
        cable_length=42.0,
        n_segments=10,
        EA=2.5e9,
        damping=0.7,
        diameter=0.05,
        gravity=9.81,
        max_equilibrium_steps=200000,
        equilibrium_tol=0.005,
        snapshot_interval=2000,
    )
    spec = CableSpec.from_cable_params(legacy_in)
    legacy_out = spec.to_cable_params(
        gravity=legacy_in.gravity,
        mode=legacy_in.mode,
        max_equilibrium_steps=legacy_in.max_equilibrium_steps,
        equilibrium_tol=legacy_in.equilibrium_tol,
        snapshot_interval=legacy_in.snapshot_interval,
    )
    assert legacy_in == legacy_out


# ---------------------------------------------------------------------------
# LumpedCableSystemParams: legacy single-line round-trip
# ---------------------------------------------------------------------------

def test_system_from_single_legacy_dict_promotes_to_one_cable():
    """Reading a legacy single-line JSON dict yields a 1-cable system."""
    legacy_dict = {
        "point_a": [500.0, 0.0, -58.0],
        "point_b": [0.0, 0.0, 0.0],
        "cable_length": 522.0,
        "n_segments": 40,
        "line_density": 348.5,
        "EA": 1.4e9,
        "damping": 0.5,
        "diameter": 0.132,
        "gravity": 9.81,
        "mode": "equilibrium",
        "max_equilibrium_steps": 500000,
        "equilibrium_tol": 0.01,
        "snapshot_interval": 10000,
    }
    sys = LumpedCableSystemParams.from_dict(legacy_dict)
    assert len(sys.cables) == 1
    assert sys.cables[0].name == DEFAULT_SINGLE_CABLE_NAME
    assert sys.cables[0].cable_length == 522.0
    assert sys.cables[0].n_points == 41
    assert sys.gravity == 9.81
    assert sys.is_single_legacy()


def test_system_single_legacy_writes_legacy_schema(tmp_path: Path):
    """A 1-cable default-name system writes flat ``point_a``/``point_b`` keys."""
    sys = LumpedCableSystemParams(
        cables=[CableSpec(name=DEFAULT_SINGLE_CABLE_NAME,
                          point_a=(1.0, 2.0, 3.0),
                          point_b=(4.0, 5.0, 6.0),
                          cable_length=10.0,
                          n_points=11)]
    )
    p = tmp_path / "single.json"
    sys.write_json(p)
    raw = json.loads(p.read_text())
    assert "point_a" in raw
    assert raw["cable_length"] == 10.0
    assert raw["n_segments"] == 10
    assert "mooring_cable" not in raw


def test_system_legacy_round_trip_via_disk(tmp_path: Path):
    sys_in = LumpedCableSystemParams.from_cable_params(CableParams())
    p = tmp_path / "legacy.json"
    sys_in.write_json(p)
    sys_out = LumpedCableSystemParams.read_json(p)
    assert len(sys_out.cables) == 1
    assert sys_out.cables[0] == sys_in.cables[0]
    assert sys_out.gravity == sys_in.gravity


# ---------------------------------------------------------------------------
# LumpedCableSystemParams: multi-line round-trip
# ---------------------------------------------------------------------------

def test_system_multi_writes_mooring_keys(tmp_path: Path):
    """A 2+ cable system writes ``mooring_<name>`` flat arrays."""
    sys = LumpedCableSystemParams(
        cables=[
            CableSpec(name="A", point_a=(0, 0, 0), point_b=(10, 0, 5), cable_length=11.0, n_points=11),
            CableSpec(name="B", point_a=(0, 0, 0), point_b=(0, 10, 5), cable_length=11.5, n_points=11),
        ],
        gravity=9.81,
    )
    p = tmp_path / "multi.json"
    sys.write_json(p)
    raw = json.loads(p.read_text())
    assert "mooring_A" in raw
    assert "mooring_B" in raw
    assert isinstance(raw["mooring_A"], list)
    assert len(raw["mooring_A"]) == 13
    assert raw["mooring_A"][0] == "A"
    assert raw["gravity"] == 9.81
    # Make sure the legacy single-line keys are NOT present:
    assert "point_a" not in raw
    assert "cable_length" not in raw


def test_system_multi_round_trip_via_disk(tmp_path: Path):
    sys_in = LumpedCableSystemParams(
        cables=[
            CableSpec(name="C01", point_a=(-104.1, 0, 0), point_b=(0, 0, 45),
                      cable_length=113.41, n_points=31, line_density=63.2514,
                      EA=1.530925e9, damping=0.5, diameter=0.10129),
            CableSpec(name="C02", point_a=(-89.904, 0, 0), point_b=(0, 0, 40.6),
                      cable_length=98.646, n_points=31, line_density=63.2514,
                      EA=1.530925e9, damping=0.5, diameter=0.10129),
        ],
        gravity=9.81,
        max_equilibrium_steps=500000,
        equilibrium_tol=0.01,
        snapshot_interval=5000,
    )
    p = tmp_path / "multi.json"
    sys_in.write_json(p)
    sys_out = LumpedCableSystemParams.read_json(p)
    assert len(sys_out.cables) == 2
    assert {c.name for c in sys_out.cables} == {"C01", "C02"}
    assert sys_out.snapshot_interval == 5000


def test_system_from_dict_detects_mooring_keys():
    """A dict with `mooring_*` keys is detected as multi-line even if it
    also has stray legacy fields."""
    d = {
        "mooring_A": ["A", 0, 0, 0, 10, 0, 5, 11.0, 11, 50.0, 1e9, 0.5, 0.05],
        "mooring_B": ["B", 0, 0, 0, 0, 10, 5, 11.5, 11, 50.0, 1e9, 0.5, 0.05],
        "gravity": 9.81,
    }
    sys = LumpedCableSystemParams.from_dict(d)
    assert len(sys.cables) == 2
    assert not sys.is_single_legacy()


def test_system_yuri_bridge_all_round_trip():
    """The actual yuri_bridge_all.json (12 cables) round-trips through
    LumpedCableSystemParams without losing any cable."""
    p = Path.home() / "Library" / "CloudStorage" / "Dropbox" / "code" / "cpp" \
        / "cable" / "gui" / "examples" / "yuri_bridge" / "inputs" / "yuri_bridge_all.json"
    if not p.exists():
        pytest.skip(f"yuri_bridge_all.json not found at {p}")
    sys = LumpedCableSystemParams.read_json(p)
    assert len(sys.cables) == 12
    names = {c.name for c in sys.cables}
    assert names == {f"C{i:02d}" for i in range(1, 13)}
    # Sanity-check one cable
    c01 = next(c for c in sys.cables if c.name == "C01")
    assert c01.cable_length == 113.41
    assert c01.n_points == 31
