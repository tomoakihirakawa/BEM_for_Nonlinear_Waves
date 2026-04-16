"""End-to-end test that runs the real cable_solver binary.

Skipped if the binary is not built. Checks that the Python package can
drive the solver and that the output JSON has the expected shape.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

import pytest

from pycable.params import CableParams
from pycable.solver_discovery import find_cable_solver, CableSolverNotFound


@pytest.fixture(scope="module")
def solver_path() -> Path:
    try:
        return find_cable_solver(must_exist=True)
    except CableSolverNotFound as e:
        pytest.skip(f"cable_solver not built: {e}")


def test_cable_solver_runs_on_catenary_example(solver_path: Path, tmp_path: Path):
    """Run the binary on examples/synthetic/catenary_500m.json and inspect result.json."""
    # tests/test_integration.py -> parents[1] = gui/, so examples/ is right next door.
    gui_dir = Path(__file__).resolve().parents[1]
    example = gui_dir / "examples" / "synthetic" / "catenary_500m.json"
    assert example.is_file(), f"Example input missing: {example}"

    out_dir = tmp_path / "run_out"
    out_dir.mkdir()

    result = subprocess.run(
        [str(solver_path), str(example), str(out_dir)],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"cable_solver failed (exit {result.returncode})\n"
        f"--- stdout tail ---\n{result.stdout[-2000:]}\n"
        f"--- stderr ---\n{result.stderr}"
    )

    result_json = out_dir / "result.json"
    assert result_json.is_file(), f"result.json not produced in {out_dir}"

    data = json.loads(result_json.read_text())

    # Schema sanity
    assert data["n_nodes"] == 41, "expected 41 nodes for n_segments=40"
    positions = data["positions"]
    assert len(positions) == 41
    assert all(len(p) == 3 for p in positions)

    # Endpoints should be pinned to the input.
    params = CableParams.read_json(example)
    for got, want in zip(positions[0], params.point_a):
        assert abs(got - want) < 1e-6, f"first node drifted from point_a: {got} vs {want}"
    for got, want in zip(positions[-1], params.point_b):
        assert abs(got - want) < 1e-6, f"last node drifted from point_b: {got} vs {want}"

    # Catenary should sag — interior nodes must dip below the straight-line chord.
    zs = [p[2] for p in positions]
    min_z = min(zs)
    chord_min_z = min(params.point_a[2], params.point_b[2])
    assert min_z < chord_min_z, (
        f"no sag: min(z)={min_z} not below chord min z={chord_min_z}"
    )

    # Tensions all non-negative; top > bottom for a hanging chain.
    tensions = data["tensions"]
    assert all(t >= 0 for t in tensions)
    assert data["top_tension"] >= data["bottom_tension"]


def test_cable_solver_stdout_emits_snapshots(solver_path: Path, tmp_path: Path):
    """The SNAPSHOT contract the GUI relies on must be intact."""
    example = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "synthetic"
        / "bridge_cable_small.json"
    )
    out_dir = tmp_path / "run_out"
    out_dir.mkdir()

    result = subprocess.run(
        [str(solver_path), str(example), str(out_dir)],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0

    # At least one SNAPSHOT line must be present and parse as JSON.
    snapshot_lines = [l for l in result.stdout.splitlines() if l.startswith("SNAPSHOT ")]
    assert snapshot_lines, "no SNAPSHOT lines in stdout — GUI live update would break"
    sample = json.loads(snapshot_lines[0][len("SNAPSHOT "):])
    assert "iter" in sample
    assert "norm_v" in sample
    assert "positions" in sample
    assert isinstance(sample["positions"], list)
    assert len(sample["positions"]) >= 2
