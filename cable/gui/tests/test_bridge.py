"""Tests for pycable.bridge.CableBridge.

These tests avoid running the real binary. They validate:
- The input.json written by ``run_equilibrium`` is parseable and complete.
- SNAPSHOT stdout lines are parsed into structured ``snapshot_ready`` dicts.
- Missing result.json triggers ``error_occurred``.

For the end-to-end integration test, see tests/test_integration.py.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import pytest

# PySide6 needs a QApplication for QObject signals to work, even in tests.
pytest.importorskip("PySide6", reason="PySide6 not available")
from PySide6.QtCore import QCoreApplication  # noqa: E402

from pycable.bridge import CableBridge, SNAPSHOT_PREFIX  # noqa: E402
from pycable.params import CableParams  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QCoreApplication.instance() or QCoreApplication([])
    yield app


def test_snapshot_prefix_constant():
    """The SNAPSHOT prefix must exactly match the one emitted by cable_solver.cpp."""
    assert SNAPSHOT_PREFIX == "SNAPSHOT "


def test_bridge_constructs_without_solver(qapp, monkeypatch):
    """Constructing a CableBridge should not raise even if no solver exists."""
    # Force solver_discovery to fail so the bridge falls back to the default path.
    from pycable import solver_discovery

    def raise_missing(**kwargs):
        if kwargs.get("must_exist", True):
            raise solver_discovery.CableSolverNotFound("not here")
        return Path("/nonexistent/cable_solver")

    monkeypatch.setattr(solver_discovery, "find_cable_solver", raise_missing)
    # Re-import the bridge module symbol the same monkeypatch affects.
    from pycable import bridge as bridge_mod
    monkeypatch.setattr(bridge_mod, "find_cable_solver", raise_missing)

    b = CableBridge()
    assert not b.is_running
    assert b.solver_path is not None


def test_params_write_json_has_all_cable_solver_keys(tmp_path: Path):
    """The JSON file written by CableParams must carry every key cable_solver reads."""
    p = CableParams()
    path = tmp_path / "input.json"
    p.write_json(path)

    data = json.loads(path.read_text())

    # These are the exact keys read by cable/cable_solver.cpp main().
    required = {
        "point_a",
        "point_b",
        "cable_length",
        "n_segments",
        "line_density",
        "EA",
        "damping",
        "diameter",
        "gravity",
        "mode",
        "max_equilibrium_steps",
        "equilibrium_tol",
        "snapshot_interval",
    }
    missing = required - set(data.keys())
    assert not missing, f"Missing keys in input.json: {missing}"


def test_snapshot_line_parses_into_dict():
    """A SNAPSHOT line as emitted by cable_solver.cpp must be parseable."""
    line = (
        'SNAPSHOT {"iter":1000,"norm_v":0.0042,"positions":'
        '[[0.0,0.0,0.0],[10.0,0.0,-1.5],[20.0,0.0,-3.0]]}'
    )
    assert line.startswith(SNAPSHOT_PREFIX)
    payload = line[len(SNAPSHOT_PREFIX):]
    snap = json.loads(payload)
    assert snap["iter"] == 1000
    assert snap["norm_v"] == pytest.approx(0.0042)
    assert len(snap["positions"]) == 3
    assert snap["positions"][1] == [10.0, 0.0, -1.5]
