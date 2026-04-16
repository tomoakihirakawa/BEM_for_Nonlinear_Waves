"""Tests for pycable.solver_discovery.find_cable_solver."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from pycable import solver_discovery
from pycable.solver_discovery import (
    CableSolverNotFound,
    find_cable_solver,
)


def _make_fake_binary(path: Path) -> Path:
    """Create an executable stub file at ``path``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\necho fake\n")
    os.chmod(path, 0o755)
    return path


def test_pycable_solver_path_honored_when_valid(monkeypatch, tmp_path: Path):
    fake = _make_fake_binary(tmp_path / "cable_solver")
    monkeypatch.setenv("PYCABLE_SOLVER_PATH", str(fake))
    monkeypatch.delenv("CABLE_DYNAMICS_ROOT", raising=False)

    result = find_cable_solver(must_exist=True)
    assert result == fake.resolve()


def test_pycable_solver_path_strict_when_invalid(monkeypatch, tmp_path: Path):
    """Explicit override must NOT silently fall through if the path is bad."""
    bad = tmp_path / "does_not_exist" / "cable_solver"
    monkeypatch.setenv("PYCABLE_SOLVER_PATH", str(bad))
    monkeypatch.delenv("CABLE_DYNAMICS_ROOT", raising=False)

    with pytest.raises(CableSolverNotFound) as excinfo:
        find_cable_solver(must_exist=True)
    assert "PYCABLE_SOLVER_PATH" in str(excinfo.value)


def test_pycable_solver_path_strict_not_raised_when_must_exist_false(
    monkeypatch, tmp_path: Path
):
    """Display mode (``must_exist=False``) never raises."""
    bad = tmp_path / "does_not_exist" / "cable_solver"
    monkeypatch.setenv("PYCABLE_SOLVER_PATH", str(bad))
    monkeypatch.delenv("CABLE_DYNAMICS_ROOT", raising=False)

    # Should not raise.
    result = find_cable_solver(must_exist=False)
    assert result == bad.expanduser().resolve()


def test_cable_dynamics_root_fallback(monkeypatch, tmp_path: Path):
    """CABLE_DYNAMICS_ROOT points at a custom repo layout."""
    monkeypatch.delenv("PYCABLE_SOLVER_PATH", raising=False)
    fake_root = tmp_path / "custom_cable_repo"
    fake = _make_fake_binary(fake_root / "build" / "cable_solver")
    monkeypatch.setenv("CABLE_DYNAMICS_ROOT", str(fake_root))

    result = find_cable_solver(must_exist=True)
    assert result == fake.resolve()


def test_all_missing_raises_with_helpful_message(monkeypatch, tmp_path: Path):
    """When nothing is found, the error mentions how to build the binary."""
    monkeypatch.delenv("PYCABLE_SOLVER_PATH", raising=False)
    monkeypatch.setenv("CABLE_DYNAMICS_ROOT", str(tmp_path / "nope"))

    # Redirect the fallback candidates to non-existent locations by
    # temporarily overriding _fallback_candidates.
    original = solver_discovery._fallback_candidates

    def empty_candidates():
        return [tmp_path / "nope" / "build" / "cable_solver"]

    monkeypatch.setattr(solver_discovery, "_fallback_candidates", empty_candidates)

    # Also stub shutil.which to ensure no system-wide cable_solver confuses us.
    import shutil as _shutil
    monkeypatch.setattr(_shutil, "which", lambda name: None)

    with pytest.raises(CableSolverNotFound) as excinfo:
        find_cable_solver(must_exist=True)

    msg = str(excinfo.value)
    assert "cable_solver binary not found" in msg
    assert "cmake" in msg
