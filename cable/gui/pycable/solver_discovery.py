"""Locate the compiled cable_solver binary.

This module works in two repository layouts:

1. **Dev layout** (``cpp/cable/gui/pycable/solver_discovery.py``):
   the binary is built at ``cpp/cable/build_solver/cable_solver``.

2. **Public layout** (``~/CableDynamics/cable/gui/pycable/solver_discovery.py``
   after ``sync_all_public.sh``): the binary is built at
   ``~/CableDynamics/build/cable_solver``.

Search order:
    1. ``PYCABLE_SOLVER_PATH`` — **strict override**. If set but invalid, raises.
    2. ``$CABLE_DYNAMICS_ROOT/build/cable_solver``
    3. Dev layout: ``<cable>/build_solver/cable_solver``
    4. Public layout: ``<cable>/../build/cable_solver``
    5. ``~/CableDynamics/build/cable_solver`` (home convention)
    6. ``shutil.which("cable_solver")``
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Optional


class CableSolverNotFound(FileNotFoundError):
    """Raised when the cable_solver binary cannot be located."""


_BUILD_INSTRUCTIONS = """\
Build the cable_solver binary first.

Dev tree (cpp/cable/):
    cd <cpp>/cable/build_solver
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j$(sysctl -n hw.logicalcpu)

Public tree (~/CableDynamics/):
    cd ~/CableDynamics && mkdir -p build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j$(sysctl -n hw.logicalcpu)

Or set PYCABLE_SOLVER_PATH to the absolute binary path.
"""


def _fallback_candidates() -> list[Path]:
    """Candidates tried when PYCABLE_SOLVER_PATH is not set.

    Layout assumptions (the path of *this* file is):
        .../cable/gui/pycable/solver_discovery.py
                            parents[0] = pycable/
                            parents[1] = gui/
                            parents[2] = cable/         (cpp/cable or ~/CableDynamics/cable)
                            parents[3] = cable parent   (cpp/ or ~/CableDynamics/)
    """
    here = Path(__file__).resolve()
    cable_dir = here.parents[2]
    cable_parent = here.parents[3]

    candidates: list[Path] = []

    cable_root_env = os.environ.get("CABLE_DYNAMICS_ROOT")
    if cable_root_env:
        candidates.append(Path(cable_root_env).expanduser() / "build" / "cable_solver")

    # Dev layout: cpp/cable/build_solver/cable_solver
    candidates.append(cable_dir / "build_solver" / "cable_solver")
    # Public layout: ~/CableDynamics/build/cable_solver
    candidates.append(cable_parent / "build" / "cable_solver")
    # Home convention fallback
    candidates.append(Path.home() / "CableDynamics" / "build" / "cable_solver")

    return candidates


def _is_executable(p: Path) -> bool:
    return p.is_file() and os.access(p, os.X_OK)


def find_cable_solver(must_exist: bool = True) -> Path:
    """Return the path to the cable_solver binary.

    Parameters
    ----------
    must_exist : bool
        If True (default), raise ``CableSolverNotFound`` when the binary
        is not present. If False, return a best-guess path even if it
        doesn't exist (useful for display in GUI status / errors).
    """
    # 1. Strict override
    env = os.environ.get("PYCABLE_SOLVER_PATH")
    if env:
        p = Path(env).expanduser()
        if _is_executable(p):
            return p.resolve()
        if must_exist:
            raise CableSolverNotFound(
                f"PYCABLE_SOLVER_PATH is set to {p!s} but no executable "
                f"file exists there.\n\n{_BUILD_INSTRUCTIONS}"
            )
        return p.resolve()

    # 2. Fallback candidate locations
    for p in _fallback_candidates():
        if _is_executable(p):
            return p.resolve()

    # 3. PATH
    which = shutil.which("cable_solver")
    if which:
        return Path(which).resolve()

    if not must_exist:
        # Return the conventional build path for display purposes.
        return (Path.home() / "CableDynamics" / "build" / "cable_solver").resolve()

    searched = "\n".join(f"  - {p}" for p in _fallback_candidates())
    raise CableSolverNotFound(
        f"cable_solver binary not found.\n\nSearched:\n"
        f"  - $PYCABLE_SOLVER_PATH (unset)\n"
        f"{searched}\n  - $PATH (shutil.which)\n\n{_BUILD_INSTRUCTIONS}"
    )
