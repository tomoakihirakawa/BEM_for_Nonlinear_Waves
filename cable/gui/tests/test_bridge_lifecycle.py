"""Lifecycle tests for CableBridge that spin up a real QCoreApplication.

These complement ``test_bridge.py`` (pure unit) and ``test_integration.py``
(happy-path end-to-end) by exercising the error / stop paths that are
easy to regress — QProcess fires both ``errorOccurred`` AND ``finished``
for terminated processes, and we have to coalesce those into exactly one
``finished`` emit per run.

Skipped if the real cable_solver binary isn't built.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

pytest.importorskip("PySide6", reason="PySide6 not available")
from PySide6.QtCore import QCoreApplication, QTimer  # noqa: E402

from pycable.bridge import CableBridge  # noqa: E402
from pycable.params import CableParams  # noqa: E402
from pycable.solver_discovery import find_cable_solver, CableSolverNotFound  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QCoreApplication.instance() or QCoreApplication([])
    yield app


@pytest.fixture(scope="module")
def solver_path() -> Path:
    try:
        return find_cable_solver(must_exist=True)
    except CableSolverNotFound as e:
        pytest.skip(f"cable_solver not built: {e}")


def _run_event_loop(app: QCoreApplication, max_ms: int) -> None:
    """Spin the event loop until app.quit() is called or max_ms elapses."""
    timer = QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(app.quit)
    timer.start(max_ms)
    app.exec()


def test_finished_fires_exactly_once_on_happy_path(qapp, solver_path):
    """Normal run: started fires once, finished fires once, no errors."""
    bridge = CableBridge()
    counts = {"started": 0, "finished": 0, "error": 0}

    bridge.started.connect(lambda: counts.__setitem__("started", counts["started"] + 1))
    bridge.finished.connect(lambda: counts.__setitem__("finished", counts["finished"] + 1))
    bridge.error_occurred.connect(lambda _: counts.__setitem__("error", counts["error"] + 1))
    bridge.finished.connect(qapp.quit)

    example = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "synthetic"
        / "catenary_500m.json"
    )
    params = CableParams.read_json(example)
    QTimer.singleShot(0, lambda: bridge.run_equilibrium(params))

    _run_event_loop(qapp, max_ms=30_000)

    assert counts["started"] == 1, f"started fired {counts['started']} times"
    assert counts["finished"] == 1, f"finished fired {counts['finished']} times"
    assert counts["error"] == 0, f"unexpected {counts['error']} error(s)"


def test_finished_fires_exactly_once_on_user_stop(qapp, solver_path):
    """stop() mid-run: finished fires once, zero errors (not treated as a crash)."""
    bridge = CableBridge()
    counts = {"started": 0, "finished": 0, "error": 0, "snapshots": 0, "result": 0}

    bridge.started.connect(lambda: counts.__setitem__("started", counts["started"] + 1))
    bridge.finished.connect(lambda: counts.__setitem__("finished", counts["finished"] + 1))
    bridge.error_occurred.connect(lambda _: counts.__setitem__("error", counts["error"] + 1))
    bridge.snapshot_ready.connect(lambda _: counts.__setitem__("snapshots", counts["snapshots"] + 1))
    bridge.result_ready.connect(lambda _: counts.__setitem__("result", counts["result"] + 1))
    bridge.finished.connect(qapp.quit)

    # Config that will not converge in the allowed time budget
    params = CableParams(
        point_a=(0.0, 0.0, 0.0),
        point_b=(1000.0, 0.0, 0.0),
        cable_length=1010.0,
        n_segments=100,
        equilibrium_tol=1e-10,
        snapshot_interval=500,
        max_equilibrium_steps=2_000_000,
    )
    QTimer.singleShot(0, lambda: bridge.run_equilibrium(params))
    QTimer.singleShot(200, bridge.stop)

    _run_event_loop(qapp, max_ms=15_000)

    assert counts["started"] == 1
    assert counts["finished"] == 1, f"finished fired {counts['finished']} times, expected 1"
    assert counts["error"] == 0, (
        f"user-requested stop should not emit error_occurred; got {counts['error']}"
    )
    assert counts["result"] == 0, "result_ready should NOT fire on a stopped run"
    assert counts["snapshots"] >= 1, "at least one SNAPSHOT should have been received before stop"


def test_second_run_resets_state(qapp, solver_path):
    """A second run after a first should reset run_ended / stopped flags."""
    bridge = CableBridge()
    done_count = [0]

    def on_finished():
        done_count[0] += 1
        if done_count[0] == 1:
            # Kick off second run once the first finishes
            QTimer.singleShot(50, lambda: bridge.run_equilibrium(params))
        elif done_count[0] == 2:
            qapp.quit()

    bridge.finished.connect(on_finished)

    example = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "synthetic"
        / "bridge_cable_small.json"
    )
    params = CableParams.read_json(example)
    QTimer.singleShot(0, lambda: bridge.run_equilibrium(params))

    _run_event_loop(qapp, max_ms=30_000)

    assert done_count[0] == 2, f"expected 2 finished events across 2 runs, got {done_count[0]}"
