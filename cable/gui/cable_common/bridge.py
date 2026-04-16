"""CableBridge — QProcess wrapper around the C++ cable_solver binary.

The bridge:

1. Writes a ``CableParams`` (single-line) or ``LumpedCableSystemParams``
   (multi-line) instance as ``input.json`` into a temp dir.
2. Launches ``cable_solver input.json tmp_dir/`` via ``QProcess``.
3. Streams stdout line by line. SNAPSHOT lines are parsed and re-emitted
   as structured ``snapshot_ready`` signals for live GUI animation.
   Multi-line snapshots include a ``"cable"`` field naming which line is
   being updated, so a multi-cable view can demux them.
4. On finish, reads ``tmp_dir/result.json`` and emits ``result_ready``.

Signals
-------
started : ()
finished : ()
log_received : (str)
snapshot_ready : (dict)   # single: {"iter", "norm_v", "positions"}
                          # multi:  {"iter", "norm_v", "positions", "cable"}
result_ready : (dict)     # full result.json content as a dict (single or multi)
error_occurred : (str)
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Optional, Union

from PySide6.QtCore import QObject, QProcess, QProcessEnvironment, Signal

from cable_common.params import CableParams, LumpedCableSystemParams
from cable_common.solver_discovery import find_cable_solver, CableSolverNotFound


SNAPSHOT_PREFIX = "SNAPSHOT "


class CableBridge(QObject):
    """Bridges the GUI to the C++ cable_solver binary via QProcess."""

    started = Signal()
    finished = Signal()
    log_received = Signal(str)
    snapshot_ready = Signal(dict)
    result_ready = Signal(dict)
    error_occurred = Signal(str)

    def __init__(self, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self._process: Optional[QProcess] = None
        self._tmp_dir: Optional[str] = None
        self._stdout_log: list[str] = []
        # True if the current run was terminated by stop(). Suppresses the
        # otherwise-noisy "exited with code 15" error on clean cancellation.
        self._stopped_by_user: bool = False
        # True if we've already fired finished for the current run. QProcess
        # emits both errorOccurred and finished for terminated processes, so
        # we have to coalesce them into a single lifecycle event.
        self._run_ended: bool = False
        try:
            self._solver_path: Path = find_cable_solver(must_exist=False)
        except CableSolverNotFound:
            self._solver_path = Path.home() / "CableDynamics" / "build" / "cable_solver"

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def solver_path(self) -> Path:
        return self._solver_path

    @solver_path.setter
    def solver_path(self, path: str | Path) -> None:
        self._solver_path = Path(os.path.expanduser(str(path)))

    @property
    def is_running(self) -> bool:
        return (
            self._process is not None
            and self._process.state() == QProcess.ProcessState.Running
        )

    # ------------------------------------------------------------------
    # Public control
    # ------------------------------------------------------------------

    def run_equilibrium(self, params: CableParams) -> None:
        """Start a solver process for a single cable (legacy form).

        Equivalent to ``run_system(LumpedCableSystemParams.from_cable_params(params))``
        but kept as a separate entry point for backward compatibility with
        callers that still use the legacy ``CableParams`` API directly.
        """
        sys_params = LumpedCableSystemParams.from_cable_params(params)
        self.run_system(sys_params)

    def run_system(self, params: "LumpedCableSystemParams") -> None:
        """Start a solver process for a multi-cable system.

        Single-cable systems with the default name "cable" round-trip
        through the legacy single-line schema (so cable_solver behaves
        bit-for-bit identically to the old `run_equilibrium` path). Larger
        systems use the multi-line ``mooring_<name>`` schema.

        Emits error_occurred on setup failure.
        """
        if self.is_running:
            self.error_occurred.emit("Solver is already running")
            return

        # Refresh the solver path in case the user built the binary after import.
        try:
            self._solver_path = find_cable_solver(must_exist=True)
        except CableSolverNotFound as e:
            self.error_occurred.emit(str(e))
            return

        # Set up a fresh temp directory for this run.
        self._cleanup_tmp_dir()
        self._tmp_dir = tempfile.mkdtemp(prefix="pycable_")
        input_path = Path(self._tmp_dir) / "input.json"
        params.write_json(input_path)

        self._stdout_log = []
        self._stopped_by_user = False
        self._run_ended = False

        self._process = QProcess(self)
        self._process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        env = QProcessEnvironment.systemEnvironment()
        self._process.setProcessEnvironment(env)

        self._process.readyReadStandardOutput.connect(self._on_stdout)
        self._process.readyReadStandardError.connect(self._on_stderr)
        self._process.finished.connect(self._on_finished)
        self._process.errorOccurred.connect(self._on_error)

        args = [str(input_path), str(self._tmp_dir)]
        self._process.start(str(self._solver_path), args)
        self.started.emit()

    def stop(self) -> None:
        """Terminate a running solver, killing if it doesn't exit in 3 s."""
        if self._process and self.is_running:
            self._stopped_by_user = True
            self._process.terminate()
            if not self._process.waitForFinished(3000):
                self._process.kill()

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _on_stdout(self) -> None:
        if not self._process:
            return
        data = bytes(self._process.readAllStandardOutput()).decode("utf-8", errors="replace")
        for line in data.splitlines():
            self._stdout_log.append(line)

            if line.startswith(SNAPSHOT_PREFIX):
                payload = line[len(SNAPSHOT_PREFIX):]
                try:
                    snap = json.loads(payload)
                except json.JSONDecodeError:
                    self.log_received.emit(f"[snapshot parse error] {line}")
                    continue
                self.snapshot_ready.emit(snap)
                iter_val = snap.get("iter", "?")
                norm_v = snap.get("norm_v", float("nan"))
                self.log_received.emit(f"iter {iter_val}  |v|_max={norm_v:.6g}")
            else:
                self.log_received.emit(line)

    def _on_stderr(self) -> None:
        if not self._process:
            return
        data = bytes(self._process.readAllStandardError()).decode("utf-8", errors="replace")
        for line in data.splitlines():
            self._stdout_log.append(f"[stderr] {line}")
            self.log_received.emit(f"[stderr] {line}")

    def _on_finished(self, exit_code: int, exit_status: QProcess.ExitStatus) -> None:
        if self._run_ended:
            return
        self._run_ended = True
        self._process = None
        self.finished.emit()

        # User pressed Stop — exit code is irrelevant, we just cleanly end.
        if self._stopped_by_user:
            return

        if exit_code != 0:
            tail = "\n".join(self._stdout_log[-20:])
            self.error_occurred.emit(
                f"cable_solver exited with code {exit_code}\n\n--- last log ---\n{tail}"
            )
            return

        if not self._tmp_dir:
            self.error_occurred.emit("No output directory tracked for this run")
            return

        result_path = Path(self._tmp_dir) / "result.json"
        if not result_path.is_file():
            self.error_occurred.emit(f"Result file not produced: {result_path}")
            return

        try:
            with result_path.open() as f:
                data = json.load(f)
        except Exception as e:
            self.error_occurred.emit(f"Failed to parse result.json: {e}")
            return

        self.result_ready.emit(data)

    def _on_error(self, err: QProcess.ProcessError) -> None:
        # QProcess fires both errorOccurred AND finished for terminated
        # processes. We let _on_finished be the single lifecycle exit point
        # and only escalate here when the process can't start at all — in
        # that case finished never fires.
        if err == QProcess.ProcessError.FailedToStart and not self._run_ended:
            self._run_ended = True
            self.error_occurred.emit(
                f"Failed to start cable_solver at {self._solver_path}"
            )
            self._process = None
            self.finished.emit()

    def _cleanup_tmp_dir(self) -> None:
        if self._tmp_dir and os.path.isdir(self._tmp_dir):
            shutil.rmtree(self._tmp_dir, ignore_errors=True)
        self._tmp_dir = None

    def __del__(self) -> None:
        try:
            self._cleanup_tmp_dir()
        except Exception:
            pass
