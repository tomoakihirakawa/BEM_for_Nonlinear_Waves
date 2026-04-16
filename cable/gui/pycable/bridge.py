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

from .params import CableParams, LumpedCableSystemParams
from .solver_discovery import find_cable_solver, CableSolverNotFound


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
        self._output_dir: Optional[str] = None  # user-specified persistent output dir
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
    def output_dir(self) -> Optional[str]:
        """User-specified persistent output directory. If set, results are
        saved here instead of a temporary directory."""
        return self._output_dir

    @output_dir.setter
    def output_dir(self, path: Optional[str | Path]) -> None:
        self._output_dir = str(path) if path else None

    @property
    def last_result_dir(self) -> Optional[str]:
        """Directory where the last run's results were written."""
        return self._tmp_dir

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

    def run_system(self, params: "LumpedCableSystemParams",
                   extra_json: dict | None = None,
                   source_path: "Path | None" = None) -> None:
        """Start a solver process for a multi-cable system.

        Parameters
        ----------
        params : LumpedCableSystemParams
            Cable system to solve.
        extra_json : dict, optional
            Additional keys to merge into the input JSON (e.g.
            ``initial_condition``, ``tension_top``, ``tension_bottom``).
            These override the params-generated keys.

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

        # Set up the output directory for this run.
        # If user specified a persistent output_dir, use a timestamped
        # subdirectory there. Otherwise use ~/.cache/pycable/runs/<timestamp>/.
        # Results always persist — no temp directories.
        self._cleanup_tmp_dir()
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        if self._output_dir:
            self._tmp_dir = os.path.join(self._output_dir, timestamp)
        else:
            cache_root = os.path.join(Path.home(), ".cache", "pycable", "runs")
            self._tmp_dir = os.path.join(cache_root, timestamp)
        os.makedirs(self._tmp_dir, exist_ok=True)

        # If source_path is already a solver-readable settings/per-cable file,
        # preserve that schema. Regenerating through LumpedCableSystemParams can
        # turn per-cable dynamic input into legacy multi-line input, and that
        # path intentionally only supports static equilibrium.
        source_data: dict | None = None
        if source_path is not None and Path(source_path).exists():
            try:
                with Path(source_path).open() as f:
                    sd = json.load(f)
                if "input_files" in sd:
                    src_dir = Path(source_path).parent
                    source_data = dict(sd)
                    # The copied settings.json lives in the run directory, so
                    # make relative input_files absolute before launching.
                    abs_input_files = []
                    for fname in sd["input_files"]:
                        fp = Path(fname).expanduser()
                        if not fp.is_absolute():
                            fp = src_dir / fp
                        abs_input_files.append(str(fp.resolve()))
                    source_data["input_files"] = abs_input_files
                elif "end_a_position" in sd:
                    source_data = dict(sd)
            except Exception:
                source_data = None

        if source_data is not None:
            input_path = Path(self._tmp_dir) / "input.json"
            if extra_json:
                source_data.update(extra_json)
            with input_path.open("w") as f:
                json.dump(source_data, f, indent=2)
        else:
            input_path = Path(self._tmp_dir) / "input.json"
            params.write_json(input_path)

            # Merge extra keys (initial_condition, tension, etc.) into the JSON
            if extra_json:
                with input_path.open() as f:
                    d = json.load(f)
                d.update(extra_json)
                with input_path.open("w") as f:
                    json.dump(d, f, indent=2)

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
                if "t" in snap:
                    # Dynamic mode: show time instead of |v|_max
                    self.log_received.emit(f"t={snap['t']:.3f}s  iter {iter_val}")
                else:
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
            # Settings-mode (per-cable output): aggregate *_result.json files
            # into the consolidated schema the GUI expects.
            per_cable = sorted(Path(self._tmp_dir).glob("*_result.json"))
            if not per_cable:
                self.error_occurred.emit(f"Result file not produced: {result_path}")
                return
            agg = {"n_cables": len(per_cable), "converged": True,
                   "computation_time_ms": 0.0, "cables": {}}
            for pc in per_cable:
                try:
                    with pc.open() as f:
                        d = json.load(f)
                except Exception:
                    continue
                name = d.get("name", pc.stem.replace("_result", ""))
                # Normalize: dynamic results have positions_final; static have positions
                positions = d.get("positions") or d.get("positions_final") or []
                # Dynamic per-cable result lacks per-node tensions; fall back to
                # top/bottom tensions or empty list.
                tensions = d.get("tensions", [])
                top_t = d.get("top_tension")
                bot_t = d.get("bottom_tension")
                if isinstance(top_t, list): top_t = top_t[-1] if top_t else 0
                if isinstance(bot_t, list): bot_t = bot_t[-1] if bot_t else 0
                # If per-node tensions missing, synthesize a linear ramp from
                # bottom_tension to top_tension so the color map has data.
                if (not tensions) and positions and top_t is not None and bot_t is not None:
                    n = len(positions)
                    if n >= 2:
                        tensions = [bot_t + (top_t - bot_t) * i / (n - 1)
                                    for i in range(n)]
                agg["cables"][name] = {
                    "n_nodes": len(positions),
                    "positions": positions,
                    "tensions": tensions,
                    "top_tension": top_t or 0,
                    "bottom_tension": bot_t or 0,
                    "converged": d.get("converged", True),
                }
                agg["computation_time_ms"] += d.get("computation_time_ms", 0)
            self.result_ready.emit(agg)
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
        # Results are always persistent — never auto-delete.
        self._tmp_dir = None

    def __del__(self) -> None:
        pass
