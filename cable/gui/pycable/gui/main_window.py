"""MainWindow — the single pycable GUI window."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QDockWidget,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from ..bridge import CableBridge
from ..params import CableParams, CableSpec, LumpedCableSystemParams, PerCableParams
from .bodies_list import BodiesListWidget, BodySpec
from .body_editor import BodyEditorDialog
from .lines_list import LinesListWidget
from .time_player import TimePlayerWidget
from .log_panel import LogPanel
from .run_history import RunHistoryWidget
from .setup_panel import SetupPanel
from .view_3d import View3D, cable_color_by_index


class MainWindow(QMainWindow):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("pycable — Cable Dynamics GUI")
        screen = QApplication.primaryScreen().availableGeometry()
        self.resize(min(1280, screen.width() - 50), min(800, screen.height() - 50))

        # ---- Central 3D view with time player below ----
        self._view_3d = View3D()
        self._time_player = TimePlayerWidget()
        self._time_player.time_changed.connect(self._on_time_changed)
        self._time_player.show_static_changed.connect(
            self._view_3d.set_static_cables_visible
        )
        self._time_player.contour_settings_changed.connect(self._on_contour_settings_changed)
        # Per-run animation frames: {cable_name: [(t, positions_array), ...]}
        self._animation_frames: dict[str, list] = {}

        central = QWidget()
        central_layout = QVBoxLayout(central)
        central_layout.setContentsMargins(0, 0, 0, 0)
        central_layout.setSpacing(0)
        central_layout.addWidget(self._view_3d, stretch=1)
        central_layout.addWidget(self._time_player)
        self.setCentralWidget(central)

        # ---- Left dock: lines list (top) + parameter forms (below) ----
        self._lines_list = LinesListWidget()
        self._setup_panel = SetupPanel()

        # Run All button — will be placed next to Run/Stop in SetupPanel
        self._btn_run_all = QPushButton("Run All")
        self._btn_run_all.setToolTip("Run full system — solve all cables (Ctrl+R)")
        self._btn_run_all.clicked.connect(self._on_run_system_triggered)
        # Insert Run All into SetupPanel's Run button row
        self._setup_panel.add_run_all_button(self._btn_run_all)

        self._bodies_list = BodiesListWidget()
        self._loaded_bodies: list[BodySpec] = []
        self._loaded_body_paths: list[Path] = []
        self._bodies_list.body_double_clicked.connect(self._on_body_double_clicked)

        left_container = QWidget()
        left_layout = QVBoxLayout(left_container)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(4)
        left_layout.addWidget(QLabel("Bodies in system"))
        left_layout.addWidget(self._bodies_list)
        left_layout.addWidget(QLabel("Cables in system"))
        left_layout.addWidget(self._lines_list)
        left_layout.addWidget(self._setup_panel, stretch=1)

        left_scroll = QScrollArea()
        left_scroll.setWidget(left_container)
        left_scroll.setWidgetResizable(True)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        left_dock = QDockWidget("Parameters", self)
        left_dock.setWidget(left_scroll)
        left_dock.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetMovable
            | QDockWidget.DockWidgetFeature.DockWidgetFloatable
        )
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, left_dock)

        self._lines_list.cable_selected.connect(self._on_lines_list_selection)

        # ---- Bottom dock: log + history tabs ----
        self._log_panel = LogPanel()
        self._run_history = RunHistoryWidget()
        self._run_history.load_result.connect(self._on_load_history_result)

        from PySide6.QtWidgets import QTabWidget, QSizePolicy
        bottom_tabs = QTabWidget()
        bottom_tabs.addTab(self._log_panel, "Log")
        bottom_tabs.addTab(self._run_history, "History")
        bottom_tabs.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        bottom_tabs.setMinimumHeight(200)

        bottom_dock = QDockWidget("Log / History", self)
        bottom_dock.setWidget(bottom_tabs)
        bottom_dock.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetMovable
            | QDockWidget.DockWidgetFeature.DockWidgetFloatable
        )
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, bottom_dock)

        # ---- Solver bridge ----
        self._bridge = CableBridge(self)
        self._bridge.started.connect(self._on_started)
        self._bridge.finished.connect(self._on_finished)
        self._bridge.log_received.connect(self._log_panel.append_line)
        self._bridge.snapshot_ready.connect(self._on_snapshot)
        self._bridge.result_ready.connect(self._on_result)
        self._bridge.error_occurred.connect(self._on_error)

        # ---- Wire panel signals ----
        self._setup_panel.run_requested.connect(self._on_run_requested)
        self._setup_panel.stop_requested.connect(self._bridge.stop)
        self._setup_panel.params_changed.connect(self._on_params_changed)

        # ---- Menu bar ----
        self._build_menus()

        # Restore persisted output directory
        saved_output = QSettings("pycable", "pycable").value("output_dir", "")
        if saved_output and Path(saved_output).is_dir():
            self._bridge.output_dir = saved_output
            self._setup_panel.set_output_dir(saved_output)

        # Track the last loaded file for window title and status
        self._last_loaded_path: Path | None = None

        # ---- Loaded system (Phase 4) ----
        # Single source of truth for "what cables are loaded". For legacy
        # single-line input this is a 1-cable system. For multi-line input
        # this contains all the cables in the file (or merged from many
        # files via Open Multiple).
        self._loaded_system: LumpedCableSystemParams | None = None
        # Per-cable params for files loaded via settings/per-cable format.
        # Indexed in the same order as _loaded_system.cables.
        # Contains initial_condition, tension_top/bottom etc. that CableSpec lacks.
        self._loaded_per_cable: list["PerCableParams"] = []

        # ---- System run state (Run full system) ----
        # True while a ``bridge.run_system`` invocation is active. Used to
        # decide whether to demux SNAPSHOTs by ``cable`` field into
        # View3D.append_cable calls vs the single-cable working curve.
        self._system_run_active: bool = False

        # ---- Sequential multi-file batch state (legacy load_input_files) ----
        # Used by `File → Open multiple input JSONs…` when the user selects
        # multiple separate single-line JSON files. Each file is run as a
        # separate single-line invocation and accumulated into the 3D view.
        # The new "Run system" path is independent of this queue.
        self._multi_queue: list[tuple[Path, CableParams]] = []
        self._multi_total: int = 0
        self._multi_done: int = 0

        # Initialize with a default 1-cable system so the form has something
        # sensible to display on startup.
        default_sys = LumpedCableSystemParams.from_cable_params(self._setup_panel.collect_params())
        self._set_loaded_system(default_sys, source_path=None)
        self._on_params_changed(self._setup_panel.collect_params())

    # ------------------------------------------------------------------
    # Multi-cable helpers
    # ------------------------------------------------------------------

    @property
    def _in_multi_mode(self) -> bool:
        return self._multi_total > 0 and self._multi_done < self._multi_total

    def _exit_multi_mode(self) -> None:
        self._multi_queue = []
        self._multi_total = 0
        self._multi_done = 0

    # ------------------------------------------------------------------
    # Loaded system helpers (Phase 4)
    # ------------------------------------------------------------------

    @property
    def loaded_system(self) -> LumpedCableSystemParams | None:
        return self._loaded_system

    def _set_loaded_system(self,
                           system: LumpedCableSystemParams,
                           source_path: Path | None) -> None:
        """Replace the loaded system, refresh the lines list and form."""
        self._loaded_system = system
        self._lines_list.populate_from_system(system)

        # Show the first cable's params in the form (or do nothing if empty).
        if system.cables:
            cp = system.cables[0].to_cable_params(
                gravity=system.gravity,
                mode=system.mode,
                max_equilibrium_steps=system.max_equilibrium_steps,
                equilibrium_tol=system.equilibrium_tol,
                snapshot_interval=system.snapshot_interval,
            )
            # Propagate system-level fluid/wind to the displayed CableParams
            cp.fluid = system.fluid
            cp.fluid_density = system.fluid_density
            cp.drag_Cd = system.drag_Cd
            cp.wind_type = system.wind_type
            cp.wind_U_mean = system.wind_U_mean
            cp.wind_turbulence_intensity = system.wind_turbulence_intensity
            cp.wind_integral_time_scale = system.wind_integral_time_scale
            cp.wind_seed = system.wind_seed
            cp.dt = system.dt
            cp.t_end = system.t_end
            cp.output_interval = system.output_interval
            self._setup_panel.set_params(cp)
            # Apply per-cable initial condition info if available
            if self._loaded_per_cable:
                pc = self._loaded_per_cable[0]
                self._setup_panel.set_initial_condition(
                    pc.initial_condition, pc.tension_top, pc.tension_bottom
                )
        if source_path is not None:
            self._last_loaded_path = source_path

    def _on_lines_list_selection(self, idx: int) -> None:
        """Selecting a cable in the list updates the form and 3D preview."""
        if not self._loaded_system or idx < 0 or idx >= len(self._loaded_system.cables):
            return
        spec = self._loaded_system.cables[idx]
        params = spec.to_cable_params(
            gravity=self._loaded_system.gravity,
            mode=self._loaded_system.mode,
            max_equilibrium_steps=self._loaded_system.max_equilibrium_steps,
            equilibrium_tol=self._loaded_system.equilibrium_tol,
            snapshot_interval=self._loaded_system.snapshot_interval,
        )
        # Propagate system-level fluid/wind to the displayed CableParams
        sys = self._loaded_system
        params.fluid = sys.fluid
        params.fluid_density = sys.fluid_density
        params.drag_Cd = sys.drag_Cd
        params.wind_type = sys.wind_type
        params.wind_U_mean = sys.wind_U_mean
        params.wind_turbulence_intensity = sys.wind_turbulence_intensity
        params.wind_integral_time_scale = sys.wind_integral_time_scale
        params.wind_seed = sys.wind_seed
        params.dt = sys.dt
        params.t_end = sys.t_end
        params.output_interval = sys.output_interval
        self._setup_panel.set_params(params)
        # Apply per-cable initial condition info if available
        if idx < len(self._loaded_per_cable):
            pc = self._loaded_per_cable[idx]
            self._setup_panel.set_initial_condition(
                pc.initial_condition, pc.tension_top, pc.tension_bottom
            )
        # Form's params_changed signal will fire and update the preview
        # (suppressed if there are accumulated result cables in the view).

    # ------------------------------------------------------------------
    # Menu
    # ------------------------------------------------------------------

    def _build_menus(self) -> None:
        menubar = self.menuBar()
        file_menu = menubar.addMenu("&File")

        open_action = QAction("&Open input JSON…", self)
        open_action.setShortcut(QKeySequence.StandardKey.Open)
        open_action.triggered.connect(self._on_open_triggered)
        file_menu.addAction(open_action)

        open_multi_action = QAction("Open &multiple input JSONs…", self)
        open_multi_action.setShortcut(QKeySequence("Shift+Ctrl+O"))
        open_multi_action.triggered.connect(self._on_open_multiple_triggered)
        file_menu.addAction(open_multi_action)

        clear_multi_action = QAction("&Clear all cables from view", self)
        clear_multi_action.triggered.connect(self._on_clear_cables_triggered)
        file_menu.addAction(clear_multi_action)

        file_menu.addSeparator()

        set_output_action = QAction("Set &output directory…", self)
        set_output_action.triggered.connect(self._on_set_output_dir)
        file_menu.addAction(set_output_action)

        file_menu.addSeparator()

        run_menu = menubar.addMenu("&Run")

        run_system_action = QAction("Run &full system (all cables)", self)
        run_system_action.setShortcut(QKeySequence("Ctrl+R"))
        run_system_action.triggered.connect(self._on_run_system_triggered)
        run_menu.addAction(run_system_action)

        save_action = QAction("&Save current params as JSON…", self)
        save_action.setShortcut(QKeySequence.StandardKey.Save)
        save_action.triggered.connect(self._on_save_triggered)
        file_menu.addAction(save_action)

        file_menu.addSeparator()

        # --- Recent Files submenu ---
        self._recent_menu = file_menu.addMenu("Recent files")
        self._rebuild_recent_menu()

        file_menu.addSeparator()

        quit_action = QAction("&Quit", self)
        quit_action.setShortcut(QKeySequence.StandardKey.Quit)
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

    def _on_open_triggered(self) -> None:
        """Open a CableParams JSON file and populate the form."""
        start_dir = ""
        if self._last_loaded_path is not None:
            start_dir = str(self._last_loaded_path.parent)
        path_str, _ = QFileDialog.getOpenFileName(
            self,
            "Open cable input JSON",
            start_dir,
            "Cable input JSON (*.json);;All files (*)",
        )
        if not path_str:
            return
        self.load_input_file(path_str)

    def _on_open_multiple_triggered(self) -> None:
        """Open multiple CableParams JSON files and queue them for batch run.

        Uses the instance-based ``QFileDialog`` API with
        ``FileMode.ExistingFiles`` explicitly set, plus
        ``DontUseNativeDialog`` so multi-select works reliably on macOS.
        (The native NSOpenPanel path via ``getOpenFileNames`` has been
        observed to only return one file on some macOS versions.)
        """
        start_dir = ""
        if self._last_loaded_path is not None:
            start_dir = str(self._last_loaded_path.parent)

        dialog = QFileDialog(self, "Open multiple cable input JSONs", start_dir)
        dialog.setFileMode(QFileDialog.FileMode.ExistingFiles)
        dialog.setNameFilter("Cable input JSON (*.json);;All files (*)")
        dialog.setOption(QFileDialog.Option.DontUseNativeDialog, True)
        if dialog.exec() != QFileDialog.DialogCode.Accepted:
            return
        path_strs = dialog.selectedFiles()
        if not path_strs:
            return
        # If a single file was selected, route through load_input_file
        # which handles all formats (settings.json, per-cable, legacy).
        if len(path_strs) == 1:
            self.load_input_file(path_strs[0])
            return
        self.load_input_files([Path(p) for p in path_strs])

    # ------------------------------------------------------------------
    # Recent files
    # ------------------------------------------------------------------

    _MAX_RECENT = 10

    def _recent_paths(self) -> list[str]:
        settings = QSettings("pycable", "pycable")
        return settings.value("recent_files", [], type=list)

    def _add_recent(self, path: Path) -> None:
        paths = self._recent_paths()
        s = str(path)
        if s in paths:
            paths.remove(s)
        paths.insert(0, s)
        paths = paths[: self._MAX_RECENT]
        QSettings("pycable", "pycable").setValue("recent_files", paths)
        self._rebuild_recent_menu()

    def _rebuild_recent_menu(self) -> None:
        self._recent_menu.clear()
        paths = self._recent_paths()
        if not paths:
            a = self._recent_menu.addAction("(no recent files)")
            a.setEnabled(False)
            return
        for p_str in paths:
            p = Path(p_str)
            action = self._recent_menu.addAction(p.name)
            action.setToolTip(p_str)
            action.setData(p_str)
            action.triggered.connect(self._on_recent_file_triggered)
        self._recent_menu.addSeparator()
        clear_action = self._recent_menu.addAction("Clear recent files")
        clear_action.triggered.connect(self._on_clear_recent)

    def _on_recent_file_triggered(self) -> None:
        action = self.sender()
        if action:
            self.load_input_file(action.data())

    def _on_clear_recent(self) -> None:
        QSettings("pycable", "pycable").setValue("recent_files", [])
        self._rebuild_recent_menu()

    def _on_set_output_dir(self) -> None:
        """Let the user pick a persistent output directory for results."""
        current = self._bridge.output_dir or ""
        d = QFileDialog.getExistingDirectory(
            self, "Set output directory", current
        )
        if not d:
            return
        self._bridge.output_dir = d
        self._setup_panel.set_output_dir(d)
        self._log_panel.append_line(f"[pycable] output directory set to: {d}")
        # Persist in QSettings
        QSettings("pycable", "pycable").setValue("output_dir", d)

    def _on_clear_cables_triggered(self) -> None:
        """Remove all accumulated result cables from the 3D view."""
        self._view_3d.clear_all_result_cables()
        self._view_3d.reset_camera()
        self._log_panel.append_line("[pycable] cleared all result cables from view")

    def _on_save_triggered(self) -> None:
        """Write the current form contents to a JSON file."""
        start_dir = ""
        if self._last_loaded_path is not None:
            start_dir = str(self._last_loaded_path)
        path_str, _ = QFileDialog.getSaveFileName(
            self,
            "Save cable input JSON",
            start_dir,
            "Cable input JSON (*.json);;All files (*)",
        )
        if not path_str:
            return
        params = self._setup_panel.collect_params()
        try:
            params.write_json(path_str)
        except Exception as e:
            QMessageBox.critical(self, "pycable — save failed", str(e))
            return
        self._setup_panel.set_status(f"Saved: {Path(path_str).name}")
        self._log_panel.append_line(f"[pycable] saved {path_str}")

    # ------------------------------------------------------------------
    # Public: load a JSON file (also used by CLI --load flag)
    # ------------------------------------------------------------------

    def load_input_file(self, path: str | Path) -> None:
        """Load a cable input JSON file into the form and refresh the preview.

        Auto-detects single-line (legacy ``point_a``/``point_b``/...) vs
        multi-line (``mooring_<name>`` flat arrays) schemas via
        ``LumpedCableSystemParams.read_json``. In either case the result is
        stored as the active ``_loaded_system`` and the Lines list is
        populated accordingly.
        """
        p = Path(path).expanduser().resolve()
        try:
            system = LumpedCableSystemParams.read_json(p)
        except Exception as e:
            QMessageBox.critical(
                self, "pycable — load failed", f"Failed to parse {p}:\n{e}"
            )
            return

        # Also load per-cable params to preserve initial_condition/tension info,
        # and extract output_directory if present.
        self._loaded_per_cable = []
        try:
            import json as _json
            with p.open() as _f:
                _top = _json.load(_f)

            # Extract output_directory from the loaded JSON
            raw_output_dir = _top.get("output_directory", "")
            if raw_output_dir:
                out_dir = (p.parent / raw_output_dir).resolve()
                self._bridge.output_dir = str(out_dir)
                self._setup_panel.set_output_dir(str(out_dir))
                QSettings("pycable", "pycable").setValue("output_dir", str(out_dir))
                self._log_panel.append_line(f"[pycable] output_directory: {out_dir}")

            self._loaded_bodies = []
            self._loaded_body_paths = []
            if "input_files" in _top:
                for fname in _top["input_files"]:
                    cable_path = p.parent / fname
                    if cable_path.exists():
                        with cable_path.open() as _fb:
                            _sub = _json.load(_fb)
                        if "RigidBody" in str(_sub.get("type", "")):
                            self._loaded_bodies.append(BodySpec(
                                name=str(_sub.get("name", cable_path.stem)),
                                type=str(_sub.get("type", "RigidBody")),
                                velocity=list(_sub.get("velocity", [])) if isinstance(
                                    _sub.get("velocity"), list) else [],
                            ))
                            self._loaded_body_paths.append(cable_path)
                            continue  # body file, not a cable
                        self._loaded_per_cable.append(PerCableParams.read_json(cable_path))
            elif "end_a_position" in _top:
                self._loaded_per_cable.append(PerCableParams.from_dict(_top))
            elif "point_a" in _top:
                legacy = CableParams.from_dict(_top)
                pc = PerCableParams.from_cable_params(legacy, name=str(_top.get("name", p.stem)))
                pc.initial_condition = str(_top.get("initial_condition", "length"))
                pc.tension = float(_top.get("tension", 0.0))
                pc.tension_top = float(_top.get("tension_top", 0.0))
                pc.tension_bottom = float(_top.get("tension_bottom", 0.0))
                self._loaded_per_cable.append(pc)
            self._bodies_list.populate(self._loaded_bodies)
        except Exception:
            pass  # fallback: no per-cable info

        # Leaving multi-mode (sequential batch) if we were in it
        self._exit_multi_mode()
        # Clear any accumulated result cables from a previous batch so the
        # new file's cables start on a clean canvas.
        self._view_3d.clear_all_result_cables()
        self._set_loaded_system(system, source_path=p)
        self.setWindowTitle(f"pycable — {p.name}")
        n = len(system.cables)
        self._setup_panel.set_status(
            f"Loaded: {p.name} ({n} cable{'s' if n != 1 else ''})"
        )
        self._log_panel.append_line(f"[pycable] loaded {p}  ({n} cables)")
        self._add_recent(p)

    def load_input_files(self, paths: list[Path]) -> None:
        """Queue a list of CableParams JSON files for batch execution.

        Each file is parsed, validated, and then run sequentially through
        the cable_solver. Results are appended to the 3D view as a
        colored cable under a stable name (the file stem, e.g. "C01").
        """
        if self._bridge.is_running:
            QMessageBox.warning(
                self,
                "pycable — busy",
                "Solver is currently running. Wait for it to finish first.",
            )
            return
        if not paths:
            return

        # Parse all files up front so we fail fast on bad JSON
        queued: list[tuple[Path, CableParams]] = []
        for p in paths:
            pr = Path(p).expanduser().resolve()
            try:
                params = CableParams.read_json(pr)
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "pycable — load failed",
                    f"Failed to parse {pr}:\n{e}",
                )
                return
            queued.append((pr, params))

        # Clear previous accumulated cables and set up the queue
        self._view_3d.clear_all_result_cables()
        self._view_3d._clear_working_actors()
        self._multi_queue = queued
        self._multi_total = len(queued)
        self._multi_done = 0
        self._last_loaded_path = queued[0][0]

        self._log_panel.clear()
        self._log_panel.append_line(
            f"[pycable] batch load: {self._multi_total} cables queued"
        )
        for p, _ in queued:
            self._log_panel.append_line(f"  queued: {p.name}")

        self._start_next_in_queue()

    def _start_next_in_queue(self) -> None:
        """Kick off the next cable in the multi-queue, or finish the batch."""
        if self._multi_done >= self._multi_total:
            # Batch done
            names = ", ".join(self._view_3d.result_cable_names)
            self._log_panel.append_line(
                f"[pycable] batch complete: {self._multi_total} cables — {names}"
            )
            self._setup_panel.set_status(
                f"Batch done — {self._multi_total} cables in view"
            )
            self.setWindowTitle(
                f"pycable — {self._multi_total} cables loaded"
            )
            self._exit_multi_mode()
            self._view_3d.reset_camera()
            return

        path, params = self._multi_queue[self._multi_done]
        # Drive the form to show which cable is currently running
        self._setup_panel.set_params(params)
        self.setWindowTitle(f"pycable — {path.name} ({self._multi_done + 1}/{self._multi_total})")
        self._setup_panel.set_status(
            f"Running {path.name} ({self._multi_done + 1}/{self._multi_total})…"
        )
        self._log_panel.append_line(
            f"[pycable] ({self._multi_done + 1}/{self._multi_total}) {path.name}  "
            f"length={params.cable_length:g}  n_segments={params.n_segments}"
        )
        # Start the solver for this cable — do NOT clear prior result cables
        self._bridge.run_equilibrium(params)

    # ------------------------------------------------------------------
    # Bridge callbacks
    # ------------------------------------------------------------------

    def _on_run_requested(self, params: CableParams) -> None:
        """Run only the cable currently shown in the setup form ("Run selected").

        Mirrors the form edits back into ``_loaded_system`` so that a user
        who edits the selected cable and hits Run gets the modified
        parameters saved into the loaded system as well.

        If the loaded JSON is a settings.json with body files (cantilever,
        sinusoidal, etc.), the selected cable is run inside a filtered
        settings.json that still references the body files — so body-driven
        dynamic analysis of a single cable works.
        """
        self._sync_form_into_loaded_system(params)
        self._system_run_active = False
        self._log_panel.clear()
        self._log_panel.append_line(f"[pycable] solver = {self._bridge.solver_path}")
        self._log_panel.append_line(
            f"[pycable] run selected cable  length={params.cable_length:g}  "
            f"n_segments={params.n_segments}  EA={params.EA:g}  mode={params.mode}"
        )
        self._view_3d.clear_for_new_run()
        self._view_3d.show_endpoints(params.point_a, params.point_b)

        # Path A: the loaded source has body files → generate a filtered
        # settings.json that references [bodies + selected cable file] and
        # pass it as source_path so bridge.run_system takes the
        # ``use_source_directly`` branch and preserves body driving.
        selected_name = None
        if (self._loaded_system is not None
                and self._loaded_system.cables):
            idx = max(0, self._lines_list.selected_index())
            if idx < len(self._loaded_system.cables):
                selected_name = self._loaded_system.cables[idx].name
        tmp_settings = self._maybe_make_filtered_settings(
            params, cables_filter=selected_name,
        )
        if tmp_settings is not None:
            self._log_panel.append_line(
                f"[pycable] body-preserving run via {tmp_settings.name}"
            )
            # LumpedCableSystemParams is still required as first arg but
            # unused when source_path triggers use_source_directly.
            sys_params = LumpedCableSystemParams.from_cable_params(params)
            self._bridge.run_system(sys_params, source_path=tmp_settings)
            return

        # Path B: no bodies and the loaded cable came from canonical
        # per-cable JSON. Keep that schema so dynamic mode and tension-based
        # initial conditions are preserved.
        per_cable_input = self._write_current_per_cable_input(
            params, selected_name=selected_name,
        )
        if per_cable_input is not None:
            self._log_panel.append_line(
                f"[pycable] per-cable run via {per_cable_input.name}"
            )
            sys_params = LumpedCableSystemParams.from_cable_params(
                params, name=self._current_per_cable_name(selected_name)
            )
            self._bridge.run_system(sys_params, source_path=per_cable_input)
            return

        # Path C: no loaded per-cable source — fall back to the historical
        # standalone single-cable JSON.
        # JSON as before.
        extra = self._collect_extra_json()
        sys_params = LumpedCableSystemParams.from_cable_params(params)
        self._bridge.run_system(sys_params, extra_json=extra)

    def _current_per_cable_name(self, selected_name: str | None = None) -> str:
        """Return the stable name for the currently selected cable."""
        if selected_name:
            return selected_name
        idx = max(0, self._lines_list.selected_index())
        if self._loaded_system is not None and idx < len(self._loaded_system.cables):
            return self._loaded_system.cables[idx].name
        if idx < len(self._loaded_per_cable):
            return self._loaded_per_cable[idx].name
        return "cable"

    def _per_cable_from_form(
        self,
        params: CableParams,
        selected_name: str | None = None,
    ) -> PerCableParams | None:
        """Build canonical per-cable input from the visible form.

        Existing PerCableParams are used as a base so body labels and endpoint
        motion keys survive form edits. The visible form then overrides the
        solver/geometry/material fields.
        """
        if not self._loaded_per_cable:
            return None
        idx = max(0, self._lines_list.selected_index())
        base = self._loaded_per_cable[idx] if idx < len(self._loaded_per_cable) else None
        if base is None:
            return None

        pc = PerCableParams.from_dict(base.to_dict())
        pc.name = self._current_per_cable_name(selected_name)
        pc.end_a_position = tuple(params.point_a)
        pc.end_b_position = tuple(params.point_b)
        pc.cable_length = params.cable_length
        pc.n_points = params.n_segments + 1
        pc.line_density = params.line_density
        pc.EA = params.EA
        pc.damping = params.damping
        pc.diameter = params.diameter
        pc.gravity = params.gravity
        pc.mode = params.mode
        pc.max_equilibrium_steps = params.max_equilibrium_steps
        pc.equilibrium_tol = params.equilibrium_tol
        pc.snapshot_interval = params.snapshot_interval
        pc.dt = params.dt
        pc.t_end = params.t_end
        pc.output_interval = params.output_interval
        pc.fluid = params.fluid
        pc.fluid_density = params.fluid_density
        pc.drag_Cd = params.drag_Cd
        pc.wind_type = params.wind_type
        pc.wind_U_mean = params.wind_U_mean
        pc.wind_turbulence_intensity = params.wind_turbulence_intensity
        pc.wind_integral_time_scale = params.wind_integral_time_scale
        pc.wind_seed = params.wind_seed

        ic = self._setup_panel.combo_initial_cond.currentText()
        pc.initial_condition = ic
        if ic == "tension":
            pc.tension_top = self._setup_panel.spin_tension_top.value()
            pc.tension_bottom = self._setup_panel.spin_tension_bot.value()
            pc.tension = 0.0
        else:
            pc.tension = 0.0
            pc.tension_top = 0.0
            pc.tension_bottom = 0.0
        return pc

    def _write_current_per_cable_input(
        self,
        params: CableParams,
        selected_name: str | None = None,
    ) -> Path | None:
        """Write a canonical per-cable input for the next solver run."""
        pc = self._per_cable_from_form(params, selected_name=selected_name)
        if pc is None:
            return None

        from datetime import datetime
        import re

        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", pc.name).strip("_") or "cable"
        base = Path(self._bridge.output_dir) if self._bridge.output_dir else (
            Path.home() / ".cache" / "pycable" / "runs"
        )
        tmp_dir = base / f"{ts}_input_{safe_name}"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        input_path = tmp_dir / "input.json"
        pc.write_json(input_path)
        return input_path

    def _maybe_make_filtered_settings(
        self,
        params: CableParams,
        cables_filter: str | None = None,
    ) -> Path | None:
        """Write a settings.json derived from the currently-loaded one, so
        that (a) mode/dt/t_end from the form override the on-disk values
        and (b) an optional cables_filter keeps only one cable.

        Returns the written path, or None if the current load state does
        not support this (no bodies, no source, etc.).

        ``cables_filter``:
        - None → include all non-body cables
        - str  → include only the cable whose name matches
        """
        if not self._loaded_bodies or not self._loaded_body_paths:
            return None
        if self._last_loaded_path is None or not self._last_loaded_path.exists():
            return None
        try:
            import json as _json
            with self._last_loaded_path.open() as f:
                orig = _json.load(f)
        except Exception:
            return None
        if "input_files" not in orig:
            return None

        src_dir = self._last_loaded_path.parent
        body_file_names: list[str] = []
        cable_file_names: list[str] = []
        for fname in orig["input_files"]:
            fpath = src_dir / fname
            if not fpath.exists():
                continue
            try:
                with fpath.open() as ff:
                    sub = _json.load(ff)
            except Exception:
                continue
            if "RigidBody" in str(sub.get("type", "")):
                body_file_names.append(fname)
                continue
            sub_name = str(sub.get("name", Path(fname).stem))
            if cables_filter is None or sub_name == cables_filter:
                cable_file_names.append(fname)

        if not cable_file_names or not body_file_names:
            return None

        # Destination tag: "selected_C01" or "all" for traceability.
        tag = f"selected_{cables_filter}" if cables_filter else "all"

        from datetime import datetime
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        base = Path(self._bridge.output_dir) if self._bridge.output_dir else (
            Path.home() / ".cache" / "pycable" / "runs"
        )
        tmp_dir = base / f"{ts}_{tag}"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        tmp_settings = tmp_dir / "settings.json"

        # Reference by absolute path so the filtered file can live anywhere.
        abs_input_files = [
            str((src_dir / n).resolve())
            for n in (body_file_names + cable_file_names)
        ]

        filtered: dict = {
            "input_files": abs_input_files,
            "gravity": float(params.gravity),
            "mode": str(params.mode),
            "max_equilibrium_steps": int(params.max_equilibrium_steps),
            "equilibrium_tol": float(params.equilibrium_tol),
            "snapshot_interval": int(params.snapshot_interval),
        }
        if params.mode == "dynamic":
            filtered["dt"] = float(params.dt)
            filtered["t_end"] = float(params.t_end)
            if params.output_interval > 0:
                filtered["output_interval"] = float(params.output_interval)
        # Propagate fluid/wind settings — otherwise the derived file falls
        # back to solver default "water" which breaks bridge runs that need
        # fluid=air.
        from ..params import _fluid_wind_to_dict
        filtered.update(_fluid_wind_to_dict(params))
        extra = self._collect_extra_json()
        if extra:
            filtered.update(extra)

        import json as _json
        with tmp_settings.open("w") as f:
            _json.dump(filtered, f, indent=2)
        return tmp_settings

    def _on_run_system_triggered(self) -> None:
        """Run the whole loaded ``LumpedCableSystemParams`` in one invocation."""
        if self._bridge.is_running:
            QMessageBox.warning(
                self, "pycable — busy",
                "Solver is currently running. Wait for it to finish first.",
            )
            return
        if self._loaded_system is None or not self._loaded_system.cables:
            QMessageBox.information(
                self, "pycable", "No cables loaded. Open a JSON file first.",
            )
            return
        # Mirror any pending form edits back into the selected cable so
        # Run-system reflects the user's latest tweaks.
        self._sync_form_into_loaded_system(self._setup_panel.collect_params())

        self._system_run_active = True
        self._view_3d.clear_all_result_cables()
        self._view_3d.clear_for_new_run()
        self._log_panel.clear()
        self._log_panel.append_line(
            f"[pycable] solver = {self._bridge.solver_path}"
        )
        n = len(self._loaded_system.cables)
        self._log_panel.append_line(
            f"[pycable] run full system  n_cables={n}"
        )
        for ci, cable in enumerate(self._loaded_system.cables):
            self._log_panel.append_line(
                f"  [{ci + 1}/{n}] {cable.name}  length={cable.cable_length:g}  "
                f"n_points={cable.n_points}"
            )
        form_params = self._setup_panel.collect_params()
        # If bodies are loaded, generate a derived settings.json so form
        # overrides (mode/dt/t_end) take effect even when running full
        # system.  _maybe_make_filtered_settings(cables_filter=None)
        # includes every cable.
        tmp_settings = self._maybe_make_filtered_settings(
            form_params, cables_filter=None,
        )
        if tmp_settings is not None:
            self._log_panel.append_line(
                f"[pycable] body-preserving run via {tmp_settings.name}"
            )
            self._bridge.run_system(self._loaded_system,
                                    source_path=tmp_settings)
            return

        if len(self._loaded_per_cable) == 1 and not self._loaded_bodies:
            per_cable_input = self._write_current_per_cable_input(form_params)
            if per_cable_input is not None:
                self._log_panel.append_line(
                    f"[pycable] per-cable run via {per_cable_input.name}"
                )
                sys_params = LumpedCableSystemParams.from_cable_params(
                    form_params, name=self._current_per_cable_name()
                )
                self._bridge.run_system(sys_params, source_path=per_cable_input)
                return

        extra = self._collect_extra_json()
        self._bridge.run_system(self._loaded_system, extra_json=extra,
                                source_path=self._last_loaded_path)

    def _collect_extra_json(self) -> dict | None:
        """Collect initial_condition / tension from the GUI form into a dict
        that will be merged into the solver input JSON."""
        ic = self._setup_panel.combo_initial_cond.currentText()
        if ic != "tension":
            return None
        t_top = self._setup_panel.spin_tension_top.value()
        t_bot = self._setup_panel.spin_tension_bot.value()
        if t_top <= 0 and t_bot <= 0:
            return None
        return {
            "initial_condition": "tension",
            "tension_top": t_top,
            "tension_bottom": t_bot,
        }

    def _sync_form_into_loaded_system(self, params: CableParams) -> None:
        """Write the current form's CableParams back into the selected cable.

        Called before Run so that a form edit is picked up by both
        "Run selected" and "Run full system".
        """
        if self._loaded_system is None or not self._loaded_system.cables:
            return
        idx = max(0, self._lines_list.selected_index())
        if idx >= len(self._loaded_system.cables):
            return
        old = self._loaded_system.cables[idx]
        self._loaded_system.cables[idx] = CableSpec.from_cable_params(
            params, name=old.name
        )
        # Shared top-level scalars too, since the form edits them too.
        self._loaded_system.gravity = params.gravity
        self._loaded_system.mode = params.mode
        self._loaded_system.max_equilibrium_steps = params.max_equilibrium_steps
        self._loaded_system.equilibrium_tol = params.equilibrium_tol
        self._loaded_system.snapshot_interval = params.snapshot_interval
        self._loaded_system.dt = params.dt
        self._loaded_system.t_end = params.t_end
        self._loaded_system.output_interval = params.output_interval
        # Fluid/wind are system-level — sync from form into loaded system.
        self._loaded_system.fluid = params.fluid
        self._loaded_system.fluid_density = params.fluid_density
        self._loaded_system.drag_Cd = params.drag_Cd
        self._loaded_system.wind_type = params.wind_type
        self._loaded_system.wind_U_mean = params.wind_U_mean
        self._loaded_system.wind_turbulence_intensity = params.wind_turbulence_intensity
        self._loaded_system.wind_integral_time_scale = params.wind_integral_time_scale
        self._loaded_system.wind_seed = params.wind_seed

    def _on_started(self) -> None:
        self._setup_panel.set_running(True)
        self._btn_run_all.setEnabled(False)
        self._setup_panel.set_status("Running…")
        # Reset playback frames and cached strain tables for the new run.
        self._animation_frames = {}
        self._anim_ref_lengths = {}
        self._anim_strain_clim = None
        self._time_player.clear()
        self._view_3d.clear_fast_cables()
        self._view_3d.clear_static_cables()

    def _on_finished(self) -> None:
        self._setup_panel.set_running(False)
        self._btn_run_all.setEnabled(True)
        # Build the union of all times across cables for the time slider.
        all_times = set()
        for frames in self._animation_frames.values():
            for t, _pos in frames:
                all_times.add(t)
        if all_times:
            self._time_player.set_times(sorted(all_times))

    def _on_body_double_clicked(self, idx: int) -> None:
        """Open BodyEditorDialog for the selected body."""
        if idx < 0 or idx >= len(self._loaded_bodies):
            return
        body = self._loaded_bodies[idx]
        path = self._loaded_body_paths[idx] if idx < len(self._loaded_body_paths) else None
        if path is None:
            return
        dlg = BodyEditorDialog(body, path, parent=self)
        if dlg.exec() == BodyEditorDialog.DialogCode.Accepted:
            # Refresh list label (motion may have changed)
            self._bodies_list.populate(self._loaded_bodies)
            self._log_panel.append_line(
                f"[pycable] body '{body.name}' updated → {path.name}"
            )

    def _on_time_changed(self, t: float) -> None:
        """Redraw every loaded cable at the frame closest to time ``t`` with
        per-node tension-proxy contour (segment strain relative to the
        equilibrium catenary in the first frame)."""
        from .view_3d import cable_color_by_index
        if not self._animation_frames:
            return
        # When the scrubber drives the view, remove the static final-tension
        # actors so the animation is the single authoritative display layer.
        self._view_3d.clear_all_result_cables()
        loaded = self._loaded_system
        names = [c.name for c in loaded.cables] if loaded else list(self._animation_frames)

        # Cache the reference (first-frame) segment lengths and the global
        # strain range once per run so the colormap is stable across frames.
        if not hasattr(self, "_anim_ref_lengths") or not self._anim_ref_lengths:
            self._anim_ref_lengths = {}
            for name, frames in self._animation_frames.items():
                if not frames:
                    continue
                p0 = frames[0][1]
                seg = np.linalg.norm(np.diff(p0, axis=0), axis=1)
                self._anim_ref_lengths[name] = np.maximum(seg, 1e-12)
            # Global strain range for a shared colormap
            smax = 0.0
            for name, frames in self._animation_frames.items():
                ref = self._anim_ref_lengths.get(name)
                if ref is None:
                    continue
                for _tt, pp in frames:
                    seg = np.linalg.norm(np.diff(pp, axis=0), axis=1)
                    strain = (seg - ref) / ref
                    smax = max(smax, float(np.max(np.abs(strain))))
            self._anim_strain_clim = (-smax, smax) if smax > 0 else (-1e-6, 1e-6)

        for idx, name in enumerate(names):
            frames = self._animation_frames.get(name)
            if not frames:
                continue
            best = min(frames, key=lambda tp: abs(tp[0] - t))
            positions = best[1]
            ref = self._anim_ref_lengths.get(name)
            scalars = None
            if ref is not None and len(ref) == positions.shape[0] - 1:
                seg = np.linalg.norm(np.diff(positions, axis=0), axis=1)
                strain_seg = (seg - ref) / ref  # per segment
                # Convert to per-node by averaging neighboring segments
                strain_node = np.empty(positions.shape[0], dtype=float)
                strain_node[0] = strain_seg[0]
                strain_node[-1] = strain_seg[-1]
                strain_node[1:-1] = 0.5 * (strain_seg[:-1] + strain_seg[1:])
                scalars = strain_node
            cs = self._time_player.contour_settings
            clim = self._anim_strain_clim
            if not cs["auto_range"]:
                clim = (cs["vmin"], cs["vmax"])
            self._view_3d.update_cable_fast(
                name, positions,
                color=cable_color_by_index(idx),
                scalars=scalars,
                clim=clim,
                cmap=cs["cmap"],
            )
        self._view_3d.render_now()

    def _on_contour_settings_changed(self, _settings: dict) -> None:
        """Re-render the current animation frame with the new contour settings.

        Forces recreation of PyVista actors so the new colormap / clim takes
        effect immediately (otherwise the cached actors keep the old cmap).
        """
        self._view_3d.clear_fast_cables()
        if self._animation_frames and self._time_player._times:
            idx = self._time_player._slider.value()
            if 0 <= idx < len(self._time_player._times):
                self._on_time_changed(self._time_player._times[idx])

    def _on_snapshot(self, snap: dict) -> None:
        if self._in_multi_mode:
            # In sequential multi-file mode, live curve updates would destroy
            # the already-accumulated result cables. Skip per-snapshot preview
            # and let the user see the final shape on result_ready instead.
            iter_val = snap.get("iter", "?")
            self._setup_panel.set_status(
                f"[{self._multi_done + 1}/{self._multi_total}] iter {iter_val}  "
                f"|v|={snap.get('norm_v', 0.0):.4g}"
            )
            return

        positions = snap.get("positions")
        if positions is None:
            return
        positions_arr = np.asarray(positions, dtype=float)
        iter_val = snap.get("iter", "?")

        # Collect time-series frame for post-run playback (dynamic mode only,
        # when snapshot carries a "t" field).
        t_snap = snap.get("t")
        cable_name = snap.get("cable")
        if t_snap is not None and cable_name is not None:
            # First frame for this cable → capture as static equilibrium reference.
            if cable_name not in self._animation_frames:
                self._view_3d.set_static_cable(cable_name, positions_arr)
            self._animation_frames.setdefault(cable_name, []).append(
                (float(t_snap), positions_arr.copy())
            )

        # Multi-line SNAPSHOT demux: when cable_solver is running in
        # multi-line mode each snapshot is tagged with the cable name.
        # We append each into the accumulated result-cable actors so the
        # user sees all N cables moving together during equilibration.
        if cable_name is not None:
            loaded = self._loaded_system
            idx = 0
            if loaded is not None:
                for i, c in enumerate(loaded.cables):
                    if c.name == cable_name:
                        idx = i
                        break
            self._view_3d.append_cable(
                cable_name, positions_arr, color=cable_color_by_index(idx)
            )
            self._setup_panel.set_status(
                f"{cable_name}  iter {iter_val}  |v|={snap.get('norm_v', 0.0):.4g}"
            )
            return

        # Legacy single-line SNAPSHOT path.
        self._view_3d.update_curve(positions_arr)
        self._setup_panel.set_status(f"iter {iter_val}  |v|={snap.get('norm_v', 0.0):.4g}")

    def _on_result(self, data: dict) -> None:
        # Multi-line result.json schema: {"n_cables", "converged",
        # "computation_time_ms", "cables": {name: {positions, tensions,
        # top_tension, bottom_tension, ...}}}. Detect by presence of
        # "cables" dict.
        if isinstance(data.get("cables"), dict):
            self._on_result_multi(data)
            return

        positions = np.asarray(data.get("positions", []), dtype=float)
        tensions = np.asarray(data.get("tensions", []), dtype=float)
        converged = data.get("converged", False)
        elapsed = data.get("computation_time_ms", 0.0)
        top_t = data.get("top_tension", float("nan"))
        bot_t = data.get("bottom_tension", float("nan"))

        if self._in_multi_mode:
            # Append this cable to the accumulated scene and advance the queue
            path, _ = self._multi_queue[self._multi_done]
            name = path.stem  # "C01"
            color = cable_color_by_index(self._multi_done)
            if positions.size:
                self._view_3d.append_cable(name, positions, color=color)
            self._log_panel.append_line(
                f"[pycable]   → {name} converged={converged}  elapsed={elapsed:.1f} ms  "
                f"top_T={top_t/1e3:.2f} kN  bot_T={bot_t/1e3:.2f} kN  color={color}"
            )
            self._multi_done += 1
            # Kick off the next one
            self._start_next_in_queue()
            return

        # Single-cable mode — replace the working cable with its tension-colored result
        if positions.size and tensions.size:
            self._view_3d.show_tensions(positions, tensions)
            self._view_3d.reset_camera()

        self._setup_panel.set_status(
            f"{'Converged' if converged else 'Did NOT converge'} in {elapsed:.0f} ms  |  "
            f"top T={top_t:.2f} N  |  bottom T={bot_t:.2f} N"
        )
        self._log_panel.append_line(
            f"[pycable] converged={converged}  elapsed={elapsed:.1f} ms  "
            f"top_tension={top_t:.2f}  bottom_tension={bot_t:.2f}"
        )

        # Record in run history
        if self._bridge._tmp_dir:
            result_path = str(Path(self._bridge._tmp_dir) / "result.json")
            cable_name = "cable"
            if self._loaded_system and self._loaded_system.cables:
                idx = max(0, self._lines_list.selected_index())
                if idx < len(self._loaded_system.cables):
                    cable_name = self._loaded_system.cables[idx].name
            self._run_history.add_entry(
                result_path, name=cable_name,
                top_tension=top_t, bottom_tension=bot_t,
                converged=converged,
            )

    def _on_result_multi(self, data: dict) -> None:
        """Handle a multi-line ``result.json`` (``cables`` dict schema).

        All cables are rendered with a shared tension colormap so that
        the color scale is consistent across the entire system.
        """
        self._system_run_active = False
        cables_dict = data.get("cables", {})
        converged = data.get("converged", False)
        elapsed = data.get("computation_time_ms", 0.0)
        n_cables = data.get("n_cables", len(cables_dict))

        loaded = self._loaded_system
        order: list[str] = []
        if loaded is not None:
            order = [c.name for c in loaded.cables if c.name in cables_dict]
            for name in cables_dict.keys():
                if name not in order:
                    order.append(name)
        else:
            order = list(cables_dict.keys())

        # Build cable data list for shared-colormap rendering
        cable_data_list: list[dict] = []
        for name in order:
            cable = cables_dict[name]
            positions = np.asarray(cable.get("positions", []), dtype=float)
            tensions = np.asarray(cable.get("tensions", []), dtype=float)
            if positions.size == 0:
                continue
            cable_data_list.append({
                "name": name,
                "positions": positions,
                "tensions": tensions,
            })
            top_t = cable.get("top_tension", float("nan"))
            bot_t = cable.get("bottom_tension", float("nan"))
            self._log_panel.append_line(
                f"[pycable]   {name} top_T={top_t:.2f} N  bot_T={bot_t:.2f} N"
            )

        # Render all cables with a shared tension color scale
        self._view_3d.show_multi_cable_tensions(cable_data_list)
        self._view_3d.reset_camera()

        self._setup_panel.set_status(
            f"{'Converged' if converged else 'Did NOT converge'} in {elapsed:.0f} ms  |  "
            f"{n_cables} cables"
        )
        self._log_panel.append_line(
            f"[pycable] system run complete  converged={converged}  "
            f"elapsed={elapsed:.1f} ms  n_cables={n_cables}"
        )

        # Record the whole system run as a single history entry.
        if self._bridge._tmp_dir:
            result_path = str(Path(self._bridge._tmp_dir) / "result.json")
            case_name = "system"
            if self._last_loaded_path:
                case_name = self._last_loaded_path.stem
            self._run_history.add_entry(
                result_path=result_path,
                name=f"{case_name} ({n_cables} cables)",
                top_tension=0,
                bottom_tension=0,
                converged=converged,
                n_cables=n_cables,
                elapsed_ms=elapsed,
            )

    def _on_load_history_result(self, result_path: str) -> None:
        """Reload a result JSON from history into the 3D view."""
        # Strip "#cable_name" fragment if present (used for multi-cable uniqueness)
        actual_path = result_path.split("#")[0]
        import json as _json
        try:
            with open(actual_path) as f:
                data = _json.load(f)
        except Exception as e:
            self._log_panel.append_line(f"[pycable] failed to load history: {e}")
            return

        # Multi-cable result — restore full system from input.json if available
        if isinstance(data.get("cables"), dict):
            # Try to load input.json from the same directory to restore
            # full cable parameters (cable_length, EA, etc.)
            input_json_path = Path(actual_path).parent / "input.json"
            if input_json_path.is_file():
                try:
                    system = LumpedCableSystemParams.read_json(input_json_path)
                    self._set_loaded_system(system, source_path=input_json_path)
                    self._log_panel.append_line(
                        f"[pycable] restored system from {input_json_path}")
                except Exception:
                    pass
            else:
                # Fallback: rebuild minimal system from result cable names
                cables_dict = data["cables"]
                specs = []
                for cname, cdata in cables_dict.items():
                    positions = cdata.get("positions", [])
                    if len(positions) >= 2:
                        specs.append(CableSpec(
                            name=cname,
                            point_a=tuple(positions[0]),
                            point_b=tuple(positions[-1]),
                        ))
                if specs:
                    system = LumpedCableSystemParams(cables=specs)
                    self._set_loaded_system(system, source_path=None)
            self._on_result_multi(data)
            self._log_panel.append_line(f"[pycable] reloaded from history: {actual_path}")
            return

        # Single-cable result
        positions = np.asarray(data.get("positions", []), dtype=float)
        tensions = np.asarray(data.get("tensions", []), dtype=float)
        if positions.size and tensions.size:
            self._view_3d.show_tensions(positions, tensions)
            self._view_3d.reset_camera()
        top_t = data.get("top_tension", 0)
        bot_t = data.get("bottom_tension", 0)
        self._setup_panel.set_status(
            f"History: top T={top_t:.2f} N  |  bottom T={bot_t:.2f} N"
        )
        self._log_panel.append_line(f"[pycable] reloaded from history: {result_path}")

    def _on_error(self, msg: str) -> None:
        self._log_panel.append_line(f"[pycable] ERROR: {msg}")
        self._setup_panel.set_status("Error — see log")
        QMessageBox.critical(self, "pycable — solver error", msg)

    def _on_params_changed(self, params: CableParams) -> None:
        if self._bridge.is_running:
            return
        # In multi-cable mode we don't draw a single-cable preview — the
        # accumulated result cables already own the scene, and a stray
        # preview would visually clutter the view. Same goes for a loaded
        # multi-line system (even before its first run).
        if self._in_multi_mode or self._view_3d.result_cable_names:
            return
        if self._loaded_system is not None and len(self._loaded_system.cables) > 1:
            return
        self._view_3d.show_endpoints(params.point_a, params.point_b)
        self._view_3d.show_initial_line(params.point_a, params.point_b, params.n_segments)
