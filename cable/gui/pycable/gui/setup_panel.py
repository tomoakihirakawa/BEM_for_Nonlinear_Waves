"""SetupPanel — parameter entry forms for CableParams.

Three collapsible groups: Geometry, Advanced, Solver. Fields match
pycable.params.CableParams, which in turn matches the cable_solver
input JSON schema 1:1.
"""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from ..params import CableParams


# Preset defaults — mirrors cable_solver.cpp::resolveFluidConfig.
_FLUID_PRESETS = {
    "water": {"density": 1000.0, "Cd": 2.5},
    "air":   {"density": 1.225,  "Cd": 1.2},
}


def _make_double_spin(
    lo: float, hi: float, default: float, decimals: int, suffix: str = ""
) -> QDoubleSpinBox:
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(decimals)
    s.setValue(default)
    if suffix:
        s.setSuffix(suffix)
    s.setAlignment(s.alignment())
    return s


class SetupPanel(QWidget):
    run_requested = Signal(object)      # CableParams
    stop_requested = Signal()
    params_changed = Signal(object)     # CableParams (live preview)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(12)

        defaults = CableParams()

        # ----------------- Geometry -----------------
        geo = QGroupBox("Geometry")
        geo_form = QFormLayout(geo)

        self.spin_ax = _make_double_spin(-1e6, 1e6, defaults.point_a[0], 3)
        self.spin_ay = _make_double_spin(-1e6, 1e6, defaults.point_a[1], 3)
        self.spin_az = _make_double_spin(-1e6, 1e6, defaults.point_a[2], 3)
        geo_form.addRow("Point A  X [m]", self.spin_ax)
        geo_form.addRow("Point A  Y [m]", self.spin_ay)
        geo_form.addRow("Point A  Z [m]", self.spin_az)
        geo_form.addRow(QLabel(""))

        self.spin_bx = _make_double_spin(-1e6, 1e6, defaults.point_b[0], 3)
        self.spin_by = _make_double_spin(-1e6, 1e6, defaults.point_b[1], 3)
        self.spin_bz = _make_double_spin(-1e6, 1e6, defaults.point_b[2], 3)
        geo_form.addRow("Point B  X [m]", self.spin_bx)
        geo_form.addRow("Point B  Y [m]", self.spin_by)
        geo_form.addRow("Point B  Z [m]", self.spin_bz)
        geo_form.addRow(QLabel(""))

        self.spin_length = _make_double_spin(0.1, 1e6, defaults.cable_length, 3)
        geo_form.addRow("Cable length [m]", self.spin_length)

        self.spin_segments = QSpinBox()
        self.spin_segments.setRange(2, 10_000)
        self.spin_segments.setValue(defaults.n_segments)
        geo_form.addRow("Segments", self.spin_segments)

        geo_form.addRow(QLabel(""))

        self.combo_initial_cond = QComboBox()
        self.combo_initial_cond.addItems(["length", "tension"])
        self.combo_initial_cond.setCurrentText("length")
        geo_form.addRow("Initial condition", self.combo_initial_cond)

        self.spin_tension_top = _make_double_spin(0.0, 1e12, 0.0, 1)
        self.spin_tension_top.setSuffix(" N")
        self.spin_tension_top.setSingleStep(1000)
        self.spin_tension_top.setEnabled(False)
        self._tension_top_label = QLabel("Tension top [N]")
        self._tension_top_label.setEnabled(False)
        geo_form.addRow(self._tension_top_label, self.spin_tension_top)

        self.spin_tension_bot = _make_double_spin(0.0, 1e12, 0.0, 1)
        self.spin_tension_bot.setSuffix(" N")
        self.spin_tension_bot.setSingleStep(1000)
        self.spin_tension_bot.setEnabled(False)
        self._tension_bot_label = QLabel("Tension bottom [N]")
        self._tension_bot_label.setEnabled(False)
        geo_form.addRow(self._tension_bot_label, self.spin_tension_bot)

        self.combo_initial_cond.currentTextChanged.connect(self._on_initial_cond_changed)

        layout.addWidget(geo)

        # ----------------- Advanced (material) -----------------
        adv = QGroupBox("Material (advanced)")
        adv.setCheckable(True)
        adv.setChecked(False)
        adv_form = QFormLayout(adv)

        self.spin_density = _make_double_spin(0.01, 1e5, defaults.line_density, 4)
        adv_form.addRow("Line density [kg/m]", self.spin_density)

        self.spin_ea = QDoubleSpinBox()
        self.spin_ea.setRange(1.0, 1e15)
        self.spin_ea.setDecimals(0)
        self.spin_ea.setSingleStep(1e8)
        self.spin_ea.setValue(defaults.EA)
        adv_form.addRow("EA [N]", self.spin_ea)

        self.spin_damping = _make_double_spin(0.0, 1e4, defaults.damping, 3)
        adv_form.addRow("Damping", self.spin_damping)

        self.spin_diameter = _make_double_spin(0.0001, 10.0, defaults.diameter, 5)
        adv_form.addRow("Diameter [m]", self.spin_diameter)

        self.spin_gravity = _make_double_spin(0.0, 100.0, defaults.gravity, 3)
        adv_form.addRow("Gravity [m/s²]", self.spin_gravity)

        layout.addWidget(adv)

        # ----------------- Fluid & Wind (advanced) -----------------
        fw = QGroupBox("Fluid & Wind (advanced)")
        fw.setCheckable(True)
        fw.setChecked(False)
        fw_form = QFormLayout(fw)

        self.cbo_fluid = QComboBox()
        self.cbo_fluid.addItems(["water", "air"])
        self.cbo_fluid.setCurrentText(defaults.fluid)
        fw_form.addRow("Fluid preset", self.cbo_fluid)

        fluid_preset = _FLUID_PRESETS[defaults.fluid]
        self.spin_fluid_density = _make_double_spin(
            0.01, 1e5,
            defaults.fluid_density if defaults.fluid_density is not None else fluid_preset["density"],
            4,
        )
        fw_form.addRow("Fluid density [kg/m³]", self.spin_fluid_density)

        self.spin_drag_cd = _make_double_spin(
            0.0, 100.0,
            defaults.drag_Cd if defaults.drag_Cd is not None else fluid_preset["Cd"],
            3,
        )
        fw_form.addRow("Drag Cd", self.spin_drag_cd)

        self.cbo_wind = QComboBox()
        self.cbo_wind.addItems(["none", "uniform", "AR1"])
        self.cbo_wind.setCurrentText(defaults.wind_type)
        fw_form.addRow("Wind type", self.cbo_wind)

        self.spin_wind_ux = _make_double_spin(-100.0, 100.0, defaults.wind_U_mean[0], 3)
        self.spin_wind_uy = _make_double_spin(-100.0, 100.0, defaults.wind_U_mean[1], 3)
        self.spin_wind_uz = _make_double_spin(-100.0, 100.0, defaults.wind_U_mean[2], 3)
        fw_form.addRow("Wind U_mean X [m/s]", self.spin_wind_ux)
        fw_form.addRow("Wind U_mean Y [m/s]", self.spin_wind_uy)
        fw_form.addRow("Wind U_mean Z [m/s]", self.spin_wind_uz)

        self.spin_wind_ti = _make_double_spin(0.0, 2.0, defaults.wind_turbulence_intensity, 3)
        fw_form.addRow("Turbulence intensity (AR1)", self.spin_wind_ti)

        self.spin_wind_tl = _make_double_spin(0.01, 1000.0, defaults.wind_integral_time_scale, 3)
        fw_form.addRow("Integral time scale [s] (AR1)", self.spin_wind_tl)

        self.spin_wind_seed = QSpinBox()
        self.spin_wind_seed.setRange(0, 2_000_000_000)
        self.spin_wind_seed.setValue(defaults.wind_seed if defaults.wind_seed is not None else 0)
        self.spin_wind_seed.setSpecialValueText("auto")
        fw_form.addRow("Wind RNG seed (0=auto)", self.spin_wind_seed)

        self.cbo_fluid.currentTextChanged.connect(self._on_fluid_preset_changed)
        self.cbo_wind.currentTextChanged.connect(self._update_wind_row_enabled)
        self._update_wind_row_enabled(defaults.wind_type)

        layout.addWidget(fw)

        # ----------------- Mode & dynamic parameters -----------------
        mode_group = QGroupBox("Mode")
        mode_form = QFormLayout(mode_group)

        self.combo_mode = QComboBox()
        self.combo_mode.addItems(["equilibrium", "dynamic"])
        self.combo_mode.setCurrentText(defaults.mode)
        mode_form.addRow("Mode", self.combo_mode)

        self.spin_dt = _make_double_spin(1e-6, 10.0, 0.005, 6, " s")
        mode_form.addRow("dt", self.spin_dt)

        self.spin_t_end = _make_double_spin(1e-3, 1e6, 5.0, 4, " s")
        mode_form.addRow("t_end", self.spin_t_end)

        self.spin_out_int = _make_double_spin(0.0, 1e6, 0.05, 6, " s")
        self.spin_out_int.setSpecialValueText("auto")
        mode_form.addRow("output_interval", self.spin_out_int)

        self._dynamic_rows = (self.spin_dt, self.spin_t_end, self.spin_out_int)
        self.combo_mode.currentTextChanged.connect(self._on_mode_changed)
        self._on_mode_changed(defaults.mode)

        layout.addWidget(mode_group)

        # ----------------- Solver -----------------
        sv = QGroupBox("Solver (advanced)")
        sv.setCheckable(True)
        sv.setChecked(False)
        sv_form = QFormLayout(sv)

        self.spin_max_steps = QSpinBox()
        self.spin_max_steps.setRange(100, 10_000_000)
        self.spin_max_steps.setValue(defaults.max_equilibrium_steps)
        self.spin_max_steps.setSingleStep(1000)
        sv_form.addRow("Max steps", self.spin_max_steps)

        self.spin_tol = _make_double_spin(1e-8, 10.0, defaults.equilibrium_tol, 6)
        sv_form.addRow("Equilibrium tol.", self.spin_tol)

        self.spin_snap = QSpinBox()
        self.spin_snap.setRange(1, 10_000_000)
        self.spin_snap.setValue(defaults.snapshot_interval)
        self.spin_snap.setSingleStep(1000)
        sv_form.addRow("Snapshot interval", self.spin_snap)

        layout.addWidget(sv)

        # ----------------- Run / Stop -----------------
        run_group = QGroupBox("Run")
        run_box = QVBoxLayout(run_group)
        self._btn_row = QHBoxLayout()
        self.btn_run = QPushButton("Run")
        self.btn_stop = QPushButton("Stop")
        self.btn_stop.setEnabled(False)
        self._btn_row.addWidget(self.btn_run)
        self._btn_row.addWidget(self.btn_stop)
        run_box.addLayout(self._btn_row)

        self._status_label = QLabel("Ready")
        self._status_label.setWordWrap(True)
        run_box.addWidget(self._status_label)

        self._output_dir_label = QLabel("Output: ~/.cache/pycable/runs/")
        self._output_dir_label.setWordWrap(True)
        self._output_dir_label.setStyleSheet("color: gray; font-size: 11px;")
        run_box.addWidget(self._output_dir_label)

        layout.addWidget(run_group)
        layout.addStretch(1)

        # ----------------- Signals -----------------
        self.btn_run.clicked.connect(self._on_run_clicked)
        self.btn_stop.clicked.connect(self.stop_requested.emit)

        # Live preview: only geometry fields drive the 3D view preview before a run.
        for spin in (
            self.spin_ax, self.spin_ay, self.spin_az,
            self.spin_bx, self.spin_by, self.spin_bz,
            self.spin_length,
        ):
            spin.valueChanged.connect(self._emit_params_changed)
        self.spin_segments.valueChanged.connect(self._emit_params_changed)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def collect_params(self) -> CableParams:
        fluid = self.cbo_fluid.currentText()
        preset = _FLUID_PRESETS[fluid]
        fd_val = self.spin_fluid_density.value()
        cd_val = self.spin_drag_cd.value()
        fd_override = fd_val if abs(fd_val - preset["density"]) > 1e-9 else None
        cd_override = cd_val if abs(cd_val - preset["Cd"]) > 1e-9 else None
        seed_val = int(self.spin_wind_seed.value())
        return CableParams(
            point_a=(self.spin_ax.value(), self.spin_ay.value(), self.spin_az.value()),
            point_b=(self.spin_bx.value(), self.spin_by.value(), self.spin_bz.value()),
            cable_length=self.spin_length.value(),
            n_segments=int(self.spin_segments.value()),
            line_density=self.spin_density.value(),
            EA=self.spin_ea.value(),
            damping=self.spin_damping.value(),
            diameter=self.spin_diameter.value(),
            gravity=self.spin_gravity.value(),
            mode=self.combo_mode.currentText(),
            dt=self.spin_dt.value(),
            t_end=self.spin_t_end.value(),
            output_interval=self.spin_out_int.value(),
            max_equilibrium_steps=int(self.spin_max_steps.value()),
            equilibrium_tol=self.spin_tol.value(),
            snapshot_interval=int(self.spin_snap.value()),
            fluid=fluid,
            fluid_density=fd_override,
            drag_Cd=cd_override,
            wind_type=self.cbo_wind.currentText(),
            wind_U_mean=(
                self.spin_wind_ux.value(),
                self.spin_wind_uy.value(),
                self.spin_wind_uz.value(),
            ),
            wind_turbulence_intensity=self.spin_wind_ti.value(),
            wind_integral_time_scale=self.spin_wind_tl.value(),
            wind_seed=seed_val if seed_val > 0 else None,
        )

    def set_params(self, params: CableParams) -> None:
        """Populate the form fields from a CableParams.

        Temporarily blocks signals while spinboxes are updated so
        ``params_changed`` fires once at the end instead of after every
        individual setValue() call.
        """
        # Block signals to avoid firing params_changed N times
        all_spins = (
            self.spin_ax, self.spin_ay, self.spin_az,
            self.spin_bx, self.spin_by, self.spin_bz,
            self.spin_length, self.spin_segments,
            self.spin_density, self.spin_ea, self.spin_damping,
            self.spin_diameter, self.spin_gravity,
            self.spin_max_steps, self.spin_tol, self.spin_snap,
            self.spin_dt, self.spin_t_end, self.spin_out_int,
            self.spin_fluid_density, self.spin_drag_cd,
            self.spin_wind_ux, self.spin_wind_uy, self.spin_wind_uz,
            self.spin_wind_ti, self.spin_wind_tl, self.spin_wind_seed,
        )
        all_combos = (self.cbo_fluid, self.cbo_wind, self.combo_mode)
        for s in all_spins:
            s.blockSignals(True)
        for c in all_combos:
            c.blockSignals(True)
        try:
            self.spin_ax.setValue(float(params.point_a[0]))
            self.spin_ay.setValue(float(params.point_a[1]))
            self.spin_az.setValue(float(params.point_a[2]))
            self.spin_bx.setValue(float(params.point_b[0]))
            self.spin_by.setValue(float(params.point_b[1]))
            self.spin_bz.setValue(float(params.point_b[2]))
            self.spin_length.setValue(float(params.cable_length))
            self.spin_segments.setValue(int(params.n_segments))
            self.spin_density.setValue(float(params.line_density))
            self.spin_ea.setValue(float(params.EA))
            self.spin_damping.setValue(float(params.damping))
            self.spin_diameter.setValue(float(params.diameter))
            self.spin_gravity.setValue(float(params.gravity))
            self.spin_max_steps.setValue(int(params.max_equilibrium_steps))
            self.spin_tol.setValue(float(params.equilibrium_tol))
            self.spin_snap.setValue(int(params.snapshot_interval))
            # Mode & dynamic fields
            self.combo_mode.setCurrentText(str(params.mode))
            self.spin_dt.setValue(float(params.dt) if params.dt > 0 else 0.005)
            self.spin_t_end.setValue(float(params.t_end) if params.t_end > 0 else 5.0)
            self.spin_out_int.setValue(float(params.output_interval))
            # Fluid & wind
            self.cbo_fluid.setCurrentText(str(params.fluid))
            preset = _FLUID_PRESETS[params.fluid]
            self.spin_fluid_density.setValue(
                float(params.fluid_density) if params.fluid_density is not None else preset["density"]
            )
            self.spin_drag_cd.setValue(
                float(params.drag_Cd) if params.drag_Cd is not None else preset["Cd"]
            )
            self.cbo_wind.setCurrentText(str(params.wind_type))
            self.spin_wind_ux.setValue(float(params.wind_U_mean[0]))
            self.spin_wind_uy.setValue(float(params.wind_U_mean[1]))
            self.spin_wind_uz.setValue(float(params.wind_U_mean[2]))
            self.spin_wind_ti.setValue(float(params.wind_turbulence_intensity))
            self.spin_wind_tl.setValue(float(params.wind_integral_time_scale))
            self.spin_wind_seed.setValue(int(params.wind_seed) if params.wind_seed is not None else 0)
        finally:
            for s in all_spins:
                s.blockSignals(False)
            for c in all_combos:
                c.blockSignals(False)
        self._update_wind_row_enabled(params.wind_type)
        self._on_mode_changed(str(params.mode))

        # Fire once after all fields are updated
        self.params_changed.emit(self.collect_params())

    def set_output_dir(self, path: str) -> None:
        """Update the output directory display."""
        if path:
            self._output_dir_label.setText(f"Output: {path}")
        else:
            self._output_dir_label.setText("Output: ~/.cache/pycable/runs/")

    def set_initial_condition(self, ic: str, tension_top: float = 0.0,
                              tension_bot: float = 0.0) -> None:
        """Set the initial condition dropdown and tension fields."""
        self.combo_initial_cond.blockSignals(True)
        self.combo_initial_cond.setCurrentText(ic)
        self.combo_initial_cond.blockSignals(False)
        is_tension = (ic == "tension")
        self.spin_tension_top.setEnabled(is_tension)
        self._tension_top_label.setEnabled(is_tension)
        self.spin_tension_bot.setEnabled(is_tension)
        self._tension_bot_label.setEnabled(is_tension)
        self.spin_tension_top.setValue(tension_top)
        self.spin_tension_bot.setValue(tension_bot)

    def set_running(self, running: bool) -> None:
        self.btn_run.setEnabled(not running)
        self.btn_stop.setEnabled(running)
        if running:
            self._status_label.setText("Running cable_solver…")

    def set_status(self, text: str) -> None:
        self._status_label.setText(text)

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _on_run_clicked(self) -> None:
        self.run_requested.emit(self.collect_params())

    def _emit_params_changed(self) -> None:
        self.params_changed.emit(self.collect_params())

    def add_run_all_button(self, btn: QPushButton) -> None:
        """Insert an external 'Run All' button into the Run button row."""
        self._btn_row.addWidget(btn)

    def _on_initial_cond_changed(self, text: str) -> None:
        is_tension = (text == "tension")
        self.spin_tension_top.setEnabled(is_tension)
        self._tension_top_label.setEnabled(is_tension)
        self.spin_tension_bot.setEnabled(is_tension)
        self._tension_bot_label.setEnabled(is_tension)

    def _on_fluid_preset_changed(self, preset_name: str) -> None:
        preset = _FLUID_PRESETS.get(preset_name)
        if not preset:
            return
        self.spin_fluid_density.setValue(float(preset["density"]))
        self.spin_drag_cd.setValue(float(preset["Cd"]))

    def _update_wind_row_enabled(self, wind_type: str) -> None:
        any_wind = wind_type in ("uniform", "AR1")
        ar1 = wind_type == "AR1"
        for s in (self.spin_wind_ux, self.spin_wind_uy, self.spin_wind_uz):
            s.setEnabled(any_wind)
        for s in (self.spin_wind_ti, self.spin_wind_tl, self.spin_wind_seed):
            s.setEnabled(ar1)

    def _on_mode_changed(self, mode: str) -> None:
        is_dynamic = (mode == "dynamic")
        for s in self._dynamic_rows:
            s.setEnabled(is_dynamic)
