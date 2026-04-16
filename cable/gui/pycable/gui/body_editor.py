"""BodyEditorDialog — edit a body's motion prescription.

Opens when the user double-clicks a body in BodiesListWidget. Writes the
updated ``velocity`` array back to the body's JSON file on Save.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import math

from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QFormLayout, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget,
)

from .bodies_list import BodySpec


_AXIS_OPTIONS = {
    "surge (+x)":  [1.0, 0.0, 0.0],
    "sway (+y)":   [0.0, 1.0, 0.0],
    "heave (+z)":  [0.0, 0.0, 1.0],
}


def _axis_label(axis: list) -> str:
    """Return the closest cardinal axis label for a 3-vector."""
    if len(axis) < 3:
        return "heave (+z)"
    ax, ay, az = axis[0], axis[1], axis[2]
    absx, absy, absz = abs(ax), abs(ay), abs(az)
    if absx >= absy and absx >= absz:
        return "surge (+x)"
    if absy >= absz:
        return "sway (+y)"
    return "heave (+z)"


class BodyEditorDialog(QDialog):
    """Simple form for editing a single BodySpec's motion."""

    def __init__(self, body: BodySpec, file_path: Path,
                 parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle(f"Edit body — {body.name}")
        self._body = body
        self._file_path = file_path

        self._name_label = QLabel(body.name)

        self._kind_combo = QComboBox()
        self._kind_combo.addItems(["fixed", "sinusoidal", "cantilever"])

        self._start_spin = QDoubleSpinBox()
        self._start_spin.setRange(0.0, 1e6)
        self._start_spin.setDecimals(4)

        self._amp_spin = QDoubleSpinBox()
        self._amp_spin.setRange(0.0, 1e6)
        self._amp_spin.setDecimals(4)

        self._period_spin = QDoubleSpinBox()
        self._period_spin.setRange(0.0001, 1e6)
        self._period_spin.setDecimals(6)

        self._axis_combo = QComboBox()
        self._axis_combo.addItems(list(_AXIS_OPTIONS.keys()))

        # --- Cantilever-specific fields ---
        self._fix_x = QDoubleSpinBox(); self._fix_x.setRange(-1e6, 1e6); self._fix_x.setDecimals(3)
        self._fix_y = QDoubleSpinBox(); self._fix_y.setRange(-1e6, 1e6); self._fix_y.setDecimals(3)
        self._fix_z = QDoubleSpinBox(); self._fix_z.setRange(-1e6, 1e6); self._fix_z.setDecimals(3)

        self._axis_beam_combo = QComboBox()
        self._axis_beam_combo.addItems(list(_AXIS_OPTIONS.keys()))
        self._axis_beam_combo.setCurrentText("heave (+z)")  # typical: tower extends upward

        self._bend_combo = QComboBox()
        self._bend_combo.addItems(list(_AXIS_OPTIONS.keys()))
        self._bend_combo.setCurrentText("surge (+x)")  # typical: bends in x

        self._length_spin = QDoubleSpinBox()
        self._length_spin.setRange(0.001, 1e6)
        self._length_spin.setDecimals(3)
        self._length_spin.setValue(1.0)

        # Mode 2 fields (cantilever only)
        self._mode2_check = QCheckBox("Enable mode 2")
        self._mode2_check.setChecked(False)

        self._amp2_spin = QDoubleSpinBox()
        self._amp2_spin.setRange(0.0, 1e6)
        self._amp2_spin.setDecimals(4)

        self._period2_spin = QDoubleSpinBox()
        self._period2_spin.setRange(0.0001, 1e6)
        self._period2_spin.setDecimals(6)

        self._period2_auto_btn = QPushButton("Auto (T₁/6.267)")
        self._period2_auto_btn.clicked.connect(
            lambda: self._period2_spin.setValue(self._period_spin.value() / 6.267))

        self._phase2_spin = QDoubleSpinBox()
        self._phase2_spin.setRange(-360.0, 360.0)
        self._phase2_spin.setDecimals(2)
        self._phase2_spin.setSuffix(" deg")

        # Populate from body.velocity
        v = body.velocity
        if v and str(v[0]) in ("sinusoidal", "sin", "cos"):
            self._kind_combo.setCurrentText("sinusoidal")
            if len(v) >= 2: self._start_spin.setValue(float(v[1]))
            if len(v) >= 3: self._amp_spin.setValue(float(v[2]))
            if len(v) >= 4: self._period_spin.setValue(float(v[3]))
            if len(v) >= 7: self._axis_combo.setCurrentText(
                _axis_label([float(v[4]), float(v[5]), float(v[6])])
            )
        elif v and str(v[0]) == "cantilever":
            self._kind_combo.setCurrentText("cantilever")
            if len(v) >= 2: self._start_spin.setValue(float(v[1]))
            if len(v) >= 3: self._amp_spin.setValue(float(v[2]))
            if len(v) >= 4: self._period_spin.setValue(float(v[3]))
            if len(v) >= 7:
                self._fix_x.setValue(float(v[4]))
                self._fix_y.setValue(float(v[5]))
                self._fix_z.setValue(float(v[6]))
            if len(v) >= 10:
                self._axis_beam_combo.setCurrentText(
                    _axis_label([float(v[7]), float(v[8]), float(v[9])])
                )
            if len(v) >= 13:
                self._bend_combo.setCurrentText(
                    _axis_label([float(v[10]), float(v[11]), float(v[12])])
                )
            if len(v) >= 14:
                self._length_spin.setValue(float(v[13]))
            if len(v) >= 17:
                self._mode2_check.setChecked(True)
                self._amp2_spin.setValue(float(v[14]))
                self._period2_spin.setValue(float(v[15]))
                self._phase2_spin.setValue(math.degrees(float(v[16])))
        else:
            self._kind_combo.setCurrentText("fixed")

        form = QFormLayout()
        form.addRow("Name", self._name_label)
        form.addRow("Motion kind", self._kind_combo)
        form.addRow("Start time [s]", self._start_spin)
        form.addRow("Amplitude [m]", self._amp_spin)
        form.addRow("Period [s]", self._period_spin)
        form.addRow("Axis (sinusoidal)", self._axis_combo)
        form.addRow("Fix point x", self._fix_x)
        form.addRow("Fix point y", self._fix_y)
        form.addRow("Fix point z", self._fix_z)
        form.addRow("Beam axis (upward)", self._axis_beam_combo)
        form.addRow("Bend direction", self._bend_combo)
        form.addRow("Beam length [m]", self._length_spin)
        form.addRow("", self._mode2_check)
        form.addRow("Amplitude mode 2 [m]", self._amp2_spin)
        period2_row = QHBoxLayout()
        period2_row.addWidget(self._period2_spin)
        period2_row.addWidget(self._period2_auto_btn)
        form.addRow("Period mode 2 [s]", period2_row)
        form.addRow("Phase mode 2", self._phase2_spin)

        self._kind_combo.currentTextChanged.connect(self._refresh_enabled)
        self._mode2_check.toggled.connect(self._refresh_enabled)
        self._refresh_enabled()

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save
            | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._save_and_accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def _refresh_enabled(self) -> None:
        kind = self._kind_combo.currentText()
        is_sin = kind == "sinusoidal"
        is_cant = kind == "cantilever"
        is_moving = is_sin or is_cant
        for w in (self._start_spin, self._amp_spin, self._period_spin):
            w.setEnabled(is_moving)
        self._axis_combo.setEnabled(is_sin)
        for w in (self._fix_x, self._fix_y, self._fix_z,
                  self._axis_beam_combo, self._bend_combo, self._length_spin):
            w.setEnabled(is_cant)
        self._mode2_check.setEnabled(is_cant)
        mode2_on = is_cant and self._mode2_check.isChecked()
        for w in (self._amp2_spin, self._period2_spin,
                  self._period2_auto_btn, self._phase2_spin):
            w.setEnabled(mode2_on)

    def _save_and_accept(self) -> None:
        kind = self._kind_combo.currentText()
        if kind == "sinusoidal":
            axis_vec = _AXIS_OPTIONS[self._axis_combo.currentText()]
            velocity = [
                "sinusoidal",
                float(self._start_spin.value()),
                float(self._amp_spin.value()),
                float(self._period_spin.value()),
                axis_vec[0], axis_vec[1], axis_vec[2],
            ]
        elif kind == "cantilever":
            axis_vec = _AXIS_OPTIONS[self._axis_beam_combo.currentText()]
            bend_vec = _AXIS_OPTIONS[self._bend_combo.currentText()]
            velocity = [
                "cantilever",
                float(self._start_spin.value()),
                float(self._amp_spin.value()),
                float(self._period_spin.value()),
                float(self._fix_x.value()),
                float(self._fix_y.value()),
                float(self._fix_z.value()),
                axis_vec[0], axis_vec[1], axis_vec[2],
                bend_vec[0], bend_vec[1], bend_vec[2],
                float(self._length_spin.value()),
            ]
            if self._mode2_check.isChecked() and self._amp2_spin.value() > 0:
                velocity.extend([
                    float(self._amp2_spin.value()),
                    float(self._period2_spin.value()),
                    math.radians(float(self._phase2_spin.value())),
                ])
        else:
            velocity = ["fixed"]

        # Update in-memory spec
        self._body.velocity = velocity

        # Write to file (preserve other keys)
        try:
            with self._file_path.open() as f:
                d = json.load(f)
        except Exception:
            d = {"name": self._body.name, "type": self._body.type}
        d["velocity"] = velocity
        d.setdefault("name", self._body.name)
        d.setdefault("type", self._body.type)
        with self._file_path.open("w") as f:
            json.dump(d, f, indent=2)

        self.accept()
