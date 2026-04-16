"""TimePlayerWidget — Paraview-like time slider + play/pause button.

Appears below the 3D view. Emits ``time_changed(float)`` when the user
scrubs the slider or the auto-play timer advances. Disabled until the
main window tells it to ``set_times(times: list[float])``.
"""

from __future__ import annotations

from typing import Optional

from PySide6.QtCore import QTimer, Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QFormLayout, QHBoxLayout, QLabel, QPushButton, QSlider, QSpinBox,
    QVBoxLayout, QWidget,
)


class ContourSettingsDialog(QDialog):
    """Popup for contour colormap and range settings."""

    def __init__(self, parent=None, *,
                 cmap: str = "coolwarm",
                 auto_range: bool = True,
                 vmin: float = 0.0,
                 vmax: float = 1.0) -> None:
        super().__init__(parent)
        self.setWindowTitle("Contour Settings")

        self._cmap_combo = QComboBox()
        self._cmap_combo.addItems([
            "coolwarm", "viridis", "plasma", "jet", "RdBu", "seismic",
        ])
        self._cmap_combo.setCurrentText(cmap)

        self._auto_check = QCheckBox("Auto range")
        self._auto_check.setChecked(auto_range)
        self._auto_check.toggled.connect(self._refresh_enabled)

        self._vmin_spin = QDoubleSpinBox()
        self._vmin_spin.setRange(-1e12, 1e12)
        self._vmin_spin.setDecimals(4)
        self._vmin_spin.setValue(vmin)

        self._vmax_spin = QDoubleSpinBox()
        self._vmax_spin.setRange(-1e12, 1e12)
        self._vmax_spin.setDecimals(4)
        self._vmax_spin.setValue(vmax)

        self._refresh_enabled(auto_range)

        form = QFormLayout()
        form.addRow("Colormap", self._cmap_combo)
        form.addRow("", self._auto_check)
        form.addRow("Min", self._vmin_spin)
        form.addRow("Max", self._vmax_spin)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def _refresh_enabled(self, auto: bool) -> None:
        self._vmin_spin.setEnabled(not auto)
        self._vmax_spin.setEnabled(not auto)

    def result_settings(self) -> dict:
        return {
            "cmap": self._cmap_combo.currentText(),
            "auto_range": self._auto_check.isChecked(),
            "vmin": self._vmin_spin.value(),
            "vmax": self._vmax_spin.value(),
        }


class TimePlayerWidget(QWidget):

    time_changed = Signal(float)
    show_static_changed = Signal(bool)
    contour_settings_changed = Signal(dict)

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)

        self._times: list[float] = []
        self._contour_settings = {
            "cmap": "coolwarm", "auto_range": True, "vmin": 0.0, "vmax": 1.0,
        }

        self._slider = QSlider(Qt.Orientation.Horizontal)
        self._slider.setMinimum(0)
        self._slider.setMaximum(0)
        self._slider.setEnabled(False)
        self._slider.valueChanged.connect(self._on_slider_changed)

        self._btn_play = QPushButton("▶")
        self._btn_play.setFixedWidth(40)
        self._btn_play.setEnabled(False)
        self._btn_play.clicked.connect(self._on_play_toggled)

        self._btn_first = QPushButton("⏮")
        self._btn_first.setFixedWidth(40)
        self._btn_first.setEnabled(False)
        self._btn_first.clicked.connect(lambda: self._slider.setValue(0))

        self._btn_last = QPushButton("⏭")
        self._btn_last.setFixedWidth(40)
        self._btn_last.setEnabled(False)
        self._btn_last.clicked.connect(
            lambda: self._slider.setValue(self._slider.maximum())
        )

        self._label = QLabel("t = —")
        self._label.setMinimumWidth(120)

        self._step_label = QLabel("step:")
        self._step_spin = QSpinBox()
        self._step_spin.setRange(1, 1000)
        self._step_spin.setValue(1)
        self._step_spin.setFixedWidth(60)
        self._step_spin.setToolTip("Advance by N frames per play tick")

        self._fps_label = QLabel("fps:")
        self._fps_spin = QSpinBox()
        self._fps_spin.setRange(1, 120)
        self._fps_spin.setValue(20)
        self._fps_spin.setFixedWidth(60)
        self._fps_spin.valueChanged.connect(self._on_fps_changed)

        self._timer = QTimer(self)
        self._timer.setInterval(50)  # 20 fps default
        self._timer.timeout.connect(self._advance_frame)

        self._chk_static = QCheckBox("static ref")
        self._chk_static.setChecked(True)
        self._chk_static.setToolTip("Show the static equilibrium shape as a gray overlay")
        self._chk_static.toggled.connect(self.show_static_changed.emit)

        self._btn_contour = QPushButton("Contour")
        self._btn_contour.setFixedWidth(70)
        self._btn_contour.setToolTip("Configure contour colormap and range")
        self._btn_contour.clicked.connect(self._open_contour_dialog)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(4, 2, 4, 2)
        layout.addWidget(self._btn_first)
        layout.addWidget(self._btn_play)
        layout.addWidget(self._btn_last)
        layout.addWidget(self._slider, stretch=1)
        layout.addWidget(self._step_label)
        layout.addWidget(self._step_spin)
        layout.addWidget(self._fps_label)
        layout.addWidget(self._fps_spin)
        layout.addWidget(self._chk_static)
        layout.addWidget(self._btn_contour)
        layout.addWidget(self._label)

    def set_times(self, times: list[float]) -> None:
        """Initialize the slider with a sorted list of available times."""
        self._times = sorted(set(times))
        self._timer.stop()
        self._btn_play.setText("▶")
        has_frames = len(self._times) >= 2
        self._slider.setEnabled(has_frames)
        self._btn_play.setEnabled(has_frames)
        self._btn_first.setEnabled(has_frames)
        self._btn_last.setEnabled(has_frames)
        if has_frames:
            self._slider.blockSignals(True)
            self._slider.setMinimum(0)
            self._slider.setMaximum(len(self._times) - 1)
            self._slider.setValue(0)
            self._slider.blockSignals(False)
            self._update_label(0)
            # Emit initial time so viewer syncs to t0
            self.time_changed.emit(self._times[0])
        else:
            self._slider.setMinimum(0)
            self._slider.setMaximum(0)
            self._label.setText("t = —")

    def clear(self) -> None:
        self.set_times([])

    def current_time(self) -> float:
        if not self._times:
            return 0.0
        i = self._slider.value()
        i = max(0, min(i, len(self._times) - 1))
        return self._times[i]

    def _on_slider_changed(self, value: int) -> None:
        if not self._times:
            return
        self._update_label(value)
        self.time_changed.emit(self._times[value])

    def _update_label(self, value: int) -> None:
        if 0 <= value < len(self._times):
            self._label.setText(f"t = {self._times[value]:.3f} s "
                                f"({value + 1}/{len(self._times)})")

    def _on_play_toggled(self) -> None:
        if self._timer.isActive():
            self._timer.stop()
            self._btn_play.setText("▶")
        else:
            if self._slider.value() >= self._slider.maximum():
                self._slider.setValue(0)
            self._timer.start()
            self._btn_play.setText("⏸")

    def _advance_frame(self) -> None:
        step = max(1, self._step_spin.value())
        v = self._slider.value() + step
        if v > self._slider.maximum():
            self._timer.stop()
            self._btn_play.setText("▶")
            self._slider.setValue(self._slider.maximum())
            return
        self._slider.setValue(v)

    # ------------------------------------------------------------------
    # Contour settings
    # ------------------------------------------------------------------

    @property
    def contour_settings(self) -> dict:
        return dict(self._contour_settings)

    def _open_contour_dialog(self) -> None:
        dlg = ContourSettingsDialog(
            self,
            cmap=self._contour_settings["cmap"],
            auto_range=self._contour_settings["auto_range"],
            vmin=self._contour_settings["vmin"],
            vmax=self._contour_settings["vmax"],
        )
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._contour_settings = dlg.result_settings()
            self.contour_settings_changed.emit(self._contour_settings)

    def _on_fps_changed(self, fps: int) -> None:
        self._timer.setInterval(max(1, int(1000 / max(1, fps))))
