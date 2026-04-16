"""BodiesListWidget — sidebar list of rigid bodies loaded alongside cables.

A "body" is a separate JSON file listed in settings.json's ``input_files``
with ``type: "RigidBody"`` (BEM convention). Each body may define a
prescribed motion via the ``velocity`` array (e.g. sinusoidal heave).
Cables reference bodies through ``end_a_body`` / ``end_b_body``.

This widget mirrors the loaded bodies so the user can see what's in the
system and which cable is attached to which body.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QListWidget, QListWidgetItem, QWidget


@dataclass
class BodySpec:
    name: str = ""
    type: str = "RigidBody"
    velocity: list[Any] = field(default_factory=list)  # e.g. ["sinusoidal", 0, 2, 20, 0, 0, 1]

    @property
    def motion_label(self) -> str:
        if not self.velocity:
            return "fixed"
        kind = str(self.velocity[0]) if self.velocity else "fixed"
        if kind in ("sinusoidal", "sin", "cos") and len(self.velocity) >= 7:
            amp, period = self.velocity[2], self.velocity[3]
            return f"{kind}  A={amp}  T={period}s"
        if kind == "cantilever" and len(self.velocity) >= 14:
            amp, period, length = self.velocity[2], self.velocity[3], self.velocity[13]
            base = f"cantilever  A_tip={amp}  T={period}s  L={length}m"
            if len(self.velocity) >= 17:
                amp2, period2 = self.velocity[14], self.velocity[15]
                return f"{base}  [+mode2 A={amp2} T={period2}s]"
            return base
        return kind


class BodiesListWidget(QListWidget):
    """List widget for bodies. Emits body_selected(int) on user selection."""

    body_selected = Signal(int)
    body_double_clicked = Signal(int)

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(60)
        self.setMaximumHeight(120)
        self.itemSelectionChanged.connect(self._on_selection_changed)
        self.itemDoubleClicked.connect(self._on_double_clicked)

    def populate(self, bodies: list[BodySpec]) -> None:
        self.blockSignals(True)
        self.clear()
        for b in bodies:
            QListWidgetItem(f"{b.name}  [{b.motion_label}]", self)
        self.blockSignals(False)

    def selected_index(self) -> int:
        row = self.currentRow()
        return row if row >= 0 else -1

    def _on_selection_changed(self) -> None:
        idx = self.selected_index()
        if idx >= 0:
            self.body_selected.emit(idx)

    def _on_double_clicked(self, _item) -> None:
        idx = self.selected_index()
        if idx >= 0:
            self.body_double_clicked.emit(idx)
