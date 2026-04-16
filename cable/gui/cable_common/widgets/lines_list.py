"""LinesListWidget — sidebar list of cables in the loaded system.

Phase 4 (2026-04-12) addition. When the loaded JSON is a multi-cable
system (`mooring_*` flat keys, or one of the new System samples), this
widget displays the list of cable names. Selecting a cable populates the
SetupPanel form below. The widget is hidden / shows a single placeholder
entry when the system has only one cable (legacy single-line mode).
"""

from __future__ import annotations

from typing import Optional

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QListWidget, QListWidgetItem, QWidget

from cable_common.params import LumpedCableSystemParams


class LinesListWidget(QListWidget):
    """A simple QListWidget that mirrors the cables of a LumpedCableSystemParams."""

    cable_selected = Signal(int)  # int = cable index (0-based)

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(80)
        self.setMaximumHeight(220)
        self.itemSelectionChanged.connect(self._on_selection_changed)

    def populate_from_system(self, system: LumpedCableSystemParams) -> None:
        """Replace the list with the cables of the given system."""
        self.blockSignals(True)
        self.clear()
        for cable in system.cables:
            QListWidgetItem(cable.name, self)
        if self.count() > 0:
            self.setCurrentRow(0)
        self.blockSignals(False)

    def selected_index(self) -> int:
        """Return the currently selected cable index (-1 if none)."""
        row = self.currentRow()
        return row if row >= 0 else -1

    def _on_selection_changed(self) -> None:
        idx = self.selected_index()
        if idx >= 0:
            self.cable_selected.emit(idx)
