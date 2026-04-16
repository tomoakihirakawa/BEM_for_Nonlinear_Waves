"""RunHistoryWidget — tracks past cable solver runs with results.

Displays a list of past runs with cable name, tension summary,
timestamp, and a Load button to re-display the result.
Persisted via QSettings.
"""

from __future__ import annotations

import json
import os
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional

from PySide6.QtCore import Qt, QSettings, Signal
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QApplication,
    QHeaderView,
    QMenu,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)


MAX_HISTORY = 50
SETTINGS_KEY = "run_history"


class RunHistoryWidget(QWidget):
    """Widget displaying past cable solver runs."""

    load_result = Signal(str)  # Emitted with result JSON path to reload

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._qsettings = QSettings("pycable", "pycable")
        self._setup_ui()

    def _setup_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(0)

        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["", "Name", "Top T [kN]", "Bot T [kN]", "Time"])
        self.tree.setRootIsDecorated(False)
        self.tree.setAlternatingRowColors(True)
        self.tree.setSortingEnabled(True)
        self.tree.setSelectionMode(QTreeWidget.ExtendedSelection)
        self.tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._on_context_menu)

        # Delete key
        for key in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            sc = QShortcut(QKeySequence(key), self.tree)
            sc.activated.connect(self._delete_selected)

        header = self.tree.header()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Fixed)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        self.tree.setColumnWidth(0, 56)

        self.tree.sortByColumn(4, Qt.SortOrder.DescendingOrder)

        layout.addWidget(self.tree)
        self.refresh()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def add_entry(self, result_path: str, name: str = "",
                  top_tension: float = 0.0, bottom_tension: float = 0.0,
                  converged: bool = False,
                  n_cables: int = 1, elapsed_ms: float = 0.0) -> None:
        """Record a completed solver run."""
        history = self._load()
        history = [h for h in history if h.get("path") != result_path]
        history.insert(0, {
            "path": result_path,
            "name": name,
            "top_tension": top_tension,
            "bottom_tension": bottom_tension,
            "converged": converged,
            "n_cables": n_cables,
            "elapsed_ms": elapsed_ms,
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        })
        history = history[:MAX_HISTORY]
        self._save(history)
        self.refresh()

    def add_multi_entry(self, output_dir: str, cables: dict,
                        converged: bool = False) -> None:
        """Record a multi-cable system run."""
        for cable_name, cable_data in cables.items():
            result_path = os.path.join(output_dir, f"{cable_name}_result.json")
            if not os.path.isfile(result_path):
                result_path = os.path.join(output_dir, "result.json")
            self.add_entry(
                result_path=result_path,
                name=cable_name,
                top_tension=cable_data.get("top_tension", 0),
                bottom_tension=cable_data.get("bottom_tension", 0),
                converged=converged,
            )

    def refresh(self) -> None:
        """Rebuild the tree from persisted history."""
        self.tree.clear()
        for entry in self._load():
            path = entry.get("path", "")
            # Strip "#cable_name" fragment for file existence check
            actual_path = path.split("#")[0]
            exists = os.path.isfile(actual_path)

            item = QTreeWidgetItem()
            name = entry.get("name", Path(path).stem if path else "?")
            item.setText(1, name if exists else f"{name}  (missing)")
            top_t = entry.get("top_tension", 0)
            bot_t = entry.get("bottom_tension", 0)
            n_cab = entry.get("n_cables", 1)
            if n_cab > 1:
                item.setText(2, f"{n_cab} cables")
                elapsed = entry.get("elapsed_ms", 0)
                item.setText(3, f"{elapsed:.0f} ms" if elapsed else "")
            else:
                item.setText(2, f"{top_t / 1e3:.1f}" if top_t else "")
                item.setText(3, f"{bot_t / 1e3:.1f}" if bot_t else "")
            item.setTextAlignment(2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            item.setTextAlignment(3, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            item.setText(4, entry.get("timestamp", ""))
            item.setData(0, Qt.ItemDataRole.UserRole, entry)

            # Rich tooltip with full details
            converged_str = "Yes" if entry.get("converged") else "No"
            n_cab = entry.get("n_cables", 1)
            elapsed = entry.get("elapsed_ms", 0)
            if n_cab > 1:
                tooltip = (
                    f"Case: {name}\n"
                    f"Cables: {n_cab}\n"
                    f"Converged: {converged_str}\n"
                    f"Elapsed: {elapsed:.0f} ms\n"
                    f"Time: {entry.get('timestamp', '')}\n"
                    f"Path: {actual_path}"
                )
            else:
                tooltip = (
                    f"Name: {name}\n"
                    f"Top tension: {top_t / 1e3:.2f} kN\n"
                    f"Bottom tension: {bot_t / 1e3:.2f} kN\n"
                    f"Converged: {converged_str}\n"
                    f"Time: {entry.get('timestamp', '')}\n"
                    f"Path: {actual_path}"
                )
            for col in range(self.tree.columnCount()):
                item.setToolTip(col, tooltip)

            self.tree.addTopLevelItem(item)

            btn = QPushButton("Load")
            btn.setFixedSize(50, 24)
            btn.setEnabled(exists)
            btn.clicked.connect(lambda checked, p=actual_path: self.load_result.emit(p))
            self.tree.setItemWidget(item, 0, btn)

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _load(self) -> list:
        try:
            data = self._qsettings.value(SETTINGS_KEY, "[]")
            if isinstance(data, list):
                return data
            return json.loads(data) if data else []
        except (json.JSONDecodeError, TypeError):
            return []

    def _save(self, history: list) -> None:
        self._qsettings.setValue(SETTINGS_KEY, json.dumps(history))

    def _delete_selected(self) -> None:
        selected = self.tree.selectedItems()
        if not selected:
            return
        paths_to_remove = set()
        for item in selected:
            entry = item.data(0, Qt.ItemDataRole.UserRole)
            if entry:
                paths_to_remove.add(entry.get("path", ""))
        if paths_to_remove:
            history = self._load()
            history = [h for h in history if h.get("path") not in paths_to_remove]
            self._save(history)
            self.refresh()

    def _on_context_menu(self, pos) -> None:
        item = self.tree.itemAt(pos)
        if not item:
            return
        entry = item.data(0, Qt.ItemDataRole.UserRole)
        if not entry:
            return
        path = entry.get("path", "")

        menu = QMenu(self)
        remove_action = menu.addAction("Remove from history")
        copy_action = menu.addAction("Copy path")
        open_action = None
        parent_dir = os.path.dirname(path)
        if os.path.isdir(parent_dir):
            open_action = menu.addAction("Open in Finder")

        action = menu.exec(self.tree.mapToGlobal(pos))
        if action == remove_action:
            self._delete_selected()
        elif action == copy_action:
            QApplication.clipboard().setText(path)
        elif open_action and action == open_action:
            subprocess.Popen(["open", parent_dir])
