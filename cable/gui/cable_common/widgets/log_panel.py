"""LogPanel — a scrolling read-only text view for solver stdout."""

from __future__ import annotations

from PySide6.QtGui import QFont, QTextCursor
from PySide6.QtWidgets import QTextEdit, QWidget


class LogPanel(QTextEdit):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setReadOnly(True)
        self.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)

        font = QFont("Menlo")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(11)
        self.setFont(font)

        self._max_lines = 2000

    def append_line(self, text: str) -> None:
        """Append one line and auto-scroll to the bottom."""
        self.append(text)
        # Cheap cap so long runs don't balloon memory.
        doc = self.document()
        while doc.blockCount() > self._max_lines:
            cursor = QTextCursor(doc)
            cursor.movePosition(QTextCursor.MoveOperation.Start)
            cursor.select(QTextCursor.SelectionType.BlockUnderCursor)
            cursor.removeSelectedText()
            cursor.deleteChar()
        self.moveCursor(QTextCursor.MoveOperation.End)
