"""QApplication bootstrap for the pycable GUI.

Usage:
    python -m pycable                     # launch with default catenary params
    python -m pycable path/to/input.json  # launch and pre-load a CableParams JSON
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from PySide6.QtWidgets import QApplication


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="pycable",
        description="Python GUI wrapper around the CableDynamics cable_solver.",
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default=None,
        help="Optional CableParams JSON file to pre-load into the form.",
    )
    # Ignore Qt's own args passed to QApplication (e.g. -style, -platform)
    return parser.parse_known_args(argv[1:])[0]


def main() -> int:
    """Entry point for ``python -m pycable`` and the ``pycable`` console script."""
    args = _parse_args(sys.argv)

    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("pycable")
    app.setOrganizationName("CableDynamics")

    # Late import so that --help on broken installs doesn't need pyvista yet.
    from .main_window import MainWindow

    window = MainWindow()
    if args.input_file:
        window.load_input_file(Path(args.input_file))
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
