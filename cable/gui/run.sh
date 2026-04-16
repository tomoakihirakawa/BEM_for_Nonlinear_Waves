#!/usr/bin/env bash
# run.sh — bootstrap venv if needed, then launch the pycable GUI.
#
# Usage (from either dev or public tree):
#   cd <cpp>/cable/gui           # dev tree
#   cd ~/CableDynamics/cable/gui # public tree
#   ./run.sh
#
# Prerequisites: Python 3.12 (`brew install python@3.12`) and the cable_solver
# binary built either at <cpp>/cable/build_solver/cable_solver (dev) or
# ~/CableDynamics/build/cable_solver (public). solver_discovery.py finds
# either automatically.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── venv location ────────────────────────────────────────────────
# Place the venv on local disk, NOT inside Dropbox/iCloud.
# Cloud-synced folders corrupt Qt framework bundles (.framework/
# inside PySide6), causing "Could not load cocoa" errors.
VENV_DIR="${PYCABLE_VENV:-$HOME/.cache/pycable_venv}"
STAMP_FILE="$VENV_DIR/.requirements_stamp"

# ── find Python ──────────────────────────────────────────────────
PYTHON_BIN="${PYTHON_BIN:-}"
if [ -z "$PYTHON_BIN" ]; then
  for py in python3.12 python3.11 python3.10 python3.13 python3 python; do
    if command -v "$py" >/dev/null 2>&1; then
      if "$py" -c "import sys; raise SystemExit(0 if (3,9)<=sys.version_info[:2]<=(3,13) else 1)" 2>/dev/null; then
        PYTHON_BIN="$py"
        break
      fi
    fi
  done
fi

if [ -z "$PYTHON_BIN" ]; then
  echo "error: Python 3.9–3.13 not found on PATH." >&2
  echo "install via:  brew install python@3.12" >&2
  exit 1
fi

PY_VER=$("$PYTHON_BIN" -c "import sys; print(f'{sys.version_info[0]}.{sys.version_info[1]}')")
echo "[run.sh] Python $PY_VER ($PYTHON_BIN)"

# ── create / rebuild venv ────────────────────────────────────────
REBUILD=0
if [ -d "$VENV_DIR" ]; then
  VENV_PYTHON="$VENV_DIR/bin/python"
  if [ ! -x "$VENV_PYTHON" ]; then
    REBUILD=1
  else
    VENV_VER=$("$VENV_PYTHON" -c "import sys; print(f'{sys.version_info[0]}.{sys.version_info[1]}')" 2>/dev/null || echo "")
    if [ "$VENV_VER" != "$PY_VER" ]; then
      REBUILD=1
    fi
  fi
else
  REBUILD=1
fi

if [ "$REBUILD" -eq 1 ]; then
  echo "[run.sh] creating venv at $VENV_DIR"
  rm -rf "$VENV_DIR"
  "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

# ── install / update packages ────────────────────────────────────
REQ_HASH=$("$PYTHON_BIN" -c "
import hashlib, pathlib, sys
p = pathlib.Path('$SCRIPT_DIR/pyproject.toml')
print(hashlib.sha256(p.read_bytes()).hexdigest())
")
STAMP_CONTENT="python=$PY_VER
pyproject_sha256=$REQ_HASH"

NEEDS_INSTALL=1
if [ -f "$STAMP_FILE" ]; then
  EXISTING_STAMP=$(cat "$STAMP_FILE")
  if [ "$EXISTING_STAMP" = "$STAMP_CONTENT" ]; then
    NEEDS_INSTALL=0
  fi
fi

if [ "$NEEDS_INSTALL" -eq 1 ]; then
  echo "[run.sh] installing pycable"
  pip install --upgrade pip
  pip install -e .
  printf '%s\n' "$STAMP_CONTENT" > "$STAMP_FILE"
fi

# ── sanity-check solver binary ───────────────────────────────────
python - <<'PY' || true
from pycable.solver_discovery import find_cable_solver, CableSolverNotFound
try:
    path = find_cable_solver(must_exist=True)
    print(f"[run.sh] cable_solver = {path}")
except CableSolverNotFound as e:
    print(f"[run.sh] WARNING: {e}")
PY

# ── launch ───────────────────────────────────────────────────────
exec python -m pycable
