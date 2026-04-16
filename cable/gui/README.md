# pycable — Python GUI for CableDynamics

A thin PySide6 / PyVista GUI that drives the `cable_solver` C++ binary.
Enter cable parameters in a form, run the equilibrium solver, and watch
the cable converge to its catenary shape in a live 3D view.

**This package does not contain physics.** All of the mechanics (tension,
drag, gravity, RK4 integration) lives in `cable/cable_solver.cpp` and is
built as the `cable_solver` C++ binary. The GUI writes an input JSON,
launches the binary, streams its stdout for live updates, and reads the
final `result.json`.

**Target OS: macOS.** Tested on Apple Silicon with Homebrew Python 3.12.

## Where this lives

This directory exists in two places:

- **Dev tree**: `<cpp>/cable/gui/` — the source of truth. Edit here.
- **Public tree**: `~/CableDynamics/cable/gui/` — synced from dev via
  `<cpp>/sync_all_public.sh`. Do not edit directly.

`solver_discovery.py` auto-detects the enclosing tree and looks for the
binary in the right place:

| Layout | cable_solver path |
|---|---|
| Dev | `<cpp>/cable/build_solver/cable_solver` |
| Public | `~/CableDynamics/build/cable_solver` |

## Prerequisites

1. **Xcode Command Line Tools**:
   ```bash
   xcode-select --install
   ```

2. **Homebrew** build deps:
   ```bash
   brew install gcc cmake lapack openblas libomp
   ```

3. **Python 3.12**:
   ```bash
   brew install python@3.12
   ```

## Build the C++ solver

### Dev tree

The dev `cpp/cable/build_solver/` directory already has a CMake setup.
Rebuild with:

```bash
cd <cpp>/cable/build_solver
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(sysctl -n hw.logicalcpu)
```

### Public tree (for handing to a friend)

```bash
cd ~/CableDynamics
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(sysctl -n hw.logicalcpu)
```

Verify either build with:

```bash
./cable_solver <cable/gui>/examples/synthetic/catenary_500m.json /tmp/cable_out/
python -m json.tool /tmp/cable_out/result.json | head -15
```

You should see `"n_nodes": 41` and a `"positions"` array with the sagged
catenary shape.

## Install and run the GUI

```bash
cd <cpp>/cable/gui         # or ~/CableDynamics/cable/gui
./run.sh
```

The first run creates a `.venv` next to `run.sh`, pip-installs pycable in
editable mode, and launches the GUI. Subsequent runs reuse the venv.

Manual path if you prefer:

```bash
cd <cpp>/cable/gui
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e .
python -m pycable
```

### GUI layout

- **Left dock — Parameters**: three collapsible groups.
  - *Geometry*: Point A (x,y,z), Point B (x,y,z), cable length, segment count.
  - *Material (advanced)*: line density, EA, damping, diameter, gravity.
  - *Solver (advanced)*: max steps, convergence tolerance, snapshot interval.
- **Central view — 3D**: PyVista canvas. Red sphere = Point A, green sphere =
  Point B. A faint gray straight line previews the initial geometry. During
  a run, a blue curve updates on every snapshot. On convergence, the cable
  is recolored by tension magnitude (viridis scale).
- **Bottom dock — Log**: live solver stdout. SNAPSHOT lines are summarized
  as `iter N  |v|_max=...`; anything else from the solver is shown verbatim.

## Run the C++ CLI directly (no GUI)

```bash
./cable_solver path/to/input.json /tmp/cable_out/
python -m json.tool /tmp/cable_out/result.json
```

The result JSON matches exactly what the GUI reads on finish — so anything
you can do in the GUI you can do via a shell script by writing your own
input JSON.

## Example input files

See [`examples/README.md`](examples/README.md) for descriptions of the
bundled example cases and parameter reference.

## Unit tests

```bash
cd <cpp>/cable/gui
source .venv/bin/activate
pip install -e '.[test]'
pytest tests/ -v
```

`test_integration.py` and `test_bridge_lifecycle.py` are skipped unless
`cable_solver` is built.

## Troubleshooting

### `cable_solver not found`

`solver_discovery.py` looks in this order:

1. `$PYCABLE_SOLVER_PATH` (strict override — must exist if set)
2. `$CABLE_DYNAMICS_ROOT/build/cable_solver`
3. `<cable>/build_solver/cable_solver` (dev layout)
4. `<cable>/../build/cable_solver` (public layout)
5. `~/CableDynamics/build/cable_solver`
6. `which cable_solver`

If the binary is somewhere else:

```bash
export PYCABLE_SOLVER_PATH=/my/custom/path/cable_solver
```

### PyVista / VTK import hang on launch

VTK 9.6+ has a known import hang under Python 3.12 on macOS. This package
pins `vtk==9.3.1`. Recreate the venv if you get an old environment:

```bash
rm -rf .venv
./run.sh
```

### PySide6 splash / event filter crash

PySide6 6.10.x has a `QSplashScreen.show()` hang and an `eventFilter`
recursion bug. Pinned below 6.10; force the pin if pip picks newer:

```bash
pip install 'PySide6<6.10'
```

### Linker can't find `lapack` / `openmp`

```bash
export PATH="/opt/homebrew/bin:$PATH"
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
```

Then clean and rebuild the C++ tree.

## License

LGPL-3.0-or-later. Same as the enclosing CableDynamics / BEM repositories.
