---
layout: default
title: GUI User Guide
---

[日本語](ja/gui.html)

# GUI User Guide

The GUI application provides a visual interface for creating, editing, and exporting simulation input files without manually writing JSON.

![GUI Screenshot](images/gui_screenshot.png)

## Requirements

| Requirement | Version | Notes |
|------------|---------|-------|
| Python | 3.11.x | 3.12+ not supported (VTK compatibility) |
| PySide6 | 6.9.3 | 6.10+ causes hang on macOS |
| macOS | 12+ | Apple Silicon recommended |

## Installation and Launch

```bash
cd gui
sh run.sh
```

The `run.sh` script automatically creates a Python virtual environment, installs dependencies, and launches the application.

## Workflow: From Case Selection to Simulation

The basic workflow is three steps: **Gallery** → **Settings** → **Export**.

### Step 1 — Select a Case from the Gallery

The left panel starts with the **Gallery** tab, which shows available simulation cases as visual cards.

1. Browse the cards and find a case closest to your target simulation (e.g., solitary wave in a tank, floating body in waves)
2. Click the card to load it

The selected case's configuration is automatically loaded into the **Settings** panel.

### Step 2 — Edit Parameters in Settings

After selecting a case, switch to the **Settings** tab (next to Gallery in the left panel). This panel contains all simulation parameters:

- **Simulation settings** — time step (`dt`), end time, solver options
- **ALE settings** — mesh relocation method, smoothing period
- **Element type** — linear or quadratic elements
- **Object Tree** — lists all objects (water, tank, wavemaker, floating bodies) with visibility toggles

You can switch between **Form** view (structured editor) and **JSON** view (raw text) at the top of the panel.

To edit individual object properties (mesh file, position, mass, boundary conditions), select an object in the Object Tree. Its details appear in the **Inspector** panel on the right side, with both Form and JSON editing modes.

### Step 3 — Export Input Files

Once settings are configured, move to the **Export** section:

1. Set the output directory and file name
2. Click **Export** to generate the JSON input files (`settings.json`, `water.json`, `tank.json`, etc.)

### Step 4 — Run the Solver

Run the time-domain solver from the terminal, passing the output directory:

```bash
./main_time_domain /path/to/exported/output_directory/
```

See [Getting Started](getting-started.html) for build instructions.

## Window Layout

### Center — 3D Viewport

- Interactive 3D visualization of the mesh geometry
- Hover over elements to see object name and type
- **Floating Toolbar** (bottom center) — zoom, reset camera, axis controls
- **BC Panel** (top right) — toggle boundary condition display

### Bottom Panel

| Tab | Description |
|-----|-------------|
| **Console** | Command execution terminal |
| **Logs** | Timestamped log messages with severity levels |
| **Output Files** | Browse generated output files |
| **History** | Export job history |
| **Checkpoints** | Saved simulation snapshots |
| **Branches** | Experiment branching and comparison |

### Status Bar

- Progress indicator (visible during export)
- Current case name, object count, solver type
- **Toggle Theme** (Ctrl+T) — switch between dark and light themes
- **EN/JP** (Ctrl+L) — switch interface language
