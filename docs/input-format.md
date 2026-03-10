---
layout: default
title: Input Format
---

[🇯🇵 日本語](ja/input-format.html)

# Input File Format

All input files use JSON format. A simulation requires a `settings.json` file and one or more body definition files.

## Directory Structure

```
my_case/
├── settings.json      # Main simulation settings
├── water.json         # Free-surface fluid definition
├── tank.json          # Wave tank walls
├── wavemaker.json     # Wavemaker boundary (optional)
├── absorber.json      # Wave absorber (optional)
└── float.json         # Floating body (optional)
```

## settings.json

```json
{
  "output_directory": "./output",
  "input_files": ["water.json", "tank.json", "wavemaker.json"],
  "dt": 0.05,
  "end_time_step": 1000,
  "ALE": "pseudo_quad",
  "ALE_period": 1
}
```

### Key Parameters

| Key | Type | Description |
|-----|------|-------------|
| `output_directory` | string | Path for output files (relative to input directory) |
| `input_files` | array | List of body definition files to load |
| `dt` | number | Time step size |
| `end_time_step` | integer | Maximum number of time steps |
| `ALE` | string | ALE method: `"linear"`, `"pseudo_quad"` |
| `ALE_period` | integer | Apply ALE every N steps |

## Body Definition Files

Each body file defines a surface mesh and its boundary conditions.

```json
{
  "name": "water",
  "type": "Fluid",
  "objfile": ["../../obj/Goring1979/water.obj"]
}
```

### Required Keys

| Key | Type | Description |
|-----|------|-------------|
| `name` | string | Unique identifier for this body |
| `type` | string | Body type (see below) |
| `objfile` | array | Path to OBJ mesh file |

### Body Types

| Type | Description |
|------|-------------|
| `Fluid` | Free-surface fluid domain |
| `RigidBody` | Rigid body with 6-DOF dynamics |
| `SoftBody` | Deformable body |
| `Absorber` | Wave-absorbing boundary |

### Optional Keys

| Key | Type | Description |
|-----|------|-------------|
| `ignore` | bool | Skip this body if `true` |
| `velocity` | array | Initial velocity `[vx, vy, vz]` |
| `angular_velocity` | array | Initial angular velocity |
| `center_of_mass` | array | Center of mass position |
| `mass` | number | Body mass |
| `MOI` | array | Moments of inertia |
