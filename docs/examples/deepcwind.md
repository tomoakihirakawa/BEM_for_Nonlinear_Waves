---
layout: default
title: "Example: DeepCWind Floating Platform"
---

[Japanese / 日本語](../ja/examples/deepcwind.html)

# DeepCWind Floating Platform

This example demonstrates a wave-body interaction simulation of the OC4
DeepCWind semi-submersible floating platform. The floating body undergoes
six-degree-of-freedom (6-DOF) rigid body motion coupled with the nonlinear
free-surface flow.

## Reference

Robertson, A., Jonkman, J., Masciola, M., Song, H., Goupee, A., Coulling, A.,
and Luan, C. (2014). *Definition of the Semisubmersible Floating System for
Phase II of OC4*. National Renewable Energy Laboratory, Technical Report
NREL/TP-5000-60601.

## Description

The DeepCWind semi-submersible is a three-column floating platform designed to
support offshore wind turbines. The platform consists of:

- Three offset columns arranged in a triangular pattern
- A central column (tower base)
- Pontoons and cross braces connecting the columns
- Heave plates at the base of each offset column

In this BEM simulation, the platform is modeled as a single rigid body floating
in a wave tank. The solver computes the fully nonlinear wave-body interaction
including:

- Nonlinear free-surface boundary conditions
- 6-DOF rigid body dynamics (surge, sway, heave, roll, pitch, yaw)
- Hydrostatic and hydrodynamic pressure integration on the body surface
- Wave generation and absorption at the tank boundaries

## Input File Structure

The input files for DeepCWind cases are located in `bem/input_files/` under
directories named `DeepCWind_*`. A typical case includes five input files.

### settings.json

```json
{
    "max_dt": 0.5,
    "end_time_step": 10000000,
    "end_time": 100,
    "element": "linear",
    "ALE": "pseudo_quad",
    "ALEPERIOD": "1",
    "output_directory": "./output",
    "input_files": [
        "tank.json",
        "water.json",
        "float.json",
        "absorber.json"
    ],
    "meshing_options": [
        "surface_flip"
    ]
}
```

| Parameter | Description |
|-----------|-------------|
| `max_dt` | Maximum time step size (seconds). For floating body problems, 0.5 s is typical. |
| `end_time` | Simulation end time. Set long enough for the body to reach steady oscillation. |
| `ALE` | `pseudo_quad` is recommended for floating body simulations. |

### water.json

Defines the fluid domain mesh. The mesh must conform to the floating body
surface at the initial waterline.

```json
{
    "name": "water",
    "type": "Fluid",
    "objfile": "path/to/water.obj"
}
```

### tank.json

Defines the rigid tank walls (bottom and far-field boundaries).

```json
{
    "name": "tank",
    "type": "RigidBody",
    "isFixed": true,
    "objfile": "path/to/tank.obj"
}
```

### float.json

Defines the floating body with its mass properties and 6-DOF dynamics.

```json
{
    "name": "float",
    "type": "RigidBody",
    "velocity": "floating",
    "mass": 13556760.999999998,
    "COM": [0, 0, 171.93],
    "MOI": [13947000000.0, 15552000000.0, 13692000000.0],
    "objfile": "path/to/float.obj"
}
```

| Parameter | Description |
|-----------|-------------|
| `velocity` | `"floating"` enables 6-DOF free-floating dynamics |
| `mass` | Total mass of the floating body (kg) |
| `COM` | Center of mass coordinates [x, y, z] (meters) |
| `MOI` | Moments of inertia [Ixx, Iyy, Izz] (kg m^2) |
| `objfile` | Surface mesh of the floating body |

The mass, center of mass, and moments of inertia must be consistent with the
platform design. The mesh file defines the wetted surface geometry.

### absorber.json

Defines the wave absorbing boundary to prevent reflection from the tank walls.

```json
{
    "name": "absorber",
    "type": "Absorber",
    "isFixed": true,
    "objfile": "path/to/absorber.obj",
    "wave_theory_L": [3.7, 150.0, 180.0, 0, 180.0]
}
```

| Parameter | Description |
|-----------|-------------|
| `type` | `"Absorber"` activates the wave absorption algorithm |
| `wave_theory_L` | Wave absorption parameters: [wavelength, ...geometry parameters] |

## Key Simulation Parameters

When setting up a DeepCWind simulation, pay attention to the following:

- **Time step**: Use `max_dt` of 0.25--0.5 s. Smaller values improve accuracy
  but increase computation time.
- **ALE scheme**: `pseudo_quad` provides better mesh quality for large body
  motions than `linear`.
- **Mesh resolution**: The fluid mesh must be fine enough near the waterline to
  resolve the body geometry and wave-body interaction.
- **End time**: Allow sufficient time for transients to decay and steady-state
  oscillation to develop (typically 50--100 s or more).

## Build and Run

```bash
cd bem/build
cmake -DSOURCE_FILE=../main_time_domain.cpp ../..
make -j$(sysctl -n hw.logicalcpu)
./main ../../bem/input_files/DeepCWind_DT0d5_ELEMlinear_ALEpseudo_quad_ALEPERIOD1/
```

Replace the input directory path with the specific DeepCWind case you wish to
run.

## Expected Results

- The floating body settles to its hydrostatic equilibrium position
- Under wave loading, the platform exhibits 6-DOF oscillatory motion
- Heave, pitch, and surge are the dominant response modes
- Hydrodynamic forces on the body surface can be extracted from the output
- The wave absorber effectively damps outgoing waves at the tank boundaries

## See Also

- [Floating Body Theory](../theory/floating-body.html)
