---
layout: default
title: Home
---

[日本語](ja/)

# BEM for Nonlinear Waves

A boundary element method (BEM) solver for nonlinear free-surface wave problems.

![GUI Screenshot](images/gui_screenshot.png)

## Overview

This software solves potential flow problems with nonlinear free-surface boundary conditions using the boundary element method. It implements the Mixed Eulerian-Lagrangian (MEL) approach for time-domain simulations and supports frequency-domain analysis.

### Key Capabilities

- **Nonlinear wave generation** -- piston/flap wavemakers, solitary waves, irregular waves
- **Wave-body interaction** -- floating bodies with 6-DOF rigid body dynamics
- **Fast Multipole Method** -- O(N) evaluation of boundary integrals
- **ALE mesh management** -- adaptive remeshing with Laplacian smoothing
- **Mooring line dynamics** -- lumped-mass cable model
- **Metal GPU acceleration** -- M2L transformations on Apple Silicon

## Documentation

### Getting Started

- [Getting Started](getting-started.html) -- Build and run your first simulation

### Theory

- [Theory Overview](theory/) -- Mathematical formulation
  - [Boundary Integral Equation](theory/boundary-integral.html) -- BIE derivation, linear elements, coefficient matrices
  - [Floating Body Dynamics](theory/floating-body.html) -- 6-DOF motion, pressure computation, $\phi_t$ problem
  - [Wave Generation](theory/wave-generation.html) -- Piston/flap wavemakers, solitary waves
  - [ALE Mesh](theory/ale-mesh.html) -- Arbitrary Lagrangian-Eulerian mesh management
  - [Fast Multipole Method](theory/fmm.html) -- FMM acceleration, M2L methods

### Reference

- [Input Format](input-format.html) -- JSON input file reference

### Examples

- [Goring (1979)](examples/goring1979.html) -- Solitary wave generation and propagation
- [DeepCWind](examples/deepcwind.html) -- Floating wind turbine platform

## License

LGPL-3.0-or-later. See [LICENSE](https://github.com/tomoakihirakawa/BEM_for_Nonlinear_Waves/blob/main/LICENSE).
