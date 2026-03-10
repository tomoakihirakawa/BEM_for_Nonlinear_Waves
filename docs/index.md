---
layout: default
title: Home
---

[🇯🇵 日本語](ja/)

# BEM for Nonlinear Waves

A boundary element method (BEM) solver for nonlinear free-surface wave problems.

![GUI Screenshot](images/gui_screenshot.png)

## Overview

This software solves potential flow problems with nonlinear free-surface boundary conditions using the boundary element method. It implements the Mixed Eulerian-Lagrangian (MEL) approach for time-domain simulations and supports frequency-domain analysis.

### Key Capabilities

- **Nonlinear wave generation** — piston/flap wavemakers, solitary waves, irregular waves
- **Wave-body interaction** — floating bodies with 6-DOF rigid body dynamics
- **Fast Multipole Method** — O(N) evaluation of boundary integrals
- **ALE mesh management** — adaptive remeshing with Laplacian smoothing
- **Mooring line dynamics** — lumped-mass cable model
- **Metal GPU acceleration** — M2L transformations on Apple Silicon

## Getting Started

See the [Getting Started](getting-started.html) guide for build instructions.

## Documentation

- [Getting Started](getting-started.html) — Build and run your first simulation
- [Theory](theory.html) — Mathematical formulation
- [Input Format](input-format.html) — JSON input file reference

## Examples

- [Goring (1979)](examples/goring1979.html) — Solitary wave generation and propagation

## License

LGPL-3.0-or-later. See [LICENSE](https://github.com/tomoakihirakawa/BEM_for_Nonlinear_Waves/blob/main/LICENSE).
