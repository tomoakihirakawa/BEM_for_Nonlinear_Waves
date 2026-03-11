---
layout: default
title: Getting Started
---

[日本語](ja/getting-started.html)

# Getting Started

## Prerequisites

| Requirement | Minimum Version | Notes |
|------------|----------------|-------|
| C++ compiler | GCC 12+ or Clang 15+ | GCC recommended for `__float128` support |
| CMake | 3.16+ | |
| LAPACK/BLAS | — | Required for linear algebra |
| OpenMP | — | Optional but recommended |

### macOS (Homebrew)

```bash
brew install gcc cmake lapack
```

### Ubuntu/Debian

```bash
sudo apt install g++ cmake liblapack-dev libblas-dev libgomp1
```

## Building

### Clone the repository

```bash
git clone https://github.com/tomoakihirakawa/BEM_for_Nonlinear_Waves.git
cd BEM_for_Nonlinear_Waves
```

### Time-domain solver

```bash
mkdir -p bem/build && cd bem/build
cmake -DSOURCE_FILE=../main_time_domain.cpp ../..
cmake --build . -j$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)
```

### Frequency-domain solver

```bash
mkdir -p build-freq && cd build-freq
cmake -DSOURCE_FILE=bem/main_freq_domain.cpp ..
cmake --build . -j$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)
```

## Running a Simulation

```bash
./main_time_domain /path/to/input_files/
```

The input directory should contain:
- `settings.json` — simulation parameters (time step, duration, etc.)
- One or more body/fluid definition files (e.g., `water.json`, `tank.json`)

Each body file specifies:
- `name` — identifier
- `type` — `Fluid`, `RigidBody`, `SoftBody`, etc.
- `objfile` — path to the surface mesh (OBJ format)

## CMake Options

```bash
cmake -DSOURCE_FILE=../main_time_domain.cpp \
      -DBEM_COMPILER=gcc \
      -DUSE_TETGEN=ON \
      -DFMM_M2L_METHOD=SimpleM2L \
      ../..
```

See the [README](https://github.com/tomoakihirakawa/BEM_for_Nonlinear_Waves#cmake-options) for all available options.

## Building without TetGen

TetGen is licensed under AGPL-3.0. To build without it:

```bash
cmake -DUSE_TETGEN=OFF -DSOURCE_FILE=../main_time_domain.cpp ../..
```

Tetrahedralization features will be disabled, but the BEM solver remains fully functional for surface-only computations.
