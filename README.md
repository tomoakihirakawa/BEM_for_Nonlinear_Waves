# BEM for Nonlinear Waves

A boundary element method (BEM) solver for nonlinear free-surface wave problems, with applications to wave-body interaction, coastal engineering, and floating offshore structures.

![GUI Screenshot](docs/images/gui_screenshot.png)

## Features

- **Fully nonlinear BEM** with mixed Eulerian-Lagrangian formulation
- **Fast Multipole Method (FMM)** for O(N) boundary integral evaluation
- **Arbitrary Lagrangian-Eulerian (ALE)** mesh motion for free-surface tracking
- **Time-domain and frequency-domain** solvers
- **Mooring line dynamics** (cable module)
- **Metal GPU acceleration** for M2L transformations (Apple Silicon, optional)
- **OpenMP parallelization**
- **Quad-precision arithmetic** (__float128) for FMM translation operators

## Repository Structure

```
BEM_for_Nonlinear_Waves/
├── lib/                  # Shared library (mesh, geometry, FMM, ODE, etc.)
│   ├── include/          # Header files
│   └── src/              # Implementation files
├── bem/                  # BEM solver (time-domain & frequency-domain)
│   ├── main_time_domain.cpp
│   ├── main_freq_domain.cpp
│   ├── BEM_*.hpp         # BEM modules
│   └── metal_m2l/        # Metal GPU acceleration (optional)
├── cable/                # Cable/mooring dynamics module
├── third_party/tetgen/   # TetGen (AGPL-3.0, header + pre-built library)
├── obj/                  # Benchmark mesh files
└── CMakeLists.txt        # Build system
```

## Quick Start

### Prerequisites

- **C++ compiler**: GCC 12+ (recommended) or Clang 15+
- **CMake** 3.16+
- **LAPACK/BLAS**
- **OpenMP** (optional but recommended)

On macOS with Homebrew:
```bash
brew install gcc cmake lapack
```

### Build

```bash
git clone https://github.com/tomoakihirakawa/BEM_for_Nonlinear_Waves.git
cd BEM_for_Nonlinear_Waves

# Time-domain solver
mkdir -p bem/build && cd bem/build
cmake -DSOURCE_FILE=../main_time_domain.cpp ../..
cmake --build . -j$(sysctl -n hw.logicalcpu)

# Frequency-domain solver (from project root)
mkdir -p build-freq && cd build-freq
cmake -DSOURCE_FILE=bem/main_freq_domain.cpp ..
cmake --build . -j$(sysctl -n hw.logicalcpu)
```

### Run

```bash
./main_time_domain /path/to/input_files/
```

Input files are JSON-based. See `bem/input_files/` for examples.

## CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `BEM_COMPILER` | `gcc` | Compiler toolchain (`gcc` or `clang`) |
| `BEM_ENABLE_OPENMP` | `ON` | Enable OpenMP parallelization |
| `BEM_DISABLE_FLOAT128` | `OFF` | Disable quad-precision (for Clang/libc++) |
| `FMM_M2L_METHOD` | `SimpleM2L` | FMM M2L translation method |
| `USE_TETGEN` | `ON` | Enable TetGen tetrahedralization |
| `USE_METAL_M2L` | `OFF` | Enable Metal GPU for M2L (Apple Silicon) |

## TetGen Notice

This project optionally uses [TetGen](https://wias-berlin.de/software/tetgen/) for tetrahedral mesh generation. TetGen is licensed under **AGPL-3.0**. Only the header file and a pre-built static library are included; the TetGen source code is not distributed with this project. Set `USE_TETGEN=OFF` to build without TetGen.

TetGen will be replaced by a custom advancing-front method in future releases.

## License

This project is licensed under the [LGPL-3.0-or-later](LICENSE).

Third-party components have their own licenses:
- TetGen: AGPL-3.0 (see `third_party/tetgen/LICENSE`)
- ankerl::unordered_dense: MIT
- pdqsort: zlib

## Citation

If you use this software in your research, please cite:

```bibtex
@software{hirakawa2025bem,
  author = {Hirakawa, Tomoaki},
  title = {BEM for Nonlinear Waves},
  year = {2025},
  url = {https://github.com/tomoakihirakawa/BEM_for_Nonlinear_Waves}
}
```
