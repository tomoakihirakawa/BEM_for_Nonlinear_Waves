# TetGen

This directory contains a pre-built TetGen library for macOS (arm64).

**TetGen** is a program for generating tetrahedral meshes. It is developed by
Hang Si at the Weierstrass Institute for Applied Analysis and Stochastics (WIAS),
Berlin.

- Website: https://wias-berlin.de/software/tetgen/
- License: **AGPL-3.0** (see `LICENSE` file)

## Contents

- `tetgen.h` — TetGen C++ header (API declarations)
- `libtet.a` — Pre-built static library (macOS arm64, Apple Silicon)
- `LICENSE` — AGPL-3.0 license text

## Building for Other Platforms

The pre-built `libtet.a` is for macOS arm64 only. To build for other platforms:

1. Download TetGen 1.6.0 source from https://wias-berlin.de/software/tetgen/
2. Build:
   ```bash
   g++ -O3 -DTETLIBRARY -c tetgen.cxx predicates.cxx
   ar rcs libtet.a tetgen.o predicates.o
   ```
3. Replace `libtet.a` in this directory with your build.

## Note

The TetGen **source code** (`tetgen.cxx`, `predicates.cxx`) is **not included**
in this repository to avoid AGPL license propagation. Only the header file and
pre-built library are distributed.

This dependency is optional (`-DUSE_TETGEN=OFF`) and will be replaced by a
custom advancing-front mesher in future releases.
