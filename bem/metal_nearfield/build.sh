#!/bin/bash
# Build script for metal_nearfield library
#
# This library must be built separately with Apple Clang.
# The resulting dylib can then be linked with the GCC-compiled main program.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"

echo "=== Building metal_nearfield library ==="
echo "Source directory: ${SCRIPT_DIR}"
echo "Build directory: ${BUILD_DIR}"

# Create build directory
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
make -j$(sysctl -n hw.logicalcpu)

# Copy results to parent directory for easy access
cp libmetal_nearfield.dylib "${SCRIPT_DIR}/"
cp metal_nearfield.metallib "${SCRIPT_DIR}/" 2>/dev/null || true

echo ""
echo "=== Build complete ==="
echo "Library: ${SCRIPT_DIR}/libmetal_nearfield.dylib"
echo ""
echo "To use in the main BEM program:"
echo "  1. Add to CMakeLists.txt:"
echo "     target_link_libraries(\${BASE_NAME} PRIVATE \"\${CMAKE_CURRENT_SOURCE_DIR}/metal_nearfield/libmetal_nearfield.dylib\")"
echo "     target_include_directories(\${BASE_NAME} PRIVATE \"\${CMAKE_CURRENT_SOURCE_DIR}/metal_nearfield\")"
echo ""
echo "  2. Include in your code:"
echo "     #include \"metal_nearfield_wrapper.hpp\""
echo ""
