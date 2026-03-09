# ==========================================================================
# CMake Project: BEM for Nonlinear Waves
#
# Supports GCC (g++) and Clang on macOS / Linux.
# Features: quad-precision (__float128), OpenMP, Metal GPU (Apple Silicon).
# ==========================================================================

cmake_minimum_required(VERSION 3.16)

# ==========================================================================
# Compiler selection
#
# NOTE: CMake caches the compiler chosen at first configure.
#       To switch compilers, use a fresh build directory.
#       If -DCMAKE_CXX_COMPILER=... is given, it takes precedence.
# ==========================================================================
set(BEM_COMPILER "gcc" CACHE STRING "Select compiler toolchain: gcc or clang")
set_property(CACHE BEM_COMPILER PROPERTY STRINGS "gcc" "clang")
option(BEM_ENABLE_OPENMP "Enable OpenMP parallelization" ON)

set(_BEM_DISABLE_FLOAT128_DEFAULT OFF)
if(BEM_COMPILER STREQUAL "clang")
  set(_BEM_DISABLE_FLOAT128_DEFAULT ON)
endif()
option(BEM_DISABLE_FLOAT128 "Disable std::float128_t / quadmath codepaths (for AppleClang/libc++)" ${_BEM_DISABLE_FLOAT128_DEFAULT})

# ==========================================================================
# FMM M2L Method Selection
# ==========================================================================
set(FMM_M2L_METHOD "SimpleM2L" CACHE STRING "FMM M2L translation method")
set_property(CACHE FMM_M2L_METHOD PROPERTY STRINGS
    "SimpleM2L"
    "FourierM2L_double"
    "FourierM2L_DoubleDouble"
    "FourierM2L_BlockDecomposition"
    "PlaneWaveM2L"
)

# ==========================================================================
# Metal GPU acceleration for M2L (Apple Silicon only)
# ==========================================================================
option(USE_METAL_M2L "Use Metal GPU acceleration for M2L transformation" OFF)

# ==========================================================================
# TetGen (optional, AGPL-3.0 licensed)
# ==========================================================================
option(USE_TETGEN "Use TetGen for tetrahedralization (AGPL-3.0 licensed)" ON)

if(NOT DEFINED CMAKE_CXX_COMPILER OR CMAKE_CXX_COMPILER STREQUAL "")
  if(BEM_COMPILER STREQUAL "gcc")
    # Prefer Homebrew g++-15/14/13/12
    set(_BEM_GXX_VERSIONS 15 14 13 12)
    foreach(ver IN LISTS _BEM_GXX_VERSIONS)
      find_program(_BEM_GXX NAMES "g++-${ver}")
      if(_BEM_GXX)
        message(STATUS "Found g++-${ver}: ${_BEM_GXX}")
        set(CMAKE_CXX_COMPILER "${_BEM_GXX}" CACHE FILEPATH "C++ compiler" FORCE)
        break()
      endif()
    endforeach()
    if(NOT DEFINED CMAKE_CXX_COMPILER OR CMAKE_CXX_COMPILER STREQUAL "")
      message(FATAL_ERROR "No suitable g++ version found (tried g++-15/14/13/12).")
    endif()
  elseif(BEM_COMPILER STREQUAL "clang")
    # Prefer Homebrew LLVM clang++ (needed for OpenMP on macOS)
    set(_BEM_LLVM_HINTS "/opt/homebrew/opt/llvm/bin" "/usr/local/opt/llvm/bin")
    find_program(_BEM_CLANGXX NAMES "clang++" PATHS ${_BEM_LLVM_HINTS} NO_DEFAULT_PATH)
    if(NOT _BEM_CLANGXX)
      find_program(_BEM_CLANGXX NAMES "clang++")
    endif()
    if(NOT _BEM_CLANGXX)
      message(FATAL_ERROR "clang++ not found. Install LLVM (e.g. Homebrew) or pass -DCMAKE_CXX_COMPILER=/path/to/clang++")
    endif()
    message(STATUS "Found clang++: ${_BEM_CLANGXX}")
    set(CMAKE_CXX_COMPILER "${_BEM_CLANGXX}" CACHE FILEPATH "C++ compiler" FORCE)
  else()
    message(FATAL_ERROR "BEM_COMPILER must be 'gcc' or 'clang' (got: ${BEM_COMPILER})")
  endif()
endif()

project(BEM_Project CXX)

# ==========================================================================
# Build type (default: Release)
# ==========================================================================
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type selected, defaulting to Release")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Terminal color codes
if(NOT WIN32 AND CMAKE_COLOR_MAKEFILE)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(Green "${Esc}[32m")
  set(Magenta "${Esc}[35m")
  set(Yellow "${Esc}[33m")
  set(Red "${Esc}[31m")
  set(BoldCyan "${Esc}[1;36m")
endif()

# ==========================================================================
# Project-wide settings
# ==========================================================================
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_UNITY_BUILD OFF)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ccache
option(CCACHE_ENABLE "Use ccache if available" ON)
if(CCACHE_ENABLE)
    find_program(CCACHE_EXE ccache)
    if(CCACHE_EXE)
        set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_EXE}")
        message(STATUS "ccache enabled")
    endif()
endif()

# ==========================================================================
# Executable definition
# ==========================================================================
set(DEFAULT_SOURCE_FILE "main.cpp")

if(CMAKE_ARGS)
  list(FIND CMAKE_ARGS "-DSOURCE_FILE" SOURCE_FILE_IDX)
  if(NOT SOURCE_FILE_IDX EQUAL -1)
    list(GET CMAKE_ARGS ${SOURCE_FILE_IDX} SOURCE_FILE_ARG)
    string(REPLACE "-DSOURCE_FILE=" "" SOURCE_FILE ${SOURCE_FILE_ARG})
  endif()
endif()

if(NOT DEFINED SOURCE_FILE)
  set(SOURCE_FILE ${DEFAULT_SOURCE_FILE})
endif()

if(DEFINED OUTPUT_NAME)
  set(BASE_NAME ${OUTPUT_NAME})
else()
  get_filename_component(BASE_NAME ${SOURCE_FILE} NAME_WE)
endif()

# ==========================================================================
# Source files (split from header-only for parallel compilation)
# ==========================================================================
set(EXTRA_SOURCES
    "${PROJECT_SOURCE_DIR}/lib/src/Network.cpp"
    "${PROJECT_SOURCE_DIR}/lib/src/networkPoint.cpp"
    "${PROJECT_SOURCE_DIR}/lib/src/networkLine.cpp"
    "${PROJECT_SOURCE_DIR}/lib/src/networkFace.cpp"
    "${PROJECT_SOURCE_DIR}/lib/src/networkTetra.cpp"
)

set(VALID_SOURCES ${SOURCE_FILE})
foreach(SRC ${EXTRA_SOURCES})
    if(EXISTS "${SRC}")
        list(APPEND VALID_SOURCES "${SRC}")
        message(STATUS "Added source: ${SRC}")
    else()
        message(WARNING "Source file not found: ${SRC}")
    endif()
endforeach()

add_executable(${BASE_NAME} ${VALID_SOURCES})

# ==========================================================================
# Target properties
# ==========================================================================

# --- Include directories ---
target_include_directories(${BASE_NAME} PUBLIC
  "${PROJECT_SOURCE_DIR}/lib/include"
)

# --- Compile options ---
target_compile_options(${BASE_NAME} PRIVATE
    -march=native
    -funroll-loops
    $<$<CXX_COMPILER_ID:GNU>:-fconcepts>
    -pthread
    $<$<CONFIG:DEBUG>:-Og>
    $<$<CONFIG:RELEASE>:-O3>
)

# --- OpenMP ---
if(BEM_ENABLE_OPENMP)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${BASE_NAME} PRIVATE -fopenmp)
    target_link_libraries(${BASE_NAME} PRIVATE gomp)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
      target_link_libraries(${BASE_NAME} PRIVATE OpenMP::OpenMP_CXX)
    else()
      find_path(BEM_OMP_INCLUDE_DIR omp.h PATHS "/opt/homebrew/opt/libomp/include" "/usr/local/opt/libomp/include")
      find_library(BEM_OMP_LIBRARY NAMES omp PATHS "/opt/homebrew/opt/libomp/lib" "/usr/local/opt/libomp/lib")
      if(BEM_OMP_INCLUDE_DIR AND BEM_OMP_LIBRARY)
        target_compile_options(${BASE_NAME} PRIVATE -Xpreprocessor -fopenmp)
        target_include_directories(${BASE_NAME} PRIVATE "${BEM_OMP_INCLUDE_DIR}")
        target_link_libraries(${BASE_NAME} PRIVATE "${BEM_OMP_LIBRARY}")
      else()
        message(FATAL_ERROR "OpenMP requested but not found for Clang. Install libomp or use Homebrew LLVM.")
      endif()
    endif()
  endif()
else()
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(${BASE_NAME} PRIVATE -fopenmp-simd)
  endif()
endif()

# --- Libraries ---
target_link_libraries(${BASE_NAME} PRIVATE
    $<$<CXX_COMPILER_ID:GNU>:quadmath>
    lapack
    blas
)

# --- TetGen (conditional) ---
if(USE_TETGEN)
  set(LIBTET_PATH "${PROJECT_SOURCE_DIR}/third_party/tetgen/libtet.a")
  target_compile_definitions(${BASE_NAME} PRIVATE USE_TETGEN)
  target_include_directories(${BASE_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/third_party/tetgen")
  target_link_libraries(${BASE_NAME} PRIVATE "${LIBTET_PATH}")
endif()

# ==========================================================================
# Metal GPU Nearfield Integration (optional)
# ==========================================================================
option(USE_METAL_NEARFIELD "Use Metal GPU acceleration for nearfield direct integration" OFF)

if(USE_METAL_NEARFIELD)
  set(METAL_NF_DIR "${PROJECT_SOURCE_DIR}/bem/metal_nearfield")
  set(METAL_NF_LIB "${METAL_NF_DIR}/build/libmetal_nearfield.dylib")

  if(EXISTS "${METAL_NF_LIB}")
    target_compile_definitions(${BASE_NAME} PRIVATE USE_METAL_NEARFIELD)
    target_include_directories(${BASE_NAME} PRIVATE "${METAL_NF_DIR}")
    target_link_libraries(${BASE_NAME} PRIVATE "${METAL_NF_LIB}")

    add_custom_command(TARGET ${BASE_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${METAL_NF_LIB}"
        "$<TARGET_FILE_DIR:${BASE_NAME}>/libmetal_nearfield.dylib"
      COMMENT "Copying Metal Nearfield library"
    )
    if(EXISTS "${METAL_NF_DIR}/build/metal_nearfield.metallib")
      add_custom_command(TARGET ${BASE_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${METAL_NF_DIR}/build/metal_nearfield.metallib"
          "$<TARGET_FILE_DIR:${BASE_NAME}>/metal_nearfield.metallib"
        COMMENT "Copying Metal Nearfield shader library"
      )
    endif()
  else()
    message(WARNING "Metal Nearfield library not found at ${METAL_NF_LIB}. Build it first with: cd bem/metal_nearfield && mkdir build && cd build && cmake .. && make")
    set(USE_METAL_NEARFIELD OFF)
  endif()
endif()

# ==========================================================================
# Metal GPU M2L Integration (optional)
# ==========================================================================
if(USE_METAL_M2L)
  set(METAL_M2L_DIR "${PROJECT_SOURCE_DIR}/bem/metal_m2l")
  set(METAL_M2L_LIB "${METAL_M2L_DIR}/build/libmetal_m2l.dylib")

  if(EXISTS "${METAL_M2L_LIB}")
    target_compile_definitions(${BASE_NAME} PRIVATE USE_METAL_M2L)
    target_include_directories(${BASE_NAME} PRIVATE "${METAL_M2L_DIR}")
    target_link_libraries(${BASE_NAME} PRIVATE "${METAL_M2L_LIB}")

    add_custom_command(TARGET ${BASE_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${METAL_M2L_LIB}"
        "$<TARGET_FILE_DIR:${BASE_NAME}>/libmetal_m2l.dylib"
      COMMENT "Copying Metal M2L library"
    )
    if(EXISTS "${METAL_M2L_DIR}/build/metal_m2l.metallib")
      add_custom_command(TARGET ${BASE_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${METAL_M2L_DIR}/build/metal_m2l.metallib"
          "$<TARGET_FILE_DIR:${BASE_NAME}>/metal_m2l.metallib"
        COMMENT "Copying Metal M2L shader library"
      )
    endif()
  else()
    message(WARNING "Metal M2L library not found at ${METAL_M2L_LIB}. Build it first with: cd bem/metal_m2l && mkdir build && cd build && cmake .. && make")
    set(USE_METAL_M2L OFF)
  endif()
endif()

# ==========================================================================
# Configuration summary
# ==========================================================================
string(REPEAT "=" 70 DIVIDER)
message(STATUS "")
message(STATUS "${Yellow}${DIVIDER}${ColourReset}")
message(STATUS "${BoldCyan}  BEM for Nonlinear Waves - Configuration Summary${ColourReset}")
message(STATUS "${Yellow}${DIVIDER}${ColourReset}")
message(STATUS "${Green}Target Name       : ${Magenta}${BASE_NAME}${ColourReset}")
message(STATUS "${Green}BEM_COMPILER      : ${Magenta}${BEM_COMPILER}${ColourReset}")
message(STATUS "${Green}BEM_ENABLE_OPENMP : ${Magenta}${BEM_ENABLE_OPENMP}${ColourReset}")
message(STATUS "${Green}BEM_DISABLE_FLOAT128 : ${Magenta}${BEM_DISABLE_FLOAT128}${ColourReset}")
message(STATUS "${Green}FMM_M2L_METHOD    : ${Magenta}${FMM_M2L_METHOD}${ColourReset}")
message(STATUS "${Green}USE_TETGEN        : ${Magenta}${USE_TETGEN}${ColourReset}")
message(STATUS "${Green}USE_METAL_M2L     : ${Magenta}${USE_METAL_M2L}${ColourReset}")
message(STATUS "${Green}USE_METAL_NEARFIELD : ${Magenta}${USE_METAL_NEARFIELD}${ColourReset}")
message(STATUS "${Green}C++ Compiler      : ${Magenta}${CMAKE_CXX_COMPILER}${ColourReset}")
message(STATUS "${Green}C++ Standard      : ${Magenta}gnu++${CMAKE_CXX_STANDARD}${ColourReset}")
message(STATUS "${Green}Build Type        : ${Magenta}${CMAKE_BUILD_TYPE}${ColourReset}")
if(BEM_COMPILER STREQUAL "gcc" AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  message(STATUS "${Yellow}NOTE: BEM_COMPILER=gcc but compiler ID is '${CMAKE_CXX_COMPILER_ID}'. Use a fresh build dir to switch.${ColourReset}")
elseif(BEM_COMPILER STREQUAL "clang" AND NOT (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
  message(STATUS "${Yellow}NOTE: BEM_COMPILER=clang but compiler ID is '${CMAKE_CXX_COMPILER_ID}'. Use a fresh build dir to switch.${ColourReset}")
endif()
message(STATUS "${Yellow}To build: ${BoldCyan}cmake --build . -j$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)${ColourReset}")
message(STATUS "${Yellow}${DIVIDER}${ColourReset}")

if(BEM_DISABLE_FLOAT128)
  target_compile_definitions(${BASE_NAME} PRIVATE BEM_DISABLE_FLOAT128=1)
endif()

# FMM M2L method compile definition
target_compile_definitions(${BASE_NAME} PRIVATE USE_${FMM_M2L_METHOD})
