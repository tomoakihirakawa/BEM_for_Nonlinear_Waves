# ==========================================================================
# CMake Project: wave_spectrum Python module
#
# pybind11 bindings for RandomWaterWaveTheory (JONSWAP, Bretschneider)
# ==========================================================================

cmake_minimum_required(VERSION 3.16)

# ==========================================================================
# Compiler settings — gcc for C++23 / GNU extensions
# ==========================================================================
if(NOT DEFINED CMAKE_CXX_COMPILER OR CMAKE_CXX_COMPILER STREQUAL "")
    foreach(ver 15 14 13 12)
        find_program(_GXX NAMES "g++-${ver}")
        if(_GXX)
            set(CMAKE_CXX_COMPILER "${_GXX}")
            message(STATUS "Using g++-${ver}: ${_GXX}")
            break()
        endif()
    endforeach()
    if(NOT _GXX)
        message(FATAL_ERROR "No suitable g++ found (tried g++-15/14/13/12). Install with: brew install gcc")
    endif()
endif()

project(wave_spectrum CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# ==========================================================================
# Python and pybind11
# ==========================================================================
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG QUIET)

if(NOT pybind11_FOUND)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
        OUTPUT_VARIABLE pybind11_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(pybind11_DIR)
        list(APPEND CMAKE_PREFIX_PATH "${pybind11_DIR}")
        find_package(pybind11 CONFIG REQUIRED)
    else()
        message(FATAL_ERROR "pybind11 not found. Install with: pip install pybind11")
    endif()
endif()

# ==========================================================================
# Paths
# ==========================================================================
get_filename_component(CPP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(INCLUDE_DIR "${CPP_ROOT}/lib/include")

# ==========================================================================
# pybind11 module (header-only, no .cpp sources needed)
# ==========================================================================
pybind11_add_module(wave_spectrum bindings.cpp)

target_include_directories(wave_spectrum PRIVATE "${INCLUDE_DIR}")

target_compile_options(wave_spectrum PRIVATE
    -march=native
    -O3
    -Wno-deprecated-declarations
    -Wno-unused-variable
    -Wno-unused-but-set-variable
)

# quadmath library (for __float128 support in gcc)
find_library(QUADMATH_LIBRARY quadmath
    PATHS /opt/homebrew/lib/gcc/current /opt/homebrew/lib/gcc/15
    NO_DEFAULT_PATH
)
if(NOT QUADMATH_LIBRARY)
    find_library(QUADMATH_LIBRARY quadmath)
endif()
if(QUADMATH_LIBRARY)
    target_link_libraries(wave_spectrum PRIVATE ${QUADMATH_LIBRARY})
endif()

# ==========================================================================
# Install
# ==========================================================================
install(TARGETS wave_spectrum LIBRARY DESTINATION .)
