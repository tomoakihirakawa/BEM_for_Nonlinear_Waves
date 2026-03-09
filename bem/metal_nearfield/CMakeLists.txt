# ==========================================================================
# CMakeLists.txt for metal_nearfield library
#
# This library must be compiled with Apple Clang (not GCC) because:
# - Metal requires Objective-C++
# - Metal.framework is only accessible via Apple toolchains
#
# Build separately:
#   cd metal_nearfield
#   mkdir build && cd build
#   cmake .. -DCMAKE_BUILD_TYPE=Release
#   make
#
# This produces: libmetal_nearfield.dylib
# ==========================================================================

cmake_minimum_required(VERSION 3.16)

# Force Apple Clang
set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")

project(metal_nearfield LANGUAGES CXX OBJCXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_OBJCXX_STANDARD 17)
set(CMAKE_OBJCXX_STANDARD_REQUIRED ON)

# Enable ARC (Automatic Reference Counting)
set(CMAKE_OBJCXX_FLAGS "${CMAKE_OBJCXX_FLAGS} -fobjc-arc")

# Find required frameworks
find_library(METAL_FRAMEWORK Metal REQUIRED)
find_library(FOUNDATION_FRAMEWORK Foundation REQUIRED)

# Create the shared library
add_library(metal_nearfield SHARED
    metal_nearfield.mm
)

target_include_directories(metal_nearfield PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(metal_nearfield PRIVATE
    ${METAL_FRAMEWORK}
    ${FOUNDATION_FRAMEWORK}
)

# Set library properties
set_target_properties(metal_nearfield PROPERTIES
    VERSION 1.0.0
    SOVERSION 1
    PUBLIC_HEADER "metal_nearfield.h"
    INSTALL_RPATH "@loader_path"
)

# Compile Metal shader to metallib (optional - requires Xcode)
find_program(METAL_COMPILER metal HINTS /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin)
if(METAL_COMPILER)
    message(STATUS "Found Metal compiler: ${METAL_COMPILER}")
    set(METAL_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/metal_nearfield.metal")
    set(METAL_AIR "${CMAKE_CURRENT_BINARY_DIR}/metal_nearfield.air")
    set(METAL_LIB "${CMAKE_CURRENT_BINARY_DIR}/metal_nearfield.metallib")

    add_custom_command(
        OUTPUT ${METAL_AIR}
        COMMAND ${METAL_COMPILER} -c ${METAL_SOURCE} -o ${METAL_AIR} -fno-fast-math
        DEPENDS ${METAL_SOURCE}
        COMMENT "Compiling Metal shader to AIR (with -fno-fast-math for Kahan summation)"
    )

    find_program(METALLIB_LINKER metallib HINTS /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin)
    if(METALLIB_LINKER)
        add_custom_command(
            OUTPUT ${METAL_LIB}
            COMMAND ${METALLIB_LINKER} ${METAL_AIR} -o ${METAL_LIB}
            DEPENDS ${METAL_AIR}
            COMMENT "Linking Metal shader library"
        )

        add_custom_target(metal_shader ALL DEPENDS ${METAL_LIB})
        add_dependencies(metal_nearfield metal_shader)

        # Copy metallib to output directory
        add_custom_command(TARGET metal_nearfield POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy ${METAL_LIB} $<TARGET_FILE_DIR:metal_nearfield>/metal_nearfield.metallib
            COMMENT "Copying metallib to output directory"
        )
    endif()
else()
    message(STATUS "Metal compiler not found - shader will be compiled at runtime")
    message(STATUS "  (Install Xcode.app for precompiled metallib)")
endif()

# Always copy .metal source for runtime compilation fallback
add_custom_command(TARGET metal_nearfield POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/metal_nearfield.metal" $<TARGET_FILE_DIR:metal_nearfield>/metal_nearfield.metal
    COMMENT "Copying Metal source for runtime compilation"
)

# Installation
install(TARGETS metal_nearfield
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include
)

if(METALLIB_LINKER)
    install(FILES ${METAL_LIB} DESTINATION lib)
endif()

# Print summary
message(STATUS "")
message(STATUS "=== metal_nearfield library configuration ===")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Metal framework: ${METAL_FRAMEWORK}")
message(STATUS "Foundation framework: ${FOUNDATION_FRAMEWORK}")
message(STATUS "")
