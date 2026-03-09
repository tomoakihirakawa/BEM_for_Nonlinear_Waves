/*

/opt/homebrew/bin/g++-13 -std=gnu++2b -x c++-header -fopenmp -fconcepts -pthread
-Ofast -march=native
-I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds
-I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/../include
-I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/../../include -o
./pch.hpp.gch ./pch.hpp

-isysroot
/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX15.2.sdk

 -arch arm64
-mmacosx-version-min=15.1

 */

#ifndef PCH_HPP
#define PCH_HPP

// Standard Library Headers
#include <algorithm>
#include <array>
#include <chrono> // IWYU pragma: export
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <execution>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <ranges>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/utsname.h>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/* -------------------------------------------------------------------------- */

#include "RigidBodyDynamics.hpp"
#include "XML.hpp"
#include "basic.hpp"
#include "basic_IO.hpp"
#include "basic_alias.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "basic_arithmetic_vector_operations.hpp"
#include "basic_constants.hpp"
#include "basic_exception.hpp"
#include "basic_geometry.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "basic_statistics.hpp"
#include "basic_vectors.hpp"
#include "integrationOfODE.hpp"
#include "interpolations.hpp"
#include "kernelFunctions.hpp"
#include "lib_measurement.hpp"
#include "lib_multipole_expansion.hpp"
#include "lib_multipole_expansion_constants.hpp"
#include "lib_spatial_partitioning.hpp"
#include "minMaxOfFunctions.hpp"
#include "rootFinding.hpp"
#include "vtkWriter.hpp"

//

#ifdef USE_TETGEN
#include "tetgen.h"
#endif

#endif // PCH_HPP