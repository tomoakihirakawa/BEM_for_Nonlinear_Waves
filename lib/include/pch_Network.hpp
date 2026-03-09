/*

/opt/homebrew/bin/g++-14 -std=gnu++2b -x c++-header -fopenmp -fconcepts -pthread -Ofast -march=native -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/../include -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/../../include -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX15.2.sdk -o ./pch_Network.hpp.gch ./pch_Network.hpp

 */

#ifndef PCH_HPP
#define PCH_HPP

// Standard Library Headers
#include <sys/utsname.h>
#include <unistd.h>
#include <algorithm>
#include <array>
#include <chrono>
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
#include <tuple>
#include <type_traits>
#include <typeinfo>
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

#include "tetgen1.6.0/tetgen.h"

//

#include "Network.hpp"

#endif  // PCH_HPP