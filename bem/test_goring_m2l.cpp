/*
 * Test M2L with actual BEM data from Goring1979 case
 * Read two time steps and perform M2L translation
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <chrono>
#include "../../include/lib_plane_wave_expansion_full.hpp"
#include "../../include/lib_plane_wave_m2l.hpp"

using cmplx = std::complex<double>;

// Simple VTU reader for point coordinates
struct Point3D {
    double x, y, z;
};

std::vector<Point3D> readVTUPoints(const std::string& filename) {
    std::vector<Point3D> points;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return points;
    }

    std::string line;
    bool in_points = false;

    while (std::getline(file, line)) {
        if (line.find("<Points>") != std::string::npos) {
            in_points = true;
            continue;
        }

        if (line.find("</Points>") != std::string::npos) {
            break;
        }

        if (in_points && line.find("NumberOfComponents=\"3\"") != std::string::npos) {
            // Next line contains the actual coordinates
            std::getline(file, line);
            std::istringstream iss(line);
            double x, y, z;
            while (iss >> x >> y >> z) {
                points.push_back({x, y, z});
            }
        }
    }

    return points;
}

// Create multipole expansion from point distribution
template<int N>
void createMultipoleFromPoints(
    const std::vector<Point3D>& points,
    const Point3D& center,
    std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>& MM
) {
    // Zero out
    for (auto& m : MM) {
        std::get<0>(m) = cmplx(0.0);
        std::get<1>(m) = cmplx(0.0);
    }

    // Use simple approximation: assume unit charge at each point
    // and compute monopole and dipole moments

    double total_charge = static_cast<double>(points.size());

    // M_0^0: monopole
    std::get<0>(MM[0]) = cmplx(total_charge, 0.0);

    // M_1^0, M_1^±1: dipole moments
    if (N >= 1) {
        double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
        for (const auto& p : points) {
            sum_x += (p.x - center.x);
            sum_y += (p.y - center.y);
            sum_z += (p.z - center.z);
        }

        // M_1^0 ~ z-component
        std::get<0>(MM[2]) = cmplx(sum_z / total_charge, 0.0);  // idx = 1*1 + 1 + 0 = 2

        // M_1^±1 ~ x, y components (simplified)
        // For real-valued, we approximate
        std::get<0>(MM[1]) = cmplx(sum_x / total_charge, 0.0);  // idx = 1*1 + 1 - 1 = 1
        std::get<0>(MM[3]) = cmplx(sum_y / total_charge, 0.0);  // idx = 1*1 + 1 + 1 = 3
    }
}

template<int N>
void runM2LTest(const std::vector<Point3D>& points1, const std::vector<Point3D>& points2,
                const Point3D& center1, const Point3D& center2) {

    std::cout << "\n=== M2L Test with Order N = " << N << " ===" << std::endl;

    using MM_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;
    using L_Type = std::array<std::tuple<cmplx, cmplx>, (N+1)*(N+1)>;

    // Create multipole expansion for source box
    MM_Type source_MM;
    createMultipoleFromPoints<N>(points1, center1, source_MM);

    std::cout << "Source box: " << points1.size() << " points" << std::endl;
    std::cout << "Center: (" << center1.x << ", " << center1.y << ", " << center1.z << ")" << std::endl;
    std::cout << "M_0^0 = " << std::get<0>(source_MM[0]) << std::endl;

    // Translation vector
    Tddd translation = {
        center2.x - center1.x,
        center2.y - center1.y,
        center2.z - center1.z
    };

    double dist = std::sqrt(translation[0]*translation[0] +
                           translation[1]*translation[1] +
                           translation[2]*translation[2]);

    std::cout << "Target box: " << points2.size() << " points" << std::endl;
    std::cout << "Center: (" << center2.x << ", " << center2.y << ", " << center2.z << ")" << std::endl;
    std::cout << "Translation distance: " << dist << " m" << std::endl;

    // Perform M2L
    L_Type target_L;
    for (auto& l : target_L) {
        std::get<0>(l) = cmplx(0.0);
        std::get<1>(l) = cmplx(0.0);
    }

    auto start = std::chrono::high_resolution_clock::now();

    PlaneWaveM2L::PlaneWaveM2L<N, PlaneWaveQuadrature::Precision::Digit6>::translateMultipoleToLocal(
        source_MM, target_L, translation, false);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "L_0^0 = " << std::get<0>(target_L[0]) << std::endl;
    std::cout << "Computation time: " << duration.count() << " μs" << std::endl;

    // Check for NaN
    bool has_nan = false;
    for (const auto& l : target_L) {
        if (std::isnan(std::abs(std::get<0>(l)))) {
            has_nan = true;
            break;
        }
    }

    if (has_nan) {
        std::cout << "✗ Result contains NaN" << std::endl;
    } else {
        std::cout << "✓ Result is finite" << std::endl;
    }
}

int main(int argc, char** argv) {
    std::cout << std::scientific << std::setprecision(6);

    std::cout << "╔════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  M2L Test with Goring1979 BEM Data                            ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════╝" << std::endl;

    // Read two time steps
    std::string file1 = "/Users/tomoaki/Goring1979/output/water_100.vtu";
    std::string file2 = "/Users/tomoaki/Goring1979/output/water_101.vtu";

    std::cout << "\nReading time step 100..." << std::endl;
    auto points1 = readVTUPoints(file1);

    std::cout << "Reading time step 101..." << std::endl;
    auto points2 = readVTUPoints(file2);

    if (points1.empty() || points2.empty()) {
        std::cerr << "Failed to read VTU files" << std::endl;
        std::cerr << "Points in step 100: " << points1.size() << std::endl;
        std::cerr << "Points in step 101: " << points2.size() << std::endl;

        // Create synthetic test instead
        std::cout << "\n--- Using synthetic data instead ---" << std::endl;

        std::mt19937 gen(42);
        std::uniform_real_distribution<double> dis(0.0, 1.0);

        points1.clear();
        points2.clear();

        for (int i = 0; i < 100; ++i) {
            points1.push_back({dis(gen), dis(gen), dis(gen)});
            points2.push_back({dis(gen) + 3.0, dis(gen), dis(gen)});
        }
    }

    std::cout << "\nPoints in step 100: " << points1.size() << std::endl;
    std::cout << "Points in step 101: " << points2.size() << std::endl;

    // Compute centers (simple average)
    Point3D center1 = {0.0, 0.0, 0.0};
    for (const auto& p : points1) {
        center1.x += p.x;
        center1.y += p.y;
        center1.z += p.z;
    }
    center1.x /= points1.size();
    center1.y /= points1.size();
    center1.z /= points1.size();

    Point3D center2 = {0.0, 0.0, 0.0};
    for (const auto& p : points2) {
        center2.x += p.x;
        center2.y += p.y;
        center2.z += p.z;
    }
    center2.x /= points2.size();
    center2.y /= points2.size();
    center2.z /= points2.size();

    // Run M2L tests for different orders
    runM2LTest<4>(points1, points2, center1, center2);
    runM2LTest<5>(points1, points2, center1, center2);
    runM2LTest<6>(points1, points2, center1, center2);
    runM2LTest<7>(points1, points2, center1, center2);
    runM2LTest<8>(points1, points2, center1, center2);
    runM2LTest<9>(points1, points2, center1, center2);
    runM2LTest<10>(points1, points2, center1, center2);

    std::cout << "\n╔════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  Test completed                                                ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════╝" << std::endl;

    return 0;
}
