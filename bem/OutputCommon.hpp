#pragma once
#include <chrono>
#include <ctime>
#include <filesystem>

struct OutputContext {
   double dt{};
   int time_step{};
   double simulation_time{};
   std::filesystem::path output_directory;
   std::clock_t cpu_clock_start{};
   std::chrono::high_resolution_clock::time_point wall_clock_start;
};