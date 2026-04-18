#pragma once

#include <string>

struct PVDWriter;

struct outputInfo {
  std::string pvd_file_name;
  std::string vtu_file_name;
  PVDWriter* PVD = nullptr;
  outputInfo() = default;
};
