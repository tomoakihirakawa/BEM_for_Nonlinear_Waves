#include <fstream>
#include <iostream>
#include "Hadzic2005.hpp"
#include "basic.hpp"

int main(int argc, char* argv[]) {
   double start_time = 0;
   Hadzic2005 hadzic2005(start_time);

   std::ofstream outFile("Hadzic2005Flap_check.dat");
   if (!outFile.is_open()) {
      std::cerr << "Failed to open the output file." << std::endl;
      return 1;
   }

   for (const auto t : Subdivide(0., 10., 1000)) {
      auto velocity = hadzic2005.getVelocity(t);
      outFile << t;
      for (const auto& v : velocity)
         outFile << " " << v;
      outFile << " " << hadzic2005.getAngle(t) << std::endl;
   }

   outFile.close();
   return 0;
}
