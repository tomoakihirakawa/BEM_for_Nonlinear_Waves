#define DEM
#include "Network.hpp"
#include "SPH_weightingFunctions.hpp"
#include "vtkWriter.hpp"

int main(int arg, char** argv) {

   Tddd O = {0, 0, 0};
   double dx = 1;
   auto f = [](const Tddd& X) { return std::sin(2. * M_PI * Total(X)); };
   auto grad_f = [](const Tddd& X) {
      auto v = 2. * M_PI * std::cos(2. * M_PI * Total(X));
      return Tddd{v, v, v};
   };
   //    auto integral_f = [](const Tddd& X) { return std::sin(2. * M_PI * Total(X)); };
   w_Bspline5();
};