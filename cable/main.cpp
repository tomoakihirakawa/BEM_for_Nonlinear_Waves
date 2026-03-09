#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "Network.hpp"
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "basic_vectors.hpp"
#include "integrationOfODE.hpp"
#include "vtkWriter.hpp"

// (cd builds/build_cable; python3 ../../extract_comments.py README.md ./ ../../)
// echo '(cd builds/build_cable; python3 ../../extract_comments.py README.md ./ ../../)'

/*DOC_EXTRACT 0_cable_dynamics

# ケーブルの動的解析

## 直線要素を用いたシミュレーション

NOTE: 弦の振動を支配する方程式として，波動方程式$\frac{\partial^2 u}{\partial t^2} = c^2 \frac{\partial^2 u}{\partial x^2}$よく紹介される．
この方程式は，ある固定した点$x$における弦の変位$u$の加速度が，弦の曲げ剛性$c^2$かける曲率に比例することを表している．

直線で結ばれた節点上にケーブルの自重を集中させ，その節点に働く張力や重力から，節点の運動を追っていく．

剛性は，ヤング率$E$と断面積$A$から$EA$．
張力$T$は，$T = EA \frac{\Delta L}{L}$となる．

オイラー法，Leap-Frog法，Runge-Kutta法を用いて，弾性体の動きをシミュレーション．

* 剛性$[N/m]$:$1400 \times 10^6$
* 減衰$[N/(m/s^2)]$:$0.9$
* 自然長$[m]$:$1$

![sample.gif](./sample.gif)

`const double stiffness = 10000;`の場合
![sample_2.gif](./sample_2.gif)

チェーン

## 実行方法

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

\insert{1_MooringLine}

*/

// #define USE_EULER
#define USE_LEAP_FROG
// #define USE_RK4
const std::string filename = "cable_leapfrog_K14E8_NL01d0.json";

const double natural_length = 1.0;              //! [m]
const double stiffness = 14 * std::pow(10, 8);  //! [N/m]
const double damp = .9;                         //! [N/(m/s^2)]
const double w = 348.5;                         //! [kg/m]
const double dt = 0.0001;                       //! [s]
const int max_step = 1000000;
const std::array<double, 3> gravity = {0, 0, -9.81};

struct Node {

   const double mass = natural_length * w;
   std::array<double, 3> X;
   std::array<double, 3> velocity;
   std::array<double, 3> acceleration;
   std::array<double, 3> tension;
   //
   Node(std::array<double, 3> X, std::array<double, 3> velocity, std::array<double, 3> acceleration)
       : X(X), velocity(velocity), acceleration(acceleration) {}

   double t;
   LeapFrog<decltype(X)> LPFG;

   RungeKutta<decltype(X)> RK_x;
   RungeKutta<decltype(X)> RK_v;
};
// nodes 30
std::vector<Node> nodes;
/* -------------------------------------------------------------------------- */

void simulateCableDynamics(double t, double dt) {

   auto tension = [&](const int i) {
      std::array<double, 3> force;
      force.fill(0.);
      if (i - 1 >= 0) {
         auto v = nodes[i - 1].X - nodes[i].X;
         double disp = Norm(v) - natural_length;
         double strain = disp / natural_length;
         if (disp > 0.) force += stiffness * strain * Normalize(v);
         auto relative_velocity = nodes[i - 1].velocity - nodes[i].velocity;
         // force -= damp * Projection(-relative_velocity, v) / dt;
         force += damp * relative_velocity / dt;
      }
      if (i + 1 < nodes.size()) {
         auto v = nodes[i + 1].X - nodes[i].X;
         double disp = Norm(v) - natural_length;
         double strain = disp / natural_length;
         if (disp > 0.) force += stiffness * strain * Normalize(v);
         auto relative_velocity = nodes[i + 1].velocity - nodes[i].velocity;
         // force -= damp * Projection(-relative_velocity, v) / dt;
         force += damp * relative_velocity / dt;
      }

      return force;
   };

#ifdef USE_LEAP_FROG
   const int total_steps = 2;
#elif defined USE_RK4
   const int total_steps = 4;
#else
   const int total_steps = 1;
#endif

   for (size_t step = 0; step < total_steps; ++step) {
      for (size_t i = 0; i < nodes.size(); ++i) {

         nodes[i].tension = tension(i);
         nodes[i].acceleration = nodes[i].tension / nodes[i].mass + gravity;

         //! 端点が静止状態の場合
         if (i == 0 || i == nodes.size() - 1)
            nodes[i].acceleration = {0., 0., 0.};

         // if (t < 30)
         //    if (i == nodes.size() - 1) {
         //       const double T = 3;
         //       const double w = 2 * M_PI / T;
         //       //! 最後の節点は，10まで最後の接点をぐるぐる回す
         //       // nodes[i].velocity = Cross(nodes[i].X, {0., 0., 1.});
         //       if (Between(t, {2 * T, 3 * T}) || Between(t, {4 * T, 5 * T}))
         //          nodes[i].acceleration = {0., 50 * std::cos(w * t), 50 * std::cos(w * t)};
         //       else
         //          nodes[i].acceleration = -nodes[i].velocity / dt;
         //    }
      }

#ifdef USE_EULER
      //! オイラー法
      for (size_t i = 1; i < nodes.size(); ++i) {
         nodes[i].velocity += nodes[i].acceleration * dt;
         nodes[i].X += nodes[i].velocity * dt;
      }
#elif defined USE_LEAP_FROG
      for (size_t i = 1; i < nodes.size(); ++i) {
         nodes[i].LPFG.push(nodes[i].acceleration);
         nodes[i].t = nodes[i].LPFG.get_t();
         nodes[i].X = nodes[i].LPFG.get_x();
         nodes[i].velocity = nodes[i].LPFG.get_v();
      }
#elif defined USE_RK4
      for (size_t i = 1; i < nodes.size(); ++i) {
         nodes[i].RK_x.push(nodes[i].velocity);
         nodes[i].RK_v.push(nodes[i].acceleration);
         nodes[i].t = nodes[i].RK_x.get_t();
         nodes[i].X = nodes[i].RK_x.get_x();
         nodes[i].velocity = nodes[i].RK_v.get_x();
      }
#endif
   }
}

/* -------------------------------------------------------------------------- */

int main() {
   const double L = 30;                     //! 全長 [m]
   int n = std::round(L / natural_length);  //! 節点数

   for (int i = 0; i < n + 1; ++i)
      nodes.push_back(Node({i * natural_length, 0, 0}, {0., 0., 0.}, {0., 0., 0.}));
   // Example setup
   std::vector<double> times;
   std::vector<std::vector<std::array<double, 3>>> positions;
   std::vector<std::vector<std::array<double, 3>>> tensions;

   double t = 0;
#ifdef USE_LEAP_FROG
   for (auto& node : nodes)
      node.LPFG.initialize(dt, t, node.X, node.velocity);
#endif

   /* --------------------------------------------- */

   for (int i = 0; i < max_step; ++i) {
#ifdef USE_RK4
      for (auto& node : nodes) {
         node.RK_x.initialize(dt, t, node.X, 4);
         node.RK_v.initialize(dt, t, node.velocity, 4);
      }
#endif
      std::cout << "t = " << t << std::endl;
      simulateCableDynamics(t, dt);

      if (i % 1000 == 0) {
         times.push_back(t);
         positions.push_back({});
         tensions.push_back({});
         for (const Node& node : nodes) {
            if (!isFinite(node.X))
               break;
            positions.back().push_back(node.X);
            tensions.back().push_back(node.tension);
         }
      }
      t += dt;
   }

   /* ----------------------------- OUTPUT AS JSON ----------------------------- */
#ifdef USE_EULER
   std::ofstream myfile("cable_euler.json");
#elif defined USE_LEAP_FROG
   std::ofstream myfile(filename);
#elif defined USE_RK4
   std::ofstream myfile("cable_rk4.json");
#endif
   myfile << "{";
   std::array<std::string, 2> bracket = {"[", "]"};
   myfile << "\"time\" : " << std::make_tuple(times, bracket) << ", ";
   myfile << "\"position\" : " << std::make_tuple(positions, bracket) << ", ";
   myfile << "\"tension\" : " << std::make_tuple(tensions, bracket);
   myfile << "}";
   myfile.close();

   /* ----------------------------- OUTPUT AS VTK ----------------------------- */
   PVDWriter pvd("/Users/tomoaki/Cable/line.pvd");
   auto net = new Network;
   std::vector<networkPoint*> points;
   int i = 0;
   for (auto k = 0; k < positions.size(); ++k)
      if (k % 10 == 0) {
         auto t = times[i];
         if (i == 0) {
            for (const auto& X : positions[i])
               points.push_back(new networkPoint(net, X));
            for (auto j = 0; j < points.size() - 1; ++j)
               new networkLine(net, points[j], points[j + 1]);
         } else {
            for (auto j = 0; j < points.size(); ++j)
               points[j]->setXSingle(positions[i][j]);
         }
         auto filename = "/Users/tomoaki/Cable/line" + std::to_string(i) + ".vtp";
         std::ofstream ofs(filename);
         vtkPolygonWrite(ofs, net->getLines());
         pvd.push(filename, t);
         ++i;
      }
   pvd.output();
   /* -------------------------------------------------------------------------- */
   return 0;
}
