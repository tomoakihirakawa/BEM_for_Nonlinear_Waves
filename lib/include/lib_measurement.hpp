#ifndef lib_measurement_H
#define lib_measurement_H

#include <chrono>
#include <tuple>
#include <vector>

// struct Timer {
//    // std::chrono::time_point<std::chrono::high_resolution_clock> start;
//    //! 一つ目は，前回からの経過時間，二つ目は作成した時間からの経過時間
//    std::chrono::time_point<std::chrono::system_clock> start;
//    std::chrono::time_point<std::chrono::system_clock> lasttime;
//    Timer() : start(std::chrono::system_clock::now()), lasttime(start){};
//    std::array<double, 2> operator()() {
//       auto tmp = this->lasttime;
//       this->lasttime = std::chrono::system_clock::now();
//       // return {this->lasttime - tmp, this->lasttime - this->start};
//       return {(std::chrono::duration<double>(this->lasttime - tmp)).count(),
//               (std::chrono::duration<double>(this->lasttime - this->start)).count()};
//    };
// };

/* ------------------------------------------------------ */

struct TimeWatch {
  using clock = std::chrono::steady_clock;
  static_assert(clock::is_steady, "steady clock required");
  clock::time_point start, last;

  TimeWatch() noexcept : start(clock::now()), last(start) {}

  void reset() noexcept { start = last = clock::now(); }

  //! 一つ目は，前回からの経過時間，二つ目は作成した時間からの経過時間
  // std::vector<std::chrono::duration<double>> get()
  std::array<double, 2> get() {
    auto now = clock::now();
    double dt = std::chrono::duration<double>(now - last).count();
    double tt = std::chrono::duration<double>(now - start).count();
    last = now;
    return {dt, tt};
  };

  std::array<double, 2> operator()() {
    auto now = clock::now();
    double dt = std::chrono::duration<double>(now - last).count();
    double tt = std::chrono::duration<double>(now - start).count();
    last = now;
    return {dt, tt};
  };
};

#endif