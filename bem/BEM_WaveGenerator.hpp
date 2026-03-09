#pragma once

#include "Hadzic2005.hpp"
#include "Network.hpp"

inline std::string normalizeVelocityName(std::string value) {
  auto is_space = [](unsigned char ch) { return ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || ch == '\f' || ch == '\v'; };
  while (!value.empty() && is_space(static_cast<unsigned char>(value.front())))
    value.erase(value.begin());
  while (!value.empty() && is_space(static_cast<unsigned char>(value.back())))
    value.pop_back();
  for (auto &ch : value) {
    if (ch >= 'A' && ch <= 'Z')
      ch = static_cast<char>(ch - 'A' + 'a');
  }
  return value;
}

/*DOC_EXTRACT 0_5_WAVE_GENERATION

## 陽に与えられる境界条件（造波装置など）

造波理論については，\cite{Dean1991}のp.170を参照．
造波板となるobjectに速度を与えることで，造波装置などを模擬する．

*/

inline Tddd relativeVelocity(const Network *net, const networkPoint *p, double t) {
  try {
    if (net->inputJSON.find("relative_velocity")) {
      auto strings = net->inputJSON.at("relative_velocity");
      std::string name = strings[0];

      if (name.find("swimming") != std::string::npos) {
        if (strings.size() < 5)
          throw std::runtime_error("strings size < 5");

        double start = std::stod(strings[1]);
        if (t < start)
          return {0., 0., 0.};

        double Y_tail_target = std::stod(strings[2]);
        double T = std::stod(strings[3]);
        double L = std::stod(strings[4]);

        double a_target = p->initialX[0];
        double w_offset = p->initialX[1];

        double omega = 2. * M_PI / T;
        double k = 2. * M_PI / L;

        // --- 修正1: Ramp-up (助走) 処理 ---
        // 最初の2周期を使って、振幅を 0 -> 100% に滑らかにする。
        // これにより、直線状態から急にS字の速度を与えられた際の初期ショック(片寄り)を防ぐ。
        double ramp_duration = 2.0 * T;
        double time_active = t - start;
        double ramp_factor = (time_active < ramp_duration) ? (time_active / ramp_duration) : 1.0;

        // 係数 A の決定 (Ramp-up適用)
        // Mathematicaと同じ: theta = A * a * Cos(...)
        double A_base = Y_tail_target * (k / L);
        double A = A_base * ramp_factor;

        // --- 速度の数値積分 ---
        int steps = 30; // 精度向上のため少し増やす
        double da = a_target / steps;

        double vx_spine = 0.0;
        double vy_spine = 0.0;

        double theta_current = 0.0;
        double theta_dot_current = 0.0;

        for (int i = 0; i <= steps; ++i) {
          double xi = i * da;

          // Mathematicaと完全一致させる位相: omega*t - k*a
          // ※ startを基準にする場合は t -> (t - start)
          double phase = omega * time_active - k * xi;

          // 角度関数: theta = A * a * Cos(...)
          double theta = A * xi * std::cos(phase);

          // 角速度: d(theta)/dt
          // d/dt(Cos) = -sin * omega
          // ※ Ramp-up中の A の時間微分(dA/dt)の影響は微小なので無視してOK
          double theta_dot = -A * xi * omega * std::sin(phase);

          if (i == steps) {
            theta_current = theta;
            theta_dot_current = theta_dot;
          }

          // 積分: v = integral( v_element ) da
          // v_element_x = d(cos theta)/dt = -sin(theta) * theta_dot
          // v_element_y = d(sin theta)/dt =  cos(theta) * theta_dot
          // ※ 台形則の簡易版として、区間の中間ではなく端点加算で十分近似
          if (i > 0) {
            vx_spine += (-std::sin(theta) * theta_dot) * da;
            vy_spine += (std::cos(theta) * theta_dot) * da;
          }
        }

        // --- 体表面の速度 (回転成分の付加) ---
        // v_surf = v_spine + w_offset * (dn/dt)
        // n = (-sin, cos) -> dn/dt = theta_dot * (-cos, -sin)
        double vx_surf = vx_spine - w_offset * std::cos(theta_current) * theta_dot_current;
        double vy_surf = vy_spine - w_offset * std::sin(theta_current) * theta_dot_current;

        Tddd v_local = {vx_surf, vy_surf, 0.0};

        // グローバル座標系へ回転
        return Dot(net->Q.Rv(), v_local);
      }
    }
    return {0., 0., 0.};
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw;
  }
};

// inline Tddd relativeVelocity(const Network *net, const networkPoint *p, double t) {
//   try {
//     if (net->inputJSON.find("relative_velocity")) {
//       auto strings = net->inputJSON.at("relative_velocity");
//       std::string name = strings[0];
//       if (name.contains("swimming")) {
//         // Example: ["swimming", "start", "A", "T", "lambda"]
//         if (strings.size() < 5)
//           throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "strings size < 5");
//         double start = std::stod(strings[1]);
//         if (t < start)
//           return {0., 0., 0.};

//         double A = std::stod(strings[2]);
//         double T = std::stod(strings[3]);
//         double L = std::stod(strings[4]);
//         double w = 2. * M_PI / T;
//         double k = 2. * M_PI / L;
//         A *= std::pow(std::abs(p->initialX[0]) / L, 2.); // 振幅Aは，変位ではなく，速度の振幅として与える
//         double vy = A * w * std::sin(w * (t - start) - k * std::abs(p->initialX[0]));
//         // ローカル座標系での速度ベクトルを返す
//         Tddd v_local = {0., vy, 0.};
//         // 物体の姿勢(Q)に合わせて回転させて返す
//         return Dot(net->Q.Rv(), v_local);
//       }
//     }
//     return {0., 0., 0.};
//   } catch (std::exception &e) {
//     std::cerr << e.what() << std::endl;
//     throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in relativeVelocity");
//   }
// };

inline Tddd velocityOnBody(const Network *net, const networkPoint *p, double t) { return net->velocityRigidBody(p->X) + relativeVelocity(net, p, t); }

inline T6d velocity(const std::string &name, const std::vector<std::string> strings, networkPoint *p, double t) {
  // std::cout << "velocity(point): " << name << std::endl;
  if ((name.contains("velocity") || name.contains("flow"))) {
    if (strings.size() >= 6) {
      double start = std::stod(strings[1]);
      if (t >= start) {
        double a = std::abs(std::stod(strings[2]));
        double T = std::abs(std::stod(strings[3]));
        double w = std::abs(2. * M_PI / T);
        double h = std::abs(std::stod(strings[4]));
        double z_surface = std::stod(strings[5]);

        DispersionRelation DS;
        if (name.contains("wave_length"))
          DS.set_L_h(T, h); // T is treated as L
        else
          DS.set_w_h(w, h);

        w = DS.w;
        double k = DS.k;

        auto [x, y, z] = p->X;
        z -= z_surface;
        t -= start;
        double wtkx = w * t - k * x;
        double kzh = k * (z + h);
        double kh = k * h;
        double phase_shift = (strings.size() > 6) ? std::stod(strings[6]) : 0.;

        T6d linear_stokes_wave_velocity = {a * w * std::cosh(kzh) / std::sinh(kh) * std::cos(wtkx + phase_shift), 0., -a * w * std::sinh(kzh) / std::sinh(kh) * std::sin(wtkx + phase_shift), 0., 0., 0.};
        return linear_stokes_wave_velocity;
      } else
        return {0., 0., 0., 0., 0., 0.};
    } else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "string must be == 6");
  } else if (name.contains("weird")) {
    if (strings.size() == 6) {
      double a = std::abs(std::stod(strings[2]));
      double w = std::abs(2 * M_PI / std::stod(strings[3]));
      double h = std::abs(std::stod(strings[4]));
      auto [x, y, z] = p->X;
      DispersionRelation DS(w, h);
      double k = std::abs(DS.k);
      return {a * w * std::cos(w * t - k * z), 0., 0., 0., 0., 0.};
    } else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "string must be == 6");
  }
  return {0., 0., 0., 0., 0., 0.};
};

inline T6d velocity(const std::string &name, const std::vector<std::string> strings, double t) {
  // std::cout << "velocity(global): " << name << std::endl;
  auto g = _GRAVITY_;
  try {
    if (name.contains("Goring1979")) {
      if (strings.size() < 4)
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "string must be == 4");
      const double start = std::stod(strings[1]);
      const double H = std::abs(std::stod(strings[2]));
      const double h = std::abs(std::stod(strings[3]));
      const double c = std::sqrt(g * (H + h));
      const double k = std::sqrt(3. * H / (4. * h * h * h));
      double eta = H * std::pow(std::cosh(k * (-c * (t - start))), -2);
      return {c * eta / (h + eta), 0., 0., 0, 0, 0};

    } else if (name.contains("Retzler2000")) {
      const std::vector<Tdd> sample = {{-0.15, 0.},      {-0.11, 0.},      {-0.1, 0.},      {-0.09, 0.},      {-0.08, 0.},      {-0.07, 0.},      {-0.06, 0.},     {-0.0504, 0.0079}, {-0.024, 0.0693}, {-0.0006, 0.1376}, {0.0283, 0.2485},  {0.0504, 0.3594},  {0.059, 0.403},    {0.0695, 0.4426},  {0.083, 0.4792},   {0.1002, 0.5168}, {0.1107, 0.5416}, {0.1211, 0.5663}, {0.1322, 0.5911}, {0.15, 0.6198}, {0.1617, 0.6307},
                                       {0.1764, 0.6446}, {0.1887, 0.6604}, {0.201, 0.6792}, {0.2182, 0.6802}, {0.2305, 0.6446}, {0.2508, 0.5743}, {0.2785, 0.403}, {0.3, 0.2356},     {0.3197, 0.104},  {0.332, 0.0277},   {0.3504, -0.0396}, {0.3707, -0.0851}, {0.3891, -0.0832}, {0.4002, -0.0683}, {0.4254, -0.0307}, {0.45, -0.0099},  {0.5, 0.},        {0.55, 0.},       {0.6, 0.},        {0.65, 0.},     {0.7, 0.}};
      double start = std::stod(strings[1]);
      auto [time, value] = Transpose(sample);
      const auto intp = InterpolationBspline(3, time, value);
      return {intp(t - start), 0., 0., 0., 0., 0.};

    } else if (name.contains("Chaplin2000")) {
      double start = std::stod(strings[1]);
      if (t < start)
        return {0., 0., 0., 0., 0., 0.};
      double A, w;
      if (strings.size() > 3) {
        A = std::stod(strings[2]);
        w = std::stod(strings[3]);
      } else
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "string must be > 3");
      auto v = A * w * std::sin(w * (t - start));
      return {0., v, 0., 0., 0., 0.};

    } else if (name.contains("flap")) {
      double start = std::stod(strings[1]);
      if (strings.size() > 7) {
        double A = std::abs(std::stod(strings[2]));
        double T = std::abs(std::stod(strings[3]));
        double h = std::abs(std::stod(strings[4]));
        double l = std::abs(std::stod(strings[5]));

        DispersionRelation DS;
        if (name.contains("wave_length"))
          DS.set_L_h(T, h);
        else
          DS.set_w_h(2 * M_PI / T, h);

        double w = DS.w;
        double k = DS.k;
        double L = DS.L;

        if (name.contains("wave_steepness")) {
          auto eps = A;
          auto H = eps * L;
          A = H / 2.;
        }

        Tddd axis = {std::stod(strings[6]), std::stod(strings[7]), std::stod(strings[8])};
        const double kh = k * h;
        const double H = 2 * A;
        const double S = H * kh / (4. * std::sinh(kh)) * (std::sinh(2. * kh) + 2. * kh) / (kh * std::sinh(kh) - std::cosh(kh) + 1.);
        const double a = S / 2.;
        double dthetadt = -((a * l * w * std::cos(t * w)) / (std::pow(l, 2) + std::pow(a, 2) * std::pow(std::sin(t * w), 2)));
        if (name.contains("negative"))
          dthetadt *= -1.;
        auto [wx, wy, wz] = -dthetadt * Normalize(axis);
        return {0., 0., 0., wx, wy, wz};
      } else
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "string must be > 7");

    } else if (name.contains("piston")) {
      double start = std::stod(strings[1]);
      if (t < start)
        return {0., 0., 0., 0., 0., 0.};
      if (strings.size() > 7) {
        const double A = std::abs(std::stod(strings[2]));
        const double h = std::abs(std::stod(strings[4]));
        double T = std::abs(std::stod(strings[3]));

        DispersionRelation DS;
        if (name.contains("wave_length"))
          DS.set_L_h(T, h);
        else
          DS.set_w_h(2 * M_PI / T, h);

        double w = DS.w;
        double k = DS.k;
        const double kh2 = 2. * k * h;
        const double F = 2. * (std::cosh(kh2) - 1.) / (kh2 + std::sinh(kh2));
        const double H = 2. * A;
        const double S = H / F;
        const double smoothing_function = 1.;
        const double shift = 0.;
        double dsdt = smoothing_function * S / 2. * w * std::sin(w * (t - start) + shift);
        Tddd axis = {std::stod(strings[5]), std::stod(strings[6]), std::stod(strings[7])};
        if (name.contains("negative"))
          dsdt *= -1.;
        return {dsdt * axis[0], dsdt * axis[1], dsdt * axis[2], 0., 0., 0.};
      } else
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "string must be > 7");

    } else if (name.contains("sinusoidal") || name.contains("sin") || name.contains("cos")) {
      if (strings.size() < 5)
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "invalid string size");

      double start = std::stod(strings[1]);
      if (t < start)
        return {0, 0, 0, 0, 0, 0};
      else {
        double a = std::stod(strings[2]);
        double w = std::abs(2 * M_PI / std::stod(strings[3]));
        double A = (name.contains("cos") ? a * w * std::sin(w * (t - start)) : a * w * std::cos(w * (t - start)));
        T6d axis{};
        for (size_t i = 4; i < strings.size() && i < 10; ++i)
          axis[i - 4] = A * std::stod(strings[i]);
        return axis;
      }
    } else if (name.contains("constant") || name.contains("const")) {
      std::cout << "constant velocity" << std::endl;
      double start = std::stod(strings[1]);
      if (t >= start)
        return {std::stod(strings[2]), std::stod(strings[3]), std::stod(strings[4]), std::stod(strings[5]), std::stod(strings[6]), std::stod(strings[7])};
    } else if (name.contains("file")) {
      double a = std::stod(strings[2]);
      T6d axis = {std::stod(strings[3]), std::stod(strings[4]), std::stod(strings[5]), std::stod(strings[6]), std::stod(strings[7]), std::stod(strings[8])};
      return a * axis;
    } else if (name.contains("Hadzic2005")) {
      double start = std::stod(strings[1]);
      Hadzic2005 hadzic2005(start);
      return hadzic2005.getVelocity(t);
    }
    std::cout << "velocity: " << name << " does not exist" << std::endl;
    return {0., 0., 0., 0., 0., 0.};
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "error in velocity");
  }
};

inline T6d acceleration(const std::string &name, const std::vector<std::string> strings, const double t) {
  if (name.contains("Hadzic2005")) {
    double start = std::stod(strings[1]);
    Hadzic2005 hadzic2005(start);
    return hadzic2005.getAccel(t);
  }
  return {0., 0., 0., 0., 0., 0.};
};
