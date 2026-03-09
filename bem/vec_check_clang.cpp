#include "lib_multipole_expansion.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

volatile double g_veccheck_sink = 0.0;

__attribute__((noinline)) static void veccheck_run() {
  constexpr int32_t N = 256;

  target4FMM t({0.0, 0.0, 0.0});

  t.near_indices.resize(static_cast<std::size_t>(N));
  t.near_weights_phi_D.resize(static_cast<std::size_t>(N));
  t.near_weights_phin_D.resize(static_cast<std::size_t>(N));
  t.near_weights_phi_F.resize(static_cast<std::size_t>(N));
  t.near_weights_phin_F.resize(static_cast<std::size_t>(N));

  for (int32_t i = 0; i < N; ++i) {
    const std::size_t ui = static_cast<std::size_t>(i);
    t.near_indices[ui] = i;
    t.near_weights_phi_D[ui] = 1.0 / (1.0 + static_cast<double>(i));
    t.near_weights_phin_D[ui] = 2.0 / (1.0 + static_cast<double>(i));
    t.near_weights_phi_F[ui] = static_cast<float>(t.near_weights_phi_D[ui]);
    t.near_weights_phin_F[ui] = static_cast<float>(t.near_weights_phin_D[ui]);
  }

  t.near_run_base_idx = {0};
  t.near_run_pos = {0};
  t.near_run_len = {N};

  std::vector<double> phi(static_cast<std::size_t>(N));
  std::vector<double> phin(static_cast<std::size_t>(N));
  std::vector<float> phi_f(static_cast<std::size_t>(N));
  std::vector<float> phin_f(static_cast<std::size_t>(N));

  std::vector<const double *> phi_ptr(static_cast<std::size_t>(N));
  std::vector<const double *> phin_ptr(static_cast<std::size_t>(N));

  for (int32_t i = 0; i < N; ++i) {
    const std::size_t ui = static_cast<std::size_t>(i);
    phi[ui] = 0.1 + static_cast<double>(i) * 0.01;
    phin[ui] = 0.2 + static_cast<double>(i) * 0.02;
    phi_f[ui] = static_cast<float>(phi[ui]);
    phin_f[ui] = static_cast<float>(phin[ui]);
    phi_ptr[ui] = &phi[ui];
    phin_ptr[ui] = &phin[ui];
  }

  const auto r_dense = t.integrateNearFieldDense(phi.data(), phin.data());
  const auto r_float = t.integrateNearFieldFloat(phi_f.data(), phin_f.data());
  const auto r_ptr = t.integrateNearField(phi_ptr.data(), phin_ptr.data());

  g_veccheck_sink += r_dense[0] + r_dense[1] + r_float[0] + r_float[1] + r_ptr[0] + r_ptr[1];
}

int main() {
  veccheck_run();
  return g_veccheck_sink != 0.0 ? 0 : 1;
}

