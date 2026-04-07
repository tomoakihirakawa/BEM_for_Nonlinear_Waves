/**
 * @file bindings.cpp
 * @brief pybind11 bindings for wave spectrum classes
 *
 * Wraps RandomWaterWaveTheory, WaterWaveTheory, and DispersionRelation
 * from rootFinding.hpp for Python-based verification and visualization.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "rootFinding.hpp"

namespace py = pybind11;

PYBIND11_MODULE(wave_spectrum, m) {
  m.doc() = "Wave spectrum module: JONSWAP, Bretschneider-Mitsuyasu, and wave component discretization";

  // ---- Enums ----

  py::enum_<SpectrumType>(m, "SpectrumType")
      .value("BRETSCHNEIDER_MITSUYA", SpectrumType::BRETSCHNEIDER_MITSUYA)
      .value("JONSWAP", SpectrumType::JONSWAP);

  py::enum_<WaveParamMode>(m, "WaveParamMode")
      .value("H13_T13", WaveParamMode::H13_T13)
      .value("H13_TP", WaveParamMode::H13_TP)
      .value("HM0_T13", WaveParamMode::HM0_T13)
      .value("HM0_TP", WaveParamMode::HM0_TP);

  // ---- DispersionRelation ----

  py::class_<DispersionRelation>(m, "DispersionRelation")
      .def(py::init<>())
      .def("set_T_h", &DispersionRelation::set_T_h)
      .def("set_w_h", &DispersionRelation::set_w_h)
      .def_readwrite("L", &DispersionRelation::L)
      .def_readwrite("T", &DispersionRelation::T)
      .def_readwrite("w", &DispersionRelation::w)
      .def_readwrite("k", &DispersionRelation::k);

  // ---- WaterWaveTheory ----

  py::class_<WaterWaveTheory>(m, "WaterWaveTheory")
      .def(py::init<>())
      .def("set_T_h", &WaterWaveTheory::set_T_h)
      .def("set_w_h", &WaterWaveTheory::set_w_h)
      .def("set_L_h", &WaterWaveTheory::set_L_h)
      .def("phi", py::overload_cast<const std::array<double, 3>&, const double>(&WaterWaveTheory::phi, py::const_))
      .def("eta", py::overload_cast<const std::array<double, 3>&, const double>(&WaterWaveTheory::eta, py::const_))
      .def("etaZeroAtRest", &WaterWaveTheory::etaZeroAtRest)
      .def("gradPhi", &WaterWaveTheory::gradPhi)
      .def("gradPhi_t", py::overload_cast<const std::array<double, 3>&, const double>(&WaterWaveTheory::gradPhi_t, py::const_))
      .def_readwrite("h", &WaterWaveTheory::h)
      .def_readwrite("L", &WaterWaveTheory::L)
      .def_readwrite("T", &WaterWaveTheory::T)
      .def_readwrite("w", &WaterWaveTheory::w)
      .def_readwrite("k", &WaterWaveTheory::k)
      .def_readwrite("c", &WaterWaveTheory::c)
      .def_readwrite("A", &WaterWaveTheory::A)
      .def_readwrite("bottom_z", &WaterWaveTheory::bottom_z)
      .def_readwrite("phase_shift", &WaterWaveTheory::phase_shift)
      .def_readwrite("theta", &WaterWaveTheory::theta);

  // ---- RandomWaterWaveTheory ----

  py::class_<RandomWaterWaveTheory>(m, "RandomWaterWaveTheory")
      // Constructors / factories
      .def(py::init<>())
      .def(py::init<double, double, double, double>(),
           py::arg("H13"), py::arg("T13"), py::arg("h"), py::arg("bottom_z"),
           "Backward-compatible constructor (Bretschneider, H13_T13)")
      .def_static("create", &RandomWaterWaveTheory::create,
           py::arg("spectrum"), py::arg("mode"),
           py::arg("height"), py::arg("period"),
           py::arg("gamma"), py::arg("h"), py::arg("bottom_z"))
      .def_static("Bretschneider", &RandomWaterWaveTheory::Bretschneider,
           py::arg("H13"), py::arg("T13"), py::arg("h"), py::arg("bottom_z"))
      .def_static("JONSWAP_H13_T13", &RandomWaterWaveTheory::JONSWAP_H13_T13,
           py::arg("H13"), py::arg("T13"), py::arg("gamma"),
           py::arg("h"), py::arg("bottom_z"))
      .def_static("JONSWAP_Hm0_Tp", &RandomWaterWaveTheory::JONSWAP_Hm0_Tp,
           py::arg("Hm0"), py::arg("Tp"), py::arg("gamma"),
           py::arg("h"), py::arg("bottom_z"))

      // Properties
      .def_readwrite("h", &RandomWaterWaveTheory::h)
      .def_readwrite("bottom_z", &RandomWaterWaveTheory::bottom_z)
      .def_readwrite("H13", &RandomWaterWaveTheory::H13)
      .def_readwrite("Hm0", &RandomWaterWaveTheory::Hm0)
      .def_readwrite("T13", &RandomWaterWaveTheory::T13)
      .def_readwrite("Tp", &RandomWaterWaveTheory::Tp)
      .def_readwrite("gamma", &RandomWaterWaveTheory::gamma)
      .def_readwrite("betaJ", &RandomWaterWaveTheory::betaJ)
      .def_readwrite("f_min", &RandomWaterWaveTheory::f_min)
      .def_readwrite("f_max", &RandomWaterWaveTheory::f_max)
      .def_readwrite("df", &RandomWaterWaveTheory::df)
      .def_readwrite("reference_wavelength", &RandomWaterWaveTheory::reference_wavelength)
      .def_readonly_static("N", &RandomWaterWaveTheory::N)
      .def_readwrite("mode", &RandomWaterWaveTheory::mode)
      .def_readwrite("spectrum_type", &RandomWaterWaveTheory::spectrum_type)

      // Spectrum evaluation
      .def("spectrum", &RandomWaterWaveTheory::spectrum, py::arg("f"),
           "Evaluate S(f) at frequency f [Hz]")

      // Physical field evaluation
      .def("eta", &RandomWaterWaveTheory::eta, py::arg("X"), py::arg("t"))
      .def("phi", &RandomWaterWaveTheory::phi, py::arg("X"), py::arg("t"))
      .def("gradPhi", &RandomWaterWaveTheory::gradPhi, py::arg("X"), py::arg("t"))
      .def("gradPhi_t", &RandomWaterWaveTheory::gradPhi_t, py::arg("X"), py::arg("t"))

      // Access to discretized wave components
      .def("get_components", [](const RandomWaterWaveTheory &self) {
        std::vector<double> freqs, amps, phases, wavenumbers;
        freqs.reserve(self.N);
        amps.reserve(self.N);
        phases.reserve(self.N);
        wavenumbers.reserve(self.N);
        for (std::size_t i = 0; i < self.N; i++) {
          const auto &w = self.waves[i];
          freqs.push_back(w->w / (2.0 * M_PI));  // omega -> f [Hz]
          amps.push_back(w->A);
          phases.push_back(w->phase_shift);
          wavenumbers.push_back(w->k);
        }
        return py::dict(
            py::arg("f") = freqs,
            py::arg("A") = amps,
            py::arg("phase") = phases,
            py::arg("k") = wavenumbers
        );
      }, "Get discretized wave components: f [Hz], A [m], phase [rad], k [1/m]")

      // Vectorized spectrum evaluation
      .def("spectrum_array", [](const RandomWaterWaveTheory &self, py::array_t<double> f_array) {
        auto buf = f_array.request();
        auto *ptr = static_cast<double *>(buf.ptr);
        std::size_t n = buf.size;
        py::array_t<double> result(n);
        auto *out = static_cast<double *>(result.request().ptr);
        for (std::size_t i = 0; i < n; i++)
          out[i] = self.spectrum(ptr[i]);
        return result;
      }, py::arg("f"), "Evaluate S(f) for a numpy array of frequencies")

      // Vectorized eta evaluation
      .def("eta_array", [](const RandomWaterWaveTheory &self,
                           std::array<double, 3> X, py::array_t<double> t_array) {
        auto buf = t_array.request();
        auto *ptr = static_cast<double *>(buf.ptr);
        std::size_t n = buf.size;
        py::array_t<double> result(n);
        auto *out = static_cast<double *>(result.request().ptr);
        for (std::size_t i = 0; i < n; i++)
          out[i] = self.eta(X, ptr[i]);
        return result;
      }, py::arg("X"), py::arg("t"), "Evaluate eta for a numpy array of times")

      .def("__repr__", [](const RandomWaterWaveTheory &self) {
        std::ostringstream os;
        os << self;
        return os.str();
      });
}
