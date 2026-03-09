#pragma once
#include "BEM.hpp"
#include "BEM_BoundaryValues.hpp"
#include "OutputCommon.hpp"

namespace OutputJSON {

namespace {

struct HydroWrench {
  Tddd force{0., 0., 0.};
  Tddd torque{0., 0., 0.};
  double area = 0.0;
};

inline bool faceActsOnBody(const networkFace *f, const Network *body) {
  if (!f || !body)
    return false;
  if (!f->Neumann)
    return false;
  const auto [p0, p1, p2] = f->getPoints();
  for (const auto *p : {p0, p1, p2}) {
    const auto effectiveFaces = getEffectiveContactFaces(p);
    const bool touches_body = std::ranges::any_of(effectiveFaces, [&](const auto *F) { return F && (F->getNetwork() == body); });
    if (!touches_body)
      return false;
  }
  return true;
}

inline HydroWrench integratePressureOnBodyFaces(const std::unordered_set<networkFace *> &allFaces, const Network *body) {
  HydroWrench out;
  if (!body)
    return out;

  // Fast, low-order surface integral:
  // - pressure is piecewise linear on each triangle -> use vertex average
  // - geometry uses precomputed face normal/centroid/area
  for (const auto *f : allFaces) {
    if (!faceActsOnBody(f, body))
      continue;

    const auto [p0, p1, p2] = f->getPoints();
    const double p = (p0->pressure + p1->pressure + p2->pressure) / 3.0;
    const Tddd df = p * f->area * f->normal;
    out.force += df;
    out.torque += Cross(f->centroid - body->COM, df);
    out.area += f->area;
  }

  return out;
}

} // namespace

void write_step(const OutputContext &ctx,
                JSONoutput &jsonout,
                const std::vector<Network *> &FluidObject,
                const std::vector<Network *> &RigidBodyObject,
                const std::vector<Network *> &SoftBodyObject,
                const std::vector<JSON> &MeasurementJSONs,
                const std::unordered_set<networkFace *> &allFaces,
                const std::vector<double> &eq_of_motion,
                double unknownsize,
                double time_setup,
                double time_solve,
                double ilu_build_time,
                double ilu_apply_time_sum,
                double gmres_iter_time_sum,
                double A_sparse_nnz,
                double A_sparse_avg_nnz) {
  // timing
  const double cpu_time_in_seconds = static_cast<double>(std::clock() - ctx.cpu_clock_start) / CLOCKS_PER_SEC;
  const double wall_clock_time_in_seconds = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - ctx.wall_clock_start).count();

  // basic
  jsonout.push("simulation_time", ctx.simulation_time);
  jsonout.push("time_step", ctx.time_step);
  jsonout.push("dt", ctx.dt);
  jsonout.push("eq_of_motion", eq_of_motion);
  jsonout.push("unknownsize", unknownsize);
  jsonout.push("time_setup", time_setup);
  jsonout.push("time_solve", time_solve);
  jsonout.push("ilu_build_time", ilu_build_time);
  jsonout.push("ilu_apply_time_sum", ilu_apply_time_sum);
  jsonout.push("gmres_iter_time_sum", gmres_iter_time_sum);
  jsonout.push("A_sparse_nnz", A_sparse_nnz);
  jsonout.push("A_sparse_avg_nnz", A_sparse_avg_nnz);
  jsonout.push("cpu_time", cpu_time_in_seconds);
  jsonout.push("wall_clock_time", wall_clock_time_in_seconds);

  // fluid mesh size
  for (const auto *water : FluidObject) {
    jsonout.push(water->getName() + "_face_size", static_cast<double>(water->getBoundaryFaces().size()));
    jsonout.push(water->getName() + "_point_size", static_cast<double>(water->getPoints().size()));
  }

  // rigid/soft basic states (safe subset)
  for (const auto *net : RigidBodyObject) {
    jsonout.push(net->getName() + "_COM", net->COM);
    jsonout.push(net->getName() + "_velocity", net->velocity);
    jsonout.push(net->getName() + "_accel", net->acceleration);
    jsonout.push(net->getName() + "_euler", Tddd{net->Q.roll(), net->Q.pitch(), net->Q.yaw()});
  }
  for (const auto *net : SoftBodyObject) {
    jsonout.push(net->getName() + "_point_size", static_cast<double>(net->getPoints().size()));
    jsonout.push(net->getName() + "_face_size", static_cast<double>(net->getBoundaryFaces().size()));
  }

  // measurements: line intersections (optional)
  for (const auto &js : MeasurementJSONs) {
    js.find("position", [&](const std::vector<std::string> &value) {
      if (value.size() != 6)
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "position must have 6 elements");
      const Tddd X0 = {std::stod(value[0]), std::stod(value[1]), std::stod(value[2])};
      const Tddd X1 = {std::stod(value[3]), std::stod(value[4]), std::stod(value[5])};
      const auto name_key = js.at("name")[0] + "_intersection";

      Tddd nearest_intersection = {1E+20, 1E+20, 1E+20};
      const Tddd mid = (X0 + X1) * 0.5;

      int intersect_fail_count = 0;
      for (auto *f : allFaces) {
        try {
          auto [isintersect, X, _] = IntersectQ_(T2Tddd{X0, X1}, (T3Tddd)(*f));
          if (isintersect && Norm(X - mid) < Norm(nearest_intersection - mid))
            nearest_intersection = X;
        } catch (const std::exception &) {
          ++intersect_fail_count; // skip degenerate/ill-conditioned faces
        }
      }
      if (intersect_fail_count > 0) {
        std::cerr << Yellow << "[OutputJSON] skipped " << intersect_fail_count
                  << " failed IntersectQ_ calls for measurement " << js.at("name")[0]
                  << colorReset << std::endl;
      }
      jsonout.push(name_key, nearest_intersection);
    });
  }

  // hydrodynamic force/moment on each rigid body (from BEM pressure on contacting fluid faces)
  for (const auto *net : RigidBodyObject) {
    if (!net)
      continue;
    const auto wrench = integratePressureOnBodyFaces(allFaces, net);
    jsonout.push(net->getName() + "_hydro_force", wrench.force);
    jsonout.push(net->getName() + "_hydro_torque", wrench.torque);
    jsonout.push(net->getName() + "_hydro_area", wrench.area);
  }

  // write result.json
  std::ofstream os(ctx.output_directory / "result.json");
  jsonout.output(os);
  os.close();
}

} // namespace OutputJSON
