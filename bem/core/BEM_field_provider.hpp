#pragma once

#include <optional>
#include <type_traits>

#include "BEM_reference_state.hpp"

namespace BEMMeshPipeline {

enum class FieldStatus {
  ok,
  not_applicable,
  reference_missing,
  fallback_failure
};

struct FieldQuery {
  const Network* water = nullptr;
  const networkPoint* p = nullptr;
  const networkLine* l = nullptr;
  const networkFace* face = nullptr;
  Tddd X = {0., 0., 0.};
  Tddd normal = {0., 0., 0.};
  NodeFaceBoundaryType bc_type = NodeFaceBoundaryType::Dirichlet;
  NodeRelocationSurface surface = NodeRelocationSurface::linear;
};

struct FieldResult {
  FieldStatus status = FieldStatus::reference_missing;
  std::optional<double> phi;
  std::optional<double> phi_t;
  std::optional<double> phin;
  std::optional<Tddd> projected_X;
};

inline double interpolateSnapshotScalar(const std::array<double, 6>& values,
                                        const Tdd& param,
                                        NodeRelocationSurface scheme) {
  if (scheme == NodeRelocationSurface::linear) {
    const double w0 = param[0];
    const double w1 = param[1];
    const double w2 = 1.0 - w0 - w1;
    return w0 * values[0] + w1 * values[1] + w2 * values[2];
  }

  const auto N = TriShape<6>(param[0], param[1]);
  double value = 0.;
  for (std::size_t i = 0; i < values.size(); ++i)
    value += N[i] * values[i];
  return value;
}

inline FieldResult evaluateSnapshotField(const SnapshotReferenceState& reference, const FieldQuery& q) {
  if (reference.snap.faces.empty())
    return {FieldStatus::reference_missing};

  auto projection = findNearestSnapshotFace(q.X, reference.snap, q.surface, reference.index.get());
  if (!projection.face)
    return {FieldStatus::reference_missing};

  return {
      FieldStatus::ok,
      interpolateSnapshotScalar(projection.face->field.phi6, projection.param, q.surface),
      interpolateSnapshotScalar(projection.face->field.phi_t6, projection.param, q.surface),
      std::nullopt,
      projection.nearest};
}

inline FieldResult evaluateInitialConditionField(const InitialConditionReferenceState& reference, const FieldQuery& q) {
  if (!q.water || !q.water->ic_phi)
    return {FieldStatus::reference_missing};
  if (q.bc_type != NodeFaceBoundaryType::Dirichlet)
    return {FieldStatus::not_applicable};

  return {
      FieldStatus::ok,
      q.water->ic_phi(q.X, reference.simulation_time),
      std::nullopt,
      std::nullopt,
      std::nullopt};
}

inline FieldResult evaluate(const ReferenceState& reference, const FieldQuery& q) {
  return std::visit([&](const auto& concrete_reference) -> FieldResult {
    using ConcreteReference = std::decay_t<decltype(concrete_reference)>;
    if constexpr (std::is_same_v<ConcreteReference, SnapshotReferenceState>)
      return evaluateSnapshotField(concrete_reference, q);
    else
      return evaluateInitialConditionField(concrete_reference, q);
  },
                    reference);
}

inline const char* fieldStatusName(FieldStatus status) {
  switch (status) {
  case FieldStatus::ok:
    return "ok";
  case FieldStatus::not_applicable:
    return "not_applicable";
  case FieldStatus::reference_missing:
    return "reference_missing";
  case FieldStatus::fallback_failure:
    return "fallback_failure";
  }
  return "unknown";
}

} // namespace BEMMeshPipeline
