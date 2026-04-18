#pragma once

#include "OutputCommon.hpp"
#include "Network.hpp"

namespace OutputJSON {

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
                double A_sparse_avg_nnz);

} // namespace OutputJSON
