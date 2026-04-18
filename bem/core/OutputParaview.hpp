#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include "BEM_output_info.hpp"
#include "OutputCommon.hpp"
#include "Network.hpp"
#include "VPM.hpp"

namespace OutputParaView {

void mk_vtu_quadratic(const std::string &filename,
                      const V_netFp &Faces,
                      const VV_VarForOutput &VV_name_comp_mapPVd = {});

void write_shell_step(const OutputContext &ctx,
                      const std::map<std::string, outputInfo> &NetOutputInfo,
                      const std::vector<Network *> &FluidObject);

void write_aft_candidates_step(const OutputContext &ctx,
                               const std::map<std::string, outputInfo> &NetOutputInfo,
                               const std::vector<Network *> &FluidObject);

void write_step(const OutputContext &ctx,
                const std::map<std::string, outputInfo> &NetOutputInfo,
                const std::vector<Network *> &FluidObject,
                const std::vector<Network *> &RigidBodyObject,
                const std::vector<Network *> &SoftBodyObject,
                const std::unordered_set<networkFace *> &allFaces);

void write_vpm_vtp(const std::filesystem::path &path, const VortexMethod &vpm);
void write_vpm(const OutputContext &ctx, const VortexMethod &vpm, PVDWriter &pvd);

} // namespace OutputParaView
