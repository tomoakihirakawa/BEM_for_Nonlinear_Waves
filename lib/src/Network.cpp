#include "Network.hpp"
#include "basic_surface_geometry.hpp"
#include "pch.hpp"

#include <cmath>
#include <cstring>
#include <limits>

namespace {

using ContactFaceCandidate = std::tuple<networkFace*, Tddd, double>;

constexpr int max_contact_faces = 10;
constexpr double normal_angle_diff_detection_range = 60 * M_PI / 180;

template <int N>
double directedHausdorffDistanceImpl(const V_netFp& query_faces, const V_netFp& target_faces) {
  constexpr auto t0t1_samples = UniformPointsOnTriangle_00_10_01<N>();

  double max_distance = 0.;
  bool has_query_face = false;
  bool has_valid_target = false;

  for (auto* f : query_faces) {
    if (!f)
      continue;
    has_query_face = true;
    V_netFp filtered_target_faces;
    filtered_target_faces.reserve(target_faces.size());
    for (auto* target_f : target_faces)
      if (target_f && target_f->getNetwork() != f->getNetwork())
        filtered_target_faces.emplace_back(target_f);
    if (filtered_target_faces.empty())
      continue;
    has_valid_target = true;
    const auto X012 = ToX(f);
    for (const auto& [t0, t1] : t0t1_samples) {
      const auto x = X012[0] * t0 + X012[1] * t1 + X012[2] * (1. - t0 - t1);
      const auto x_nearest = ::Nearest(x, filtered_target_faces);
      max_distance = std::max(max_distance, Norm(x - x_nearest));
    }
  }

  if (!has_query_face)
    return 0.;
  if (!has_valid_target)
    return std::numeric_limits<double>::infinity();
  return max_distance;
}

V_netFp uniqueValidFaces(const V_netFp& faces) {
  V_netFp ret;
  ret.reserve(faces.size());
  std::unordered_set<networkFace*> seen;
  seen.reserve(faces.size());
  for (auto* f : faces) {
    if (f && seen.emplace(f).second)
      ret.emplace_back(f);
  }
  return ret;
}

V_netFp boundaryFacesFromNetworks(const V_Netp& networks) {
  V_netFp ret;
  std::unordered_set<networkFace*> seen;
  for (auto* net : networks) {
    if (!net)
      continue;
    for (auto* f : net->getBoundaryFaces()) {
      if (f && seen.emplace(f).second)
        ret.emplace_back(f);
    }
  }
  return ret;
}

// contactAcceptanceAngle is defined in Network.hpp (inline)

std::vector<ContactFaceCandidate> collectBroadPhaseCandidates(const ContactDetectable& detectable,
                                                              const std::vector<Network*>& objects,
                                                              const bool include_self_network,
                                                              const std::vector<networkFace*>& check_faces) {
  const auto& pos = detectable.getPosition();
  std::vector<ContactFaceCandidate> candidates;
  for (const auto& object : objects) {
    object->BucketSurfaces.apply(pos, detectable.contact_range, [&](networkFace* const& f) {
      if (!include_self_network && f->getNetwork() == detectable.getNetwork())
        return;
      if (std::ranges::any_of(candidates, [&](const auto& candidate) { return f == std::get<0>(candidate); }))
        return;
      const auto f_vertices = ToX(f);
      const auto Y = Nearest(pos, f_vertices);
      const auto distance = Norm(Y - pos);
      if (detectable.contact_range < distance)
        return;
      const bool any_close_normal = std::ranges::any_of(check_faces, [&](const auto& F) {
        return isFacing(F->normal, f_vertices, M_PI / 2);
      });
      if (any_close_normal)
        candidates.emplace_back(f, Y, distance);
    });
  }
  std::stable_sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return std::get<2>(a) < std::get<2>(b); });
  return candidates;
}

std::vector<ContactFaceCandidate> selectNarrowPhaseCandidates(const ContactDetectable& detectable,
                                                              const std::vector<networkFace*>& check_faces,
                                                              const std::vector<ContactFaceCandidate>& candidates) {
  std::vector<bool> keep(candidates.size(), true);

  for (size_t i = 0; i < candidates.size(); ++i) {
    if (!keep[i])
      continue;

    const auto [fi, _, distance] = candidates[i];
    const auto angle = contactAcceptanceAngle(distance, detectable.contact_range, 60., 30.);
    if (std::ranges::none_of(check_faces, [&](const auto& f) { return isFlat(fi->normal, -f->normal, angle); })) {
      keep[i] = false;
      continue;
    }

    for (size_t j = i + 1; j < candidates.size(); ++j) {
      if (!keep[j])
        continue;
      const auto [fj, __, ___] = candidates[j];
      if (isFlat(fi->normal, fj->normal, M_PI / 180))
        keep[j] = false;
    }
  }

  std::vector<ContactFaceCandidate> filtered;
  filtered.reserve(std::min(candidates.size(), static_cast<size_t>(max_contact_faces)));
  for (size_t i = 0; i < candidates.size(); ++i) {
    if (!keep[i])
      continue;
    filtered.emplace_back(candidates[i]);
    if (filtered.size() >= max_contact_faces)
      break;
  }
  return filtered;
}

void mapBoundaryFacesToNearestContactFaces(ContactDetectable& detectable,
                                           const std::vector<networkFace*>& surfaces,
                                           const std::vector<ContactFaceCandidate>& contact_faces) {
  for (const auto& f : surfaces) {
    auto& d = detectable.dof(f);
    d.contact_opponent_faces.clear();
    for (const auto& F_X : contact_faces) {
      auto* contact_face = std::get<0>(F_X);
      if (isFlat(f->normal, contact_face->normal, normal_angle_diff_detection_range) ||
          isFlat(f->normal, -contact_face->normal, normal_angle_diff_detection_range)) {
        d.contact_opponent_faces.push_back(contact_face);
      }
    }
  }
}

} // namespace

/* -------------------------------------------------------------------------- */
/*                  ContactDetectable::addContactFaces                        */
/* -------------------------------------------------------------------------- */

void ContactDetectable::addContactFaces(const std::vector<Network*>& objects, bool include_self_network) {
  auto check_faces = this->getFacesForContactCheck();
  auto surfaces = this->getBoundaryFaces();
  const auto broad_phase = collectBroadPhaseCandidates(*this, objects, include_self_network, check_faces);
  const auto narrow_phase = selectNarrowPhaseCandidates(*this, check_faces, broad_phase);
  mapBoundaryFacesToNearestContactFaces(*this, surfaces, narrow_phase);
}

/* -------------------------------------------------------------------------- */

Network::Network(const std::string& filename, const std::string& name_IN)
    : CoordinateBounds(Tddd{{0., 0., 0.}}), RigidBodyDynamics(), name(name_IN), filename(filename), BucketFaces(CoordinateBounds(Tddd{{0., 0., 0.}}), 1.), BucketSurfaces(CoordinateBounds(Tddd{{0., 0., 0.}}), 1.), BucketPoints(CoordinateBounds(Tddd{{0., 0., 0.}}), 1.), BucketTetras(CoordinateBounds(Tddd{{0., 0., 0.}}), 1.), IGNORE(false), grid_pull_depth(0), velocity_name_start({"fixed", 0.}), inputJSON(), octreeOfFaces(nullptr), octreeOfPoints(nullptr), surfaceNet(nullptr) {
  if (filename.contains(".obj") || filename.contains(".off")) {
    Load3DFile objLoader(filename);
    if (!objLoader.f_v.empty())
      this->setFaces(objLoader.f_v,
                     this->setPoints(objLoader.v)); // indexの書き換えも可能だがする必要は今のところない
    if (!objLoader.l_v.empty())
      this->setLines(objLoader.l_v,
                     this->setPoints(objLoader.v)); // indexの書き換えも可能だがする必要は今のところない

    this->displayStates();
  } else if (filename == "file_name_is_not_given") {
    // std::cout << "filename is not given" << std::endl;
  } else {
    std::stringstream ss;
    ss << "file format is not supported: " << filename;
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
  }
};

/* -------------------------------------------------------------------------- */

bool Network::MemberQ(networkPoint* const& p_IN) const { return (this->Points.find(p_IN) != this->Points.end()); };
bool Network::MemberQ(networkFace* const& f_IN) const { return (this->Faces.find(f_IN) != this->Faces.end()); };
void Network::add(netPp const p_IN) { this->Points.emplace(p_IN); };
void Network::add(netLp const l_IN) { this->Lines.emplace(l_IN); };
void Network::add(netFp const f_IN) { this->Faces.emplace(f_IN); };
void Network::add(netTp const t_IN) { this->Tetras.emplace(t_IN); };
void Network::add(const V_netPp& ps_IN) {
  for (const auto& p : ps_IN)
    this->add(p);
};
//
void Network::add(const V_netFp& fs_IN) {
  for (const auto& f : fs_IN)
    this->add(f);
};
//
bool Network::erase(networkFace* f) {
  auto it = this->Faces.find(f);
  if (it != this->Faces.end()) {
    this->Faces.erase(it);
    return true;
  } else
    return false;
};
void Network::erase(const V_netFp& fs) {
  for (auto& f : fs)
    erase(f);
};
bool Network::erase(networkPoint* p) {
  auto it = this->Points.find(p);
  if (it != this->Points.end()) {
    this->Points.erase(it);
    return true;
  } else
    return false;
};
bool Network::erase(networkLine* l) {
  auto it = this->Lines.find(l);
  if (it != this->Lines.end()) {
    this->Lines.erase(it);
    return true;
  } else
    return false;
};

// b$ =================================================================== */
// b$       ３次元空間分割に関する．　バケツ．Faces,Points,ParametricPoints
// b$ =================================================================== */

void Network::makeBuckets(const double spacing) {
  this->makeBucketPoints(spacing);
  this->makeBucketFaces(spacing);
  this->makeBucketTetras(spacing);
};

void Network::makeBucketTetras(double spacing) {
  if (spacing > 1E+10)
    spacing = this->getScale() / 10.;
  VerbosePrint(this->getName(), " makeBucketTetras(", spacing, ")");
  this->setGeometricPropertiesForce();
  this->BucketTetras.clear();
  this->BucketTetras.initialize(this->scaledBounds(expand_bounds), spacing);
  double particlize_spacing;
  std::array<double, 3> mean;
  for (const auto& t : this->getTetras()) {
    //! add extra points to make sure that the face is included in the bucket
    bool fully_inside_a_cell_found = false;
    auto [p0, p1, p2, p3] = t->Points;
    mean = (p0->X + p1->X + p2->X + p3->X) / 4.;
    auto [inserted, ijk] = this->BucketTetras.addAndGetIndices(mean, t);
    if (inserted && ::InsideQ(t->bounds, this->BucketTetras.getBounds(ijk).bounds)) {
      fully_inside_a_cell_found = true;
      break;
    }

    if (!fully_inside_a_cell_found) {
      for (auto& f : t->Faces) {
        double particlize_spacing = std::min(spacing, Min(extLength(f->getLines()))) / 20.;
        for (const auto& [xyz, t0t1] : triangleIntoPoints(f->getXVertices(), particlize_spacing))
          this->BucketTetras.add(xyz, t);
      }
    }
  }
  this->BucketTetras.setVector();
  VerbosePrint(this->getName(), " makeBucketTetras done");
};

/*面は点と違って，複数のバケツ（セル）と接することがある*/
void Network::makeBucketFaces(double spacing) {
  if (spacing > 1E+10)
    spacing = this->getScale() / 10.;
  VerbosePrint(this->getName(), " makeBucketFaces(", spacing, ")");
  this->setGeometricPropertiesForce();
  VerbosePrint("setGeometricProperties done");
  this->BucketFaces.clear();
  VerbosePrint("BucketFaces.clear done");
  this->BucketSurfaces.clear();
  VerbosePrint("BucketSurfaces.clear done");
  this->BucketFaces.initialize(this->scaledBounds(expand_bounds), spacing);
  VerbosePrint("BucketFaces.initialize done");
  this->BucketSurfaces.initialize(this->scaledBounds(expand_bounds), spacing);
  VerbosePrint("BucketSurfaces.initialize done");

  auto insert = [&spacing](networkFace* f, Buckets<networkFace*>& bucket) {
    //! add extra points to make sure that the face is included in the bucket
    auto [inserted, ijk] = bucket.addAndGetIndices(Mean(f->getXVertices()), f);
    if (inserted && ::InsideQ(f->bounds, bucket.getBounds(ijk).bounds))
      return;

    double particlize_spacing = std::min(spacing, Min(extLength(f->getLines()))) / 30.;
    for (const auto& [xyz, t0t1] : triangleIntoPoints(f->getXVertices(), particlize_spacing))
      bucket.add(xyz, f);
  };

  for (auto f : this->getFaces()) {
    insert(f, this->BucketFaces);
    if (f->BoundaryQ())
      insert(f, this->BucketSurfaces);
  }
  VerbosePrint("insert done");
  this->BucketFaces.setVector();
  VerbosePrint("BucketFaces.setVector done");
  this->BucketSurfaces.setVector();
  VerbosePrint("BucketSurfaces.setVector done");
  VerbosePrint(this->getName(), " makeBucketFaces done");
};

void Network::makeBucketPoints(double spacing) {
  TimeWatch tw;
  if (spacing > 1E+10)
    spacing = this->getScale() / 10.;
  VerbosePrint(this->getName(), " makeBucketPoints(", spacing, ")");
  this->last_makeBucketPoints_spacing = spacing;
  this->setGeometricPropertiesForce();
  this->BucketPoints.initialize_add(this->scaledBounds(expand_bounds), spacing, this->getPoints());
  VerbosePrint(green, "BucketPoints.data1D.size() = ", this->BucketPoints.data1D.size(), colorReset);
  this->BucketPoints.setVector();
  if (bem_verbose())
    std::cout << this->getName() << " makeBucketPoints done" << Blue << " Elapsed time: " << Red << tw() << colorReset << " s\n";
};

//! ------------------------------------------------------ */
//!                          接触の判別                      */
//! ------------------------------------------------------ */

std::unordered_set<networkFace*> Network::getContactFacesOfPoints() const {
  std::unordered_set<networkFace*> ret;
  for (const auto& p : this->getPoints())
    ret.insert(std::begin(p->getContactFaces()), std::end(p->getContactFaces()));
  return ret;
};

std::unordered_set<networkPoint*> Network::getContactPointsOfPoints() const {
  std::unordered_set<networkPoint*> ret;
  for (const auto& p : this->getPoints())
    ret.insert(begin(p->getContactPoints()), end(p->getContactPoints()));
  return ret;
};
std::unordered_set<networkPoint*> Network::getContactPointsOfPoints(const std::vector<Network*>& nets) const {
  /*
  getContactPointsOfPointsは，自身の保有するPointsが接した点を返す．
  Pointsの保有するmap_Net_ContactPointsから指定されたNetworkのPointsを抽出している．
  */
  std::unordered_set<networkPoint*> ret;
  for (const auto& p : this->getPoints())
    ret.insert(begin(p->getContactPoints(nets)), end(p->getContactPoints(nets)));
  return ret;
};

/* ------------------------------------------------------ */
/*                          体積の計算                      */
/* ------------------------------------------------------ */
double Network::getVolume() const {
  // ガウスの定理において，F=(x,y,z)とおいてdivF=3とすると，体積積分の結果は体積の3倍となる．
  // 面積分側は(x*nx+y*ny+z*nz)を核にした面積分と体積の3倍が等しいことになる．
  double ret = 0;
  for (const auto& f : getBoundaryFaces()) {
    auto intp = interpolationTriangleLinear0101(f->getXVertices());
    for (const auto& [x0, x1, w0w1] : __GWGW5__Tuple)
      ret += w0w1 * Dot(intp(x0, x1), intp.cross(x0, x1));
  }
  return ret / 3.;
};

double Network::getVolume(std::function<std::array<double, 3>(const networkPoint*)> func_X) const {
  double ret = 0;
  std::array<std::array<double, 3>, 3> X0X1X2;
  for (const auto& f : this->getBoundaryFaces()) {
    auto [p0, p1, p2] = f->getPoints();
    X0X1X2 = {func_X(p0), func_X(p1), func_X(p2)};
    auto intp = interpolationTriangleLinear0101(X0X1X2);
    for (const auto& [x0, x1, w0w1] : __GWGW5__Tuple)
      ret += w0w1 * Dot(intp(x0, x1), intp.cross(x0, x1));
  }
  return ret / 3.;
};

double Network::GaussIntegral2(const Tddd& X) const {
  /*
  integrate(r.n/|r^3|)

   ^^^^
   ||||
  +----+
  | 4pi|-->   outside:0
  +----+
  */
  double ret = 0;
  Tddd r;
  // for (const auto &f : this->getBoundaryFaces())
  // {
  // 	interpolationTriangleLinear0101 intp(f->getXVertices());
  // 	for (const auto &[x0, x1, w0w1] : __GWGW14__Tuple)
  // 	{
  // 		r = intp(x0, x1) - X;
  // 		ret += Dot(r / std::pow(Norm(r), 3), intp.cross(x0, x1)) * w0w1;
  // 	}
  // }

  for (const auto& f : this->getBoundaryFaces()) {

    auto [X0, X1, X2] = f->getXVertices();
    X0 -= X2;
    X1 -= X2;
    auto A = X - X2;
    double tmp = 0;
    for (const auto& [t0, t1, w0w1] : __GWGW14__Tuple)
      tmp += (1. - t0) / std::pow(Norm(X0 * t0 + X1 * t1 * (1 - t0) - A), 3) * w0w1;
    ret += Dot(A, Cross(X0, X1)) * tmp;
  }
  return ret;

  // for (const auto &f : this->getBoundaryFaces())
  // {

  // 	auto [X0, X1, X2] = f->getXVertices();
  // 	X0 -= X2;
  // 	X1 -= X2;
  // 	auto A = X - X2;
  // 	double tmp = 0;
  // 	for (const auto &[t0, t1, w0w1] : __GWGW8__Tuple)
  // 		tmp += 1. / std::pow(Norm(X0 * t0 * std::pow(1 - t0, -1 / 3.) +
  // X1 * t1 * std::pow(1 - t0, 2 / 3.) - A * std::pow(1 - t0, -1 / 3.)), 3) *
  // w0w1; 	ret += Dot(A, Cross(X0, X1)) * tmp;
  // }
  // return ret;
};

double Network::GaussIntegral(const Tddd& X) const {
  double ret = 0;
  Tddd r;
  for (const auto& f : this->getBoundaryFaces()) {
    interpolationTriangleLinear0101 intp(f->getXVertices());
    for (const auto& [x0, x1, w0w1] : __GWGW14__Tuple) {
      r = intp(x0, x1) - X;
      ret += Dot(r / std::pow(Norm(r), 3), intp.cross(x0, x1)) * w0w1;
    }
  }
  return ret;
};

double Network::windingNumber(const Tddd& X) const {
  double ret = 0;
  for (const auto& f : this->getBoundaryFaces())
    ret += SolidAngle_VanOosteromAandStrackeeJ1983(X, f->getXVertices());
  return ret / (4. * M_PI);
};

//

bool Network::InsideQ(const Tddd& X) const {
  if (!CoordinateBounds::InsideQ(X))
    return false;
  else if (this->windingNumber(X) < 0.5)
    return false;
  else
    return true;
};

std::vector<double> Network::windingNumber(const std::vector<Tddd>& Xs) const {
  std::vector<double> ret(Xs.size(), 0.);
  T3Tddd V;
  for (const auto& f : this->getBoundaryFaces()) {
    V = f->getXVertices();
    for (auto i = 0; i < Xs.size(); ++i)
      ret[i] += SolidAngle_VanOosteromAandStrackeeJ1983(Xs[i], V);
  }
  for (auto& r : ret)
    r /= (4. * M_PI);
  return ret;
};

bool Network::all_of_isInside(const std::vector<Tddd>& Xs) const {
  if (!CoordinateBounds::InsideQ(X))
    return false;
  else {
    for (auto& wn : this->windingNumber(Xs)) {
      // ひとつでも小さい値があればfalse
      if (wn < 0.75)
        return false;
    }
    return true;
  }
};

/* ------------------------------------------------------ */
/*             ネットワークの姿勢に関する変数とメソッド          */
/* ------------------------------------------------------ */
void Network::rotate(const Quaternion& Q, const Tddd& barycenter) {
  for (auto& p : this->getPoints())
    p->setXSingle(Q.Rv(p->X - barycenter) + barycenter);
  setGeometricPropertiesForce();
};
void Network::scale(const double ratio, const Tddd& center) {
  for (auto& p : this->getPoints())
    p->setXSingle(center + (p->X - center) * ratio);
  this->setGeometricPropertiesForce();
};
void Network::scale(const Tddd& ratio, const Tddd& center) {
  for (auto& p : this->getPoints())
    p->setXSingle(center + (p->X - center) * ratio);
  this->setGeometricPropertiesForce();
};
void Network::translate(const Tddd& translation) {
  for (auto& p : this->getPoints())
    p->setXSingle(p->X + translation);
  this->setGeometricPropertiesForce();
};
void Network::translateFromInitialX(const Tddd& translation) {
  for (auto& p : this->getPoints())
    p->setXSingle(p->initialX + translation);
  this->setGeometricPropertiesForce();
};
void Network::resetInitialX() {
  for (auto& p : this->getPoints())
    p->initialX = p->X;
};

/* -------------------------------------------------------------------------- */
void Network::applyTransformations(const JSON& json) {
  if (json.find("center_of_mass")) {
    std::get<0>(this->center_of_mass) = stod(json()["center_of_mass"])[0];
    std::get<1>(this->center_of_mass) = stod(json()["center_of_mass"])[1];
    std::get<2>(this->center_of_mass) = stod(json()["center_of_mass"])[2];
  }
  if (json.find("ignore")) {
    this->IGNORE = stob(json.at("ignore"))[0];
  }
  // "rotation" [x, y, z, angle] format (used by GUI and BEM_inputfile_reader)
  // Legacy "rotate" [angle, x, y, z] format also supported for backward compatibility
  if (json.find("rotation")) {
    auto rot = stod(json()["rotation"]);
    if (rot.size() >= 4)
      this->rotate(rot[3], Tddd{rot[0], rot[1], rot[2]});
  } else if (json.find("rotate")) {
    auto rotate = stod(json()["rotate"]);
    if (rotate.size() > 1)
      this->rotate(rotate[0], Tddd{rotate[1], rotate[2], rotate[3]});
  }
  if (json.find("scale")) {
    auto scale = stod(json()["scale"]);
    if (scale.size() > 1)
      this->scale({scale[0], scale[1], scale[2]});
    else
      this->scale(scale[0]);
  }
  // "translation" (used by GUI) or legacy "translate"
  {
    const char* tr_key = json.find("translation") ? "translation" : (json.find("translate") ? "translate" : nullptr);
    if (tr_key) {
      auto tr = stod(json.at(tr_key));
      if (tr.size() >= 3)
        this->translate({tr[0], tr[1], tr[2]});
      resetInitialX();
    }
  }

  if (json.find("radius"))
    for (auto& p : this->getPoints())
      p->radius = stod(json.at("radius"))[0];

  if (json.find("mass")) {
    if (json.at("mass").size() == 1) {
      std::get<2>(this->inertia) = std::get<1>(this->inertia) = std::get<0>(this->inertia) = this->mass = stod(json.at("mass"))[0];
    }
    if (json.at("mass").size() == 3) {
      std::get<0>(this->inertia) = stod(json.at("mass"))[0];
      std::get<1>(this->inertia) = stod(json.at("mass"))[1];
      std::get<2>(this->inertia) = stod(json.at("mass"))[2];
    }
  }

  if (json.find("MOI")) {
    std::get<3>(this->inertia) = stod(json.at("MOI"))[0];
    std::get<4>(this->inertia) = stod(json.at("MOI"))[1];
    std::get<5>(this->inertia) = stod(json.at("MOI"))[2];
  }
  if (json.find("COM"))
    this->COM = this->initial_center_of_mass = Tddd{stod(json.at("COM"))[0], stod(json.at("COM"))[1], stod(json.at("COM"))[2]};

  // Set isFixed and adjust inertia accordingly
  // 後方互換: "fixed" キーも受け入れる（"isFixed" を優先）
  {
    const char* fixed_key = json.find("isFixed") ? "isFixed" : (json.find("fixed") ? "fixed" : nullptr);
    if (fixed_key) {
      auto is_fixed = stob(json.at(fixed_key));
      if (is_fixed.size() == 1) {
        this->isFixed.fill(is_fixed[0]);
        if (is_fixed[0])
          this->inertia.fill(1E+20);
      } else {
        for (auto i = 0; i < is_fixed.size(); ++i) {
          this->isFixed[i] = is_fixed[i];
          if (is_fixed[i])
            this->inertia[i] = 1E+20;
        }
      }
    }
  }

  // if (json.find("remesh"))
  // {
  // 	auto minlen = stod(json()["remesh"]);
  // 	if (minlen.size() > 0)
  // 		remesh(&this, minlen[0]);
  // }
  // if (json.find("coarsen"))
  // {
  // 	auto minlen = stod(json()["coarsen"]);
  // 	if (minlen.size() > 0)
  // 		coarsen(&this, minlen[0]);
  // }
  // if (json.find("reverseNormal")) {
  //    std::string TorF = json()["reverseNormal"][0];
  //    if (TorF.compare("True") == 0 || TorF.compare("true") == 0 ||
  //    TorF.compare("1") == 0) {
  //       this->reverseNormal();
  //       std::cout << "reverse done" << std::endl;
  //    }
  // }
};

/* -------------------------------------------------------------------------- */

const std::unordered_set<networkPoint*>& Network::getPoints() const { return this->Points; };

std::vector<networkPoint*> Network::getBoundaryPoints() const {
  std::unordered_set<networkPoint*> uniquePoints;
  for (const auto& f : this->getBoundaryFaces())
    for (const auto& p : f->getPoints())
      uniquePoints.emplace(p);
  std::vector<networkPoint*> ret(uniquePoints.begin(), uniquePoints.end());
  std::sort(ret.begin(), ret.end());
  return ret;
}

/* ------------------------------ 境界条件に関する設定関数 ------------------------------ */

std::unordered_set<networkPoint*> Network::getPointsCORNER() const {
  std::unordered_set<networkPoint*> ret;
  ret.reserve(this->Points.size());
  for (const auto& p : this->Points)
    if (p->CORNER)
      ret.emplace(p);
  return ret;
};
std::unordered_set<networkPoint*> Network::getPointsNeumann() const {
  std::unordered_set<networkPoint*> ret;
  ret.reserve(this->Points.size());
  for (const auto& p : this->Points)
    if (!p->getFacesNeumann().empty())
      ret.emplace(p);
  return ret;
};
std::unordered_set<networkPoint*> Network::getPointsDirichlet() const {
  std::unordered_set<networkPoint*> ret;
  ret.reserve(this->Points.size());
  for (const auto& p : this->Points)
    if (!p->getFacesDirichlet().empty())
      ret.emplace(p);
  return ret;
};

void Network::setMinDepthFromCORNER() {
  std::vector<networkPoint*> points = ToVector(this->getPoints());
  for (const auto& p : points) {
    if (p->CORNER)
      p->minDepthFromCORNER_ = p->minDepthFromCORNER = 0;
    else
      p->minDepthFromCORNER_ = p->minDepthFromCORNER = 100000000;
    //
    if (p->isMultipleNode)
      p->minDepthFromMultipleNode_ = p->minDepthFromMultipleNode = 0;
    else
      p->minDepthFromMultipleNode_ = p->minDepthFromMultipleNode = 100000000;
  }

  for (auto i = 0; i < 100; i++) {
    // #pragma omp parallel
    for (auto& p : points)
      // #pragma omp single nowait
      for (auto& q : p->getNeighbors()) {
        p->minDepthFromCORNER_ = std::min(p->minDepthFromCORNER_, q->minDepthFromCORNER + 1);
        p->minDepthFromMultipleNode_ = std::min(p->minDepthFromMultipleNode_, q->minDepthFromMultipleNode + 1);
      }
    // apply
    // #pragma omp parallel
    for (const auto& p : points)
    // #pragma omp single nowait
    {
      p->minDepthFromCORNER = p->minDepthFromCORNER_;
      p->minDepthFromMultipleNode = p->minDepthFromMultipleNode_;
    }
  }
};

/* --------------------------------------------------------------------------*/

const std::unordered_set<networkFace*>& Network::getFaces() const { return this->Faces; };
const std::vector<networkFace*>& Network::getFacesVector() const { return this->Faces_vector; }
const std::vector<networkPoint*>& Network::getPointsVector() const { return this->Points_vector; }
// Update Faces_vector based on the contents of Faces
void Network::setFacesVector() {
  this->Faces_vector.clear();
  this->Faces_vector.reserve(this->Faces.size());
  this->Faces_vector.assign(this->Faces.begin(), this->Faces.end());
}
// Update Points_vector based on the contents of Points
void Network::setPointsVector() {
  this->Points_vector.clear();
  this->Points_vector.reserve(this->Points.size());
  this->Points_vector.assign(this->Points.begin(), this->Points.end());
}

const std::unordered_set<networkTetra*>& Network::getTetras() const { return this->Tetras; };

const std::vector<networkTetra*> Network::getIsolatedTetras() const {
  std::vector<networkTetra*> ret;
  ret.reserve(this->Tetras.size());
  for (const auto& t : this->Tetras)
    if (t->IsolatedQ())
      ret.emplace_back(t);
  return ret;
};

std::unordered_set<networkPoint*> Network::getParametricPoints() const {
  std::unordered_set<networkPoint*> ret, tmp;
  for (const auto& f : this->getFaces()) {
    tmp = f->getParametricPoints();
    ret.insert(tmp.begin(), tmp.end());
  }
  return ret;
};
Tddd Network::getMeanX() const {
  Tddd ret = {0, 0, 0};
  for (const auto& p : this->Points)
    ret += p->X;
  return ret / this->Points.size();
};

void Network::getLocations(VV_d& ret) const {
  ret.resize(Points.size());
  int i(0);
  // for (const auto &p : Points)
  // 	ret[i++] = p->xyz;
  for (const auto& p : Points)
    ret[i++] = {std::get<0>(p->X), std::get<1>(p->X), std::get<2>(p->X)};
};

std::unordered_set<networkLine*> Network::getLines() const {
  std::unordered_set<networkLine*> ret;
  ret.reserve(this->Points.size() * 3);
  for (auto p : this->Points)
    for (auto l : p->getLines())
      ret.emplace(l);
  return ret;
};

std::vector<networkFace*> Network::getBoundaryFaces() const {
  std::vector<networkFace*> surfaces;
  surfaces.reserve(this->Faces.size());
  for (const auto& f : this->Faces)
    if (f->BoundaryQ())
      surfaces.emplace_back(f);
  return surfaces;
};

std::vector<networkLine*> Network::getBoundaryLines() const {
  std::unordered_set<networkLine*> ret;
  ret.reserve(this->Faces.size() * 3);
  for (const auto& f : this->getBoundaryFaces()) {
    auto [l0, l1, l2] = f->getLines();
    ret.emplace(l0);
    ret.emplace(l1);
    ret.emplace(l2);
  }
  return std::vector<networkLine*>(ret.begin(), ret.end());
};

std::vector<networkFace*> Network::getComputationalBoundaryFaces() const {
  std::vector<networkFace*> surfaces;
  surfaces.reserve(this->Faces.size());
  for (const auto& f : this->Faces)
    if (f->ComputationalBoundaryQ())
      surfaces.emplace_back(f);
  return surfaces;
};

std::vector<networkLine*> Network::getComputationalBoundaryLines() const {
  std::unordered_set<networkLine*> ret;
  ret.reserve(this->Faces.size() * 3);
  for (const auto& f : this->getComputationalBoundaryFaces()) {
    auto [l0, l1, l2] = f->getLines();
    ret.emplace(l0);
    ret.emplace(l1);
    ret.emplace(l2);
  }
  return std::vector<networkLine*>(ret.begin(), ret.end());
};

std::vector<networkPoint*> Network::getComputationalBoundaryPoints() const {
  std::unordered_set<networkPoint*> uniquePoints;
  for (const auto& f : this->getComputationalBoundaryFaces())
    for (const auto& p : f->getPoints())
      uniquePoints.emplace(p);
  std::vector<networkPoint*> ret(uniquePoints.begin(), uniquePoints.end());
  std::sort(ret.begin(), ret.end());
  return ret;
};

std::vector<networkFace*> Network::getInteriorFaces() const {
  std::vector<networkFace*> ret;
  ret.reserve(this->Faces.size());
  for (const auto& f : this->Faces)
    if (!f->BoundaryQ())
      ret.emplace_back(f);
  return ret;
};

std::vector<networkLine*> Network::getInteriorLines() const {
  std::unordered_set<networkLine*> ret;
  ret.reserve(this->Lines.size());
  for (const auto& l : this->getLines())
    if (std::ranges::none_of(l->getFaces(), [](networkFace* f) { return f->BoundaryQ(); }))
      ret.emplace(l);
  return std::vector<networkLine*>(ret.begin(), ret.end());
};

std::vector<networkPoint*> Network::getInteriorPoints() const {
  std::vector<networkPoint*> ret;
  ret.reserve(this->Points.size());
  for (const auto& p : this->Points)
    if (std::ranges::none_of(p->getFaces(), [](networkFace* f) { return f->BoundaryQ(); }))
      ret.emplace_back(p);
  return ret;
};

std::unordered_set<networkLine*> Network::getLinesUO() const {
  std::unordered_set<networkLine*> ret;
  ret.reserve(3 * this->Faces.size());
  for (const auto& f : this->Faces) {
    auto [l0, l1, l2] = f->getLines();
    ret.emplace(l0);
    ret.emplace(l1);
    ret.emplace(l2);
  }
  return ret;
};

V_netLp Network::getLinesIntxn() const {
  V_netLp ret(0);
  for (const auto& p : this->getPoints())
    for (const auto& l : p->Lines)
      if (l->isIntxn())
        ret.emplace_back(l);
  return DeleteDuplicates(ret);
};

/* -------------------------------------------------------------------------- */

namespace {
inline void hash_combine(std::size_t& seed, std::size_t value) { seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2); }
} // namespace

std::size_t Network::computeTopologySignature() const {
  // true_quadratic では Face->Lines の対応順も幾何に効くため、
  // Face->Points だけでなく Face->Lines / Line->Points も含める。
  // サイズは3コンテナすべてチェック（追加/削除検出）。
  std::size_t h = 0;
  hash_combine(h, this->Points.size());
  hash_combine(h, this->Faces.size());
  hash_combine(h, this->Lines.size());

  for (const auto& f : this->Faces) {
    hash_combine(h, reinterpret_cast<std::uintptr_t>(f));
    auto [p0, p1, p2] = f->getPoints();
    hash_combine(h, reinterpret_cast<std::uintptr_t>(p0));
    hash_combine(h, reinterpret_cast<std::uintptr_t>(p1));
    hash_combine(h, reinterpret_cast<std::uintptr_t>(p2));
    auto [l0, l1, l2] = f->getLines();
    hash_combine(h, reinterpret_cast<std::uintptr_t>(l0));
    hash_combine(h, reinterpret_cast<std::uintptr_t>(l1));
    hash_combine(h, reinterpret_cast<std::uintptr_t>(l2));
  }

  for (const auto& l : this->Lines) {
    hash_combine(h, reinterpret_cast<std::uintptr_t>(l));
    auto [p0, p1] = l->getPoints();
    hash_combine(h, reinterpret_cast<std::uintptr_t>(p0));
    hash_combine(h, reinterpret_cast<std::uintptr_t>(p1));
  }
  return h;
}

std::size_t Network::computeGeometrySignature(double /*inv_eps*/) const {
  // ビットパターン直接ハッシュ: quantize (long double + llround) より高速。
  // 変更検出は現行より保守的（ビット完全一致、1pm丸め許容なし）。
  // -0.0 と +0.0 を同一視するため正規化する。
  std::size_t h = 0;
  std::uint64_t bits;
  for (const auto& p : this->Points) {
    const auto& X = p->X;
    auto hash_coord = [&](double v) {
      std::memcpy(&bits, &v, sizeof(double));
      if (bits == 0x8000000000000000ULL)
        bits = 0; // -0 → +0
      hash_combine(h, bits);
    };
    hash_coord(std::get<0>(X));
    hash_coord(std::get<1>(X));
    hash_coord(std::get<2>(X));
  }
  for (const auto& l : this->Lines) {
    auto hash_vec = [&](const Tddd& X) {
      auto hash_coord = [&](double v) {
        std::memcpy(&bits, &v, sizeof(double));
        if (bits == 0x8000000000000000ULL)
          bits = 0;
        hash_combine(h, bits);
      };
      hash_coord(std::get<0>(X));
      hash_coord(std::get<1>(X));
      hash_coord(std::get<2>(X));
    };
    hash_vec(l->X_mid);
    hash_vec(l->corner_midpoint_offset);
  }
  return h;
}

std::size_t Network::topologySignature() const { return computeTopologySignature(); }

std::size_t Network::geometrySignature(double inv_eps) const { return computeGeometrySignature(inv_eps); }

void Network::assignPointFaceIndices() {
  int count = 0;
  for (const auto& p : this->getPoints())
    p->index = count++;
  count = 0;
  for (const auto& f : this->getFaces())
    f->index = count++;
}

bool Network::setGeometricProperties() {
  const std::size_t topo = computeTopologySignature();
  const std::size_t geom = computeGeometrySignature();
  if (geometry_signature_valid && topo == topology_signature_cached && geom == geometry_signature_cached) {
    assignPointFaceIndices();
    return false;
  }
  this->setGeometricPropertiesImpl();
  topology_signature_cached = topo;
  geometry_signature_cached = geom;
  geometry_signature_valid = true;
  return true;
}

void Network::setGeometricPropertiesImpl() {
  try {
    if (!this->getPoints().empty()) {
      for (const auto& l : this->getLines())
        l->setBoundsSingle();

      for (const auto& p : this->getPoints())
        p->setFaces(); // % point->Lines must have been determined

      // Surfacesの法線方向は外向きに揃える：テトラとは逆向きにする
      // for (const auto &f : this->getBoundaryFaces()) {
      //   auto [t0, t1] = f->getTetras();
      //   auto [p0, p1, p2] = f->getPoints();
      //   if (t0 != nullptr) {
      //     auto X = Mean(t0->vertices);
      //     if (CrossDot(p1->X - p0->X, p2->X - p0->X, X - (p0->X + p1->X + p2->X) / 3.) > 0.) {
      //       std::cout << "reversing face normal " << f << std::endl;
      //       f->reverse();
      //     }
      //   } else if (t1 != nullptr) {
      //     auto X = Mean(t1->vertices);
      //     if (CrossDot(p1->X - p0->X, p2->X - p0->X, X - (p0->X + p1->X + p2->X) / 3.) > 0.) {
      //       std::cout << "reversing face normal " << f << std::endl;
      //       f->reverse();
      //     }
      //   }
      // }
      //

      for (const auto& f : this->getFaces())
        f->setGeometricProperties(ToX(f->getPoints())); // @ face->Lines must have been determined

      for (const auto& f : this->getBoundaryFaces())
        f->setDodecaPoints();

      for (const auto& t : this->getTetras())
        t->setProperties(ToX(t->Points));

      CoordinateBounds::setBounds(ToX(this->getPoints()));
      assignPointFaceIndices();

    } else {
      CoordinateBounds::setBounds(Tddd{{0., 0., 0.}});
    }
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Error in Network::setGeometricPropertiesForce(): " + std::string(e.what()));
  }
};

void Network::setGeometricPropertiesForce() {
  this->setGeometricPropertiesImpl();
  topology_signature_cached = computeTopologySignature();
  geometry_signature_cached = computeGeometrySignature();
  geometry_signature_valid = true;
};

std::array<bool, 4> Network::validateConectivity() {
  std::array<bool, 4> ret = {true, true, true, true};
  auto& validPoints = std::get<0>(ret);
  auto& validLines = std::get<1>(ret);
  auto& validFaces = std::get<2>(ret);
  auto& validTetras = std::get<3>(ret);
  // Check Points

  std::cout << Yellow << "validating points - lines";
  for (const auto& p : this->getPoints()) {
    for (const auto& l : p->getLines())
      if (std::find(l->getPoints().begin(), l->getPoints().end(), p) == l->getPoints().end()) {
        validPoints = false;
        break;
      }
  }

  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  std::cout << Yellow << "validating points - faces";
  for (const auto& p : this->getPoints()) {
    for (const auto& f : p->getFaces()) {
      if (std::find(f->getPoints().begin(), f->getPoints().end(), p) == f->getPoints().end()) {
        validPoints = false;
        break;
      }
    }
  }

  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  std::cout << Yellow << "validating points - tetras";
  for (const auto& p : this->getPoints()) {
    for (const auto& t : p->Tetras) {
      if (std::find(t->Points.begin(), t->Points.end(), p) == t->Points.end()) {
        validPoints = false;
        break;
      }
    }
  }

  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  /* --------------------------------------------------------------------------
   */

  std::cout << Green << "validating lines - points";
  // Check Lines
  for (const auto& l : this->getLines()) {
    for (const auto& p : l->getPoints()) {
      if (std::find(p->getLines().begin(), p->getLines().end(), l) == p->getLines().end()) {
        validLines = false;
        break;
      }
    }
  }

  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  std::cout << Green << "validating lines - faces";
  for (const auto& l : this->getLines()) {
    for (const auto& f : l->getFaces()) {
      if (std::find(f->getLines().begin(), f->getLines().end(), l) == f->getLines().end()) {
        validLines = false;
        break;
      }
    }
  }

  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  std::cout << Green << "validating lines - tetras";
  for (const auto& l : this->getLines()) {
    for (const auto& t : l->Tetras) {
      if (std::find(t->Lines.begin(), t->Lines.end(), l) == t->Lines.end()) {
        validLines = false;
        break;
      }
    }
  }

  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  /* --------------------------------------------------------------------------
   */

  // Check Faces
  std::cout << Magenta << "validating faces - points";
  for (const auto& f : this->getFaces()) {
    for (const auto& p : f->getPoints()) {
      if (std::find(p->getFaces().begin(), p->getFaces().end(), f) == p->getFaces().end()) {
        validFaces = false;
        break;
      }
    }
  }
  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  std::cout << Magenta << "validating faces - lines";
  for (const auto& f : this->getFaces()) {
    for (const auto& l : f->getLines()) {
      if (std::find(l->getFaces().begin(), l->getFaces().end(), f) == l->getFaces().end()) {
        validFaces = false;
        break;
      }
    }
  }

  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  std::cout << Magenta << "validating faces - tetras" << colorReset;
  for (const auto& f : this->getFaces()) {
    for (const auto& t : f->Tetras)
      if (t != nullptr) {
        if (std::find(t->Faces.begin(), t->Faces.end(), f) == t->Faces.end()) {
          validFaces = false;
          break;
        }
      }
  }

  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  /* --------------------------------------------------------------------------
   */

  std::cout << "validating tetras - points";
  for (const auto& t : this->Tetras) {
    for (const auto& p : t->Points) {
      if (std::find(p->Tetras.begin(), p->Tetras.end(), t) == p->Tetras.end()) {
        validTetras = false;
        break;
      }
    }
  }

  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  std::cout << "validating tetras - lines";
  for (const auto& t : this->Tetras) {
    for (const auto& l : t->Lines) {
      if (std::find(l->Tetras.begin(), l->Tetras.end(), t) == l->Tetras.end()) {
        validTetras = false;
        break;
      }
    }
  }

  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  std::cout << "validating tetras - faces";
  for (const auto& t : this->Tetras) {
    for (const auto& f : t->Faces) {
      if (std::find(f->Tetras.begin(), f->Tetras.end(), t) == f->Tetras.end()) {
        validTetras = false;
        break;
      }
    }
  }

  std::cout << Red << "{validPoints, validLines, validFaces, validTetras} = " << ret << colorReset << std::endl;

  return ret;
}

/* -------------------------------------------------------------------------- */

V_netPp Network::setPoints(const std::vector<Tddd>& v_IN) {
  std::unordered_map<std::shared_ptr<Tddd>, int> map_Tddd_Int;
  std::vector<Tddd> vecTddd(v_IN.size());
  map_Tddd_Int.reserve(v_IN.size());
  int n = 0;
  for (const auto& v : v_IN) {
    std::shared_ptr<Tddd> x(new Tddd(v));
    vecTddd[n] = *x;
    map_Tddd_Int[x] = n++;
  }
  CoordinateBounds tmp(MinMaxColumns(vecTddd));
  CoordinateBounds bound(tmp.scaledBounds(1.11111111111));

  Buckets<std::shared_ptr<Tddd>> bucket(bound, bound.getScale() / 20.);
  std::vector<networkPoint*> ret(v_IN.size());
  int overlap_index = 0, overlaps = 0;
  for (const auto& [x, n] : map_Tddd_Int) {
    auto [i, j, k] = bucket.indices(*x);
    overlap_index = 0;
    for (const auto& tdd : bucket.data[i][j][k])
      if (Norm(*tdd - *x) < 1E-10) {
        overlap_index = map_Tddd_Int[tdd];
        break;
      }
    if (!overlap_index) {
      bucket.add(*x, x);
      ret[n] = new networkPoint(this, *x);
    } else {
      ret[n] = ret[overlap_index];
      overlaps++;
    }
  }
  return ret;
};

V_netPp Network::setPoints(const VV_d& v_IN) {
  std::unordered_map<std::shared_ptr<Tddd>, int> map_Tddd_Int;
  std::vector<Tddd> vecTddd(v_IN.size());
  map_Tddd_Int.reserve(v_IN.size());
  int n = 0;
  for (const auto& v : v_IN) {
    if (v.size() != 3)
      throw std::runtime_error("v.size() != 3");
    std::shared_ptr<Tddd> x(new Tddd({v[0], v[1], v[2]}));
    vecTddd[n] = *x;
    map_Tddd_Int[x] = n++;
  }
  CoordinateBounds tmp(MinMaxColumns(vecTddd));
  CoordinateBounds bound(tmp.scaledBounds(1.111111111111));

  Buckets<std::shared_ptr<Tddd>> bucket(bound, bound.getScale() / 20.);
  std::vector<networkPoint*> ret(v_IN.size());
  int overlap_index = 0, overlaps = 0;
  for (const auto& [x, n] : map_Tddd_Int) {
    auto [i, j, k] = bucket.indices(*x);
    overlap_index = 0;
    for (const auto& tdd : bucket.data[i][j][k])
      if (Norm(*tdd - *x) < 1E-10) {
        overlap_index = map_Tddd_Int[tdd];
        break;
      }
    if (!overlap_index) {
      bucket.add(*x, x);
      ret[n] = new networkPoint(this, *x);
    } else {
      ret[n] = ret[overlap_index];
      overlaps++;
    }
  }
  // std::cout << "overlaped nodes : " << overlaps << std::endl;
  // std::cout << "setPoints elapsed time :" << timer() << std::endl;
  return ret;
};

/* -------------------------------------------------------------------------- */

void Network::setFaces(const std::vector<T3Tddd>& v_IN, const double resolution) {
  std::vector<std::array<networkPoint*, 3>> v_Points(v_IN.size());
  CoordinateBounds tmp(v_IN);
  CoordinateBounds bound(tmp.scaledBounds(1.11111111111));
  Buckets<networkPoint*, 0 /*一時的なものでFMMは使わないのでこれでいい*/> bucket;
  bucket.initialize(bound, bound.getScale() / 20.);

  auto findOverlap = [&](const Tddd& X) -> networkPoint* {
    for (const auto& p : bucket.getData(X))
      if (Norm(p->X - X) < resolution)
        return p;
    return nullptr;
  };

  networkPoint* p = nullptr;
  int i = 0;

  for (const auto& position : v_IN) {
    int j = 0;
    for (const auto& X : position) {
      auto p = findOverlap(X);
      if (p == nullptr) {
        p = v_Points[i][j] = new networkPoint(this, X);
        bucket.add(X, p);
      } else
        v_Points[i][j] = p;
      j++;
    }
    i++;
  }

  for (const auto& [p0, p1, p2] : v_Points)
    new networkFace(this, p0, p1, p2);

  this->setGeometricPropertiesForce();
};

/* -------------------------------------------------------------------------- */

void Network::setFaces(const VV_i& f_v, const V_netPp& points) {
  try {
    int count = 0;
    for (const auto& index : f_v) {
      if (points[index[0]] && points[index[1]] && points[index[2]])
        new networkFace(this, points[index[0]], link(points[index[0]], points[index[1]], this), points[index[1]], link(points[index[1]], points[index[2]], this), points[index[2]], link(points[index[2]], points[index[0]], this));
      else {
        std::stringstream ss;
        ss << "index = " << index << std::endl;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      }
    }
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
  this->setGeometricPropertiesForce();
};

/* -------------------------------------------------------------------------- */

void Network::setLines(const VV_i& l_v, const V_netPp& points) {
  std::cout << "l_v " << l_v << std::endl;
  try {
    for (const auto& index : l_v) {
      for (auto i = 0; i < index.size() - 1; i++) {
        if (points[index[i]] && points[index[i + 1]])
          new networkLine(this, points[index[i]], points[index[i + 1]]);
        else {
          std::stringstream ss;
          ss << "index = " << index << std::endl;
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
        }
      }
    }
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
  this->setGeometricPropertiesForce();
};

/* -------------------------------------------------------------------------- */

void Network::mergeObjFile(const std::string& filename, int comp_id) {
  Load3DFile loader(filename);
  if (loader.f_v.empty())
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                        "mergeObjFile: no faces in " + filename + " (face mesh required for multi-part rigid body)");

  // Record existing faces before adding new ones
  auto faces_before = this->Faces;

  // Add points (shared by both faces and lines from this obj)
  auto points = this->setPoints(loader.v);

  // Add faces
  this->setFaces(loader.f_v, points);

  // Add lines if present
  if (!loader.l_v.empty())
    this->setLines(loader.l_v, points);

  // Label newly added faces with component_id
  for (auto* f : this->Faces)
    if (!faces_before.contains(f))
      f->component_id = comp_id;

  this->displayStates();
};

/* -------------------------------------------------------------------------- */

void Network::DeleteParticles() {
  auto points = this->Points;
  for (const auto& p : points)
    if (std::get<0>(p->particlize_info))
      p->Delete();
};

void Network::DeleteIsolatedTetras() {
  std::vector<networkTetra*> tetras_to_delete = this->getIsolatedTetras();
  std::vector<networkFace*> faces_to_delete;
  for (auto& t : tetras_to_delete)
    for (const auto& f : t->Faces)
      faces_to_delete.push_back(f);
  for (const auto& t : tetras_to_delete)
    delete t;
  for (const auto& f : faces_to_delete)
    delete f;
  std::vector<networkLine*> lines_to_delete;
  for (const auto& l : this->getBoundaryLines())
    if (l->getFaces().empty())
      lines_to_delete.push_back(l);
  for (const auto& l : lines_to_delete)
    delete l;
  std::vector<networkPoint*> points_to_delete;
  for (const auto& p : this->getBoundaryPoints())
    if (p->getBoundaryLines().empty())
      points_to_delete.push_back(p);
  for (const auto& p : points_to_delete)
    delete p;
};

void Network::DeleteInteriorTetras() {
  std::vector<networkLine*> boundary_lines = this->getBoundaryLines();
  std::vector<networkFace*> inner_faces = this->getInteriorFaces();
  std::vector<networkLine*> inner_lines = this->getInteriorLines();
  std::vector<networkPoint*> inner_points = this->getInteriorPoints();
  // 4. 削除実行
  // 依存関係の順序: Face -> Line -> Point の順が良い（参照カウントなどがある場合）
  // ただし Networkクラスの delete はポインタ管理を適切に行う前提

  for (const auto& f : inner_faces)
    delete f;

  // Face削除後にLineのトポロジーが変わる可能性があるためチェック
  for (const auto& l : boundary_lines)
    if (!l->checkTopology())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "topology error after deleting inner faces");

  for (const auto& l : inner_lines)
    delete l;

  for (const auto& l : boundary_lines)
    if (!l->checkTopology())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "topology error after deleting inner lines");

  for (const auto& p : inner_points)
    delete p;

  for (const auto& l : boundary_lines)
    if (!l->checkTopology())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "topology error after deleting inner points");

  // 5. 残存確認（念のため）
  auto tetras = this->getTetras();
  if (!tetras.empty()) {
    std::cout << "tetrahedra are not empty after deleting inner lines and faces!" << std::endl;
    for (const auto& t : tetras)
      delete t;
  }
};

/* -------------------------------------------------------------------------- */

void Network::checkConnectivity() const {
  const auto& points = this->getPoints();
  const auto& faces = this->getFaces();
  const auto& tetras = this->getTetras();
  const auto& lines = this->getLines();

  for (auto* p : points) {
    if (!p)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "null point found");
    if (p->getNetwork() != this)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point-network mismatch");

    for (auto* l : p->Lines)
      if (l == nullptr || l->Point_A != p && l->Point_B != p)
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "line/point adjacency broken");
    for (auto* f : p->Faces)
      if (f == nullptr || std::ranges::none_of(f->Points, [p](auto* pf) { return pf == p; }))
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "face/point adjacency broken");
    for (auto* t : p->Tetras)
      if (t == nullptr || std::ranges::none_of(t->Points, [p](auto* pt) { return pt == p; }))
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "tetra/point adjacency broken");
  }

  for (auto* f : faces) {
    if (!f)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "null face found");
    if (f->getNetwork() != this)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "face-network mismatch");
    for (auto* l : f->Lines)
      if (l == nullptr || std::ranges::none_of(l->Faces, [f](auto* fl) { return fl == f; })) {
        std::cout << "f : " << f << std::endl;
        std::cout << "l : " << l << std::endl;
        std::cout << "l->Faces : " << l->getFaces() << std::endl;
        for (auto l : f->Lines)
          std::cout << "f->Lines : " << l << std::endl;

        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "line/face adjacency broken");
      }
    for (auto* p : f->Points)
      if (p == nullptr || std::ranges::none_of(p->Faces, [f](auto* fp) { return fp == f; }))
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point/face adjacency broken");
    for (auto* t : f->Tetras)
      if (t && std::ranges::none_of(t->Faces, [f](auto* ft) { return ft == f; }))
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "tetra/face adjacency broken");
    for (auto* t : f->Tetras)
      if (t && std::ranges::all_of(t->Faces, [f](auto* ft) { return ft == nullptr; }))
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "tetra/face adjacency broken");
  }

  for (auto* t : tetras) {
    if (!t)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "null tetra found");
    for (auto* f : t->Faces)
      if (f == nullptr || std::ranges::none_of(f->Tetras, [t](auto* tf) { return tf == t; }))
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "face/tetra adjacency broken");
    for (auto* l : t->Lines)
      if (l == nullptr || std::ranges::none_of(l->Tetras, [t](auto* tl) { return tl == t; }))
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "line/tetra adjacency broken");
    for (auto* p : t->Points)
      if (p == nullptr || std::ranges::none_of(p->Tetras, [t](auto* tp) { return tp == t; }))
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point/tetra adjacency broken");
  }

  for (auto* l : lines) {
    if (!l)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "null line found");
    if (l->getNetwork() != this)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "line-network mismatch");

    for (auto* p : l->getPoints())
      if (p == nullptr || std::ranges::none_of(p->Lines, [l](auto* pl) { return pl == l; }))
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point/line adjacency broken");
    for (auto* f : l->getFaces())
      if (f == nullptr || std::ranges::none_of(f->Lines, [l](auto* fl) { return fl == l; })) {
        std::cout << "f : " << f << std::endl;
        std::cout << "l : " << l << std::endl;
        std::cout << "l->Faces : " << l->getFaces() << std::endl;
        for (auto l : f->Lines)
          std::cout << "f->Lines : " << l << std::endl;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "face/line adjacency broken");
      }
    for (auto* t : l->Tetras)
      if (t == nullptr || std::ranges::none_of(t->Lines, [l](auto* tl) { return tl == l; }))
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "tetra/line adjacency broken");
  }
}

bool Network::erase_element(networkPoint* const p) {
  auto it = this->Points.find(p); // unordered_setの場合はfindを使用
  if (it != this->Points.end()) {
    this->Points.erase(it);
    return true;
  }
  return false;
}

bool Network::erase_element(networkLine* const l) {
  auto it = this->Lines.find(l); // unordered_setの場合はfindを使用
  if (it != this->Lines.end()) {
    this->Lines.erase(it);
    return true;
  }
  return false;
}

bool Network::erase_element(networkFace* const f) {
  auto it = this->Faces.find(f); // unordered_setの場合はfindを使用
  if (it != this->Faces.end()) {
    this->Faces.erase(it);
    return true;
  }
  return false;
}

bool Network::erase_element(networkTetra* const t) {
  auto it = this->Tetras.find(t); // unordered_setの場合はfindを使用
  if (it != this->Tetras.end()) {
    this->Tetras.erase(it);
    return true;
  }
  return false;
}

/* -------------------------------------------------------------------------- */

size_t Network::improveTetrahedraDelaunay(const double eps) {
  if (this->Tetras.size() == 0)
    return 0;
  size_t total_flips = 0;
  const int max_iter = 10;
  bool flipped_this_cycle = false; // 修正点 1: フリップ実行フラグ

  auto func1 = [](double before, double after) { return after / before - 1.; };
  auto func2 = [](double before, double after) { return after - before; };

  for (auto iter = 0; iter < max_iter; iter++) {

    flipped_this_cycle = false; // 毎サイクルリセット

    Flip23Candidate best23 = improveTetrahedraDelaunay23(eps, iter < max_iter * 0.1 ? func1 : func2);
    if (iter > max_iter * 0.2)
      best23.smaller_is_better = 10.;
    Flip32Candidate best32 = improveTetrahedraDelaunay32(eps, iter < max_iter * 0.1 ? func1 : func2);

    // 候補が「実際に改善するか」を先に判定
    bool improve23 = (best23.f && best23.smaller_is_better < -eps);
    bool improve32 = (best32.edge && best32.smaller_is_better < -eps);

    if (improve23 && !improve32) {
      improve(best23);
      flipped_this_cycle = true;
    } else if (!improve23 && improve32) {
      improve(best32);
      flipped_this_cycle = true;
    } else if (improve23 && improve32) {
      if (best23.smaller_is_better < best32.smaller_is_better)
        improve(best23);
      else
        improve(best32);
      flipped_this_cycle = true;
    } else {
      std::cout << "Delaunay converged after " << iter << " iterations." << std::endl;
      break;
    }

    if (flipped_this_cycle) {
      total_flips++;
      this->setGeometricPropertiesForce();
      this->checkConnectivity();
    } else
      break;

    if (iter == max_iter - 1)
      std::cerr << "Warning: improveTetrahedraDelaunay reached max_iter." << std::endl;
  }

  return total_flips;
}
/* --------------------------------------------------------------------------*/

Network::Flip32Candidate Network::improveTetrahedraDelaunay32(double eps, std::function<double(double, double)> quality_metric) {
  auto insphere_ = [](const networkPoint* a, const networkPoint* b, const networkPoint* c, const networkPoint* d, const networkPoint* e) { return insphere(a->X, b->X, c->X, d->X, e->X); };
  auto TetrahedronVolume_ = [](const networkPoint* a, const networkPoint* b, const networkPoint* c, const networkPoint* d) { return TetrahedronVolume(a->X, b->X, c->X, d->X); };

  Flip32Candidate best32;

  for (auto* l : this->getLines()) {
    auto tetras = l->Tetras;
    if (tetras.size() != 3)
      continue;
    auto const [p, q] = l->getPoints();
    if (!p || !q)
      continue;

    std::unordered_set<networkPoint*> others;
    for (auto f : l->getFaces())
      for (auto a : f->getPoints())
        if (a != p && a != q)
          others.insert(a);

    if (others.size() != 3)
      continue;

    auto others_ = ToVector(others);
    auto a = others_[0];
    auto b = others_[1];
    auto c = others_[2];

    double s_before, s_after, smaller_is_better;

    s_before = std::max({insphere_(p, q, a, b, c), insphere_(p, q, b, c, a), insphere_(p, q, c, a, b)});
    s_after = std::max(insphere_(a, b, c, p, q), insphere_(a, b, c, q, p));

    // 改善量
    if (quality_metric)
      smaller_is_better = quality_metric(s_before, s_after);
    else
      smaller_is_better = s_after - s_before;

    if (smaller_is_better < best32.smaller_is_better)
      best32 = {l, {tetras[0], tetras[1], tetras[2]}, p, q, a, b, c, smaller_is_better, s_before, s_after};
  }

  return best32;
}

void Network::improve(Flip32Candidate& best32) {
  if (!best32.edge)
    return;
  std::cout << Blue << "delaunay flip on edge " << best32.edge << " among tetras " << best32.tetras[0] << ", " << best32.tetras[1] << ", " << best32.tetras[2] << ", improvement = " << best32.smaller_is_better << " s_before = " << best32.s_before << " s_after = " << best32.s_after << colorReset << '\n';
  auto faces = best32.edge->getFaces();
  for (auto f : faces)
    delete f;
  delete best32.edge;
  // 3-2 flip
  ::genTetra(this, T_4P{best32.p, best32.a, best32.b, best32.c});
  std::cout << "genTetra 1;\n";
  ::genTetra(this, T_4P{best32.q, best32.a, best32.b, best32.c});
  std::cout << "genTetra 2;\n";
  // this->setGeometricPropertiesForce();
  // this->checkConnectivity();
}

/* --------------------------------------------------------------------------
 */

Network::Flip23Candidate Network::improveTetrahedraDelaunay23(double eps, std::function<double(const double, const double)> quality_metric) {
  auto insphere_ = [](const networkPoint* a, const networkPoint* b, const networkPoint* c, const networkPoint* d, const networkPoint* e) { return insphere(a->X, b->X, c->X, d->X, e->X); };
  auto TetrahedronVolume_ = [](const networkPoint* a, const networkPoint* b, const networkPoint* c, const networkPoint* d) { return TetrahedronVolume(a->X, b->X, c->X, d->X); };

  Flip23Candidate best23;
  // 適当な面を選び，その両側の四面体を調べる
  for (auto* f : this->getFaces()) {
    if (!f)
      continue;
    auto const [t0, t1] = f->getTetras();
    if (!t0 || !t1)
      continue;
    auto const [a, b, c] = f->getPoints();
    if (!a || !b || !c)
      continue;

    auto findOpp = [&](networkTetra* t) -> networkPoint* {
      for (auto* pp : t->Points)
        if (pp && pp != a && pp != b && pp != c)
          return pp;
      return nullptr;
    };
    auto* p = findOpp(t0);
    auto* q = findOpp(t1);
    if (!p || !q || p == q)
      continue;

    if (isFace(p, q, a) || isFace(p, q, b) || isFace(p, q, c))
      continue;

    //   double oa = orient3d(a->X, b->X, c->X, p->X);
    //   double ob = orient3d(a->X, b->X, c->X, q->X);
    //   if (std::abs(oa) < eps || std::abs(ob) < eps || oa * ob >= 0.)
    //     continue;

    // s > 0 （デローニー条件違反）
    // s <= 0 の場合は、p, qがa, b, cの外側にあることを意味する

    // 体積変化が大きい場合は実行しない
    double vol_before = TetrahedronVolume_(a, b, c, p) + TetrahedronVolume_(a, b, c, q);
    double vol_after = TetrahedronVolume_(p, q, a, b) + TetrahedronVolume_(p, q, b, c) + TetrahedronVolume_(p, q, c, a);
    if (vol_before < 0.9 * vol_after || vol_before > 1.1 * vol_after)
      continue;

    double s_before, s_after, smaller_is_better;
    s_before = std::max(insphere_(a, b, c, p, q), insphere_(a, b, c, q, p));
    s_after = std::max({insphere_(p, q, a, b, c), insphere_(p, q, b, c, a), insphere_(p, q, c, a, b)});

    if (quality_metric) {
      smaller_is_better = quality_metric(s_before, s_after);
    } else {
      smaller_is_better = s_after - s_before;
    }

    if (smaller_is_better < best23.smaller_is_better)
      best23 = {f, t0, t1, a, b, c, p, q, smaller_is_better, s_before, s_after};
  }
  return best23;
}

void Network::improve(Flip23Candidate& best23) {

  if (!best23.f)
    return;
  std::cout << Red << "delaunay flip on face " << best23.f << " between tetras " << best23.t0 << " and " << best23.t1 << ", improvement = " << best23.smaller_is_better << " s_before = " << best23.s_before << " s_after = " << best23.s_after << colorReset << '\n';
  delete best23.f;
  // 2-3 flip
  std::cout << isFace(best23.p, best23.q, best23.a) << std::endl;
  std::cout << isFace(best23.p, best23.q, best23.b) << std::endl;
  std::cout << isFace(best23.p, best23.q, best23.c) << std::endl;
  std::cout << isFace(best23.p, best23.a, best23.b) << std::endl;
  ::genTetra(this, T_4P{best23.p, best23.q, best23.a, best23.b});
  std::cout << "genTetra 1;\n";
  ::genTetra(this, T_4P{best23.p, best23.q, best23.b, best23.c});
  std::cout << "genTetra 2;\n";
  ::genTetra(this, T_4P{best23.p, best23.q, best23.c, best23.a});
  std::cout << "genTetra 3;\n";
  // this->setGeometricPropertiesForce();
  // this->checkConnectivity();
}

/* -------------------------------------------------------------------------- */

T6d Network::velocityPredefined() {
  double t = _current_time_;
  double s = M_PI / 2.;
  // double a = move_amplitude;
  double k = M_PI / 1.;
  /* ------------------------------------------------------ */
  // T6d move_dir = {std::cos(k * t), std::sin(k * t), 0., 0., 0., 0.};
  // T6d ddt_move_dir = {-k * std::sin(k * t), k * std::std::cos(k * t), 0.,
  // 0., 0., 0.};
  // // /* |U|*n_p . n_surface = phin <-- given
  // auto tmp = (-move_amplitude * std::exp(-t) * (sin(k * t - s) -
  // std::sin(-s)) + move_amplitude * std::exp(-t) * (std::cos(k * t - s) *
  // k)) * move_dir; tmp += move_amplitude * std::exp(-t) * (sin(k * t - s) -
  // std::sin(-s)) * ddt_move_dir; return tmp;
  /* ------------------------------------------------------ */
  Tddd tmp = Normalize(Tddd{1., 1., 0.});
  T6d move_dir = {std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp), 0., 0., 0.};
  return (-move_amplitude * std::exp(-t) * (std::sin(k * t - s) - std::sin(-s)) + move_amplitude * std::exp(-t) * (std::cos(k * t - s) * k)) * move_dir;
};

Tddd Network::translationPredefined() {
  double t = _current_time_;
  double s = M_PI / 2.;
  double k = M_PI / 1.;
  /* ------------------------------------------------------ */
  // Tddd move_dir = {std::cos(k * t), std::sin(k * t), 0.};
  // return move_amplitude * std::exp(-t) * (sin(k * t - s) - std::sin(-s)) *
  // move_dir;
  /* ------------------------------------------------------ */
  Tddd move_dir = Normalize(Tddd{1., 1., 0.});
  return move_amplitude * std::exp(-t) * (std::sin(k * t - s) - std::sin(-s)) * move_dir;
};

/* -------------------------------------------------------------------------- */
/*                                 クラスに属さない関数                                 */
/* -------------------------------------------------------------------------- */

bool isEdge(const networkFace* const f) {
  return std::ranges::any_of(f->getPoints(), [](const networkPoint* p) { return isEdge(p); });
};

/* -------------------------------------------------------------------------- */

T3Tddd ToX(const networkFace* const f) { return ToX(f->getPoints()); };
std::vector<T3Tddd> ToX(const std::vector<networkFace*>& fs) {
  std::vector<T3Tddd> ret(fs.size());
  int i = 0;
  for (const auto& f : fs)
    ret[i++] = ToX(f->getPoints());
  return ret;
};
std::vector<T3Tddd> ToX(const std::unordered_set<networkFace*>& fs) {
  std::vector<T3Tddd> ret(fs.size());
  int i = 0;
  for (const auto& f : fs)
    ret[i++] = ToX(f->getPoints());
  return ret;
};
/* -------------------------------------------------------------------------- */
Tddd TriangleNormal(const networkFace* const f) { return TriangleNormal(ToX(f)); };
//@ ------------------------ 抽出用関数など ----------------------- */
std::vector<Tddd> extX(const std::unordered_set<networkFace*>& fs) {
  std::vector<Tddd> ret;
  for (const auto& f : fs)
    ret.emplace_back(f->X);
  return ret;
};
std::vector<T3Tddd> extVertices(const std::unordered_set<networkFace*>& fs) {
  std::vector<T3Tddd> ret(fs.size());
  int i = 0;
  for (const auto& f : fs)
    ret[i++] = f->getXVertices();
  return ret;
};
std::vector<T2Tddd> extX(const std::unordered_set<networkLine*>& ls) {
  std::vector<T2Tddd> ret(ls.size());
  int i = 0;
  for (const auto& l : ls)
    ret[i++] = l->getLocationsTuple();
  return ret;
};
std::vector<Tddd> extNormals(const V_netFp& fs) {
  std::vector<Tddd> ret(fs.size());
  int i = 0;
  for (const auto& f : fs)
    ret[i++] = f->normal;
  return ret;
};
//@ ------------------------------------------------------ */

std::vector<netL*> link(const V_netPp& obj, Network* net) {
  try {
    std::vector<netL*> ret;
    int s = obj.size();
    for (int i = 0; i < s; i++) {
      auto l = link(obj[i], obj[(i + 1) % s], net);
      ret.emplace_back(l);
    }
    return ret;
  } catch (const error_message& e) {
    std::cerr << Red << "[link] diagnostic: obj=" << obj << colorReset << '\n';
    std::cerr << Red << "[link] " << (DuplicateFreeQ(obj) ? "no duplication...???" : "duplication found!") << colorReset << '\n';
    throw;
  }
};
netL* unlink(netP* obj, netP* obj_) {
  if (obj_ == obj)
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "a point is trying to unlink itself!");

  if (!obj || !obj_ /*if NULL*/)
    return nullptr;

  auto line = obj->getLineBetween(obj_);
  auto line_ = obj_->getLineBetween(obj);

  if (line != nullptr) {
    line->erase(obj);
    line->erase(obj_);
  }

  if (obj != nullptr)
    obj->erase(line);
  if (obj_ != nullptr)
    obj_->erase(line);
  return line;

  // if ((line == line_) && line) {
  //    line->erase(obj);
  //    line->erase(obj_);
  //    obj->erase(line);
  //    obj_->erase(line);
  //    return line;
  // } else {
  //    std::cout << Red << obj << " and " << obj_ << " are not linked" <<
  //    colorReset << std::endl; return nullptr;
  // }
};

//@ ------------------------------------------------------ */
//@                         extract                        */
//@ ------------------------------------------------------ */

/*
 * under scorer _ means the function returns FLATTEND list
 */

std::unordered_set<networkPoint*> extPointsCORNER_(const std::vector<networkPoint*>& ps) {
  std::unordered_set<networkPoint*> ret;
  for (const auto& p : ps)
    if (p->CORNER)
      ret.emplace(p);
  return ret;
};
std::unordered_set<networkPoint*> extPointsCORNER_(const std::unordered_set<networkPoint*>& ps) {
  std::unordered_set<networkPoint*> ret;
  for (const auto& p : ps)
    if (p->CORNER)
      ret.emplace(p);
  return ret;
};
std::unordered_set<networkPoint*> extPointsCornerOrNeumann_(const std::vector<networkPoint*>& ps) {
  std::unordered_set<networkPoint*> ret;
  for (const auto& p : ps)
    if (p->CORNER || p->Neumann)
      ret.emplace(p);
  return ret;
};
std::unordered_set<networkPoint*> extPointsCornerOrNeumann_(const std::unordered_set<networkPoint*>& ps) {
  std::unordered_set<networkPoint*> ret;
  for (const auto& p : ps)
    if (p->CORNER || p->Neumann)
      ret.emplace(p);
  return ret;
};

// b! ------------------------------------------------------ */
// b! ---------------------- extLines ---------------------- */
// b! ------------------------------------------------------ */
/* ------------------ for unordered_set ----------------- */
std::unordered_set<networkLine*> extLinesCORNER_(const std::unordered_set<networkFace*>& fs) {
  std::unordered_set<networkLine*> ret;
  for (const auto& f : fs)
    std::ranges::for_each(f->getLines(), [&](const auto& l) {
      if (l->CORNER) {
        ret.emplace(l);
      };
    });
  return ret;
};
std::unordered_set<networkLine*> extLines_(const std::unordered_set<networkFace*>& fs) {
  std::unordered_set<networkLine*> ret;
  for (const auto& f : fs)
    std::ranges::for_each(f->getLines(), [&](const auto& l) { ret.emplace(l); });
  return ret;
};
std::unordered_set<networkLine*> extLinesCORNER_(const std::unordered_set<networkPoint*>& ps) {
  std::unordered_set<networkLine*> ret;
  for (const auto& p : ps)
    for (const auto& l : p->getLines())
      if (l->CORNER)
        ret.emplace(l);
  return ret;
};
std::unordered_set<networkLine*> extLinesCORNER_(const networkPoint* p) {
  std::unordered_set<networkLine*> ret;
  for (const auto& l : p->getLines())
    if (l->CORNER)
      ret.emplace(l);
  return ret;
};
std::unordered_set<networkLine*> extLines_(const std::unordered_set<networkPoint*>& ps) {
  std::unordered_set<networkLine*> ret;
  for (const auto& p : ps)
    for (const auto& l : p->getLines())
      ret.emplace(l);
  return ret;
};
/* --------------------- for vector --------------------- */
std::unordered_set<networkLine*> extLinesCORNER_(const std::vector<networkFace*>& fs) {
  std::unordered_set<networkLine*> ret;
  for (const auto& f : fs)
    std::ranges::for_each(f->getLines(), [&](const auto& l) {
      if (l->CORNER) {
        ret.emplace(l);
      };
    });
  return ret;
};
std::unordered_set<networkLine*> extLines_(const std::vector<networkFace*>& fs) {
  std::unordered_set<networkLine*> ret;
  for (const auto& f : fs)
    std::ranges::for_each(f->getLines(), [&](const auto& l) { ret.emplace(l); });
  return ret;
};
std::unordered_set<networkLine*> extLinesCORNER_(const std::vector<networkPoint*>& ps) {
  std::unordered_set<networkLine*> ret;
  for (const auto& p : ps)
    for (const auto& l : p->getLines())
      if (l->CORNER)
        ret.emplace(l);
  return ret;
};
std::unordered_set<networkLine*> extLines_(const std::vector<networkPoint*>& ps) {
  std::unordered_set<networkLine*> ret;
  for (const auto& p : ps)
    for (const auto& l : p->getLines())
      ret.emplace(l);
  return ret;
};
//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
//@ ------------------------------------------------------ */
V_d extLength(const std::unordered_set<networkLine*>& ls) {
  V_d ret(ls.size());
  int i = 0;
  for (const auto& l : ls)
    ret[i++] = l->length();
  return ret;
};
//
V_d extLength(const V_netLp& ls) {
  V_d ret(ls.size());
  int i = 0;
  for (const auto& l : ls)
    ret[i++] = l->length();
  return ret;
};
Tddd extLength(const std::array<networkLine*, 3>& ls) { return {std::get<0>(ls)->length(), std::get<1>(ls)->length(), std::get<2>(ls)->length()}; };

V_d extAreas(const V_netFp& fs) {
  V_d ret(fs.size());
  int i = 0;
  for (const auto& f : fs)
    ret[i++] = f->area;
  return ret;
};

V_d extractAreas(const V_netFp& fs) {
  V_d ret(fs.size());
  int i = 0;
  for (const auto& f : fs)
    ret[i++] = f->area;
  return ret;
};
// 2021/09/06追加
std::vector<Tddd> extXtuple(const V_netPp& points) {
  std::vector<Tddd> ret(points.size());
  int i = 0;
  for (const auto& p : points)
    ret[i++] = p->X;
  return ret;
};
std::vector<Tddd> extXtuple(const V_netFp& points) {
  std::vector<Tddd> ret(points.size());
  int i = 0;
  for (const auto& p : points)
    ret[i++] = p->X;
  return ret;
};

T3Tddd extractXtuple(networkFace const* f) {
  auto [p0, p1, p2] = f->getPoints();
  // return std::make_tuple(p0->X, p1->X, p2->X);
  return {p0->X, p1->X, p2->X};
};
//
std::vector<std::vector<Tddd>> extractXtuple(const std::vector<networkLine*>& lines) {
  std::vector<std::vector<Tddd>> ret;
  for (const auto& l : lines) {
    ret.push_back({});
    auto [p, q] = l->getPoints();
    ret.rbegin()->emplace_back(p->X);
    ret.rbegin()->emplace_back(q->X);
  }
  return ret;
};

/* -------------------------------------------------------------------------- */
/*                                    四面体関連                                   */
/* -------------------------------------------------------------------------- */

netF* genFace(Network* const net, netL* const l0, netL* const l1, netL* const l2) {
  try {
    if (l0 && l1 && l2 && l0 != l1 && l0 != l2 && l1 != l2) {
      auto intx_faces = Intersection(l0->getFaces(), l1->getFaces(), l2->getFaces());
      auto s = intx_faces.size();
      if (s == 0)
        return new networkFace(net, Point(l2, l0), l0, Point(l0, l1), l1, Point(l1, l2), l2);
      else if (s == 1)
        return intx_faces[0];
      else {
        std::stringstream ss;
        ss << "too many inteserctions of faces : " << intx_faces;
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      }
    } else {
      std::stringstream ss;
      ss << "{l0,l1,l2} = {" << l0 << "," << l1 << "," << l2 << "}";
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
    }
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  };
};

std::tuple<bool, networkTetra*> genTetra(Network* const net, netP* const p0, netP* const p1, netP* p2, netP* p3) {
  try {
    if (!p0 || !p1 || !p2 || !p3)
      return {false, nullptr};

    if (!DuplicateFreeQ(T_4P{p0, p1, p2, p3}))
      return {false, nullptr};

    /* p1, p2, p3 are always directed in the outward normal direction
    //
    これは右手系になっていることを保証するための処理．右手系は計算幾何学の標準のようだ．
    // this orientation Dot(Cross(p1-p0,p2-p0),(p3-p0)) > 0
    //      p0
    //     /|\
    //    / | \
    //   /  |  \
    //  p1--|---p2
    //   \  |  /
    //    \ | /
    //     \|/
    //      p3
    */

    if (CrossDot(p1->X - p0->X, p2->X - p0->X, p3->X - p0->X) < 0)
      std::swap(p2, p3); // 右手系になる

    /*
    //       p0
    //     / | \
    //   l0  |  l1
    //   /   l2   \
    //  p1---|-l3-p2
    //   \   |    /
    //   l5  |  l4
    //     \ | /
    //      p3
    */

    auto [l3, l4, l5] = link(p1, p2, p3, net);
    auto f0 = genFace(net, l3, l4, l5);
    if (f0->Tetras[0] && f0->Tetras[1]) {
      std::cout << "(face:" << f0 << ")->Tetras:" << f0->Tetras << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "face already has 2 tetras");
    }
    auto [l1, l3_, l0] = link(p0, p2, p1, net);
    auto f1 = genFace(net, l1, l3_, l0);
    if (f1->Tetras[0] && f1->Tetras[1]) {
      std::cout << "(face:" << f1 << ")->Tetras:" << f1->Tetras << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "face already has 2 tetras");
    }
    auto [l2, l4_, l1_] = link(p0, p3, p2, net);
    auto f2 = genFace(net, l2, l4_, l1_);
    if (f2->Tetras[0] && f2->Tetras[1]) {
      std::cout << "(face:" << f2 << ")->Tetras:" << f2->Tetras << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "face already has 2 tetras");
    }
    auto [l0_, l5_, l2_] = link(p0, p1, p3, net);
    auto f3 = genFace(net, l0_, l5_, l2_);
    if (f3->Tetras[0] && f3->Tetras[1]) {
      std::cout << "(face:" << f3 << ")->Tetras:" << f3->Tetras << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "face already has 2 tetras");
    }
    if (l0 == l0_ || l1 == l1_ || l2 == l2_ || l3 == l3_ || l4 == l4_ || l5 == l5_) {
      auto [t0, t1] = f0->Tetras;
      if (t0 && Intersection(t0->Points, T_4P{p0, p1, p2, p3}).size() == 4)
        return {false, t0};
      else if (t1 && Intersection(t1->Points, T_4P{p0, p1, p2, p3}).size() == 4)
        return {false, t1};
      else {
        auto tet = new networkTetra(net, T_4P{p0, p1, p2, p3}, T_6L{l0, l1, l2, l3, l4, l5}, T_4F{f0, f1, f2, f3});
        for (const auto& p : {p0, p1, p2, p3})
          if (!MemberQ(p->Tetras, tet))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "contradictions");
        for (const auto& l : {l0, l1, l2, l3, l4, l5})
          if (!MemberQ(l->Tetras, tet))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "contradictions");
        for (const auto& f : {f0, f1, f2, f3})
          if (!MemberQ(f->Tetras, tet))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "contradictions");
        return {true, tet};
      }
    } else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "contradictions");
  } catch (const error_message&) {
    throw;
  } catch (const std::exception& e) {
    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, e.what());
  }
};

std::tuple<bool, networkTetra*> genTetra(Network* const net, const T_4P& abcd) {
  auto [a, b, c, d] = abcd;
  return genTetra(net, a, b, c, d);
};

Tddd Nearest(const Tddd& X, const networkFace* f) {
  return Nearest(X, ToX(f));
  // auto [a, b, c] = ToX(f);
  // Tddd m = (a + b + c) / 3.;
  // Tddd ab = (a + b) / 2., bc = (b + c) / 2., ca = (c + a) / 2.;
  // auto ret = Nearest(X, ToX(f));
  // auto X1 = Nearest(X, T3Tddd{a, ab, ca});
  // auto X2 = Nearest(X, T3Tddd{b, bc, ab});
  // auto X3 = Nearest(X, T3Tddd{c, ca, bc});
  // auto X4 = Nearest(X, T3Tddd{ab, bc, ca});
  // if (Norm(ret - X) > Norm(X1 - X))
  //    ret = X1;
  // if (Norm(ret - X) > Norm(X2 - X))
  //    ret = X2;
  // if (Norm(ret - X) > Norm(X3 - X))
  //    ret = X3;
  // if (Norm(ret - X) > Norm(X4 - X))
  //    ret = X4;
  // return ret;
};
double Distance(const Tddd& X, const networkFace* f) { return Norm(Nearest(X, ToX(f)) - X); };
// Tddd Nearest(const networkPoint *p, const networkFace *f) { return
// Nearest(ToX(p), ToX(f)); };
std::tuple<Tddd, networkFace*> Nearest_(const Tddd& X, const std::unordered_set<networkFace*>& faces) {
  double nearest_d = 1E+20, d;
  Tddd nearest_x, x;
  networkFace* nearest_f = nullptr;
  std::ranges::for_each(faces, [&](const auto& f) {
    x = Nearest(X, f);
    if (nearest_d >= (d = Norm(X - x))) {
      nearest_d = d;
      nearest_x = x;
      nearest_f = f;
    }
  });
  return {nearest_x, nearest_f};
};
std::tuple<Tddd, networkFace*> Nearest_(const Tddd& X, const std::vector<networkFace*>& faces) {
  double distance = 1E+20, tmp;
  Tddd ret, near;
  networkFace* F = nullptr;
  for (const auto& f : faces)
    if (distance > (tmp = Norm(X - (near = Nearest(X, f))))) {
      distance = tmp;
      ret = near;
      F = f;
    }
  return {ret, F};
};
std::tuple<Tddd, networkFace*> Nearest_(const Tddd& X, const double& r, const std::unordered_set<networkFace*>& faces) {
  double distance = 1E+20, tmp;
  Tddd ret, near;
  networkFace* F = nullptr;
  for (const auto& f : faces)
    if (distance > (tmp = Norm(X - (near = Nearest(X, f))))) {
      distance = tmp;
      ret = near;
      F = f;
    }
  return {ret, F};
};
Tddd Nearest(const Tddd& X, const std::vector<networkFace*>& faces) { return std::get<0>(Nearest_(X, faces)); };
Tddd Nearest(const Tddd& X, const std::unordered_set<networkFace*>& faces) { return std::get<0>(Nearest_(X, faces)); };
Tddd Nearest(const networkPoint* p, const std::unordered_set<networkFace*>& faces) { return Nearest(ToX(p), faces); };

double DirectedHausdorffDistance(const V_netFp& query_faces, const V_netFp& target_faces, const int triangle_subdivision) {
  const auto query = uniqueValidFaces(query_faces);
  const auto target = uniqueValidFaces(target_faces);

  if (query.empty())
    return 0.;
  if (target.empty())
    return std::numeric_limits<double>::infinity();

  switch (triangle_subdivision) {
  case 1:
    return directedHausdorffDistanceImpl<1>(query, target);
  case 2:
    return directedHausdorffDistanceImpl<2>(query, target);
  case 3:
    return directedHausdorffDistanceImpl<3>(query, target);
  case 4:
    return directedHausdorffDistanceImpl<4>(query, target);
  case 5:
    return directedHausdorffDistanceImpl<5>(query, target);
  case 6:
    return directedHausdorffDistanceImpl<6>(query, target);
  case 8:
    return directedHausdorffDistanceImpl<8>(query, target);
  case 10:
    return directedHausdorffDistanceImpl<10>(query, target);
  case 12:
    return directedHausdorffDistanceImpl<12>(query, target);
  case 15:
    return directedHausdorffDistanceImpl<15>(query, target);
  default:
    return directedHausdorffDistanceImpl<6>(query, target);
  }
}

double DirectedHausdorffDistance(const V_netFp& query_faces, const V_Netp& target_networks, const int triangle_subdivision) {
  return DirectedHausdorffDistance(query_faces, boundaryFacesFromNetworks(target_networks), triangle_subdivision);
}

/* -------------------------------------------------------------------------- */
/*                                    接続関係                                    */
/* -------------------------------------------------------------------------- */

networkFace* isFace(const networkPoint* const a, const networkPoint* const b, const networkPoint* const c) {
  for (const auto& p : {a, b, c})
    for (const auto& f : p->Faces)
      if (f->MemberQ(a) && f->MemberQ(b) && f->MemberQ(c))
        return f;
  return nullptr;
};

networkLine* ConnectedQ(const networkFace* const a, const networkFace* const b) {
  auto [l0, l1, l2] = a->getLines();
  if ((*l0)(a) == b)
    return l0;
  else if ((*l1)(a) == b)
    return l1;
  else if ((*l2)(a) == b)
    return l2;
  return nullptr;
};

networkLine* ConnectedQ(const networkLine* const l, const networkPoint* const a) {
  // 両方が参照し合う状態かどうか
  /*
   __line__
   |      |
   |     *|-><--* point
   --------
  */
  for (const auto& L : a->getLines())
    if (L == l) {
      auto [p0, p1] = L->getPoints();
      if (p0 == a || p1 == a)
        return L;
    }
  return nullptr;
};

networkLine* ConnectedQ(const networkPoint* const a, const networkLine* const l) { return ConnectedQ(l, a); };

networkLine* ConnectedQ(const networkPoint* const a, const networkPoint* const b) {
  // 両方が参照し合う状態かどうか
  /*
   point  __line__
   *--><--|*     |
          |     *|-><--* point
          --------
  */
  networkLine* L;
  for (const auto& l : a->getLines())
    if (L = ConnectedQ(l, a))
      if (ConnectedQ(L, b))
        return L;
  return nullptr;
};

T_3L ConnectedQ(const networkPoint* const a, const networkPoint* const b, const networkPoint* const c) {
  return T_3L{ConnectedQ(a, b), ConnectedQ(b, c), ConnectedQ(c, a)};

  // for (const auto &la : a->getLines())
  //    if ((*la)(a) == b)
  //       for (const auto &lb : b->getLines())
  //          if ((*lb)(b) == c)
  //             for (const auto &lc : c->getLines())
  //                if ((*lc)(c) == a)
  //                   return T_3L{la, lb, lc};
  // return T_3L{nullptr, nullptr, nullptr};
};

T_3L ConnectedQ(const T_3P& abc) { return ConnectedQ(std::get<0>(abc), std::get<1>(abc), std::get<2>(abc)); };

/* -------------------------------------------------------------------------- */
/*                                     置換                                     */
/* -------------------------------------------------------------------------- */

void dual_replace(netP* p, netL* line_before, netL* line_after) {
  p->replace(line_before, line_after);
  line_before->replace(p, nullptr); // lineのpointの数は固定:2
  line_after->replace(nullptr, p);  // lineのpointの数は固定:2
}

void dual_replace(netP* p, netF* face_before, netF* face_after) {
  p->replace(face_before, face_after);
  face_before->replace(p, nullptr); // faceのpointの数は固定:3
  face_after->replace(nullptr, p);  // faceのpointの数は固定:3
}

/* -------------------------------------------------------------------------- */

void dual_replace(netL* l, netP* point_before, netP* point_after) {
  l->replace(point_before, point_after);
  point_before->erase(l); //! pointのlineの数は固定ではない
  point_after->add(l);    //! pointのlineの数は固定ではない
}

// void dual_replace(netL *l, netF *face_before, netF *face_after) {
//   l->replace(face_before, face_after);
//   face_before->replace(l, nullptr); // faceのlineの数は固定:3
//   face_after->replace(nullptr, l);  // faceのlineの数は固定:3
// }

// これは上を更に発展させたもので，face_afterの持っていたl_beforeの代わりにlを入れる．
void dual_replace(netL* l, netF* face_before, netF* face_after, netL* l_face_before, netL* l_face_after) {
  l->replace(face_before, face_after);
  face_before->replace(l, l_face_before); // faceのlineの数は固定:3
  face_after->replace(l_face_after, l);   // faceのlineの数は固定:3
  if (l_face_before != nullptr)
    l_face_before->add(face_before); //! lineのfaceの数は固定ではない
  if (l_face_after != nullptr)
    l_face_after->erase(face_after); //! lineのfaceの数は固定ではない
}

/* -------------------------------------------------------------------------- */

void dual_replace(netF* f, netL* line_before, netL* line_after) {
  f->replace(line_before, line_after);
  line_before->erase(f); //! lineのfaceの数は固定ではない
  line_after->add(f);    //! lineのfaceの数は固定ではない
}

void dual_replace(netF* f, netL* line_before, netL* line_after, netF* f_line_before) {
  f->replace(line_before, line_after);
  line_before->replace(f, f_line_before); //! lineのfaceの数は固定ではない
  line_after->add(f);                     //! lineのfaceの数は固定ではない
}

void dual_replace(netF* f, netL* line_before, netL* line_after, netF* f_line_before, netF* f_line_after) {
  f->replace(line_before, line_after);
  line_before->replace(f, f_line_before); //! lineのfaceの数は固定ではない
  line_after->replace(f_line_after, f);   //! lineのfaceの数は固定ではない
}

void dual_replace(netF* f, netP* point_before, netP* point_after) {
  f->replace(point_before, point_after);
  point_before->erase(f); //! pointのfaceの数は固定ではない
  point_after->add(f);    //! pointのfaceの数は固定ではない
}

void dual_replace(netF* f, netP* point_before, netP* point_after, netF* f_point_before) {
  f->replace(point_before, point_after);
  point_before->replace(f, f_point_before); //! pointのfaceの数は固定ではない
  point_after->add(f);                      //! pointのfaceの数は固定ではない
}

// ---------------------------------------------------------------------------
// buildLocalFrame implementations
// ---------------------------------------------------------------------------

void networkLine::buildLocalFrame() {
  this->local_frame = {};
  auto [p0, p1] = this->getPoints();
  if (!p0 || !p1)
    return;
  Tddd edge_dir = p1->X - p0->X;
  Tddd n_avg = p0->getNormalAreaAveraged() + p1->getNormalAreaAveraged();
  double n_len = Norm(n_avg);
  if (!(n_len > 1e-10))
    return;
  Tddd n_axis = n_avg / n_len;
  Tddd u_proj = edge_dir - Dot(edge_dir, n_axis) * n_axis;
  double u_len = Norm(u_proj);
  if (!(u_len > 1e-10))
    return;
  Tddd u_axis = u_proj / u_len;
  Tddd v_axis = Cross(n_axis, u_axis);
  this->local_frame = {u_axis, v_axis, n_axis, true};
}

// ---------------------------------------------------------------------------
// Principal curvature computation via quadric fitting
// ---------------------------------------------------------------------------

void Network::computePrincipalCurvatures(const bool force_all) {

  // 共通処理: 近傍座標の抽出 → フィット → geom_curvature 更新
  auto fitAndStore = [](BEM_DOF_Base* dof,
                        const std::unordered_set<networkPoint*>& point_set,
                        const std::unordered_set<networkLine*>& line_set) {
    std::vector<Tddd> neighbor_pos;
    neighbor_pos.reserve(point_set.size() + line_set.size());
    for (auto* q : point_set)
      neighbor_pos.push_back(q->X);
    for (auto* ll : line_set)
      neighbor_pos.push_back(ll->X_mid);

    if (neighbor_pos.size() < 5)
      return;

    dof->buildLocalFrame();
    if (!dof->local_frame.valid)
      return;
    auto curv = surface_geometry::computePrincipalCurvatures(
        dof->getPosition(), dof->local_frame.u_axis, dof->local_frame.v_axis, dof->local_frame.n_axis, neighbor_pos);
    if (!curv.valid)
      return;

    dof->geom_curvature = curv;
  };

  // --- 頂点の曲率 ---
  for (auto* p : this->getBoundaryPoints()) {
    if (!force_all && p->geom_curvature.valid)
      continue;
    p->geom_curvature = {};

    std::unordered_set<networkPoint*> point_set;
    std::unordered_set<networkLine*> line_set;
    for (auto* q : p->getNeighbors())
      point_set.insert(q);
    for (auto* ll : p->getLines())
      if (!ll->getBoundaryFaces().empty())
        line_set.insert(ll);

    fitAndStore(p, point_set, line_set);
  }

  // --- 辺中点の曲率 ---
  for (auto* l : this->getBoundaryLines()) {
    if (!force_all && l->geom_curvature.valid)
      continue;
    l->geom_curvature = {};

    auto [p0, p1] = l->getPoints();
    if (!p0 || !p1)
      continue;

    std::unordered_set<networkPoint*> point_set;
    std::unordered_set<networkLine*> line_set;
    for (auto* q : p0->getNeighbors())
      point_set.insert(q);
    for (auto* q : p1->getNeighbors())
      point_set.insert(q);
    point_set.insert(p0);
    point_set.insert(p1);
    for (auto* q : point_set)
      for (auto* ll : q->getLines())
        if (ll != l && !ll->getBoundaryFaces().empty())
          line_set.insert(ll);

    fitAndStore(l, point_set, line_set);
  }
}

// ---------------------------------------------------------------------------
// replacePatch: patchB の内容で this の対応領域を置き換える
// ---------------------------------------------------------------------------
// patchB 外周点と this 上の点を座標比較で対応づけ、
// this 側の古い内部要素を削除し、patchB の要素を this に移植する。

bool Network::replacePatch(Network& patchB, double coord_eps) {
  // 1. patchB の外周点を特定（面が < 2 の辺の端点）
  std::unordered_set<networkPoint*> patch_boundary;
  std::unordered_set<networkLine*> patch_boundary_lines;
  for (auto* l : patchB.Lines)
    if (l->getBoundaryFaces().size() < 2) {
      patch_boundary_lines.insert(l);
      auto [p0, p1] = l->getPoints();
      if (p0) patch_boundary.insert(p0);
      if (p1) patch_boundary.insert(p1);
    }

  // 2. 座標マッチング: patchB 外周点 → this の点
  std::unordered_map<networkPoint*, networkPoint*> boundary_match; // copy → orig
  std::unordered_set<networkPoint*> used_orig_boundary;
  for (auto* bp : patch_boundary) {
    networkPoint* found = nullptr;
    for (auto* wp : this->Points) {
      if (Norm(bp->X - wp->X) < coord_eps) {
        if (used_orig_boundary.count(wp))
          return false; // patchB の複数境界点が同一点に潰れるのは reject
        found = wp;
        break;
      }
    }
    if (!found) return false;
    boundary_match[bp] = found;
    used_orig_boundary.insert(found);
  }

  // 3. this 側の古い内部要素を特定
  std::unordered_set<networkPoint*> orig_boundary;
  for (auto& [copy, orig] : boundary_match)
    orig_boundary.insert(orig);

  // patchB の外周辺に対応する this 側の外周辺を特定
  std::unordered_set<networkLine*> orig_boundary_lines;
  std::unordered_set<networkFace*> boundary_adjacent_faces;
  for (auto* copy_line : patch_boundary_lines) {
    auto [cp0, cp1] = copy_line->getPoints();
    if (!cp0 || !cp1)
      return false;
    auto it0 = boundary_match.find(cp0);
    auto it1 = boundary_match.find(cp1);
    if (it0 == boundary_match.end() || it1 == boundary_match.end())
      return false;
    auto* orig_line = Line(it0->second, it1->second);
    if (!orig_line)
      return false;
    orig_boundary_lines.insert(orig_line);
    for (auto* f : orig_line->getBoundaryFaces())
      if (f)
        boundary_adjacent_faces.insert(f);
  }
  if (orig_boundary_lines.empty() || boundary_adjacent_faces.empty())
    return false;

  // boundary line を跨がない面グラフの連結成分を作り，boundary に接する成分のうち
  // 最小のものを old patch とみなす
  std::unordered_set<networkFace*> visited;
  std::vector<std::vector<networkFace*>> boundary_components;
  for (auto* seed_face : boundary_adjacent_faces) {
    if (!seed_face || visited.count(seed_face))
      continue;
    std::vector<networkFace*> comp;
    std::deque<networkFace*> q = {seed_face};
    visited.insert(seed_face);
    while (!q.empty()) {
      auto* f = q.front();
      q.pop_front();
      comp.push_back(f);
      for (auto* l : f->getLines()) {
        if (!l || orig_boundary_lines.count(l))
          continue;
        for (auto* nf : l->getBoundaryFaces()) {
          if (!nf || visited.count(nf))
            continue;
          visited.insert(nf);
          q.push_back(nf);
        }
      }
    }
    if (!comp.empty())
      boundary_components.push_back(std::move(comp));
  }
  if (boundary_components.size() < 2)
    return false;

  auto best_it = std::min_element(boundary_components.begin(), boundary_components.end(),
                                  [](const auto& a, const auto& b) {
                                    return a.size() < b.size();
                                  });
  std::vector<networkFace*> old_faces = *best_it;
  std::unordered_set<networkFace*> old_face_set(old_faces.begin(), old_faces.end());

  // 古い辺: old_faces のみに属し，境界切断辺ではないものを削除対象にする
  std::vector<networkLine*> old_lines;
  for (auto* f : old_faces)
    for (auto* l : f->getLines()) {
      if (!l || orig_boundary_lines.count(l))
        continue;
      bool all_faces_inside = true;
      for (auto* lf : l->getBoundaryFaces())
        if (!old_face_set.count(lf)) {
          all_faces_inside = false;
          break;
        }
      if (all_faces_inside)
        old_lines.push_back(l);
    }
  std::sort(old_lines.begin(), old_lines.end());
  old_lines.erase(std::unique(old_lines.begin(), old_lines.end()), old_lines.end());
  std::unordered_set<networkLine*> old_line_set(old_lines.begin(), old_lines.end());

  // 古い内部点: old_faces のみからなる内部点だけを削除対象にする
  std::unordered_set<networkPoint*> old_internal_points;
  std::unordered_set<networkPoint*> old_face_points;
  for (auto* f : old_faces)
    for (auto* p : f->getPoints())
      if (p && !orig_boundary.count(p))
        old_face_points.insert(p);

  for (auto* p : old_face_points) {
    bool all_faces_inside = true;
    for (auto* pf : p->getBoundaryFaces())
      if (!old_face_set.count(pf)) {
        all_faces_inside = false;
        break;
      }
    if (!all_faces_inside)
      continue;

    bool all_lines_inside = true;
    for (auto* pl : p->getLines())
      if (!old_line_set.count(pl)) {
        all_lines_inside = false;
        break;
      }
    if (all_lines_inside)
      old_internal_points.insert(p);
  }

  // 4. patchB の外周コピー点を this の元の点に差し替える
  std::unordered_set<networkFace*> affected_faces;
  for (auto& [copy, orig] : boundary_match) {
    auto copy_lines = copy->getLines();
    for (auto* l : copy_lines)
      dual_replace(l, copy, orig);

    auto copy_faces = copy->getBoundaryFaces();
    for (auto* f : copy_faces) {
      dual_replace(f, copy, orig);
      affected_faces.insert(f);
    }
  }

  // 5. patchB の外周コピー辺を this の既存外周辺へ差し替える
  std::vector<networkLine*> boundary_line_copies(patch_boundary_lines.begin(), patch_boundary_lines.end());
  for (auto* copy_line : boundary_line_copies) {
    auto [cp0, cp1] = copy_line->getPoints();
    if (!cp0 || !cp1)
      return false;
    auto* orig_line = Line(cp0, cp1);
    if (!orig_line)
      return false;

    auto faces = copy_line->getBoundaryFaces();
    for (auto* f : faces) {
      dual_replace(f, copy_line, orig_line);
      affected_faces.insert(f);
    }

    auto [p0, p1] = copy_line->getPoints();
    if (p0) {
      p0->erase(copy_line);
      copy_line->replace(p0, nullptr);
    }
    if (p1) {
      p1->erase(copy_line);
      copy_line->replace(p1, nullptr);
    }
  }

  // 5b. replace で Points/Lines が変わった面の PLPLPL を再同期
  for (auto* f : affected_faces)
    f->syncPLPLPL();

  // face::replace(point/line) は Lines/Points を更新するが，
  // patch replacement では dodeca 用に参照される PLPLPL も確実に同期しておく。
  for (auto* f : patchB.Faces)
    f->syncPLPLPL();

  // 6. patchB の内部要素を this に移植（network ポインタ変更 + add + patchB から erase）
  //    古い要素を delete する前に移植し、外周点に新しい面/辺が接続された状態にする
  {
    auto pts = std::vector(patchB.Points.begin(), patchB.Points.end());
    for (auto* p : pts) {
      if (patch_boundary.count(p)) continue;
      p->network = this;
      this->add(p);
      patchB.erase(p);
    }
  }
  {
    auto lns = std::vector(patchB.Lines.begin(), patchB.Lines.end());
    for (auto* l : lns) {
      if (patch_boundary_lines.count(l)) continue;
      l->network = this;
      this->add(l);
      patchB.erase(l);
    }
  }
  {
    auto fcs = std::vector(patchB.Faces.begin(), patchB.Faces.end());
    for (auto* f : fcs) {
      f->network = this;
      this->add(f);
      patchB.erase(f);
    }
  }

  // 7. 古い要素を削除（面 → 辺 → 点の順、Delete() が参照解除 + Network::erase を自動で行う）
  for (auto* f : old_faces)
    delete f;
  for (auto* l : old_lines)
    delete l;
  for (auto* p : old_internal_points)
    delete p;

  // 8. 置換済みの patchB 外周コピー辺/点を削除
  for (auto* l : boundary_line_copies)
    delete l;
  for (auto& [copy, orig] : boundary_match)
    delete copy;

  // 9. 幾何プロパティ更新 + 接続チェック
  std::cout << Green << "[replacePatch] success: boundary=" << patch_boundary.size()
            << " old_faces=" << old_faces.size()
            << " old_lines=" << old_lines.size()
            << " old_points=" << old_internal_points.size()
            << colorReset << std::endl;
  this->setGeometricPropertiesForce();
  this->checkConnectivity();

  return true;
}

// ---------------------------------------------------------------------------
// copyLocalPatch: 辺の近傍を別の Network にコピーする
// ---------------------------------------------------------------------------

networkLine* Network::copyLocalPatch(Network& patch, networkLine* l, int ring_depth) const {
  auto [p0, p1] = l->getPoints();
  if (!p0 || !p1)
    return nullptr;

  // 1. seed を構築（辺の隣接面の全頂点 + ring_depth 拡張）
  std::unordered_set<networkPoint*> seed;
  for (auto* f : l->getFaces()) {
    auto [fp0, fp1, fp2] = f->getPoints();
    seed.insert(fp0);
    seed.insert(fp1);
    seed.insert(fp2);
  }
  for (int r = 0; r < ring_depth; ++r) {
    std::unordered_set<networkPoint*> next;
    for (auto* p : seed)
      for (auto* q : p->getNeighbors())
        if (q)
          next.insert(q);
    seed.insert(next.begin(), next.end());
  }

  // 2. seed 内の全点をコピー
  std::unordered_map<networkPoint*, networkPoint*> orig2copy;
  for (auto* p : seed) {
    auto* cp = new networkPoint(&patch, p->X);
    cp->CORNER = p->CORNER;
    cp->Neumann = p->Neumann;
    cp->Dirichlet = p->Dirichlet;
    cp->isMultipleNode = p->isMultipleNode;
    cp->phiphin = p->phiphin;
    cp->phiphin_t = p->phiphin_t;
    orig2copy[p] = cp;
  }

  // 3. 3頂点が全て seed 内にある boundary face をコピー
  // （face コンストラクタ内の link() で line が自動生成される）
  for (auto* f : this->getBoundaryFaces()) {
    auto [fp0, fp1, fp2] = f->getPoints();
    if (orig2copy.count(fp0) && orig2copy.count(fp1) && orig2copy.count(fp2))
      new networkFace(&patch, orig2copy[fp0], orig2copy[fp1], orig2copy[fp2]);
  }

  // 3b. パッチの辺に元の辺のデータをコピー（X_mid, phiphin, phiphin_t, midpoint_index 等）
  // copy → original の逆引き配列
  std::vector<std::pair<networkPoint*, networkPoint*>> copy2orig;
  copy2orig.reserve(orig2copy.size());
  for (auto& [orig, copy] : orig2copy)
    copy2orig.push_back({copy, orig});

  for (auto* copy_l : patch.getLines()) {
    auto [copy_p0, copy_p1] = copy_l->getPoints();
    networkPoint* orig_p0 = nullptr;
    networkPoint* orig_p1 = nullptr;
    for (auto& [copy_p, orig_p] : copy2orig) {
      if (copy_p == copy_p0) orig_p0 = orig_p;
      if (copy_p == copy_p1) orig_p1 = orig_p;
    }
    if (orig_p0 && orig_p1) {
      auto* orig_l = Line(orig_p0, orig_p1);
      if (orig_l) {
        copy_l->CORNER = orig_l->CORNER;
        copy_l->Neumann = orig_l->Neumann;
        copy_l->Dirichlet = orig_l->Dirichlet;
        copy_l->midpoint_index = orig_l->midpoint_index;
        copy_l->phiphin = orig_l->phiphin;
        copy_l->phiphin_t = orig_l->phiphin_t;
        copy_l->X_mid = orig_l->X_mid;
      }
    }
  }

  // 4. 幾何プロパティ更新
  patch.setGeometricPropertiesForce();

  // 5. コピー側の対象辺を返す
  return Line(orig2copy[p0], orig2copy[p1]);
}

// ---------------------------------------------------------------------------
// copyLocalPatch: 点の近傍を別の Network にコピーする
// ---------------------------------------------------------------------------

networkPoint* Network::copyLocalPatch(Network& patch, networkPoint* p, int ring_depth) const {
  if (!p)
    return nullptr;

  // 1. seed を構築（点 + ring_depth 拡張）
  std::unordered_set<networkPoint*> seed = {p};
  for (int r = 0; r < ring_depth; ++r) {
    std::unordered_set<networkPoint*> next;
    for (auto* q : seed)
      for (auto* r : q->getNeighbors())
        if (r)
          next.insert(r);
    seed.insert(next.begin(), next.end());
  }

  // 2. seed 内の全点をコピー
  std::unordered_map<networkPoint*, networkPoint*> orig2copy;
  for (auto* q : seed) {
    auto* cq = new networkPoint(&patch, q->X);
    cq->CORNER = q->CORNER;
    cq->isMultipleNode = q->isMultipleNode;
    cq->Neumann = q->Neumann;
    cq->Dirichlet = q->Dirichlet;
    cq->phiphin = q->phiphin;
    cq->phiphin_t = q->phiphin_t;
    orig2copy[q] = cq;
  }

  // 3. 3頂点が全て seed 内にある boundary face をコピー
  for (auto* f : this->Faces) {
    if (!f->BoundaryQ())
      continue;
    auto [fp0, fp1, fp2] = f->getPoints();
    if (orig2copy.count(fp0) && orig2copy.count(fp1) && orig2copy.count(fp2))
      new networkFace(&patch, orig2copy[fp0], orig2copy[fp1], orig2copy[fp2]);
  }

  // 3b. パッチの辺に元の辺のデータをコピー（X_mid, phiphin, phiphin_t, midpoint_index 等）
  for (auto* cl : patch.getLines()) {
    auto [cp0, cp1] = cl->getPoints();
    networkPoint* orig_p0 = nullptr;
    networkPoint* orig_p1 = nullptr;
    for (auto& [orig, copy] : orig2copy) {
      if (copy == cp0) orig_p0 = orig;
      if (copy == cp1) orig_p1 = orig;
    }
    if (orig_p0 && orig_p1) {
      auto* orig_line = Line(orig_p0, orig_p1);
      if (orig_line) {
        cl->X_mid = orig_line->X_mid;
        cl->phiphin = orig_line->phiphin;
        cl->phiphin_t = orig_line->phiphin_t;
        cl->midpoint_index = orig_line->midpoint_index;
        cl->CORNER = orig_line->CORNER;
        cl->Neumann = orig_line->Neumann;
        cl->Dirichlet = orig_line->Dirichlet;
      }
    }
  }

  // 4. 幾何プロパティ更新
  patch.setGeometricPropertiesForce();

  // 5. コピー側の対象点を返す
  return orig2copy[p];
}
