#pragma once

#include <fstream>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "BEM_step_failure.hpp"
#include "BEM_collision_settings.hpp"

struct networkPoint;
struct networkLine;
struct networkFace;
class Network;

// ============================================================================
// Surface collision detection and resolution
//
// Problem A: Non-adjacent surface collision (plunging breaker)
// Problem B: Adjacent face folding (degenerate mesh deformation)
// ============================================================================

struct CollisionZone {
  std::unordered_set<networkFace*> problem_faces;
  std::unordered_set<networkFace*> boundary_faces;
  std::unordered_set<networkPoint*> interior_points;
  std::unordered_set<networkPoint*> boundary_points;
  std::unordered_set<networkPoint*> all_points;
  std::unordered_set<networkLine*> all_lines;
  std::unordered_map<networkPoint*, Tdd> saved_phiphin;
  std::unordered_map<networkPoint*, Tdd> saved_phiphin_t;
};

struct CollisionDetectionResult {
  std::unordered_set<networkFace*> faces;
  std::vector<std::pair<networkFace*, networkFace*>> pairs;
};

// --- Function declarations ---

std::unordered_set<networkFace*> detectFoldedFaces(
    const Network& water,
    double normal_reversal_cos,
    double mean_len = 0.0);

bool isTopologicallyClose(
    const networkFace* f, const networkFace* g, int max_depth);

double triangleTriangleDistance(const networkFace* f, const networkFace* g);

CollisionDetectionResult detectNonAdjacentCollisions(
    Network& water,
    double proximity_threshold);

std::vector<CollisionZone> buildCollisionZones(
    const std::unordered_set<networkFace*>& problem_faces,
    const std::vector<std::pair<networkFace*, networkFace*>>& collision_pairs,
    int min_zone_faces);

#ifdef tetgenH
bool resolveCollisionZone(
    Network& water,
    CollisionZone& zone,
    int time_step,
    int zone_idx);
#endif

bool detectAndResolveCollisions(
    Network& water,
    int time_step,
    const CollisionSettings& settings,
    std::unordered_set<networkLine*>& protected_lines);
