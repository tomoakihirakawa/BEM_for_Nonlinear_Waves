#pragma once

// Lightweight remeshing settings — no Network.hpp dependency, safe to include from anywhere.

struct CollisionSettings {
   bool enabled = false;
   double proximity_factor = 1.5;
   double normal_reversal_cos = -0.3;
   int min_zone_faces = 3;
   bool detect_folding = true;
   bool detect_non_adjacent = true;
   bool resolve_with_tetgen = true;
};

struct SubsurfaceAltitudeRejectSettings {
   bool enabled = true;
   double min_face_altitude_rel = 0.05;
   double min_edge_angle_deg = 1.0;
   double max_edge_angle_deg = 120.0;
};
