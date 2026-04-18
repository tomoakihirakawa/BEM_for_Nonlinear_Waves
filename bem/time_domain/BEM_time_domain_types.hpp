#pragma once

enum class NodeRelocationMethod {
  none,
  ALE,
  interpolation
};

enum class NodeRelocationSurface {
  linear,
  pseudo_quadratic,
  true_quadratic
};

enum class InterpolationMidpointMode {
  nearest
};
