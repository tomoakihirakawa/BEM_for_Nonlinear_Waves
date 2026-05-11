# Remesh Architecture Notes

This note maps the active remeshing code path. The remesh implementation is
Trial-only: old method names are not converted automatically.

## Top-Level Flow

```text
main_time_domain.cpp
  -> BEMMeshPipeline::TopologyModifier::applyLocalRemesh()
    -> remesh_for_main_loop()
      -> remesh_local_patch_trial_multi_objective()  [Trial]
```

`remesh_method` must be `Trial` or `local_patch_trial_multi_objective`.

## File Roles

- `BEM_remesh_main.hpp/.cpp`
  - Validates `SimulationSettings::RemeshingSettings::remesh_method`.
  - Dispatches only to Trial.

- `BEM_remesh_local_patch.cpp`
  - Main Trial implementation.
  - Scans candidate boundary lines.
  - Copies local patches around trigger lines.
  - Runs scenario operations such as `P`, `Sc`, `Fvs`, and `Fds`.
  - Scores candidate patches with quality, Hausdorff, altitude, contact, and
    repair-aware geometry checks.
  - Applies accepted patches with `replacePatch`.
  - Runs the shared global flip/smooth polish when enabled.

- `BEM_remesh_global_passes.cpp`
  - Owns shared global `pass_flip()` and `pass_smooth()`.
  - Owns `valenceDeviationScore()`.

- `BEM_remesh_scenario_engine.hpp/.cpp`
  - Runs per-edge scenario patches used by Trial.
  - Provides `scenario_engine::runTrials()`.

- `BEM_remesh_common.hpp`
  - Shared predicates, smoothing helpers, quality helpers, settings-facing
    declarations, and Trial remesh declarations.

- `BEM_remesh_logging.hpp`
  - Trial profile, reject summary, geometry repair, and topology log helpers.

## Geometry Repair Dependencies

- `BEM_geometry_repair.hpp`
  - Repairs current water geometry against reference free-surface/body geometry.
  - Repairs BCInterface waterline midpoints.
  - Summarizes repair quality and damaged topology.

- `BEM_geometry_projector.hpp`
  - Chooses projection targets.
  - Sharp and BCInterface regions use constrained projection.
  - Flat regions use nearest-style projection.

- `Network.cpp`, `networkPoint.cpp`, `networkLine.cpp`
  - Contact face detection, contact ranges, and Hausdorff queries.

## Settings And GUI

Core settings live in `BEM_inputfile_reader.hpp`.

Trial keeps:

- `remesh_method`
- `remesh_split_scenarios`
- `remesh_collapse_scenarios`
- `remesh_flip_scenarios`
- `remesh_flip_mode`
- `remesh_smooth_mode`
- `remesh_odt_smooth_iterations`
- `remesh_global_flip_enabled`
- `remesh_global_smoothing_enabled`
- `remesh_global_smooth_iterations`
- split/collapse/flip operation limits
- `remesh_quality_*`
- `remesh_debug_*`

GUI settings are maintained in:

- `gui/app/ui/settings_editor.py`
- `gui/app/ui/remesh_workspace.py`
- `gui/app/ui/tooltips.py`
- `gui/app/integrations/modules/ui_constants.py`

## Trial Responsibilities

`BEM_remesh_local_patch.cpp` currently contains these logical blocks:

1. Boundary/contact refresh helpers.
2. Debug recorder setup.
3. Target length and quality thresholds.
4. Candidate-line predicates and priority ordering.
5. Patch-copy and old-surface projection helpers.
6. Virtual geometry repair for trial patches.
7. Trial execution and scoring.
8. Split/collapse/flip loops.
9. Global flip/smoothing polish.
10. Post-remesh quality checks.
11. Profile and optional debug summaries.

Future refactors should split only one logical block at a time and must keep
accepted patches identical unless the change is intentionally behavioral.

## Behavior-Sensitive Areas

Do not change these during mechanical cleanup:

- Candidate priority order.
- Scenario list and scenario execution order.
- Hausdorff limit and early-exit semantics.
- `replacePatch` timing.
- Boundary/contact refresh timing after `replacePatch`.
- Geometry repair target selection.
- BCInterface no-contact rejection.
- Subsurface altitude rejection/rescue.
- Smoothing step limit and iteration count.
