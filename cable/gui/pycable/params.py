"""Cable input parameter types for the cable_solver binary.

Two dataclasses live here:

- ``CableParams`` — a single cable, the original schema. Maps 1:1 to the
  legacy single-line JSON form (`point_a`, `point_b`, `cable_length`, ...)
  read by ``cable/cable_solver.cpp::main``. Backward compatible — existing
  pycable callers and the original test suite still use this.

- ``LumpedCableSystemParams`` — Phase 4 (2026-04-12) addition. Holds a
  list of ``CableParams`` plus shared top-level scalars (gravity, mode,
  solver settings). Round-trips both the legacy single-line JSON
  (when the list has exactly one cable, with the original schema) and
  the new multi-line ``mooring_<name>`` flat format used by the BEM input
  reader and the multi-line `cable_solver.cpp` mode. The two formats are
  auto-detected on read: if the JSON has any `mooring_*` keys, we treat
  it as multi-line; otherwise it's the legacy single-line form.

Both dataclasses emit *natural* JSON (numbers unquoted, arrays mixed),
which the C++ ``parseJSON`` in lib/include/basic.hpp accepts.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Tuple


Vec3 = Tuple[float, float, float]


# ---------------------------------------------------------------------------
# Fluid & Wind helpers (schema mirrors cable_solver.cpp::resolveFluidConfig).
# See cable_common/params.py for the canonical documentation.
# ---------------------------------------------------------------------------

def _fluid_wind_to_dict(p: Any) -> dict[str, Any]:
    d: dict[str, Any] = {}
    if p.fluid != "water":
        d["fluid"] = p.fluid
    if p.fluid_density is not None:
        d["fluid_density"] = float(p.fluid_density)
    if p.drag_Cd is not None:
        d["drag_Cd"] = float(p.drag_Cd)
    if p.wind_type != "none":
        d["wind_type"] = p.wind_type
        d["wind_U_mean"] = list(p.wind_U_mean)
        if p.wind_type == "AR1":
            d["wind_turbulence_intensity"] = float(p.wind_turbulence_intensity)
            d["wind_integral_time_scale"] = float(p.wind_integral_time_scale)
            if p.wind_seed is not None:
                d["wind_seed"] = int(p.wind_seed)
    return d


def _fluid_wind_from_dict(d: dict[str, Any], target: Any) -> None:
    target.fluid = str(d.get("fluid", "water"))
    target.fluid_density = (
        float(d["fluid_density"]) if "fluid_density" in d else None
    )
    target.drag_Cd = float(d["drag_Cd"]) if "drag_Cd" in d else None
    target.wind_type = str(d.get("wind_type", "none"))
    wu = d.get("wind_U_mean", [0.0, 0.0, 0.0])
    target.wind_U_mean = (float(wu[0]), float(wu[1]), float(wu[2]))
    target.wind_turbulence_intensity = float(d.get("wind_turbulence_intensity", 0.15))
    target.wind_integral_time_scale = float(d.get("wind_integral_time_scale", 5.0))
    target.wind_seed = int(d["wind_seed"]) if "wind_seed" in d else None


@dataclass
class CableParams:
    """Parameters for a single cable_solver equilibrium run."""

    # --- Geometry ---
    point_a: Vec3 = (500.0, 0.0, -58.0)
    point_b: Vec3 = (0.0, 0.0, 0.0)
    cable_length: float = 522.0
    n_segments: int = 40

    # --- Material ---
    line_density: float = 348.5      # [kg/m]
    EA: float = 1.4e9                # [N]
    damping: float = 0.5             # artificial relaxation coefficient
    diameter: float = 0.132          # [m]
    gravity: float = 9.81            # [m/s^2]

    # --- Solver ---
    mode: str = "equilibrium"
    max_equilibrium_steps: int = 500_000
    equilibrium_tol: float = 0.01
    snapshot_interval: int = 10_000

    # --- Dynamic mode ---
    dt: float = 0.005
    t_end: float = 5.0
    output_interval: float = 0.0   # 0 = auto (use dt)

    # --- Fluid & Wind (system-wide, applied to all cables) ---
    fluid: str = "water"
    fluid_density: float | None = None
    drag_Cd: float | None = None
    wind_type: str = "none"
    wind_U_mean: Vec3 = (0.0, 0.0, 0.0)
    wind_turbulence_intensity: float = 0.15
    wind_integral_time_scale: float = 5.0
    wind_seed: int | None = None

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """Return a dict matching the cable_solver input JSON schema.

        Numbers are unquoted, arrays are natural. The cable_solver
        parseJSON walks this in string form internally, but on-disk the
        file is human-readable standard JSON.
        """
        d = {
            "point_a": [float(x) for x in self.point_a],
            "point_b": [float(x) for x in self.point_b],
            "cable_length": float(self.cable_length),
            "n_segments": int(self.n_segments),
            "line_density": float(self.line_density),
            "EA": float(self.EA),
            "damping": float(self.damping),
            "diameter": float(self.diameter),
            "gravity": float(self.gravity),
            "mode": str(self.mode),
            "max_equilibrium_steps": int(self.max_equilibrium_steps),
            "equilibrium_tol": float(self.equilibrium_tol),
            "snapshot_interval": int(self.snapshot_interval),
        }
        if self.mode == "dynamic":
            d["dt"] = float(self.dt)
            d["t_end"] = float(self.t_end)
            if self.output_interval > 0:
                d["output_interval"] = float(self.output_interval)
        d.update(_fluid_wind_to_dict(self))
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "CableParams":
        """Build CableParams from a dict, filling missing keys with defaults."""
        defaults = cls()
        obj = cls(
            point_a=tuple(d.get("point_a", defaults.point_a)),
            point_b=tuple(d.get("point_b", defaults.point_b)),
            cable_length=float(d.get("cable_length", defaults.cable_length)),
            n_segments=int(d.get("n_segments", defaults.n_segments)),
            line_density=float(d.get("line_density", defaults.line_density)),
            EA=float(d.get("EA", defaults.EA)),
            damping=float(d.get("damping", defaults.damping)),
            diameter=float(d.get("diameter", defaults.diameter)),
            gravity=float(d.get("gravity", defaults.gravity)),
            mode=str(d.get("mode", defaults.mode)),
            max_equilibrium_steps=int(
                d.get("max_equilibrium_steps", defaults.max_equilibrium_steps)
            ),
            equilibrium_tol=float(d.get("equilibrium_tol", defaults.equilibrium_tol)),
            snapshot_interval=int(d.get("snapshot_interval", defaults.snapshot_interval)),
            dt=float(d.get("dt", defaults.dt)),
            t_end=float(d.get("t_end", defaults.t_end)),
            output_interval=float(d.get("output_interval", defaults.output_interval)),
        )
        _fluid_wind_from_dict(d, obj)
        return obj

    def write_json(self, path: Path | str) -> None:
        """Write the params as a cable_solver-compatible input.json."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def read_json(cls, path: Path | str) -> "CableParams":
        """Load CableParams from a cable_solver input.json."""
        with Path(path).open() as f:
            return cls.from_dict(json.load(f))


# ---------------------------------------------------------------------------
# Phase 4 (2026-04-12): Multi-cable system params
# ---------------------------------------------------------------------------

# Default name used when promoting a single CableParams into a one-element
# system (e.g. when reading a legacy single-line JSON).
DEFAULT_SINGLE_CABLE_NAME = "cable"


@dataclass
class CableSpec:
    """One cable inside a ``LumpedCableSystemParams``.

    Mirrors the BEM/cable_solver multi-line flat 13-element schema:
    ``[name, ax, ay, az, bx, by, bz, total_length, n_points, line_density,
    EA, damping, diameter]``.

    Note ``n_points`` (= ``n_segments + 1``), matching the BEM JSON form.
    """

    name: str = DEFAULT_SINGLE_CABLE_NAME
    point_a: Vec3 = (500.0, 0.0, -58.0)
    point_b: Vec3 = (0.0, 0.0, 0.0)
    cable_length: float = 522.0
    n_points: int = 41  # = n_segments + 1
    line_density: float = 348.5
    EA: float = 1.4e9
    damping: float = 0.5
    diameter: float = 0.132

    @property
    def n_segments(self) -> int:
        return max(0, int(self.n_points) - 1)

    def to_flat_array(self) -> list[Any]:
        """Encode as the 13-element BEM mooring flat array."""
        return [
            str(self.name),
            float(self.point_a[0]), float(self.point_a[1]), float(self.point_a[2]),
            float(self.point_b[0]), float(self.point_b[1]), float(self.point_b[2]),
            float(self.cable_length),
            int(self.n_points),
            float(self.line_density),
            float(self.EA),
            float(self.damping),
            float(self.diameter),
        ]

    @classmethod
    def from_flat_array(cls, value: Iterable[Any]) -> "CableSpec":
        """Decode a BEM mooring flat array (13 elements)."""
        v = list(value)
        if len(v) != 13:
            raise ValueError(f"mooring flat array must have 13 elements, got {len(v)}")
        return cls(
            name=str(v[0]),
            point_a=(float(v[1]), float(v[2]), float(v[3])),
            point_b=(float(v[4]), float(v[5]), float(v[6])),
            cable_length=float(v[7]),
            n_points=int(v[8]),
            line_density=float(v[9]),
            EA=float(v[10]),
            damping=float(v[11]),
            diameter=float(v[12]),
        )

    @classmethod
    def from_cable_params(cls, params: "CableParams",
                          name: str = DEFAULT_SINGLE_CABLE_NAME) -> "CableSpec":
        """Promote a CableParams (legacy single-line) to a CableSpec."""
        return cls(
            name=name,
            point_a=tuple(params.point_a),
            point_b=tuple(params.point_b),
            cable_length=float(params.cable_length),
            n_points=int(params.n_segments) + 1,
            line_density=float(params.line_density),
            EA=float(params.EA),
            damping=float(params.damping),
            diameter=float(params.diameter),
        )

    def to_cable_params(self,
                        gravity: float = 9.81,
                        mode: str = "equilibrium",
                        max_equilibrium_steps: int = 500_000,
                        equilibrium_tol: float = 0.01,
                        snapshot_interval: int = 10_000) -> "CableParams":
        """Demote to a CableParams for the legacy single-line GUI / form."""
        return CableParams(
            point_a=tuple(self.point_a),
            point_b=tuple(self.point_b),
            cable_length=float(self.cable_length),
            n_segments=int(self.n_segments),
            line_density=float(self.line_density),
            EA=float(self.EA),
            damping=float(self.damping),
            diameter=float(self.diameter),
            gravity=float(gravity),
            mode=str(mode),
            max_equilibrium_steps=int(max_equilibrium_steps),
            equilibrium_tol=float(equilibrium_tol),
            snapshot_interval=int(snapshot_interval),
        )


@dataclass
class LumpedCableSystemParams:
    """A collection of cables sharing top-level solver scalars.

    Round-trips two on-disk JSON formats:

    1. **Legacy single-line** — when ``cables`` has exactly one entry and
       the cable's name is the default ``"cable"``. Output uses the
       historical flat keys (``point_a``, ``point_b``, ``cable_length``, ...).
       Loaded by ``cable_solver.cpp`` legacy path.

    2. **Multi-line BEM-compatible** — when ``cables`` has more than one
       entry, or when the user explicitly requested the multi-line format.
       Output uses ``mooring_<name>`` flat 13-element arrays. Loaded by
       ``cable_solver.cpp``'s multi-line auto-detected path.
    """

    cables: list[CableSpec] = field(default_factory=list)

    # Top-level shared scalars
    gravity: float = 9.81
    mode: str = "equilibrium"
    max_equilibrium_steps: int = 500_000
    equilibrium_tol: float = 0.01
    snapshot_interval: int = 10_000

    # Dynamic mode
    dt: float = 0.005
    t_end: float = 5.0
    output_interval: float = 0.0  # 0 = auto

    # Fluid & Wind (system-wide, applied to all cables by cable_solver.cpp)
    fluid: str = "water"
    fluid_density: float | None = None
    drag_Cd: float | None = None
    wind_type: str = "none"
    wind_U_mean: Vec3 = (0.0, 0.0, 0.0)
    wind_turbulence_intensity: float = 0.15
    wind_integral_time_scale: float = 5.0
    wind_seed: int | None = None

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_cable_params(cls, params: "CableParams",
                          name: str = DEFAULT_SINGLE_CABLE_NAME) -> "LumpedCableSystemParams":
        """Wrap a single legacy CableParams as a one-cable system."""
        return cls(
            cables=[CableSpec.from_cable_params(params, name=name)],
            gravity=float(params.gravity),
            mode=str(params.mode),
            max_equilibrium_steps=int(params.max_equilibrium_steps),
            equilibrium_tol=float(params.equilibrium_tol),
            snapshot_interval=int(params.snapshot_interval),
            dt=float(params.dt),
            t_end=float(params.t_end),
            output_interval=float(params.output_interval),
            fluid=params.fluid,
            fluid_density=params.fluid_density,
            drag_Cd=params.drag_Cd,
            wind_type=params.wind_type,
            wind_U_mean=params.wind_U_mean,
            wind_turbulence_intensity=params.wind_turbulence_intensity,
            wind_integral_time_scale=params.wind_integral_time_scale,
            wind_seed=params.wind_seed,
        )

    def is_single_legacy(self) -> bool:
        """Return True if this system can be safely written as legacy single-line."""
        return len(self.cables) == 1 and self.cables[0].name == DEFAULT_SINGLE_CABLE_NAME

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """Render to a dict in either single-line or multi-line schema.

        Routing rule:
        - Empty system: top-level scalars only (degenerate)
        - Single cable with default name "cable": legacy single-line schema
        - Otherwise: multi-line ``mooring_<name>`` schema
        """
        if not self.cables:
            d = {
                "gravity": float(self.gravity),
                "mode": str(self.mode),
                "max_equilibrium_steps": int(self.max_equilibrium_steps),
                "equilibrium_tol": float(self.equilibrium_tol),
                "snapshot_interval": int(self.snapshot_interval),
            }
            d.update(_fluid_wind_to_dict(self))
            return d

        if self.is_single_legacy():
            # Legacy single-line schema (matches CableParams.to_dict()).
            cp = self.cables[0].to_cable_params(
                gravity=self.gravity,
                mode=self.mode,
                max_equilibrium_steps=self.max_equilibrium_steps,
                equilibrium_tol=self.equilibrium_tol,
                snapshot_interval=self.snapshot_interval,
            )
            cp.dt = self.dt
            cp.t_end = self.t_end
            cp.output_interval = self.output_interval
            cp.fluid = self.fluid
            cp.fluid_density = self.fluid_density
            cp.drag_Cd = self.drag_Cd
            cp.wind_type = self.wind_type
            cp.wind_U_mean = self.wind_U_mean
            cp.wind_turbulence_intensity = self.wind_turbulence_intensity
            cp.wind_integral_time_scale = self.wind_integral_time_scale
            cp.wind_seed = self.wind_seed
            return cp.to_dict()

        # Multi-line BEM-compatible schema.
        d: dict[str, Any] = {}
        for cable in self.cables:
            d[f"mooring_{cable.name}"] = cable.to_flat_array()
        d["gravity"] = float(self.gravity)
        d["mode"] = str(self.mode)
        d["max_equilibrium_steps"] = int(self.max_equilibrium_steps)
        d["equilibrium_tol"] = float(self.equilibrium_tol)
        d["snapshot_interval"] = int(self.snapshot_interval)
        if self.mode == "dynamic":
            d["dt"] = float(self.dt)
            d["t_end"] = float(self.t_end)
            if self.output_interval > 0:
                d["output_interval"] = float(self.output_interval)
        d.update(_fluid_wind_to_dict(self))
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "LumpedCableSystemParams":
        """Auto-detect single-line vs multi-line vs per-cable and build.

        Detection order:
        1. ``end_a_position`` present → per-cable format (single cable)
        2. Any ``mooring_*`` / ``cable_*`` key with 13-element list → multi-line
        3. Otherwise → legacy single-line (``point_a`` / ``point_b``)
        """
        defaults = cls()

        # --- Per-cable format ---
        if "end_a_position" in d:
            pc = PerCableParams.from_dict(d)
            spec = CableSpec(
                name=pc.name,
                point_a=tuple(pc.end_a_position),
                point_b=tuple(pc.end_b_position),
                cable_length=pc.cable_length,
                n_points=pc.n_points,
                line_density=pc.line_density,
                EA=pc.EA,
                damping=pc.damping,
                diameter=pc.diameter,
            )
            obj = cls(
                cables=[spec],
                gravity=pc.gravity,
                mode=pc.mode,
                max_equilibrium_steps=pc.max_equilibrium_steps,
                equilibrium_tol=pc.equilibrium_tol,
                snapshot_interval=pc.snapshot_interval,
                dt=pc.dt,
                t_end=pc.t_end,
                output_interval=pc.output_interval,
            )
            _fluid_wind_from_dict(d, obj)
            return obj

        # --- Multi-line format ---
        multi_keys: list[tuple[str, list[Any]]] = []
        for key, value in d.items():
            if not isinstance(value, list):
                continue
            if (key.startswith("mooring_") or key.startswith("cable_")) and len(value) == 13:
                multi_keys.append((key, value))

        if multi_keys:
            cables = [CableSpec.from_flat_array(v) for _, v in multi_keys]
            obj = cls(
                cables=cables,
                gravity=float(d.get("gravity", defaults.gravity)),
                mode=str(d.get("mode", defaults.mode)),
                max_equilibrium_steps=int(d.get("max_equilibrium_steps", defaults.max_equilibrium_steps)),
                equilibrium_tol=float(d.get("equilibrium_tol", defaults.equilibrium_tol)),
                snapshot_interval=int(d.get("snapshot_interval", defaults.snapshot_interval)),
                dt=float(d.get("dt", defaults.dt)),
                t_end=float(d.get("t_end", defaults.t_end)),
                output_interval=float(d.get("output_interval", defaults.output_interval)),
            )
            _fluid_wind_from_dict(d, obj)
            return obj

        # --- Legacy single-line format ---
        legacy = CableParams.from_dict(d)
        return cls.from_cable_params(
            legacy,
            name=str(d.get("name", DEFAULT_SINGLE_CABLE_NAME)),
        )

    def write_json(self, path: Path | str) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def read_json(cls, path: Path | str) -> "LumpedCableSystemParams":
        """Read a cable JSON file (any format) or a settings.json.

        Detection order:
        1. ``input_files`` key → settings mode: load each referenced file
        2. ``end_a_position`` → per-cable
        3. ``mooring_*`` 13-element → multi-line
        4. ``point_a`` → legacy single-line
        """
        path = Path(path)
        with path.open() as f:
            d = json.load(f)

        # --- Settings mode ---
        if "input_files" in d:
            input_dir = path.parent
            all_cables: list[CableSpec] = []
            for fname in d["input_files"]:
                cable_path = input_dir / fname
                # Skip body files (BEM convention: type contains "RigidBody")
                if cable_path.exists():
                    with cable_path.open() as f_body:
                        sub_d = json.load(f_body)
                    if "RigidBody" in str(sub_d.get("type", "")):
                        continue
                sub = cls.read_json(cable_path)
                all_cables.extend(sub.cables)
            defaults = cls()
            obj = cls(
                cables=all_cables,
                gravity=float(d.get("gravity", defaults.gravity)),
                mode=str(d.get("mode", defaults.mode)),
                max_equilibrium_steps=int(d.get("max_equilibrium_steps", defaults.max_equilibrium_steps)),
                equilibrium_tol=float(d.get("equilibrium_tol", defaults.equilibrium_tol)),
                snapshot_interval=int(d.get("snapshot_interval", defaults.snapshot_interval)),
                dt=float(d.get("dt", defaults.dt)),
                t_end=float(d.get("t_end", defaults.t_end)),
                output_interval=float(d.get("output_interval", defaults.output_interval)),
            )
            _fluid_wind_from_dict(d, obj)
            return obj

        return cls.from_dict(d)


# ---------------------------------------------------------------------------
# Per-Cable JSON format (new canonical form)
# ---------------------------------------------------------------------------

@dataclass
class PerCableParams:
    """Self-contained single-cable JSON — the new canonical format.

    Each cable file declares its own endpoints, material properties, and
    (optionally) body labels and motion prescriptions. Both BEM and
    cable_solver read this format.

    Keys map 1:1 to flat JSON keys accepted by ``cable_solver.cpp``'s
    per-cable detection path (``end_a_position`` presence).
    """

    name: str = "cable"

    # --- Endpoints ---
    end_a_position: Vec3 = (0.0, 0.0, 0.0)
    end_b_position: Vec3 = (0.0, 0.0, 0.0)
    end_a_body: str = ""      # BEM linkage label ("" = world / fixed)
    end_b_body: str = ""

    # --- Endpoint motion (Phase B: dynamic mode) ---
    end_a_motion: str = "fixed"
    end_b_motion: str = "fixed"
    # sinusoidal motion params (used when motion != "fixed")
    end_a_motion_dof: str = ""
    end_a_motion_amplitude: float = 0.0
    end_a_motion_frequency: float = 0.0
    end_a_motion_phase: float = 0.0
    end_b_motion_dof: str = ""
    end_b_motion_amplitude: float = 0.0
    end_b_motion_frequency: float = 0.0
    end_b_motion_phase: float = 0.0

    # --- Initial condition ---
    initial_condition: str = "length"   # "length" or "tension"
    tension: float = 0.0               # legacy: single target tension [N]
    tension_top: float = 0.0           # target top-end tension [N]
    tension_bottom: float = 0.0        # target bottom-end tension [N]

    # --- Cable properties ---
    cable_length: float = 100.0
    n_points: int = 41
    line_density: float = 348.5    # kg/m
    EA: float = 1.4e9              # N
    damping: float = 0.5
    diameter: float = 0.132        # m

    # --- Solver ---
    gravity: float = 9.81
    mode: str = "equilibrium"
    max_equilibrium_steps: int = 500_000
    equilibrium_tol: float = 0.01
    snapshot_interval: int = 10_000

    # --- Dynamic mode (Phase B) ---
    dt: float = 0.0
    t_end: float = 0.0
    output_interval: float = 0.0

    # --- Fluid & Wind (system-wide) ---
    fluid: str = "water"
    fluid_density: float | None = None
    drag_Cd: float | None = None
    wind_type: str = "none"
    wind_U_mean: Vec3 = (0.0, 0.0, 0.0)
    wind_turbulence_intensity: float = 0.15
    wind_integral_time_scale: float = 5.0
    wind_seed: int | None = None

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """Render as a flat dict matching the per-cable JSON schema."""
        d: dict[str, Any] = {}
        d["name"] = self.name
        d["end_a_position"] = list(self.end_a_position)
        d["end_b_position"] = list(self.end_b_position)
        if self.end_a_body:
            d["end_a_body"] = self.end_a_body
        if self.end_b_body:
            d["end_b_body"] = self.end_b_body
        if self.end_a_motion != "fixed":
            d["end_a_motion"] = self.end_a_motion
            if self.end_a_motion_dof:
                d["end_a_motion_dof"] = self.end_a_motion_dof
                d["end_a_motion_amplitude"] = self.end_a_motion_amplitude
                d["end_a_motion_frequency"] = self.end_a_motion_frequency
                d["end_a_motion_phase"] = self.end_a_motion_phase
        if self.end_b_motion != "fixed":
            d["end_b_motion"] = self.end_b_motion
            if self.end_b_motion_dof:
                d["end_b_motion_dof"] = self.end_b_motion_dof
                d["end_b_motion_amplitude"] = self.end_b_motion_amplitude
                d["end_b_motion_frequency"] = self.end_b_motion_frequency
                d["end_b_motion_phase"] = self.end_b_motion_phase
        if self.initial_condition != "length":
            d["initial_condition"] = self.initial_condition
        if self.initial_condition == "tension":
            if self.tension_top > 0:
                d["tension_top"] = self.tension_top
            if self.tension_bottom > 0:
                d["tension_bottom"] = self.tension_bottom
            if self.tension > 0 and self.tension_top == 0 and self.tension_bottom == 0:
                d["tension"] = self.tension
        d["cable_length"] = self.cable_length
        d["n_points"] = self.n_points
        d["line_density"] = self.line_density
        d["EA"] = self.EA
        d["damping"] = self.damping
        d["diameter"] = self.diameter
        d["gravity"] = self.gravity
        d["mode"] = self.mode
        d["max_equilibrium_steps"] = self.max_equilibrium_steps
        d["equilibrium_tol"] = self.equilibrium_tol
        d["snapshot_interval"] = self.snapshot_interval
        if self.mode == "dynamic":
            d["dt"] = self.dt
            d["t_end"] = self.t_end
            if self.output_interval > 0:
                d["output_interval"] = self.output_interval
        d.update(_fluid_wind_to_dict(self))
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "PerCableParams":
        """Build from a flat dict (per-cable JSON)."""
        def _vec3(key: str, default: Vec3) -> Vec3:
            v = d.get(key)
            if v is not None and len(v) >= 3:
                return (float(v[0]), float(v[1]), float(v[2]))
            return default

        obj = cls(
            name=str(d.get("name", "cable")),
            end_a_position=_vec3("end_a_position", (0, 0, 0)),
            end_b_position=_vec3("end_b_position", (0, 0, 0)),
            end_a_body=str(d.get("end_a_body", "")),
            end_b_body=str(d.get("end_b_body", "")),
            end_a_motion=str(d.get("end_a_motion", "fixed")),
            end_b_motion=str(d.get("end_b_motion", "fixed")),
            end_a_motion_dof=str(d.get("end_a_motion_dof", "")),
            end_a_motion_amplitude=float(d.get("end_a_motion_amplitude", 0)),
            end_a_motion_frequency=float(d.get("end_a_motion_frequency", 0)),
            end_a_motion_phase=float(d.get("end_a_motion_phase", 0)),
            end_b_motion_dof=str(d.get("end_b_motion_dof", "")),
            end_b_motion_amplitude=float(d.get("end_b_motion_amplitude", 0)),
            end_b_motion_frequency=float(d.get("end_b_motion_frequency", 0)),
            end_b_motion_phase=float(d.get("end_b_motion_phase", 0)),
            initial_condition=str(d.get("initial_condition", "length")),
            tension=float(d.get("tension", 0)),
            tension_top=float(d.get("tension_top", 0)),
            tension_bottom=float(d.get("tension_bottom", 0)),
            cable_length=float(d.get("cable_length", 100)),
            n_points=int(d.get("n_points", 41)),
            line_density=float(d.get("line_density", 348.5)),
            EA=float(d.get("EA", 1.4e9)),
            damping=float(d.get("damping", 0.5)),
            diameter=float(d.get("diameter", 0.132)),
            gravity=float(d.get("gravity", 9.81)),
            mode=str(d.get("mode", "equilibrium")),
            max_equilibrium_steps=int(d.get("max_equilibrium_steps", 500_000)),
            equilibrium_tol=float(d.get("equilibrium_tol", 0.01)),
            snapshot_interval=int(d.get("snapshot_interval", 10_000)),
            dt=float(d.get("dt", 0)),
            t_end=float(d.get("t_end", 0)),
            output_interval=float(d.get("output_interval", 0)),
        )
        _fluid_wind_from_dict(d, obj)
        return obj

    def write_json(self, path: Path | str) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def read_json(cls, path: Path | str) -> "PerCableParams":
        with Path(path).open() as f:
            return cls.from_dict(json.load(f))

    # ------------------------------------------------------------------
    # Conversion from/to legacy formats
    # ------------------------------------------------------------------

    @classmethod
    def from_cable_params(cls, p: "CableParams",
                          name: str = "cable") -> "PerCableParams":
        """Convert a legacy CableParams to per-cable format."""
        return cls(
            name=name,
            end_a_position=tuple(p.point_a),
            end_b_position=tuple(p.point_b),
            cable_length=p.cable_length,
            n_points=p.n_segments + 1,
            line_density=p.line_density,
            EA=p.EA,
            damping=p.damping,
            diameter=p.diameter,
            gravity=p.gravity,
            mode=p.mode,
            max_equilibrium_steps=p.max_equilibrium_steps,
            equilibrium_tol=p.equilibrium_tol,
            snapshot_interval=p.snapshot_interval,
            dt=p.dt,
            t_end=p.t_end,
            output_interval=p.output_interval,
            fluid=p.fluid,
            fluid_density=p.fluid_density,
            drag_Cd=p.drag_Cd,
            wind_type=p.wind_type,
            wind_U_mean=p.wind_U_mean,
            wind_turbulence_intensity=p.wind_turbulence_intensity,
            wind_integral_time_scale=p.wind_integral_time_scale,
            wind_seed=p.wind_seed,
        )

    def to_cable_params(self) -> "CableParams":
        """Convert to legacy CableParams (loses body/motion info)."""
        return CableParams(
            point_a=tuple(self.end_a_position),
            point_b=tuple(self.end_b_position),
            cable_length=self.cable_length,
            n_segments=max(0, self.n_points - 1),
            line_density=self.line_density,
            EA=self.EA,
            damping=self.damping,
            diameter=self.diameter,
            gravity=self.gravity,
            mode=self.mode,
            max_equilibrium_steps=self.max_equilibrium_steps,
            equilibrium_tol=self.equilibrium_tol,
            snapshot_interval=self.snapshot_interval,
            dt=self.dt,
            t_end=self.t_end,
            output_interval=self.output_interval,
            fluid=self.fluid,
            fluid_density=self.fluid_density,
            drag_Cd=self.drag_Cd,
            wind_type=self.wind_type,
            wind_U_mean=self.wind_U_mean,
            wind_turbulence_intensity=self.wind_turbulence_intensity,
            wind_integral_time_scale=self.wind_integral_time_scale,
            wind_seed=self.wind_seed,
        )

    @classmethod
    def from_cable_spec(cls, spec: "CableSpec",
                        gravity: float = 9.81,
                        mode: str = "equilibrium",
                        **solver_kwargs: Any) -> "PerCableParams":
        """Convert a CableSpec to per-cable format."""
        return cls(
            name=spec.name,
            end_a_position=tuple(spec.point_a),
            end_b_position=tuple(spec.point_b),
            cable_length=spec.cable_length,
            n_points=spec.n_points,
            line_density=spec.line_density,
            EA=spec.EA,
            damping=spec.damping,
            diameter=spec.diameter,
            gravity=gravity,
            mode=mode,
            **solver_kwargs,
        )
