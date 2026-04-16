"""View3D — PyVista-backed 3D visualization of the cable state.

Responsibilities:
- Show fixed endpoints A and B as colored spheres (single cable working state).
- Show the initial straight-line cable as a faint preview before a run.
- Live-update the cable curve from each SNAPSHOT during a solve.
- Show final per-segment tension as a colormap on the cable.
- Accumulate multiple solved cables as a dict of actor groups, so a
  batch of N cables can be displayed simultaneously in one scene.
"""

from __future__ import annotations

import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from PySide6.QtWidgets import QWidget


# Categorical palette used to color cables when many are shown at once
# (tab10 + 2 from tab20). Long enough for 12 cables.
_CABLE_COLOR_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#bcbd22",
    "#17becf", "#aec7e8", "#ffbb78", "#98df8a",
]


def cable_color_by_index(i: int) -> str:
    """Return a stable hex color for the i-th cable in a batch."""
    return _CABLE_COLOR_PALETTE[i % len(_CABLE_COLOR_PALETTE)]


class View3D(QWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._plotter = QtInteractor(self)
        self._plotter.set_background("white")

        # Orientation marker (XYZ triad) in the lower-left corner so the
        # user can tell which direction is up at a glance. Red = X,
        # green = Y, blue = Z. interactive=False means the widget stays
        # fixed and only the main scene responds to mouse drags.
        self._plotter.add_axes(
            color="black",
            line_width=3,
            xlabel="X",
            ylabel="Y",
            zlabel="Z (up)",
            interactive=False,
        )

        from PySide6.QtWidgets import QVBoxLayout
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._plotter.interactor)

        # Actor handles so we can update in place.
        # -- "Working cable" (single-cable mode: form ↔ active bridge run) --
        self._endpoint_actors: list = []
        self._initial_line_actor = None
        self._live_line_actor = None
        self._tension_actor = None
        self._bounds_visible = False

        # -- "Result cables" (multi-cable mode: each entry owns a fixed
        # color and a list of actors — endpoints + curve — that persist
        # across subsequent run_equilibrium calls). Keys are the cable
        # name (e.g. "C01"), values are a list of vtk actors. --
        self._result_cables: dict[str, list] = {}

        self._camera_reset_done = False

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def show_endpoints(self, point_a, point_b) -> None:
        """No-op: endpoint spheres are disabled (polyline endpoints are self-evident)."""
        return
        for a in self._endpoint_actors:
            self._plotter.remove_actor(a, render=False)
        self._endpoint_actors = []

        radius = _endpoint_radius(point_a, point_b)
        sphere_a = pv.Sphere(radius=radius, center=np.asarray(point_a, float))
        sphere_b = pv.Sphere(radius=radius, center=np.asarray(point_b, float))
        a1 = self._plotter.add_mesh(sphere_a, color="red", render=False)
        a2 = self._plotter.add_mesh(sphere_b, color="green", render=False)
        self._endpoint_actors = [a1, a2]
        self._refresh_bounds()
        self._plotter.render()

    def show_initial_line(self, point_a, point_b, n_segments: int) -> None:
        """Faint straight line between A and B as pre-run preview."""
        if self._initial_line_actor is not None:
            self._plotter.remove_actor(self._initial_line_actor, render=False)
            self._initial_line_actor = None

        A = np.asarray(point_a, float)
        B = np.asarray(point_b, float)
        n = max(2, int(n_segments) + 1)
        pts = np.linspace(A, B, n)
        poly = _chain_polyline(pts)
        self._initial_line_actor = self._plotter.add_mesh(
            poly, color="lightgray", line_width=2, render=False
        )
        self._refresh_bounds()
        self._plotter.render()

    def update_curve(self, positions: np.ndarray) -> None:
        """Replace the live cable curve with new positions. Called on every SNAPSHOT."""
        positions = np.asarray(positions, float)
        if positions.ndim != 2 or positions.shape[1] != 3 or positions.shape[0] < 2:
            return

        if self._live_line_actor is not None:
            self._plotter.remove_actor(self._live_line_actor, render=False)
            self._live_line_actor = None

        # Hide the initial-line preview once a real snapshot is drawn.
        if self._initial_line_actor is not None:
            self._plotter.remove_actor(self._initial_line_actor, render=False)
            self._initial_line_actor = None

        poly = _chain_polyline(positions)
        self._live_line_actor = self._plotter.add_mesh(
            poly, color="steelblue", line_width=4, render=False
        )
        if not self._camera_reset_done:
            self._plotter.reset_camera()
            self._camera_reset_done = True
        self._refresh_bounds()
        self._plotter.render()

    def show_tensions(self, positions: np.ndarray, tensions: np.ndarray) -> None:
        """Color-map the final cable by per-node tension magnitude."""
        positions = np.asarray(positions, float)
        tensions = np.asarray(tensions, float)
        if positions.ndim != 2 or positions.shape[1] != 3 or len(positions) < 2:
            return

        # Remove any live curve so we don't double-draw.
        for handle_name in ("_live_line_actor", "_tension_actor", "_initial_line_actor"):
            actor = getattr(self, handle_name)
            if actor is not None:
                self._plotter.remove_actor(actor, render=False)
                setattr(self, handle_name, None)

        poly = _chain_polyline(positions)
        poly.point_data["tension"] = tensions[: poly.n_points]
        self._tension_actor = self._plotter.add_mesh(
            poly,
            scalars="tension",
            cmap="viridis",
            line_width=6,
            scalar_bar_args={"title": "Tension [N]"},
            render=False,
        )
        self._refresh_bounds()
        self._plotter.render()

    def reset_camera(self) -> None:
        self._plotter.reset_camera()
        self._plotter.render()
        self._camera_reset_done = True

    def clear_for_new_run(self) -> None:
        """Drop live / tension overlays but keep endpoints + initial line."""
        for handle_name in ("_live_line_actor", "_tension_actor"):
            actor = getattr(self, handle_name)
            if actor is not None:
                self._plotter.remove_actor(actor, render=False)
                setattr(self, handle_name, None)
        self._camera_reset_done = False
        self._plotter.render()

    # ------------------------------------------------------------------
    # Multi-cable API — accumulate results from a batch of cables
    # ------------------------------------------------------------------

    def update_cable_fast(
        self,
        name: str,
        positions: np.ndarray,
        color: str,
        scalars: np.ndarray | None = None,
        clim: tuple[float, float] | None = None,
        cmap: str = "coolwarm",
    ) -> None:
        """Fast in-place update of a cable's polyline for animation scrubbing.

        When ``scalars`` is provided, the polyline is colored by that field
        using ``cmap`` colormap instead of the solid ``color``. Recreates the
        actor only when switching between solid-color and scalar modes, or
        on first use; otherwise mutates mesh.points / mesh.point_data in place.
        """
        positions = np.ascontiguousarray(np.asarray(positions, float))
        if positions.ndim != 2 or positions.shape[1] != 3 or positions.shape[0] < 2:
            return

        if not hasattr(self, "_fast_meshes"):
            # name -> (poly, actor, mode, cmap) ; mode ∈ {"solid", "scalar"}
            self._fast_meshes: dict[str, tuple] = {}

        desired_mode = "scalar" if scalars is not None else "solid"
        entry = self._fast_meshes.get(name)

        # Recreate the actor if the rendering mode or colormap changed, or first call.
        old_cmap = entry[3] if entry and len(entry) > 3 else None
        need_recreate = (entry is None or entry[2] != desired_mode
                         or (desired_mode == "scalar" and old_cmap != cmap))
        if need_recreate:
            if entry is not None:
                try:
                    self._plotter.remove_actor(entry[1], render=False)
                except Exception:
                    pass
            poly = _chain_polyline(positions)
            if desired_mode == "scalar":
                poly.point_data["tension"] = np.asarray(scalars, float)[: poly.n_points]
                actor = self._plotter.add_mesh(
                    poly, scalars="tension", cmap=cmap,
                    clim=clim, line_width=4, show_scalar_bar=False, render=False,
                )
            else:
                actor = self._plotter.add_mesh(
                    poly, color=color, line_width=4, render=False
                )
            self._fast_meshes[name] = (poly, actor, desired_mode, cmap)
        else:
            poly, _actor, _mode = entry
            poly.points = positions
            if desired_mode == "scalar":
                poly.point_data["tension"] = np.asarray(scalars, float)[: poly.n_points]
            poly.Modified()

    def clear_fast_cables(self) -> None:
        """Remove all fast-update polylines."""
        if not hasattr(self, "_fast_meshes"):
            return
        for entry in self._fast_meshes.values():
            # entry = (poly, actor, mode)
            try:
                self._plotter.remove_actor(entry[1], render=False)
            except Exception:
                pass
        self._fast_meshes = {}

    def set_static_cable(self, name: str, positions: np.ndarray) -> None:
        """Store a reference 'static equilibrium' cable shape as a light-gray
        background overlay. Call once per cable (typically the first frame
        of the dynamic run)."""
        positions = np.ascontiguousarray(np.asarray(positions, float))
        if positions.ndim != 2 or positions.shape[1] != 3 or positions.shape[0] < 2:
            return
        if not hasattr(self, "_static_meshes"):
            # name -> (positions, actor_or_None)
            self._static_meshes: dict[str, tuple] = {}
        if not hasattr(self, "_static_visible"):
            self._static_visible = True
        # Remove old actor for this name
        entry = self._static_meshes.pop(name, None)
        if entry is not None and entry[1] is not None:
            try:
                self._plotter.remove_actor(entry[1], render=False)
            except Exception:
                pass
        # Create new actor only if currently visible
        actor = None
        if self._static_visible:
            poly = _chain_polyline(positions)
            actor = self._plotter.add_mesh(
                poly, color="lightgray", line_width=2,
                opacity=0.6, render=False,
            )
        self._static_meshes[name] = (positions.copy(), actor)

    def clear_static_cables(self) -> None:
        """Remove all static-equilibrium reference overlays."""
        if not hasattr(self, "_static_meshes"):
            return
        for _pos, actor in self._static_meshes.values():
            if actor is None:
                continue
            try:
                self._plotter.remove_actor(actor, render=False)
            except Exception:
                pass
        self._static_meshes = {}

    def set_static_cables_visible(self, visible: bool) -> None:
        """Show or hide the static-equilibrium reference overlays by actually
        removing/re-creating the actors (the most reliable across VTK wrappers)."""
        if not hasattr(self, "_static_meshes"):
            return
        self._static_visible = bool(visible)
        updated: dict[str, tuple] = {}
        for name, (positions, actor) in self._static_meshes.items():
            # Remove current actor if any
            if actor is not None:
                try:
                    self._plotter.remove_actor(actor, render=False)
                except Exception:
                    pass
            if visible:
                poly = _chain_polyline(positions)
                new_actor = self._plotter.add_mesh(
                    poly, color="lightgray", line_width=2,
                    opacity=0.6, render=False,
                )
                updated[name] = (positions, new_actor)
            else:
                updated[name] = (positions, None)
        self._static_meshes = updated
        self._plotter.render()

    def render_now(self) -> None:
        """Trigger a single render (use after a batch of update_cable_fast)."""
        self._plotter.render()

    def append_cable(
        self,
        name: str,
        positions: np.ndarray,
        color: str,
    ) -> None:
        """Add one solved cable to the accumulated result scene.

        Existing result cables are left untouched. The working-cable
        actors (endpoints / initial line / live curve / tension overlay)
        are cleared so they don't clutter the multi-cable view.

        Parameters
        ----------
        name : str
            Identifier for this cable (e.g. "C01"). If a cable with the
            same name is already in the scene, it is replaced.
        positions : np.ndarray, shape (N, 3)
            Node coordinates from the solver result.
        color : str
            Hex or CSS color string for this cable. Use ``cable_color_by_index``
            to pick from the built-in categorical palette.
        """
        positions = np.asarray(positions, float)
        if positions.ndim != 2 or positions.shape[1] != 3 or positions.shape[0] < 2:
            return

        # Clear single-cable working actors to avoid double-drawing
        self._clear_working_actors()

        # Replace any previous actors under this name
        if name in self._result_cables:
            for a in self._result_cables[name]:
                self._plotter.remove_actor(a, render=False)
            del self._result_cables[name]

        # Curve as polyline (endpoint spheres disabled per user preference)
        poly = _chain_polyline(positions)
        actor_curve = self._plotter.add_mesh(
            poly, color=color, line_width=4, render=False
        )

        self._result_cables[name] = [actor_curve]

        self._refresh_bounds()
        self._plotter.render()

    def show_multi_cable_tensions(
        self,
        cables: list[dict],
    ) -> None:
        """Display multiple cables with a shared tension colormap.

        Parameters
        ----------
        cables : list of dict
            Each dict: ``{"name": str, "positions": ndarray (N,3),
            "tensions": ndarray (N,)}``.
        """
        self._clear_working_actors()
        self.clear_all_result_cables()
        # Also remove fast-update animation cables so the tension display
        # doesn't stack on top of the live scrubber polylines.
        self.clear_fast_cables()

        if not cables:
            return

        # Compute global tension range across all cables
        all_tensions = np.concatenate([c["tensions"] for c in cables])
        t_min = float(np.min(all_tensions))
        t_max = float(np.max(all_tensions))
        if t_max <= t_min:
            t_max = t_min + 1.0  # avoid degenerate range

        # Compute uniform sizes from the overall scene bounding box.
        all_pts = np.concatenate([np.asarray(c["positions"], float) for c in cables])
        scene_span = float(np.linalg.norm(all_pts.max(axis=0) - all_pts.min(axis=0)))
        endpoint_radius = max(0.05, scene_span * 0.004)
        node_point_size = max(3, scene_span * 0.012)

        first = True
        for cable_data in cables:
            name = cable_data["name"]
            positions = np.asarray(cable_data["positions"], float)
            tensions = np.asarray(cable_data["tensions"], float)
            if positions.ndim != 2 or positions.shape[1] != 3 or len(positions) < 2:
                continue

            # Node markers (lumped mass positions) — small dots at each node
            node_cloud = pv.PolyData(positions)
            node_cloud.point_data["tension"] = tensions[: len(positions)]
            actor_nodes = self._plotter.add_mesh(
                node_cloud,
                scalars="tension",
                cmap="viridis",
                clim=[t_min, t_max],
                point_size=node_point_size,
                render_points_as_spheres=True,
                show_scalar_bar=False,
                render=False,
            )

            # Cable curve with shared tension colormap
            poly = _chain_polyline(positions)
            poly.point_data["tension"] = tensions[: poly.n_points]
            actor_curve = self._plotter.add_mesh(
                poly,
                scalars="tension",
                cmap="viridis",
                clim=[t_min, t_max],
                line_width=4,
                show_scalar_bar=first,
                scalar_bar_args={"title": "Tension [N]"} if first else {},
                render=False,
            )
            first = False

            self._result_cables[name] = [actor_nodes, actor_curve]

        self._refresh_bounds()
        self._plotter.render()

    def clear_all_result_cables(self) -> None:
        """Remove every cable accumulated by ``append_cable``."""
        for actors in self._result_cables.values():
            for a in actors:
                self._plotter.remove_actor(a, render=False)
        self._result_cables.clear()
        self._camera_reset_done = False
        self._plotter.render()

    @property
    def result_cable_names(self) -> list[str]:
        """Names currently in the multi-cable result dict."""
        return list(self._result_cables.keys())

    def _clear_working_actors(self) -> None:
        """Remove single-cable working actors (endpoints / line / tension)."""
        for actor in self._endpoint_actors:
            self._plotter.remove_actor(actor, render=False)
        self._endpoint_actors = []
        for handle_name in ("_initial_line_actor", "_live_line_actor", "_tension_actor"):
            actor = getattr(self, handle_name)
            if actor is not None:
                self._plotter.remove_actor(actor, render=False)
                setattr(self, handle_name, None)

    # ------------------------------------------------------------------
    # Bounding box with axis labels and ticks
    # ------------------------------------------------------------------

    def _refresh_bounds(self) -> None:
        """(Re)draw the labeled bounding box around the current scene.

        Called whenever the scene geometry changes. ``show_bounds`` in
        PyVista replaces the previous bounds actor in place, so we can
        safely call this repeatedly.

        If the scene is empty (no actors added yet), the computed bounds
        are degenerate and show_bounds may fail — we guard against that
        with a try/except.
        """
        try:
            self._plotter.show_bounds(
                grid="back",
                location="outer",
                ticks="both",
                xtitle="X [m]",
                ytitle="Y [m]",
                ztitle="Z [m]",
                color="black",
                font_size=10,
                bold=False,
                minor_ticks=False,
                padding=0.05,
            )
        except Exception:
            # Empty scene or PyVista version mismatch — no bounds to draw.
            pass


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def _chain_polyline(points: np.ndarray) -> pv.PolyData:
    """Build a PolyData polyline through ``points`` (shape (N, 3))."""
    n = len(points)
    # PyVista line connectivity: [n, 0, 1, 2, ..., n-1]
    lines = np.empty(n + 1, dtype=np.int64)
    lines[0] = n
    lines[1:] = np.arange(n)
    poly = pv.PolyData(points)
    poly.lines = lines
    return poly


def _endpoint_radius(a, b) -> float:
    """Sphere radius proportional to endpoint separation — auto-scale so the dots are visible."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    d = float(np.linalg.norm(b - a))
    return max(0.02, d * 0.012)
