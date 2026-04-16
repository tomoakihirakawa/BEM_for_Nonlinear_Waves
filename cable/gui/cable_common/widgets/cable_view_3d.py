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
        """Red sphere at A, green sphere at B."""
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

        # Endpoints as small spheres
        radius = _endpoint_radius(positions[0], positions[-1])
        sphere_a = pv.Sphere(radius=radius, center=positions[0])
        sphere_b = pv.Sphere(radius=radius, center=positions[-1])
        actor_a = self._plotter.add_mesh(sphere_a, color=color, render=False)
        actor_b = self._plotter.add_mesh(sphere_b, color=color, render=False)

        # Curve as polyline
        poly = _chain_polyline(positions)
        actor_curve = self._plotter.add_mesh(
            poly, color=color, line_width=4, render=False
        )

        self._result_cables[name] = [actor_a, actor_b, actor_curve]

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
