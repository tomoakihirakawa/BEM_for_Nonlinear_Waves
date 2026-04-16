"""pycable — Python GUI wrapper around the CableDynamics C++ cable_solver.

Public API (stable):
    CableParams       — input parameters dataclass
    CableBridge       — QObject wrapper around the cable_solver binary
    find_cable_solver — locate the cable_solver binary on disk
"""

from .params import CableParams
from .bridge import CableBridge
from .solver_discovery import find_cable_solver

__version__ = "0.1.0"

__all__ = ["CableParams", "CableBridge", "find_cable_solver", "__version__"]
