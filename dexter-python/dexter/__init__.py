from dexter._core import Geometry, Qfactor, Currents, Bfield, Harmonic, Perturbation
from dexter._core import (
    InitialConditions,
    MappingParameters,
    Evolution,
    Frequencies,
    Particle,
)
from dexter._core import HeapInitialConditions, Heap

from importlib.util import find_spec

# gtk3agg backend needs PyGObject(gi), which needs a C compiler to be installed.
# gkt4agg spams warnings for no reason
if find_spec("gi") is not None:
    import matplotlib
    import matplotlib.pyplot
    import matplotlib

    matplotlib.use("gtk3agg")
    matplotlib.pyplot.rcParams["text.usetex"] = True

__all__ = [
    "Geometry",
    "Qfactor",
    "Currents",
    "Bfield",
    "Harmonic",
    "Perturbation",
    "InitialConditions",
    "MappingParameters",
    "Evolution",
    "Frequencies",
    "Particle",
    "HeapInitialConditions",
    "Heap",
]
