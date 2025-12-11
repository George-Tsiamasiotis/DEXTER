from dexter._core import (
    NcGeometry,
    NcQfactor,
    UnityQfactor,
    NcCurrent,
    LarCurrent,
    NcBfield,
    NcHarmonic,
    NcPerturbation,
)

from importlib.util import find_spec

_LAR_NETCDF_PATH = "./crates/equilibrium/lar_netcdf.nc"

# gtk3agg backend needs PyGObject(gi), which needs a C compiler to be installed.
# gkt4agg spams warnings for no reason
if find_spec("gi") is not None:
    import matplotlib
    import matplotlib.pyplot
    import matplotlib

    matplotlib.use("gtk3agg")
    matplotlib.pyplot.rcParams["text.usetex"] = True

__all__ = [
    "_LAR_NETCDF_PATH",
    "NcGeometry",
    "NcQfactor",
    "UnityQfactor",
    "NcCurrent",
    "LarCurrent",
    "NcBfield",
    "NcHarmonic",
    "NcPerturbation",
]
