from dexter._core import NcGeometry
from dexter._core import UnityQfactor, ParabolicQfactor, NcQfactor
from dexter._core import LarCurrent, NcCurrent
import matplotlib.pyplot

_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/test_netcdf.nc"
_TOROIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/toroidal_test_netcdf.nc"
_POLOIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/poloidal_test_netcdf.nc"

matplotlib.pyplot.rcParams["text.usetex"] = True

__all__ = [
    "_TOROIDAL_TEST_NETCDF_PATH",
    "_POLOIDAL_TEST_NETCDF_PATH",
    "_TEST_NETCDF_PATH",
    "NcGeometry",
    "UnityQfactor",
    "ParabolicQfactor",
    "NcQfactor",
    "LarCurrent",
    "NcCurrent",
]
