from typing import TypeAlias

_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/test_netcdf.nc"
_TOROIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/toroidal_test_netcdf.nc"
_POLOIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/poloidal_test_netcdf.nc"

from .objects import NcGeometry
from .objects import UnityQfactor, ParabolicQfactor, NcQfactor
from .objects import LarCurrent, NcCurrent

Geometry: TypeAlias = NcGeometry
""" Availiable 'Geometry' Objects"""

Qfactor: TypeAlias = UnityQfactor | ParabolicQfactor | NcQfactor
""" Availiable 'Qfactor' Objects"""

Current: TypeAlias = LarCurrent | NcCurrent
""" Availiable 'Current' Objects"""

__all__ = [
    "_TOROIDAL_TEST_NETCDF_PATH",
    "_POLOIDAL_TEST_NETCDF_PATH",
    "_TEST_NETCDF_PATH",
    "Geometry",
    "Qfactor",
    "Current",
    "NcGeometry",
    "UnityQfactor",
    "ParabolicQfactor",
    "NcQfactor",
    "LarCurrent",
    "NcCurrent",
]
