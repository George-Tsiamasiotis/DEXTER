from typing import TypeAlias

_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/test_netcdf.nc"
_TOROIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/toroidal_test_netcdf.nc"
_POLOIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/poloidal_test_netcdf.nc"

from .objects import LarGeometry, NcGeometry
from .objects import UnityQfactor, ParabolicQfactor, NcQfactor
from .objects import LarCurrent, NcCurrent
from .objects import LarBfield, NcBfield
from .objects import CosHarmonic, NcHarmonic
from .objects import CosPerturbation, NcPerturbation
from .objects import perturbation

Geometry: TypeAlias = LarGeometry | NcGeometry
"""Available 'Geometry' Objects"""

Qfactor: TypeAlias = UnityQfactor | ParabolicQfactor | NcQfactor
"""Available 'Qfactor' Objects"""

Current: TypeAlias = LarCurrent | NcCurrent
"""Available 'Current' Objects"""

Bfield: TypeAlias = LarBfield | NcBfield
"""Available 'Bfield' Objects"""

Harmonic: TypeAlias = CosHarmonic | NcHarmonic
"""Available 'Harmonic' Objects"""

Perturbation: TypeAlias = CosPerturbation | NcPerturbation
"""Available 'Perturbation' Objects"""

__all__ = [
    "_TOROIDAL_TEST_NETCDF_PATH",
    "_POLOIDAL_TEST_NETCDF_PATH",
    "_TEST_NETCDF_PATH",
    "Geometry",
    "Qfactor",
    "Current",
    "Bfield",
    "Harmonic",
    "Perturbation",
    "LarGeometry",
    "NcGeometry",
    "UnityQfactor",
    "ParabolicQfactor",
    "NcQfactor",
    "LarCurrent",
    "NcCurrent",
    "LarBfield",
    "NcBfield",
    "CosHarmonic",
    "NcHarmonic",
    "CosPerturbation",
    "NcPerturbation",
    "perturbation",
]
