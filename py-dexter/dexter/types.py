r"""Type Aliases used throughout the package."""

import numpy as np
from typing import Literal, TypeAlias
from dexter import NcGeometry
from dexter import UnityQfactor, ParabolicQfactor, NcQfactor
from dexter import LarCurrent, NcCurrent

Geometry: TypeAlias = NcGeometry
"""Geometry Objects.

Corresponds to objects that implement the 'Geometry' trait.
"""

Qfactor: TypeAlias = UnityQfactor | ParabolicQfactor | NcQfactor
"""Qfactor Object.

Corresponds to objects that implement the 'Qfactor' trait.
"""

Current: TypeAlias = LarCurrent | NcCurrent
"""Current Object.

Corresponds to objects that implement the 'Current' trait.
"""

FluxCoordinate: TypeAlias = Literal["Toroidal", "Poloidal"]
r"""Magnetic flux coordinates $\psi$ and $\psi_p$."""

NetCDFVersion: TypeAlias = str
"""The netCDF convention version (SemVer)."""

EquilibriumType: TypeAlias = Literal["Numerical", "Analytical"]
"""Described the type of equilibrium the object represents."""

Interp1DType: TypeAlias = Literal[
    "Linear", "Cubic", "Cubic Periodic", "Akima", "Akima Periodic", "Steffen"
]
"""Available 1D Interpolation types (case-insensitive)."""

Interp2DType: TypeAlias = Literal["Bilinear", "Bicubic"]
"""Available 2D Interpolation types (case-insensitive)."""

FluxState: TypeAlias = Literal["Good", "Bad", "None"]
"""Flux coordinate state.

    - Good: values exist and are increasing. Can be used as an x coordinate for evaluations.
    - Bad: values exist but are not monotonic. Can only be used as y-data.
    - None: Variable does not exist or is empty.
"""

ArrayShape: TypeAlias = tuple[int, ...]
"""Shape of a numpy array."""

Array1: TypeAlias = np.ndarray[tuple[int], np.dtype[np.float64]]
"""1D numpy array. """

Array2: TypeAlias = np.ndarray[tuple[int, int], np.dtype[np.float64]]
"""2D numpy array. """
