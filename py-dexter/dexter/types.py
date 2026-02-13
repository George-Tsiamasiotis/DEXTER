r"""Type Aliases used throughout the package."""

import numpy as np
from typing import Literal, TypeAlias

# =================== Common

ArrayShape: TypeAlias = tuple[int, ...]
"""Shape of a numpy array."""

Array1: TypeAlias = np.ndarray[tuple[int], np.dtype[np.float64]]
"""1D numpy array. """

Array2: TypeAlias = np.ndarray[tuple[int, int], np.dtype[np.float64]]
"""2D numpy array. """

# =================== Equilibrium

NetCDFVersion: TypeAlias = str
"""The netCDF convention version (SemVer)."""

FluxCoordinate: TypeAlias = Literal["Toroidal", "Poloidal"]
r"""Magnetic flux coordinates $\psi$ and $\psi_p$."""

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

FluxWall: TypeAlias = tuple[Literal["Toroidal"] | Literal["Poloidal"], float]
"""Helper type to define a ParabolicQfactor with respect to one of the two fluxes’ values at the wall."""

PhaseMethod: TypeAlias = (
    Literal["Zero", "Average", "Resonance", "Interpolation"]
    | tuple[Literal["Custom"], float]
)
r""" Defines the calculation method of the phase $\phi$ in a Numerical Harmonic.

    - Zero: Corresponds to $\phi = 0$.
    - Average: Corresponds to $\phi = const =$ the average of all the values of the `phase_array`.
    - Resonance: Corresponds to $\phi = const =$ the value of $\phi$ at the resonance $m/n$. In the
      case that the resonance falls outside the wall, or does not correspond to a valid q-factor
      value, it defaults to `Zero`.
    - Interpolation: Interpolation over the `phase_array`.
    - Custom(f64): Use a custom value for $\phi = const$.
"""
