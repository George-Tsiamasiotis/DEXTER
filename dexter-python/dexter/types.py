r"""Type Aliases used thoughout the package."""

import numpy as np
from typing import Literal, TypeAlias


Interp1DType: TypeAlias = Literal[
    "Linear", "Cubic", "Cubic Periodic", "Akima", "Akima Periodic", "Steffen"
]
"""Availiable 1D Interpolation types (case-insensitive)."""

Interp2DType: TypeAlias = Literal["Bilinear", "Bicubic"]
"""Availiable 2D Interpolation types (case-insensitive)."""


PhaseMethod: TypeAlias = Literal["Zero", "Average", "Resonance", "Interpolation"]
r"""
Availiable methods for the calculation of an NcHarmonic's phase $\phi(\psi_p)$
(case-insensitive).

- `Zero`: Corresponds to `φ(ψp) = 0`.
- `Average`: Corresponds to φ = const = the average of all the values of the phase array.
- `Resonance`: Corresponds to φ = const = the value of φ(ψp) at the resonance m/n.
- `Interpolation`: Interpolation over the phase array.
"""

NDArrayShape: TypeAlias = tuple[int, ...]
"""Shape of a numpy array."""

NDArray1D: TypeAlias = np.ndarray[tuple[int], np.dtype[np.float64]]
"""1D numpy array. """

NDArray2D: TypeAlias = np.ndarray[tuple[int, int], np.dtype[np.float64]]
"""2D numpy array. """

PoincareSection: TypeAlias = Literal["ConstTheta", "ConstZeta"]
r"""The $\theta=const$ or $\zeta=const$ section on which to calculate orbit intersections."""

# Enum
IntegrationStatus: TypeAlias = Literal[
    "Initialized",
    "Integrated",
    "Mapped",
    "SinglePeriodIntegrated",
    "Escaped",
    "EvaluationNaN",
    "TimedOut",
    "InvalidIntersections",
    "Failed",
]
"""A Particle's integrations status."""

# Enum
OrbitType: TypeAlias = Literal["Trapped", "Passing", "Undefined"]
"""A Particle's orbit type, calculated form its θ-span."""

CalculatedFrequency: TypeAlias = float | None
r"""A particle's calculated $\omega_\theta$, $\omega_\zeta$ or $q_{kinetic}$."""

SingePeriodIntersections: TypeAlias = tuple[int, list[int]]
r"""
The `(number of found intersections, [step of each intersection])` of the variables $\theta$
and $\psi_p$ calculated during SinglePeriod integration.
"""
