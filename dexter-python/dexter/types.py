r"""Type Aliases used thoughout the package.

Type aliases
------------
Interp1DType
    Availiable 1D Interpolation types (case-insensitive).
Interp2DType
    Availiable 2D Interpolation types (case-insensitive).
NDArrayShape
    Shape of a numpy array.
NDArray1D
    1D numpy array.
NDArray2D
    2D numpy array.
PoincareSection
    The $\theta=const$ or $\zeta=const$ section on which to calculate orbit
    intersections.
IntegrationStatus
    A Particle's integrations status.
OrbitType
    A Particle's orbit type, calculated form its θ-span.
CalculatedFrequency
    Possibly missing value for a particle's calculated $\omega_\theta$, $\omega_\zeta$
    or $q_{kinetic}$.

"""

import numpy as np
from typing import Literal, TypeAlias


Interp1DType: TypeAlias = Literal[
    "linear", "cubic", "cubic periodic", "akima", "akima periodic", "steffen"
]
Interp2DType: TypeAlias = Literal["bilinear", "bicubic"]

NDArrayShape: TypeAlias = tuple[int, ...]
NDArray1D: TypeAlias = np.ndarray[tuple[int], np.dtype[np.float64]]
NDArray2D: TypeAlias = np.ndarray[tuple[int, int], np.dtype[np.float64]]

PoincareSection: TypeAlias = Literal["ConstTheta", "ConstZeta"]

# Enum
IntegrationStatus: TypeAlias = Literal[
    "Initialized",
    "Integrated",
    "Mapped",
    "SinglePeriodIntegrated",
    "Escaped",
    "TimedOut",
    "InvalidIntersections",
    "Failed",
]

# Enum
OrbitType: TypeAlias = Literal["Trapped", "Passing", "Undefined"]

CalculatedFrequency: TypeAlias = float | None
