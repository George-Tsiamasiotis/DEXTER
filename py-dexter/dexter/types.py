r"""Type Aliases used thoughout the package."""

import numpy as np
from typing import Literal, TypeAlias


Interp1DType: TypeAlias = Literal[
    "Linear", "Cubic", "Cubic Periodic", "Akima", "Akima Periodic", "Steffen"
]
"""Availiable 1D Interpolation types (case-insensitive)."""

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
