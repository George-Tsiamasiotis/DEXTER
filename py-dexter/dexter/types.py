r"""Type Aliases used throughout the package."""

import numpy as np
from typing import Literal, TypeAlias
from collections.abc import Sequence

# =================== Common

Array1: TypeAlias = np.ndarray[tuple[int], np.dtype[np.float64]]
"""1D numpy array."""

Array2: TypeAlias = np.ndarray[tuple[int, int], np.dtype[np.float64]]
"""2D numpy array."""

ArrayShape: TypeAlias = tuple[int, ...]
"""Shape of a numpy array."""

Array: TypeAlias = np.ndarray[ArrayShape, np.dtype[np.float64]]
"""ND numpy array."""

ArrayLike: TypeAlias = float | Array | Sequence
"""Objects that can be converted to arrays, i.e. float, np.ndarray, sequences, ..."""

# =================== Equilibrium

NetCDFVersion: TypeAlias = str
"""The netCDF convention version (SemVer)."""

EquilibriumType: TypeAlias = Literal["Numerical", "Analytical"]
"""An object's equilibrium type."""

Interp1DType: TypeAlias = Literal[
    "Linear", "Cubic", "Cubic Periodic", "Akima", "Akima Periodic", "Steffen"
]
"""Available 1D Interpolation types (case-insensitive)."""

Interp2DType: TypeAlias = Literal["Bilinear", "Bicubic"]
"""Available 2D Interpolation types (case-insensitive)."""

FluxCoordinate: TypeAlias = Literal["Toroidal", "Poloidal"]
r"""Magnetic flux coordinates $\psi$ and $\psi_p$."""

FluxState: TypeAlias = Literal["Good", "Bad", "None"]
"""Flux coordinate state.

    - Good: values exist and are increasing. Can be used as an x coordinate for evaluations.
    - Bad: values exist but are not monotonic. Can only be used as y-data.
    - None: Variable does not exist or is empty.
"""

PhaseMethod: TypeAlias = (
    Literal["Zero", "Average", "Resonance", "Interpolation"]
    | tuple[Literal["Custom"], float]
)
r""" Defines the calculation method of the phase $\phi$ in a Numerical Harmonic.

    - `Zero`: Corresponds to $\phi = 0$.
    - `Average`: Corresponds to $\phi = const =$ the average of all the values of the `phase_array`.
    - `Resonance`: Corresponds to $\phi = const =$ the value of $\phi$ at the resonance $m/n$. In the
      case that the resonance falls outside the last closed flux surface, or does not correspond
      to a valid q-factor value, it defaults to `Zero`.
    - `Interpolation`: Interpolation over the `phase_array`.
    - `Custom(f64)`: Use a custom value for $\phi = const$.
"""

# =================== Simulate

CoordinateSet: TypeAlias = Literal[
    "BoozerToroidal",
    "BoozerPoloidal",
    "MixedToroidal",
    "MixedPoloidal",
]
r""" The kind of InitialConditions set.

    - `BoozerToroidal`: Initial conditions set in the $(t, \psi, \theta, \zeta, \rho, \mu)$ space.
    - `BoozerPoloidal`: Initial conditions set in the $(t, \psi_p, \theta, \zeta, \rho, \mu)$ space.
    - `MixedToroidal`: Initial conditions set in the $(t, P_\zeta, \psi, \theta, \zeta, \mu)$ space.
    - `MixedPoloidal`: Initial conditions set in the $(t, P_\zeta, \psi_p, \theta, \zeta, \mu)$ space.
"""


Intersection: TypeAlias = Literal["ConstZeta", "ConstTheta"]
r""" Defines the surface of the Poincare section.

    - `ConstTheta`: Defines a surface of $\chi_i = \theta$.
    - `ConstZeta`: Defines a surface of $\chi_i = \zeta$.
"""

IntegrationStatus: TypeAlias = Literal[
    "Initialized",
    "PartlyInitialized",
    "InvalidInitialConditions",
    "OutOfBoundsInitialization",
    "Integrated",
    "Intersected",
    "IntersectedTimedOut",
    "InvalidIntersections",
    "Escaped",
    "TimedOut(...)",
    "Failed(...)",
]
r"""The integration status of a Particle.

    - `Initialized`: Initialized by InitialConditions, not integrated.
    - `PartlyInitialized`: InitialConditions have not been fully calculated yet.
    - `InvalidInitialConditions`: Invalid InitialConditions. May occur when using Mixed variables
      with objects that cannot define them, for example Mixed Toroidal coordinates when $g(\psi)$ is
      not defined.
    - `OutOfBoundsInitialization`: InitialConditions where out of bounds.
    - `Integrated`: Reached the end of the integration successfully.
    - `Intersected`: Intersections calculation successful.
    - `ClosedPeriods(..)`: Integrated for a certain amount of θ-ψ periods.
    - `Escaped`: Escaped the last closed flux surface (LCFS).
    - `ModStateEscaped`: Escaped when performing a step on the modified system. This indicates that
      something is wrong in Hénon’s trick implementation.
    - `IntersectedTimedOut`: Calculated some intersections correctly but also timed out.
    - `InvalidIntersections`: Calculated invalid intersections.
    - `TimedOut(...)`: Timed out after a maximum number of steps.
    - `Failed(...)`: Simulation failed for unknown reasons.
"""

OrbitType: TypeAlias = Literal[
    "CoPassing",
    "CuPassing",
    "Trapped",
    "Unclassified",
    "Undefined",
]
r"""A particle's orbit type, calculated through the [`dexter.Particle.close()`] routine.

    - `CoPassing`: $\rho>0$ during the whole integration.
    - `CuPassing`: $\rho<0$ during the whole integration.
    - `Trapped`: $|\theta_{max} - \theta_{min}| < 2\pi$ The $0.999$ is needed since passing particles
      also have a span of $2\pi$ when integrated for a single period. This might cut off some trapped
      particles very close to the separatrix..
    - `Unclassified`: Failed to classify the particle’s orbit into a valid type.
    - `Undefined`: Particle has not been integrated.
"""

SteppingMethod = (
    Literal["EnergyAdaptiveStep", "ErrorAdaptiveStep"]
    | tuple[Literal["FixedStep"], float]
)
"""The stepping method of the solver.

    - `EnergyAdaptiveStep`: Forces the step size to be small enough so that the Energy difference
      from step to step is under a certain threshold. The tolerances can be adjusted with the
      energy_rel_tol and energy_abs_tol fields.
    - `ErrorAdaptiveStep`: Classic RK error estimation : Adjust the step size to minimize the
      local truncation error.
    - `FixedStep(float)`: Fixed step size.
"""

Routine: TypeAlias = Literal["None", "Integrate", "Intersect"]
"""The routine run by a Queue."""
