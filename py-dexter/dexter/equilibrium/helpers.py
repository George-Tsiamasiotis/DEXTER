"""Defines various helper functions associated with the 'dexter-equilibrium' crate."""

from dexter.types import Interp1DType, Interp2DType
from .objects import CosHarmonic, NcHarmonic, CosPerturbation, NcPerturbation
from .objects import (
    NcGeometry,
    NcQfactor,
    NcCurrent,
    NcBfield,
    NcHarmonic,
    NcPerturbation,
)


class NcEquilibrium:
    """Numerical equilibrium from an netCDF file.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp1d_type
        The type of Interpolation of the 1D quantities.
    interp2d_type
        The type of Interpolation of the 2D quantities.

    Attributes
    ----------
    geometry
        The equilibrium's geometry.
    qfactor
        The equilibrium's qfactor.
    current
        The equilibrium's current.
    bfield
        The equilibrium's bfield.
    perturbation
        todo: The equilibrium's perturbation, if they exist.

    Example
    -------

    ```python title="NcEquilibrium creation"
    >>> eq = dex.NcEquilibrium(path, "Cubic", "Bicubic")
    >>> geometry = eq.geometry
    >>> qfactor = eq.qfactor
    >>> current = eq.current
    >>> bfield = eq.bfield
    >>> perturbation = eq.perturbation

    ```
    """

    geometry: NcGeometry
    qfactor: NcQfactor
    current: NcCurrent
    bfield: NcBfield
    perturbation: None

    def __init__(
        self,
        path: str,
        interp1d_type: Interp1DType,
        interp2d_type: Interp2DType,
    ):
        self.geometry = NcGeometry(path, interp1d_type, interp2d_type)
        self.qfactor = NcQfactor(path, interp1d_type)
        self.current = NcCurrent(path, interp1d_type)
        self.bfield = NcBfield(path, interp2d_type)
        self.perturbation = None  # TODO:


def perturbation(
    harmonics: list[CosHarmonic] | list[NcHarmonic],
) -> CosPerturbation | NcPerturbation:
    """Helper function to create a `Perturbation` object.

    The type of `Perturbation` is determined by the type of the passed harmonics.

    Parameters
    ----------
    harmonics
        The harmonics that comprise the perturbation. **Must be of the same type.**

    Returns
    -------
    CosPerturbation | NcPerturbation
        The corresponding `Perturbation` type.

    Example
    -------

    ```python title="NcEquilibrium creation"
    >>> perturbation = dex.perturbation(
    ...     [
    ...         dex.CosHarmonic(alpha=1e-3, m=1, n=3, phase=0),
    ...         dex.CosHarmonic(alpha=1e-3, m=2, n=3, phase=0),
    ...     ]
    ... )

    ```
    """
    match harmonics:
        case []:
            return CosPerturbation([])
        case [*cos] if all([isinstance(harmonic, CosHarmonic) for harmonic in cos]):
            return CosPerturbation(harmonics)  # type: ignore
        case [*nc] if all([isinstance(harmonic, NcHarmonic) for harmonic in nc]):
            return NcPerturbation(harmonics)  # type: ignore
        case _:
            raise TypeError("All harmonics must be of the same type")
