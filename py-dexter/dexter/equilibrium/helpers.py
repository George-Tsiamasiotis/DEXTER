"""Defines various helper functions associated with the 'dexter-equilibrium' crate."""

from dexter.types import Interp1DType, Interp2DType
from .objects import (
    Bfield,
    Current,
    Geometry,
    NcGeometry,
    NcQfactor,
    NcCurrent,
    NcBfield,
    Qfactor,
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
    >>> (geometry, qfactor, current, bfield) = eq.objects()

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

    def objects(self) -> tuple[NcGeometry, NcQfactor, NcCurrent, NcBfield]:
        """Returns a tuple with the constructed objects."""
        return (self.geometry, self.qfactor, self.current, self.bfield)
