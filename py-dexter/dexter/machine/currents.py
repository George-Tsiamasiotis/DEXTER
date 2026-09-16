"""Defines the different plasma current objects as wrappers over `_PyCurrent`."""

import numpy as np
from semver import Version
from typing import TypeAlias

from dexter._core import _PyCurrent

from dexter.machine.flux import MagneticFlux
from dexter.machine.base import MachineObject, Current
from dexter.types import Array1, FluxCoordinateState, Interpolation1dType, NetCDFVersion


class LarCurrent(MachineObject, Current):
    r"""Analytical Large Aspect Ratio Current with $g=1$ and $I=0$.

    Notes
    -----
    No $\psi/\psi_p$ bounds checks are performed in evaluations.

    Example
    -------
    ``` py
    >>> current = dex.LarCurrent()

    ```
    """

    _r: _PyCurrent

    def __init__(self) -> None:
        self._r = _PyCurrent.build_lar()
        MachineObject.__init__(self)
        Current.__init__(self)


class NcCurrent(MachineObject, Current):
    r"""Numerical plasma current profile reconstructed from a netCDF file.

    Related quantities are computed by interpolating over the data arrays.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The 1D interpolation type.

    Attributes
    ----------
    path
        The path to the netCDF file.
    netcdf_version
        The netCDF file's version (SemVer).
    interp_type
        The 1D interpolation type.
    psi_array
        The toroidal flux' values.
    psip_array
        The poloidal flux' values.
    g_array
        The poloidal current's $g$ values.
    i_array
        The poloidal current's $I$ values.

    Example
    -------
    ``` py
    >>> current = dex.NcCurrent(path, "Cubic")

    ```
    """

    _r: _PyCurrent
    path: str
    netcdf_version: NetCDFVersion
    interp_type: Interpolation1dType

    def __init__(self, path: str, interp_type: Interpolation1dType) -> None:
        self._r = _PyCurrent.build_nc(path, interp_type)
        MachineObject.__init__(self)
        Current.__init__(self)

        self.path = self._r.path
        self.netcdf_version = Version.parse(self._r.netcdf_version)
        self.interp_type = self._r.interp_type

    @property
    def psi_array(self) -> Array1 | None:
        return self._r.get_array("psi_array")

    @property
    def psip_array(self) -> Array1 | None:
        return self._r.get_array("psip_array")

    @property
    def g_array(self) -> Array1:
        array = self._r.get_array("g_array")
        assert array is not None, "g_array always exists"
        return array

    @property
    def i_array(self) -> Array1:
        array = self._r.get_array("i_array")
        assert array is not None, "i_array always exists"
        return array


CurrentObject: TypeAlias = LarCurrent | NcCurrent
"""Available [`Current`][dexter.machine.base.Current] objects.

+ [`LarCurrent`][dexter.LarCurrent]: Analytical Large Aspect Ratio Current with $g=1$ and $I=0$.
+ [`NcCurrent`][dexter.NcCurrent]: Numerical plasma current profile reconstructed from a netCDF file.
"""
