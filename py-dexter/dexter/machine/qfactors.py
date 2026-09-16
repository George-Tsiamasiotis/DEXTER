"""Defines the different q-factor objects as wrappers over `_PyQfactor`."""

import numpy as np
from semver import Version
from typing import TypeAlias

from dexter._core import _PyQfactor

from dexter.machine.flux import MagneticFlux
from dexter.machine.base import MachineObject, Qfactor
from dexter.types import (
    Array1,
    FluxCoordinateState,
    Interpolation1dType,
    NetCDFVersion,
)


class UnityQfactor(MachineObject, Qfactor):
    r"""Analytical q-factor profile of $q=1$ and $\psi=\psi_p$.

    Parameters
    ----------
    lcfs
        The last closed flux surface. Only used for bounds checking.

    Example
    -------
    ``` py
    >>> lcfs = dex.MagneticFlux.Toroidal(0.05)
    >>> qfactor = dex.UnityQfactor(lcfs)

    ```
    """

    _r: _PyQfactor

    def __init__(self, lcfs: MagneticFlux) -> None:
        self._r = _PyQfactor.build_unity(lcfs._r)
        MachineObject.__init__(self)
        Qfactor.__init__(self)


class ParabolicQfactor(MachineObject, Qfactor):
    r"""Analytical parabolic q-factor profile.

    Parameters
    ----------
    qaxis
        The q-factor's value at the magnetic axis, $q_{axis}$.
    qlast
        The q-factor's value at the last closed flux surface, $q_{LCFS}$.

    Example
    -------
    ``` py
    >>> lcfs = dex.MagneticFlux.Toroidal(0.05)
    >>> qfactor = dex.ParabolicQfactor(1.1, 3.9, lcfs)

    ```
    """

    _r: _PyQfactor

    def __init__(self, qaxis: float, qlast: float, lcfs: MagneticFlux) -> None:
        self._r = _PyQfactor.build_parabolic(qaxis, qlast, lcfs._r)
        MachineObject.__init__(self)
        Qfactor.__init__(self)


class NcQfactor(MachineObject, Qfactor):
    r"""Numerical q-factor profile reconstructed from a netCDF file.

    Related quantities are computed by interpolating over the data arrays.

    If either psi_norm or psip_norm is missing from the netCDF file, it is calculated from the
    other by integrating $q(\psi_p)$ or $\iota(\psi)$ respectively. In the case that the calculated
    values are monotonic, the other flux can be used as a flux coordinate as well.

    Parameters
    ----------
    path
        The path to the netCDF file.
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
    q_array
        The q-factor's values.

    Example
    -------
    ``` py
    >>> qfactor = dex.NcQfactor(path, "Cubic")

    ```
    """

    _r: _PyQfactor
    path: str
    netcdf_version: NetCDFVersion
    interp_type: Interpolation1dType

    def __init__(self, path: str, interp_type: Interpolation1dType) -> None:
        self._r = _PyQfactor.build_nc(path, interp_type)
        MachineObject.__init__(self)
        Qfactor.__init__(self)

        self.path = self._r.path
        self.netcdf_version = Version.parse(self._r.netcdf_version)
        self.interp_type = self._r.interp_type

    @property
    def psi_array(self) -> Array1:
        return self._r.get_array("psi_array")

    @property
    def psip_array(self) -> Array1:
        return self._r.get_array("psip_array")

    @property
    def q_array(self) -> Array1:
        return self._r.get_array("q_array")


QfactorObject: TypeAlias = UnityQfactor | ParabolicQfactor | NcQfactor
r"""Available [`Qfactor`][dexter.machine.base.Qfactor] objects.

+ [`UnityQfactor`][dexter.UnityQfactor]: Analytical q-factor profile of $q=1$ and $\psi=\psi_p$.
+ [`ParabolicQfactor`][dexter.ParabolicQfactor]: Analytical parabolic q-factor profile.
+ [`NcQfactor`][dexter.NcQfactor]: Numerical q-factor profile reconstructed from a netCDF file.
"""
