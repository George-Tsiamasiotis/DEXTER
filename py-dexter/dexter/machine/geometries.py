"""Defines the different geometry objects as wrappers over `_PyGeometry`."""

import numpy as np
from semver import Version
from typing import TypeAlias

from dexter._core import _PyGeometry

from .utils import MagneticFlux
from .base import MachineObject, Geometry
from dexter.types import (
    Array1,
    Array2,
    ArrayShape,
    FluxCoordinateState,
    Interpolation1dType,
    Interpolation2dType,
    NetCDFVersion,
)


class LarGeometry(MachineObject, Geometry):
    r"""Analytical Large Aspect Ratio Geometry of a circular device.

    Parameters
    ----------
    baxis
        The magnetic field strength on the axis $B_0$ in $[T]$.
    raxis
        The horizontal position of the magnetic axis $R_0$ in $[m]$.
    rlast
        The $r$ coordinate's value at the last closed flux surface in $[m]$.

    Attributes
    ----------
    psi_last
        The value of the toroidal flux at the last closed flux surface, in Normalized units.

    Notes
    -----
    + No $\psi/\psi_p$ bounds checks are performed in evaluations.

    + Evaluations methods that calculate or accept $\psi_p$ as a parameter always raise an Exception, since
    $\psi_p$ is defined through the q-factor.

    + The Jacobian is not available, since it is defined through $q$, $g$, $I$ and $B$.

    + In LAR equilibria, it holds that $R_0 \equiv R_{geo}$.

    + The definitions are not very strict at the moment.

    Example
    -------
    ``` py
    >>> geometry = dex.LarGeometry(
    ...     baxis=2, # Tesla
    ...     raxis=1.75, # meters
    ...     rlast=0.5, # meters
    ... )

    ```
    """

    _r: _PyGeometry
    psi_last: MagneticFlux

    def __init__(self, baxis: float, raxis: float, rlast: float) -> None:
        self._r = _PyGeometry.build_lar(baxis, raxis, rlast)
        MachineObject.__init__(self)
        Geometry.__init__(self)

        assert self._r.psi_last is not None, "always defined"
        self.psi_last = MagneticFlux._wrap(self._r.psi_last)


class NcGeometry(MachineObject, Geometry):
    r"""Geometry of a realistic configuration.

    Stores fluxes, angles and lab variables’ data, and provides interpolation methods between them.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp1d_type
        The 1D interpolation type for the 1D quantities.
    interp2d_type
        The 2D interpolation type for the 2D quantities.

    Attributes
    ----------
    path
        The path to the netCDF file.
    netcdf_version
        The netCDF file's version (SemVer).
    interp1d_type
        The 1D interpolation type.
    interp2d_type
        The 2D interpolation type.
    shape
        The $(\psi/\psi_p,\theta)$ shape of the 2D arrays.
    psi_last
        The value of the last closed toroidal flux $\psi_{LCFS}$.
    psip_last
        The value of the last closed toroidal flux $\psi_{LCFS}$.
    psi_array
        The toroidal flux' values.
    psip_array
        The poloidal flux' values.
    theta_array
        The poloidal angle $\theta$ values.
    r_array
        The radial coordinate $r$ values, in $[m]$.
    rlab_array
        The $R$ values, in $[m]$.
    zlab_array
        The $Z$ values, in $[m]$.
    jacobian_array
        The Jacobian $J$ values, in $[m]$.

    Example
    -------
    ``` py
    >>> geometry = dex.NcGeometry(path, "Cubic", "Bicubic")

    ```
    """

    _r: _PyGeometry
    path: str
    netcdf_version: NetCDFVersion
    interp1d_type: Interpolation1dType
    interp2d_type: Interpolation2dType
    shape: ArrayShape
    psi_last: MagneticFlux | None
    psip_last: MagneticFlux | None

    def __init__(
        self,
        path: str,
        interp1d_type: Interpolation1dType,
        interp2d_type: Interpolation2dType,
    ) -> None:
        self._r = _PyGeometry.build_nc(path, interp1d_type, interp2d_type)
        MachineObject.__init__(self)
        Geometry.__init__(self)

        self.path = self._r.path
        self.netcdf_version = Version.parse(self._r.netcdf_version)
        self.interp1d_type = self._r.interp1d_type
        self.interp2d_type = self._r.interp2d_type
        self.shape = self._r.shape

        last = self._r.psi_last
        self.psi_last = MagneticFlux._wrap(last) if last is not None else None
        last = self._r.psip_last
        self.psip_last = MagneticFlux._wrap(last) if last is not None else None

    @property
    def psi_array(self) -> Array1 | None:
        return self._r.get_array("psi_array")

    @property
    def psip_array(self) -> Array1 | None:
        return self._r.get_array("psip_array")

    @property
    def theta_array(self) -> Array1:
        array = self._r.get_array("theta_array")
        assert array is not None, "theta_array always exists"
        return array

    @property
    def r_array(self) -> Array1:
        array = self._r.get_array("r_array")
        assert array is not None, "r_array always exists"
        return array

    @property
    def rlab_array(self) -> Array2:
        array = self._r.get_array2d("rlab_array")
        assert array is not None, "rlab_array always exists"
        return array

    @property
    def zlab_array(self) -> Array2:
        array = self._r.get_array2d("zlab_array")
        assert array is not None, "zlab_array always exists"
        return array

    @property
    def jacobian_array(self) -> Array2 | None:
        return self._r.get_array2d("jacobian_array")


GeometryObject: TypeAlias = LarGeometry | NcGeometry
r"""Available [`Geometry`][dexter.machine.base.Geometry] objects.

+ [`LarGeometry`][dexter.LarGeometry]: Analytical Large Aspect Ratio Geometry of a circular device.
+ [`NcGeometry`][dexter.NcGeometry]: Geometry of a realistic configuration.
"""
