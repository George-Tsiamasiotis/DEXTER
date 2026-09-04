"""Defines the different bfield objects as wrappers over `_PyBfield`."""

import numpy as np
from semver import Version
from typing import TypeAlias

from dexter._core import _PyBfield

from .utils import MagneticFlux
from .base import MachineObject, Bfield
from dexter.types import (
    Array1,
    Array2,
    ArrayShape,
    FluxCoordinateState,
    Interpolation1dType,
    Interpolation2dType,
    NetCDFVersion,
)


class LarBfield(MachineObject, Bfield):
    r"""Analytical Large Aspect Ratio magnetic field with
    $B(\psi,\theta) = 1-\sqrt{2\psi}\cos\theta$.

    No $\psi/\psi_p$ bounds checks are performed in evaluations.

    Example
    -------
    ``` py
    >>> bfield = dex.LarBfield()

    ```
    """

    _r: _PyBfield

    def __init__(self) -> None:
        self._r = _PyBfield.build_lar()
        MachineObject.__init__(self)
        Bfield.__init__(self)


class NcBfield(MachineObject, Bfield):
    r"""Numerical magnetic field profile reconstructed from a netCDF file.

    Related quantities are computed by interpolating over the data arrays.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The 2D interpolation type.
    padding
        Sets the left-right $\theta$ padding width.

        !!! note "Magnetic field $\theta$ padding"

            At the grid edges, the interpolator’s higher derivatives are not well defined. By
            left-right padding the $B$ array with extra $\psi=const$ columns, we force the
            interpolator to take $\theta$’s periodicity into account and therefore calculate the
            correct derivative values. Note that in contrast to the one-dimensional cubic spline,
            in a bicubic interpolation 3 columns are not enough to ensure periodicity, since the
            spline coefficients depend on the values of the whole array.

            According to [this](https://stackoverflow.com/a/25106574/32596387) stack overflow
            thread, the effect of the $i$-th column at the $j$-th column of the spline scales as
            $r^{|i-j|}$, where $r=\sqrt{3}-2 \approx -0.26$. Therefore, with a padding of 10, the
            effect at the $\theta=0$ boundary would be of the order of $10^{-6}$.

    Attributes
    ----------
    path
        The path to the netCDF file.
    netcdf_version
        The netCDF file's version (SemVer).
    interp_type
        The 2D interpolation type.
    baxis
        The magnetic field strength on the axis $B_0$ in $[T]$.
    padding
        The number of $\theta$ padding columns (per side).
    shape
        Returns the $(\psi/\psi_p,\theta)$ shape of the initial 1D arrays (before the padding).
    shape_padded
        Returns the $(\psi/\psi_p,\theta)$ shape of the **padded** arrays that were used to create the interpolator.
    psi_array
        The toroidal flux' values.
    psip_array
        The poloidal flux' values.
    theta_array
        The poloidal angle $\theta$ values.
    theta_array_padded
        The **padded** $\theta$ values.
    b_array
        The magnetic field $B$ values.
    b_array_padded
        The **padded** $B$ values.

    Example
    -------
    ``` py
    >>> bfield = dex.NcBfield(path, "Bicubic", padding=20)

    ```
    """

    _r: _PyBfield
    path: str
    netcdf_version: NetCDFVersion
    interp_type: Interpolation2dType
    baxis: float
    padding: float
    shape: ArrayShape
    shape_padded: ArrayShape

    def __init__(
        self,
        path: str,
        interp_type: Interpolation2dType,
        padding: int = 15,
    ) -> None:
        self._r = _PyBfield.build_nc(path, interp_type, padding)
        MachineObject.__init__(self)
        Bfield.__init__(self)

        self.path = self._r.path
        self.netcdf_version = Version.parse(self._r.netcdf_version)
        self.interp_type = self._r.interp_type
        self.baxis = self._r.baxis
        self.padding = self._r.padding
        self.shape = self._r.shape
        self.shape_padded = self._r.shape_padded

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
    def theta_array_padded(self) -> Array1:
        array = self._r.get_array("theta_array_padded")
        assert array is not None, "theta_array_padded always exists"
        return array

    @property
    def b_array(self) -> Array2:
        return self._r.get_array2d("b_array")

    @property
    def b_array_padded(self) -> Array2:
        return self._r.get_array2d("b_array_padded")


BfieldObject: TypeAlias = LarBfield | NcBfield
r"""Available [`Bfield`][dexter.machine.base.Bfield] objects.

+ [`LarBfield`][dexter.LarBfield]: Analytical Large Aspect Ratio magnetic field with $B(\psi,\theta) = 1-\sqrt{2\psi}\cos\theta$.
+ [`NcBfield`][dexter.NcBfield]: Numerical magnetic field profile reconstructed from a netCDF file.
"""
