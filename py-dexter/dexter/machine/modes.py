"""Defines the different Mode objects as wrappers over `_PyMode`."""

import numpy as np
from semver import Version
from typing import TypeAlias

from dexter._core import _PyMode

from .utils import MagneticFlux
from .base import MachineObject, Mode
from dexter.types import (
    Array1,
    FluxCoordinateState,
    Interpolation1dType,
    NetCDFVersion,
    PhaseMethod,
)


class FluteMode(MachineObject, Mode):
    r"""A simple analytical flute mode.


    A flute mode is defined as:

    $$
    m(\psi, \theta, \zeta) = \epsilon\sqrt{\dfrac{\psi}{\psi_{LCFS}}}\cos(m\theta-n\zeta+\phi)
    $$

    Parameters
    ----------
    epsilon
        The modes's "amplitude" $\epsilon$. Corresponds the value of the amplitude at the last
        closed flux surface.
    lcfs
        The Last Closed Flux Surface, with respect to which the mode is defined.
    phase
        The mode's constant phase $\phi$.

    Attributes
    ----------
    epsilon
        The modes's "amplitude" $\epsilon$. Corresponds the value of the amplitude at the last
        closed flux surface.
    lcfs
        The Last Closed Flux Surface, with respect to which the mode is defined.
    phase
        The mode's constant phase $\phi$.

    Example
    -------
    ``` py
    >>> lcfs = dex.MagneticFlux.Toroidal(0.05)
    >>> mode = dex.FluteMode(1e-4, lcfs, 3, 4, 0)

    ```
    """

    _r: _PyMode
    lcfs: MagneticFlux
    epsilon: float
    phase: float

    def __init__(
        self,
        epsilon: float,
        lcfs: MagneticFlux,
        m: int,
        n: int,
        phase: float,
    ) -> None:
        self._r = _PyMode.build_flute(epsilon, lcfs._r, m, n, phase)
        MachineObject.__init__(self)
        Mode.__init__(self)

        self.lcfs = MagneticFlux._wrap(self._r.lcfs)
        self.epsilon = self._r.epsilon
        self.phase = self._r.phase


class NcFluteMode(MachineObject, Mode):
    r"""Single perturbation flute mode from a netCDF file.

    Related quantities are computed by interpolating over the data arrays.

    A numerical flute mode is defined as:

    $$
    m(\psi, \theta, \zeta) = \sum_{m,n} \alpha_{m,n}(\psi, \theta, \zeta)\cos\big(m\theta-n\zeta+\phi(\psi)\big)
    $$

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The 1D interpolation type.
    m
        The poloidal mode number.
    n
        The toroidal mode number.
    phase_method
        The phase $\phi(\psi/\psi_p)$ calculation method.
    analytical_threshold_index
        The modes’s analytical threshold point. Defines the index of the magnetic flux values under
        which to switch to the analytica formula.

        !!! note "Numerical flute mode analytical patching"

            By definition, flute modes must behave like $\approx\sqrt\psi$ close to the axis, and
            therefore their derivative with respect to the flux must go to infinity. This is a
            behavior that splines cannot replicate, resulting to unnatural orbits close to the
            magnetic axis. To solve this, the mode switches to an analytical formula for the values
            of $\psi/\psi_p$ under a certain threshold. The threshold is defined by the flux value
            at the position index of the data array.

            #Formula
            The patch has the form $\beta\sqrt\psi + \gamma$, where $\beta$ and $\gamma$ are
            adjusted in order to ensure continuity of both $\alpha(\psi)$ and its first derivative.
            $\beta$ is calculated first by $\beta = 2\alpha' \sqrt\psi$ to ensure the correct
            value of the derivative $d\alpha/d\psi$ at the patch’s edge. Finally,
            $\gamma = \alpha - \beta\sqrt\psi$ ensures the continuity of $\alpha$ itself.

            Note that sometimes $\gamma$ may become slightly negative, resulting to $\alpha$
            becoming slightly negative extremely close to the axis. However this error should be
            negligible compared to the possible non-continuity of $\alpha$’s higher derivatives
            and/or its deviation from the actual data.

    Attributes
    ----------
    path
        The path to the netCDF file.
    netcdf_version
        The netCDF file's version (SemVer).
    interp_type
        The 1D interpolation type.
    phase_method
        The phase $\phi$ calculation method.
    analytical_threshold_index
        The analytical threshold index.
    phase_average
        The average value of the phase array, if `PhaseMethod == 'Average'`.
    psi_array
        The toroidal flux' values.
    psip_array
        The poloidal flux' values.
    alpha_array
        The amplitude $\alpha$ values.
    phase_array
        The phase $\phi$ values.

    Example
    -------
    ``` py
    >>> mode = dex.NcFluteMode(path, "Cubic", 3, 2, phase_method="Zero")
    >>> mode = dex.NcFluteMode(
    ...     path=path,
    ...     interp_type="Cubic",
    ...     m=3,
    ...     n=2,
    ...     phase_method=("Custom", 1.57),
    ...     analytical_threshold_index=0,
    ... )

    ```
    """

    _r: _PyMode
    path: str
    netcdf_version: NetCDFVersion
    interp_type: Interpolation1dType
    phase_method: PhaseMethod
    analytical_threshold_index: int
    phase_average: float | None

    def __init__(
        self,
        path: str,
        interp_type: Interpolation1dType,
        m: int,
        n: int,
        phase_method: PhaseMethod = "Interpolation",
        analytical_threshold_index: int = 3,
    ) -> None:
        self._r = _PyMode.build_nc(
            path, interp_type, m, n, phase_method, analytical_threshold_index
        )
        MachineObject.__init__(self)
        Mode.__init__(self)

        self.path = self._r.path
        self.netcdf_version = Version.parse(self._r.netcdf_version)
        self.interp_type = self._r.interp_type
        self.phase_method = self._r.phase_method
        self.analytical_threshold_index = self._r.analytical_threshold_index
        self.phase_average = (
            self._r.phase_average if self._r.phase_method == "Average" else None
        )

    @property
    def psi_array(self) -> Array1 | None:
        return self._r.get_array("psi_array")

    @property
    def psip_array(self) -> Array1 | None:
        return self._r.get_array("psip_array")

    @property
    def alpha_array(self) -> Array1:
        array = self._r.get_array("alpha_array")
        assert array is not None, "alpha_array always exists"
        return array

    @property
    def phase_array(self) -> Array1:
        array = self._r.get_array("phase_array")
        assert array is not None, "phase_array always exists"
        return array


ModeObject: TypeAlias = FluteMode | NcFluteMode
r"""Available [`Mode`][dexter.machine.base.Mode] objects.

+ [`FluteMode`][dexter.FluteMode]: A simple analytical flute mode.
+ [`NcFluteMode`][dexter.NcFluteMode]: Single perturbation flute mode from a netCDF file.
"""
