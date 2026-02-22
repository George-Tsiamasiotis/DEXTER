"""This file mirrors all the definitions made in the dexter-python Rust API.

Note
----

For the Equilibrium objects, ABCs are used for documentation purposes, to mirror the
behavior of the corresponding Traits and avoid re-writing each method. The final
object is the direct export of the wrapped pyo3 object.

Any added functionality on the final objects should be defined and documented in
their corresponding sub-package.

Attributes are documented on the wrapper for documentation building purposes.

Note
----

Types annotated as `Any` means that the actual type is not yet defined in the package
and therefore cannot be imported (circular import). The actual type is defined in the final
wrapper.
"""

from dexter.types import (
    EquilibriumType,
    FluxWall,
    Interp1DType,
    Interp2DType,
    FluxState,
    Array1,
    Array2,
    ArrayShape,
    PhaseMethod,
    NetCDFVersion,
    InitialFlux,
    IntegrationStatus,
    SteppingMethod,
)

from abc import ABC
from typing import Any, Optional

class _FluxCommuteTrait(ABC):
    """Documents the methods provided by the 'FluxCommute' trait."""

    def psip_of_psi(self, psi: float) -> float:
        r"""The $\psi_p(\psi)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """

    def psi_of_psip(self, psip: float) -> float:
        r"""The $\psi(\psi_p)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

class _GeometryTrait(ABC):
    """Documents the methods provided by the 'Geometries' trait."""

    def r_of_psi(self, psi: float) -> float:
        r"""The $r(\psi)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """

    def r_of_psip(self, psip: float) -> float:
        r"""The $r(\psi_p)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def psi_of_r(self, r: float) -> float:
        r"""The $\psi(r)$ value in Normalized Units.

        Parameters
        ----------
        r
            The radial distance $r$ in $[m]$.
        """

    def psip_of_r(self, r: float) -> float:
        r"""The $\psi_p(r)$ value in Normalized Units.

        Parameters
        ----------
        r
            The radial distance $r$ in $[m]$.
        """

    def rlab_of_psi(self, psi: float, theta: float) -> float:
        r"""The $R_{lab}(\psi, \theta)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

    def rlab_of_psip(self, psip: float, theta: float) -> float:
        r"""The $R_{lab}(\psi_p, \theta)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

    def zlab_of_psi(self, psi: float, theta: float) -> float:
        r"""The $Z_{lab}(\psi, \theta)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

    def zlab_of_psip(self, psip: float, theta: float) -> float:
        r"""The $Z_{lab}(\psi_p, \theta)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

    def jacobian_of_psi(self, psi: float, theta: float) -> float:
        r"""The $J(\psi, \theta)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

    def jacobian_of_psip(self, psip: float, theta: float) -> float:
        r"""The $J(\psi_p, \theta)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

class _QfactorTrait(ABC):
    """Documents the methods provided by the 'Qfactor' trait."""

    def q_of_psi(self, psi: float) -> float:
        r"""The $q(\psi)$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """

    def q_of_psip(self, psip: float) -> float:
        r"""The $q(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def dpsi_dpsip(self, psip: float) -> float:
        r"""The derivative $d\psi(\psi_p)/d\psi_p$ value.

        It's a good check that the values coincide with `qfactor.q_of_psip(psip)`.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def dpsip_dpsi(self, psi: float) -> float:
        r"""The derivative $d\psi_p(\psi)/d\psi$ value.

        It's a good check that the values coincide with `qfactor.iota_of_psi(psi)`.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """

    def iota_of_psi(self, psi: float) -> float:
        r"""The $\iota(\psi) = \dfrac{1}{q(\psi)}$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """

    def iota_of_psip(self, psip: float) -> float:
        r"""The $\iota(\psi_p) = \dfrac{1}{q(\psi_p)}$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

class _CurrentTrait(ABC):
    """Documents the methods provided by the 'Current' trait."""

    def g_of_psi(self, psi: float) -> float:
        r"""The $g(\psi)$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """

    def g_of_psip(self, psip: float) -> float:
        r"""The $g(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def i_of_psi(self, psi: float) -> float:
        r"""The $I(\psi)$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """

    def i_of_psip(self, psip: float) -> float:
        r"""The $I(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def dg_dpsi(self, psi: float) -> float:
        r"""The $dg(\psi)/d\psi$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """

    def dg_dpsip(self, psip: float) -> float:
        r"""The $dg(\psi_p)/d\psi_p$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi$ in Normalized Units.
        """

    def di_dpsi(self, psi: float) -> float:
        r"""The $dg(\psi)/d\psi$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """

    def di_dpsip(self, psip: float) -> float:
        r"""The $dg(\psi_p)/d\psi_p$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi$ in Normalized Units.
        """

class _BfieldTrait(ABC):
    """Documents the methods provided by the 'Bfield' trait."""

    def b_of_psi(self, psi: float, theta: float) -> float:
        r"""The $B(\psi, \theta)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

    def b_of_psip(self, psip: float, theta: float) -> float:
        r"""The $B(\psi_p, \theta)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

    def db_dpsi(self, psi: float, theta: float) -> float:
        r"""The $dB(\psi, \theta)/d\psi$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

    def db_dpsip(self, psip: float, theta: float) -> float:
        r"""The $dB(\psi_p, \theta)/d\psi_p$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

    def db_of_psi_dtheta(self, psi: float, theta: float) -> float:
        r"""The $dB(\psi, \theta)/d\theta$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

    def db_of_psip_dtheta(self, psip: float, theta: float) -> float:
        r"""The $dB(\psi_p, \theta)/d\theta$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

class _HarmonicTrait(ABC):
    """Documents the methods provided by the 'Harmonic' trait."""

    def alpha_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float:
        r"""The harmonic's amplitude $\alpha(\psi, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def alpha_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float:
        r"""The harmonic's amplitude $\alpha(\psi_p, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def phase_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float:
        r"""The harmonic's phase $\phi(\psi, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def phase_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float:
        r"""The harmonic's phase $\phi(\psi_p, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def h_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float:
        r"""The full harmonic's value $h(\psi, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def h_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float:
        r"""The full harmonic's value $h(\psi_p, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dh_dpsi(self, psi: float, theta: float, zeta: float, t: float) -> float:
        r"""The harmonic's derivative with respect to $\psi$, $\partial h(\psi, \theta, \zeta, t)/\partial\psi$
        in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dh_dpsip(self, psip: float, theta: float, zeta: float, t: float) -> float:
        r"""The harmonic's derivative with respect to $\psi_p$, $\partial h(\psi_p, \theta, \zeta, t)/\partial \psi_p$
        in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dh_of_psi_dtheta(
        self, psi: float, theta: float, zeta: float, t: float
    ) -> float:
        r"""The harmonic's derivative with respect to $\theta$, $\partial h(\psi, \theta, \zeta, t)/\partial \theta$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dh_of_psip_dtheta(
        self, psip: float, theta: float, zeta: float, t: float
    ) -> float:
        r"""The harmonic's derivative with respect to $\theta$, $\partial h(\psi_p, \theta, \zeta, t)/\partial \theta$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dh_of_psi_dzeta(self, psi: float, theta: float, zeta: float, t: float) -> float:
        r"""The harmonic's derivative with respect to $\zeta$, $\partial h(\psi, \theta, \zeta, t)/\partial \zeta$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dh_of_psip_dzeta(
        self, psip: float, theta: float, zeta: float, t: float
    ) -> float:
        r"""The harmonic's derivative with respect to $\zeta$, $\partial h(\psi_p, \theta, \zeta, t)/\partial \zeta$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dh_of_psi_dt(self, psi: float, theta: float, zeta: float, t: float) -> float:
        r"""The harmonic's derivative with respect to the time $t$, $\partial h(\psi, \theta, \zeta, t)/\partial t$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dh_of_psip_dt(self, psip: float, theta: float, zeta: float, t: float) -> float:
        r"""The harmonic's derivative with respect to the time $t$, $\partial h(\psi_p, \theta, \zeta, t)/\partial t$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

class _PyPerturbation(ABC):
    """Documents all the methods provided by rust's `Perturbation` type.

    Since `Perturbation` is a concrete type generic over `Harmonic`, there is not 'evaluation' trait and all
    methods are documented here. This class contains the full behavior of *every* wrapper type, except the
    constructors.
    """

    def p_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float:
        r"""The perturbation's value $p(\psi, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def p_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float:
        r"""The perturbation's value $p(\psi_p, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dp_dpsi(self, psi: float, theta: float, zeta: float, t: float) -> float:
        r"""The perturbation's derivative with respect to $\psi$, $\partial p(\psi, \theta, \zeta, t)/\partial\psi$
        in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dp_dpsip(self, psip: float, theta: float, zeta: float, t: float) -> float:
        r"""The perturbation's derivative with respect to $\psi_p$, $\partial p(\psi_p, \theta, \zeta, t)/\partial \psi_p$
        in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dp_of_psi_dtheta(
        self, psi: float, theta: float, zeta: float, t: float
    ) -> float:
        r"""The perturbation's derivative with respect to $\theta$, $\partial p(\psi, \theta, \zeta, t)/\partial \theta$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dp_of_psip_dtheta(
        self, psip: float, theta: float, zeta: float, t: float
    ) -> float:
        r"""The perturbation's derivative with respect to $\theta$, $\partial p(\psi_p, \theta, \zeta, t)/\partial \theta$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dp_of_psi_dzeta(self, psi: float, theta: float, zeta: float, t: float) -> float:
        r"""The perturbation's derivative with respect to $\zeta$, $\partial p(\psi, \theta, \zeta, t)/\partial \zeta$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dp_of_psip_dzeta(
        self, psip: float, theta: float, zeta: float, t: float
    ) -> float:
        r"""The perturbation's derivative with respect to $\zeta$, $\partial p(\psi_p, \theta, \zeta, t)/\partial \zeta$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dp_of_psi_dt(self, psi: float, theta: float, zeta: float, t: float) -> float:
        r"""The perturbation's derivative with respect to the time $t$, $\partial p(\psi, \theta, \zeta, t)/\partial t$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def dp_of_psip_dt(self, psip: float, theta: float, zeta: float, t: float) -> float:
        r"""The perturbation's derivative with respect to the time $t$, $\partial p(\psi_p, \theta, \zeta, t)/\partial t$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """

    def __getitem__(self, index: int):
        """Makes the object indexable."""

    def __len__(self) -> int:
        """Returns the number of the contained harmonics."""

# ================================================================================================

class _PyLarGeometry(_GeometryTrait):
    """PyO3 export of `LarGeometry`. Contains the full behavior of the wrapped object."""

    equilibrium_type: EquilibriumType
    baxis: float
    raxis: float
    rwall: float
    psi_wall: float
    rlab_wall: Array1
    zlab_wall: Array1

    def __init__(self, baxis: float, raxis: float, rwall: float) -> None: ...

class _PyNcGeometry(_FluxCommuteTrait, _GeometryTrait):
    """PyO3 export of `NcGeometry`. Contains the full behavior of the wrapped object."""

    path: str
    netcdf_version: NetCDFVersion
    equilibrium_type: EquilibriumType
    interp1d_type: Interp1DType
    interp2d_type: Interp2DType
    baxis: float
    raxis: float
    zaxis: float
    rgeo: float
    rwall: float
    shape: ArrayShape
    psi_wall: float
    psip_wall: float
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    theta_array: Array1
    r_array: Array1
    rlab_array: Array2
    zlab_array: Array2
    jacobian_array: Array2
    rlab_wall: Array1
    zlab_wall: Array1

    def __init__(
        self,
        path: str,
        interp1d_type: Interp1DType,
        interp2d_type: Interp2DType,
    ) -> None: ...

# ================================================================================================

class _PyUnityQfactor(_FluxCommuteTrait, _QfactorTrait):
    """PyO3 export of `UnityQfactor`. Contains the full behavior of the wrapped object."""

    equilibrium_type: EquilibriumType

    def __init__(self): ...

class _PyParabolicQfactor(_FluxCommuteTrait, _QfactorTrait):
    """PyO3 export of `ParabolicQfactor`. Contains the full behavior of the wrapped object."""

    equilibrium_type: EquilibriumType
    qaxis: float
    qwall: float
    psi_wall: float
    psip_wall: float

    def __init__(self, qaxis: float, qwall: float, flux_wall: FluxWall) -> None: ...

class _PyNcQfactor(_FluxCommuteTrait, _QfactorTrait):
    """PyO3 export of `NcQfactor`. Contains the full behavior of the wrapped object."""

    path: str
    netcdf_version: NetCDFVersion
    equilibrium_type: EquilibriumType
    interp_type: Interp1DType
    qaxis: float
    qwall: float
    psi_wall: float
    psip_wall: float
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    q_array: Array1

    def __init__(self, path: str, interp_type: Interp1DType) -> None: ...

# ================================================================================================

class _PyLarCurrent(_CurrentTrait):
    """PyO3 export of `LarCurrent`. Contains the full behavior of the wrapped object."""

    equilibrium_type: EquilibriumType

    def __init__(self) -> None: ...

class _PyNcCurrent(_CurrentTrait):
    """PyO3 export of `NcCurrent`. Contains the full behavior of the wrapped object."""

    path: str
    netcdf_version: NetCDFVersion
    equilibrium_type: EquilibriumType
    interp_type: Interp1DType
    psi_wall: float
    psip_wall: float
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    g_array: Array1
    i_array: Array1

    def __init__(self, path: str, interp_type: Interp1DType) -> None: ...

# ================================================================================================

class _PyLarBfield(_BfieldTrait):
    """PyO3 export of `LarBfield`. Contains the full behavior of the wrapped object."""

    equilibrium_type: EquilibriumType

    def __init__(self) -> None: ...

class _PyNcBfield(_BfieldTrait):
    """PyO3 export of `NcBfield`. Contains the full behavior of the wrapped object."""

    path: str
    netcdf_version: NetCDFVersion
    equilibrium_type: EquilibriumType
    interp_type: Interp1DType
    baxis: float
    shape: ArrayShape
    psi_wall: float
    psip_wall: float
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    theta_array: Array1
    b_array: Array2

    def __init__(
        self,
        path: str,
        interp_type: Interp2DType,
    ) -> None: ...

# ================================================================================================

class _PyCosHarmonic(_HarmonicTrait):
    """PyO3 export of `CosHarmonic`. Contains the full behavior of the wrapped object."""

    equilibrium_type: EquilibriumType
    alpha: float
    phase: float
    m: int
    n: int

    def __init__(self, alpha: float, m: int, n: int, phase: float) -> None: ...

class _PyNcHarmonic(_HarmonicTrait):
    """PyO3 export of `NcHarmonic`. Contains the full behavior of the wrapped object."""

    path: str
    netcdf_version: NetCDFVersion
    equilibrium_type: EquilibriumType
    interp_type: Interp1DType
    m: int
    n: int
    phase_method: PhaseMethod
    phase_average: float
    psi_phase_resonance: float
    psip_phase_resonance: float
    psi_wall: float
    psip_wall: float
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    alpha_array: Array1
    phase_array: Array1

    def __init__(
        self,
        path: str,
        interp_type: Interp1DType,
        m: int,
        n: int,
        phase_method: PhaseMethod = "Zero",
    ) -> None: ...

# ================================================================================================

class _PyCosPerturbation(_PyPerturbation):
    """Perturbation with H=CosHarmonic. Contains the full behavior of the wrapped object."""

    harmonics: list[Any]

    def __init__(self, harmonics: list[Any]) -> None: ...

class _PyNcPerturbation(_PyPerturbation):
    """Perturbation with H=NcHarmonic. Contains the full behavior of the wrapped object."""

    harmonics: list[Any]

    def __init__(self, harmonics: list[Any]) -> None: ...

# ================================================================================================
# ================================================================================================

class _PyInitialConditions:
    """Particle initial conditions. Contains the full behavior of the wrapped object."""

    t0: float
    flux0: InitialFlux
    theta0: float
    zeta0: float
    rho0: float
    mu0: float

    def __init__(
        self,
        t0: float,
        flux0: InitialFlux,
        theta0: float,
        zeta0: float,
        rho0: float,
        mu0: float,
    ): ...

class _PyParticle:
    """Particle. Contains the full behavior of the wrapped object."""

    initial_conditions: Any
    integration_status: IntegrationStatus
    steps_taken: int
    steps_stored: int
    initial_energy: float | None
    final_energy: float | None
    energy_var: float | None
    t_array: Array1
    psi_array: Array1
    psip_array: Array1
    theta_array: Array1
    zeta_array: Array1
    rho_array: Array1
    mu_array: Array1
    ptheta_array: Array1
    pzeta_array: Array1
    energy_array: Array1

    def __init__(self, initial_conditions: Any): ...
    def integrate(
        self,
        /,
        qfactor: Any,
        current: Any,
        bfield: Any,
        perturbation: Any,
        teval: tuple[float, float],
        *,
        method: Optional[SteppingMethod],
        max_steps: Optional[int],
        first_step: Optional[float],
        safety_factor: Optional[float],
        energy_rel_tol: Optional[float],
        energy_abs_tol: Optional[float],
        error_rel_tol: Optional[float],
        error_abs_tol: Optional[float],
    ): ...
    # # fmt: off
    # Dynamic dispatch for integrate()
    def __int_uniQ_larC_larB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_uniQ_larC_larB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_uniQ_larC_ncdB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_uniQ_larC_ncdB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_uniQ_ncdC_larB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_uniQ_ncdC_larB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_uniQ_ncdC_ncdB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_uniQ_ncdC_ncdB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_parQ_larC_larB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_parQ_larC_larB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_parQ_larC_ncdB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_parQ_larC_ncdB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_parQ_ncdC_larB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_parQ_ncdC_larB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_parQ_ncdC_ncdB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_parQ_ncdC_ncdB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_ncdQ_larC_larB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_ncdQ_larC_larB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_ncdQ_larC_ncdB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_ncdQ_larC_ncdB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_ncdQ_ncdC_larB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_ncdQ_ncdC_larB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_ncdQ_ncdC_ncdB_cosP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
    def __int_ncdQ_ncdC_ncdB_ncdP(self, qfactor: Any, current: Any, bfield: Any, perturbation: Any): ...
