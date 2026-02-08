"""This file mirrors all the definitions made in the dexter-python Rust API.

Note
----

For the Equilibrium objects, ABCs are used for documentation purposes, to mirror the
behavior of the corresponding Traits and avoid re-writing each method. The final
object is the direct export of the wrapped pyo3 object.

Any added functionality on the final objects should be defined and documented in
their corresponding sub-package.

Attributes are documented on the wrapper for documentation building purposes.
"""

from dexter.types import (
    EquilibriumType,
    Interp1DType,
    Interp2DType,
    FluxState,
    Array1,
    Array2,
    ArrayShape,
    NetCDFVersion,
)

from abc import ABC

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

# ================================================================================================

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

class _PyParabolicQfactor(_FluxCommuteTrait, _QfactorTrait):
    """PyO3 export of `ParabolicQfactor`. Contains the full behavior of the wrapped object."""

    equilibrium_type: EquilibriumType
    qaxis: float
    qwall: float
    psi_wall: float
    psip_wall: float

    def __init__(self, qaxis: float, qwall: float, psi_wall: float) -> None: ...

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
    psi_wall: float | None
    psip_wall: float | None
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    g_array: Array1
    i_array: Array1

    def __init__(self, path: str, interp_type: Interp1DType) -> None: ...
