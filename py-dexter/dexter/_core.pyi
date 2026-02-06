"""This file mirrors all the definitions made in the dexter-python Rust API.

Note
----

For the Equilibrium objects, ABCs are used for documentation purposes, to mirror the behavior of
the corresponding Traits and avoid re-writing each method. Other than that, they function as
stand-alone objects.
"""

from abc import ABC

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

from dexter.plotters.qfactor_plotter import _QfactorPlotter
from dexter.plotters.current_plotter import _CurrentPlotter

class Geometry(ABC):
    """Geometries base class.

    Defines all the evaluation methods.
    Corresponds to the 'Geometry' Trait from the Rust API.
    """

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

class Qfactor(ABC):
    """Qfactors base class.

    Defines all the evaluation methods.
    Corresponds to the 'Qfactor' Trait from the Rust API.
    """

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

    def psip_of_psi(self, psi: float) -> float:
        r"""The $\psi_p(\psi)$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """

    def psi_of_psip(self, psip: float) -> float:
        r"""The $\psi(\psi_p)$ value.

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

class Current(ABC):
    """Currents base class.

    Defines all the evaluation methods.
    Corresponds to the 'Current' Trait from the Rust API.
    """

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

class NcGeometry(Geometry):
    r"""Object describing the general geometry of a numerical equilibrium.

    Stores relates scalars and arrays, and provides interpolation methods for converting
    between different variables and coordinate systems.

    Attributes
    ----------
    path
        The path of the netCDF file.
    netcdf_version
        The netCDF convention version (SemVer).
    interp1d_type
        The 1D Interpolation type.
    interp2d_type
        The 2D Interpolation type.
    equilibrium_type
        The Equilibrium's type.
    baxis
        The magnetic field strength on the magnetic axis $B_0$ in $[T]$.
    raxis
        The horizontal position of the magnetic axis $R0$ in $[m]$.
    zaxis
        The vertical position of the magnetic axis in $[m]$.
    rgeo
        The geometrical axis (device major radius) in $[m]$.
    rwall
        The value of the $r_{wall}$ coordinate at the wall in $[m]$.
    shape
        The shape of the 2 dimensional data arrays, as in $(len(\psi/\psi_p), len(\theta))$.
    psi_state
        The state of the toroidal flux coordinate.
    psip_state
        The state of the poloidal flux coordinate.
    psi_wall
        The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units.
    psip_wall
        The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units.
    psi_array
        The NetCDF $\psi$ data.
    psip_array
        The NetCDF $\psi_p$ data.
    theta_array
        The NetCDF $\theta$ data,
    theta_array
        The NetCDF $r$ data in $[m]$,
    rlab_array
        The $R_{lab}$ data array.
    zlab_array
        The $Z_{lab}$ data array.
    jacobian_array
        The Jacobian $J$ data array.
    """

    path: str
    netcdf_version: NetCDFVersion
    interp1d_type: Interp1DType
    interp2d_type: Interp2DType
    equilibrium_type: EquilibriumType
    baxis: float
    raxis: float
    zaxis: float
    rgeo: float
    rwall: float
    shape: ArrayShape
    psi_state: FluxState
    psip_state: FluxState
    psip_wall: float
    psi_wall: float
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
    ):
        """Constructs an `NcGeometry`.

        Parameters
        ----------
        path
            The path to the NetCDF file.
        interp1d_type
            The type of 1D Interpolation.
        interp2d_type
            The type of 2D Interpolation.

        Example
        -------
        ```python title="NcGeometry creation"
        >>> geometry = dex.NcGeometry(path, "Cubic", "Bicubic")

        ```
        """

class UnityQfactor(Qfactor, _QfactorPlotter):
    r"""Analytical q-factor profile of $q=1$ and $\psi=\psi_p$.

    Provides methods for calculating the quantities:

    - $q(\psi), q(\psi_p), \psi_p(\psi), \psi(\psi_p), \iota(\psi), \iota(\psi_p)$

    - $d\psi_p(\psi)/d\psi, d\psi(\psi_p)/d\psi_p$


    Attributes
    ----------
    equilibrium_type
        The Equilibrium's type.
    """

    equilibrium_type: EquilibriumType

    def __init__(self) -> None:
        """Constructs a `UnityQfactor`.

        Example
        -------

        ```python title="UnityQfactor creation"
        >>> qfactor = dex.UnityQfactor()

        ```
        """

class ParabolicQfactor(Qfactor, _QfactorPlotter):
    # FIXME:
    r"""Analytical q-factor of parabolic q(ψ) profile.

    Provides methods for calculating the quantities,

    - $q(\psi), q(\psi_p), \psi_p(\psi), \psi(\psi_p), \iota(\psi), \iota(\psi_p)$

    - $d\psi_p(\psi)/d\psi, d\psi(\psi_p)/d\psi_p$

    Described by the following formulas:

    $$
    q(\psi) = q_{axis} + (q_{wall} - q_{axis})
        \bigg( \dfrac{\psi}{\psi_{wall}} \bigg)^2
    $$

    $$
    \psi_p(\psi) = \dfrac{\psi_{wall}}{\sqrt{q_{axis}(q_{wall} - q_{axis})}}
        \arctan\bigg[ \dfrac{\psi\sqrt{q_{wall} - q_{axis}}}{\psi_{wall}\sqrt{q_{axis}}} \bigg]
    $$

    $$
    \psi(\psi_p) = \dfrac{\sqrt{q_{axis}}}{\psi_{wall}\sqrt{q_{wall} - q_{axis}}}
        \tan\bigg[ \dfrac{\sqrt{q_{axis}(q_{wall} - q_{axis})}}{\psi_{wall}}\psi_p \bigg]
    $$

    $$
    \dfrac{d\psi_p(\psi)}{d\psi} =
        \dfrac{\psi_{wall}}{q_{axis}\psi_{wall}^2 + (q_{wall} - q_{axis})\psi^2 }
    $$

    $$
    \dfrac{d\psi(\psi_p)}{d\psi_p} =
        \dfrac{q_{axis}}{\cos^2 \bigg[
        \dfrac{\sqrt{q_{axis}(q_{wall} - q_{axis})}}{\psi_{wall}}\psi_p
        \bigg]}
    $$

    $$
    q(\psi_p) = q(\psi(\psi_p))
    $$

    Attributes
    ----------
    equilibrium_type
        The Equilibrium's type.
    qaxis
        The value of $q$ on the magnetic axis.
    qwall
        The value of $q$ on the wall.
    psi_wall
        The value of the toroidal flux at the wall (if it exists).
    """

    equilibrium_type: EquilibriumType
    qaxis: float
    qwall: float
    psi_wall: float | None

    def __init__(self, qaxis: float, qwall: float, psi_wall: float) -> None:
        """Constructs a `UnityQfactor`.

        Parameters
        ----------
        qaxis
            The value of $q$ on the magnetic axis.
        qwall
            The value of $q$ on the wall.
        psi_wall
            The value of the toroidal flux at the wall (if it exists).

        Example
        -------

        ```python title="UnityQfactor creation"
        >>> qfactor = dex.ParabolicQfactor(qaxis=1.1, qwall=3.8, psi_wall=0.45)

        ```
        """

class NcQfactor(Qfactor, _QfactorPlotter):
    r"""Numerical q-factor reconstructed from a NetCDF file.

    Provides methods for calculating the quantities:

    - $q(\psi), q(\psi_p), \psi_p(\psi), \psi(\psi_p), \iota(\psi), \iota(\psi_p)$

    - $d\psi_p(\psi)/d\psi, d\psi(\psi_p)/d\psi_p$


    Attributes
    ----------
    path
        The path to the NetCDF file.
    netcdf_version
        The netCDF convention version (SemVer).
    interp_type
        The 1D Interpolation type.
    equilibrium_type
        The Equilibrium's type.
    qaxis
        The value of $q$ on the magnetic axis.
    qwall
        The value of $q$ on the wall.
    psi_wall
        The value of the toroidal flux at the wall (if it exists).
    psip_wall
        The value of the poloidal flux at the wall (if it exists).
    psi_state
        The state of the toroidal flux coordinate.
    psip_state
        The state of the poloidal flux coordinate.
    psi_array
        The NetCDF $\psi$ data used to construct the $q(\psi)$ and $I(\psi)$ splines.
    psip_array
        The NetCDF $\psi_p$ data used to construct the $q(\psi_p)$ and $I(\psi_p)$ splines.
    q_array
        The NetCDF $q$ data used to construct the $q(\psi_p)$ spline.
    """

    path: str
    netcdf_version: NetCDFVersion
    interp_type: Interp1DType
    equilibrium_type: EquilibriumType
    qaxis: float
    qwall: float
    psi_wall: float | None
    psip_wall: float | None
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    q_array: Array1

    def __init__(self, path: str, interp_type: Interp1DType) -> None:
        """Constructs an `NcQfactor`.

        Parameters
        ----------
        path
            The path to the NetCDF file.
        interp_type
            The 1D Interpolation type.

        Example
        -------

        ```python title="NcQfactor creation"
        >>> qfactor = dex.NcQfactor(path, "Steffen")

        ```
        """

class LarCurrent(Current, _CurrentPlotter):
    """Analytical Large Aspect Ratio Current with $g=1$ and $I=0$.

    Attributes
    ----------
    equilibrium_type
        The Equilibrium's type.
    """

    equilibrium_type: EquilibriumType

    def __init__(self):
        """Constructs a `LarCurrent`.

        Example
        -------
        ```python title="LarCurrent creation"
        >>> current = dex.LarCurrent()

        ```
        """

class NcCurrent(Current, _CurrentPlotter):
    r"""Numerical plasma reconstructed from a NetCDF file.

    Provides methods for calculating the quantities:

    - $g(\psi), g(\psi_p), I(\psi), I(\psi_p)$

    - $dg(\psi)/d\psi, dg(\psi_p)/d\psi_p, dI(\psi)/d\psi, dI(\psi_p)/d\psi_p$

    Attributes
    ----------
    path
        The path to the NetCDF file.
    netcdf_version
        The netCDF convention version (SemVer).
    interp_type
        The 1D Interpolation type.
    equilibrium_type
        The Equilibrium's type.
    psi_wall
        The value of the toroidal flux at the wall (if it exists).
    psip_wall
        The value of the poloidal flux at the wall (if it exists).
    psi_state
        The state of the toroidal flux coordinate.
    psip_state
        The state of the poloidal flux coordinate.
    psi_array
        The NetCDF $\psi$ data used to construct the $q(\psi)$ and $I(\psi)$ splines.
    psip_array
        The NetCDF $\psi_p$ data used to construct the $q(\psi_p)$ and $I(\psi_p)$ splines.
    g_array
        The NetCDF $g$ data used to construct the $g(\psi_p)$ spline.
    i_array
        The NetCDF $I$ data used to construct the $I(\psi_p)$ spline.
    """

    path: str
    netcdf_version: NetCDFVersion
    interp_type: Interp1DType
    equilibrium_type: EquilibriumType
    psi_wall: float | None
    psip_wall: float | None
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    g_array: Array1
    i_array: Array1

    def __init__(self, path: str, interp_type: Interp1DType) -> None:
        """Constructs an `NcCurrent`.

        Parameters
        ----------
        path
            The path to the NetCDF file.
        interp_type
            The 1D Interpolation type.

        Example
        -------
        ```python title="NcCurrent creation"
        >>> current = dex.NcCurrent(path, "Steffen")

        ```
        """
