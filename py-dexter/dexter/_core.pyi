from dexter.types import (
    Interp1DType,
    Interp2DType,
    FluxState,
    Array1,
    Array2,
    ArrayShape,
)

class NcGeometry:
    r"""Object describing the general geometry of an equilibrium.

    Stores relates scalars and arrays, and provides interpolation methods for converting
    between different variables and coordinate systems.

    Attributes
    ----------
    path
        The path of the netCDF file.
    typ1d
        The 1D Interpolation type.
    typ2d
        The 2D Interpolation type.
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
        The shape of the 2 dimensinal data arrays, as in $(len(\psi/\psi_p), len(\theta))$.
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
    typ1d: Interp1DType
    typ2d: Interp2DType
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

    def __init__(self, path: str, typ1d: Interp1DType, typ2d: Interp2DType):
        """Constructs an `NcGeometry`.

        Parameters
        ----------
        path
            The path to the NetCDF file.
        typ1d
            The type of 1D Interpolation.
        typ2d
            The type of 2D Interpolation.

        Example
        -------


        ```python title="NcGeometry creation"
        >>> geometry = dex.NcGeometry(path, "Cubic", "Bicubic")
        >>>
        >>> # r->ψp interpolation
        >>> psip = geometry.psip_wall / 2
        >>> r_of_psip = geometry.r_of_psip(psip)

        ```
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
        r"""The $J_{lab}(\psi, \theta)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

    def jacobian_of_psip(self, psip: float, theta: float) -> float:
        r"""The $J_{lab}(\psi_p, \theta)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """

class NcCurrent:
    r"""Plasma reconstructed from a NetCDF file.

    Provides methods for calculating the quantities:

    - $g(\psi), g(\psi_p), I(\psi), I(\psi_p)$

    - $dg(\psi)/d\psi, dg(\psi_p)/d\psi_p, dI(\psi)/d\psi, dI(\psi_p)/d\psi_p$

    Attributes
    ----------
    path
        The path to the NetCDF file.
    typ
        The 1D Interpolation type.
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
    typ: Interp1DType
    psi_wall: float | None
    psip_wall: float | None
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    g_array: Array1
    i_array: Array1

    def __init__(self, path: str, typ: Interp1DType) -> None:
        """Constructs an `NcCurrent`.

        Parameters
        ----------
        path
            The path to the NetCDF file.
        typ
            The 1D Interpolation type.

        Example
        -------

        ```python title="NcCurrent creation"
        >>> current = dex.NcCurrent(path, "Steffen")
        >>>
        >>> # ψp->g interpolation
        >>> psip = 0.015
        >>> g_of_psip = current.g_of_psip(psip)

        ```
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

    def __len__(self) -> int:
        """Returns the number of ψp data points."""
