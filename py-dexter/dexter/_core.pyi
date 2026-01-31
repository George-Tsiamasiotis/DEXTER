from dexter.types import (
    Interp1DType,
    FluxState,
    Array1,
)

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
