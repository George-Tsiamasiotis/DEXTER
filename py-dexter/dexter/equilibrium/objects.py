from dexter._core import _PyNcGeometry
from dexter._core import _PyUnityQfactor, _PyParabolicQfactor, _PyNcQfactor
from dexter._core import _PyLarCurrent, _PyNcCurrent
from .plotters import _FluxPlotter, _QfactorPlotter, _CurrentPlotter, _GeometryPlotter


class NcGeometry(_PyNcGeometry, _FluxPlotter, _GeometryPlotter):
    r"""Object describing the general geometry of a numerical equilibrium.

    Stores relates scalars and arrays, and provides interpolation methods for converting
    between different variables and coordinate systems.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp1d_type
        The type of 1D Interpolation.
    interp2d_type
        The type of 2D Interpolation.

    Attributes
    ----------
    path
        The path of the netCDF file.
    netcdf_version
        The netCDF convention version (SemVer).
    equilibrium_type
        The Equilibrium's type.
    interp1d_type
        The 1D Interpolation type.
    interp2d_type
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
        The shape of the 2 dimensional data arrays, as in $(len(\psi/\psi_p), len(\theta))$.
    psi_wall
        The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units.
    psip_wall
        The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units.
    psi_state
        The state of the toroidal flux coordinate.
    psip_state
        The state of the poloidal flux coordinate.
    psi_array
        The NetCDF $\psi$ data.
    psip_array
        The NetCDF $\psi_p$ data.
    theta_array
        The NetCDF $\theta$ data,
    r_array
        The NetCDF $r$ data in $[m]$,
    rlab_array
        The $R_{lab}$ data array.
    zlab_array
        The $Z_{lab}$ data array.
    jacobian_array
        The Jacobian $J$ data array.

    Example
    -------

    ```python title="NcGeometry creation"
    >>> geometry = dex.NcGeometry(path, "Cubic", "Bicubic")

    ```
    """


# ================================================================================================


class UnityQfactor(_PyUnityQfactor, _FluxPlotter, _QfactorPlotter):
    r"""Analytical q-factor profile of $q=1$ and $\psi=\psi_p$.

    Attributes
    ----------
    equilibrium_type
        The Equilibrium's type.

    Example
    -------
    ```python title="UnityQfactor creation"
    >>> geometry = dex.UnityQfactor()

    ```
    """


class ParabolicQfactor(_PyParabolicQfactor, _FluxPlotter, _QfactorPlotter):
    # FIXME:
    r"""Analytical q-factor of parabolic q(ψ) profile.

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

    Parameters
    ----------
    qaxis
        The value of $q$ on the magnetic axis.
    qwall
        The value of $q$ on the wall.
    psi_wall
        The value of the toroidal flux at the wall.

    Attributes
    ----------
    equilibrium_type
        The Equilibrium's type.
    qaxis
        The value of $q$ on the magnetic axis.
    qwall
        The value of $q$ on the wall.
    psi_wall
        The value of the toroidal flux at the wall, $\psi_{wall}$.
    psip_wall
        The value of the poloidal flux at the wall, $\psi_{p,wall}$.

    Example
    -------

    ```python title="UnityQfactor creation"
    >>> qfactor = dex.ParabolicQfactor(qaxis=1.1, qwall=3.8, psi_wall=0.45)

    ```
    """


class NcQfactor(_PyNcQfactor, _FluxPlotter, _QfactorPlotter):
    r"""Numerical q-factor reconstructed from a NetCDF file.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The 1D Interpolation type.

    Attributes
    ----------
    path
        The path to the NetCDF file.
    netcdf_version
        The netCDF convention version (SemVer).
    equilibrium_type
        The Equilibrium's type.
    interp_type
        The 1D Interpolation type.
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

    Example
    -------

    ```python title="NcQfactor creation"
    >>> qfactor = dex.NcQfactor(path, "Steffen")

    ```
    """


# ================================================================================================


class LarCurrent(_PyLarCurrent, _CurrentPlotter):
    """Analytical Large Aspect Ratio Current with $g=1$ and $I=0$.

    Attributes
    ----------
    equilibrium_type
        The Equilibrium's type.

    Example
    -------
    ```python title="LarCurrent creation"
    >>> geometry = dex.LarCurrent()

    """


class NcCurrent(_PyNcCurrent, _CurrentPlotter):
    r"""Numerical plasma reconstructed from a NetCDF file.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The 1D Interpolation type.

    Attributes
    ----------
    path
        The path to the NetCDF file.
    netcdf_version
        The netCDF convention version (SemVer).
    equilibrium_type
        The Equilibrium's type.
    interp_type
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

    Example
    -------
    ```python title="NcCurrent creation"
    >>> current = dex.NcCurrent(path, "Steffen")

    ```
    """
