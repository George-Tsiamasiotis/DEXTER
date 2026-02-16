from dexter._core import _PyLarGeometry, _PyNcGeometry
from dexter._core import _PyUnityQfactor, _PyParabolicQfactor, _PyNcQfactor
from dexter._core import _PyLarCurrent, _PyNcCurrent
from dexter._core import _PyLarBfield, _PyNcBfield
from dexter._core import _PyCosHarmonic, _PyNcHarmonic
from dexter._core import _PyCosPerturbation, _PyNcPerturbation
from ._plotters import (
    _FluxPlotter,
    _QfactorPlotter,
    _CurrentPlotter,
    _GeometryPlotter,
    _HarmonicPlotter,
)


class LarGeometry(_PyLarGeometry, _GeometryPlotter):
    r"""Analytical Large Aspect Ratio Geometry of a circular device.

    Parameters
    ----------
    baxis
        The magnetic field strength on the magnetic axis $B_0$ in $[T]$.
    raxis
        The horizontal position of the magnetic axis $R_0$ in $[m]$.
    rwall
        The value of the $r_{wall}$ coordinate at the wall in $[m]$.

    Attributes
    ----------
    equilibrium_type
        The Equilibrium's type.
    baxis
        The magnetic field strength on the magnetic axis $B_0$ in $[T]$.
    raxis
        The horizontal position of the magnetic axis $R_0$ in $[m]$.
    rwall
        The value of the $r_{wall}$ coordinate at the wall in $[m]$.
    psi_wall
        The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units.
    rlab_wall
        Last $R_{lab}$ values that correspond to the device's walls.
    zlab_wall
        Last $Z_{lab}$ values that correspond to the device's walls.

    Example
    -------
    ```python title="LarGeometry creation"
    >>> geometry = dex.LarGeometry(baxis=2, raxis=1.75, rwall=0.5)

    ```
    """


class NcGeometry(_PyNcGeometry, _FluxPlotter, _GeometryPlotter):
    r"""Object describing the general geometry of a numerical equilibrium.

    Stores relates scalars and arrays, and provides interpolation methods for converting
    between different variables and coordinate systems.

    Related quantities are computed by interpolating over the data arrays.

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
        The horizontal position of the magnetic axis $R_0$ in $[m]$.
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
    rlab_wall
        Last $R_{lab}$ values that correspond to the device's walls.
    zlab_wall
        Last $Z_{lab}$ values that correspond to the device's walls.

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
    >>> qfactor = dex.UnityQfactor()

    ```
    """


class ParabolicQfactor(_PyParabolicQfactor, _FluxPlotter, _QfactorPlotter):
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
    \psi(\psi_p) = \dfrac{\psi_{wall}\sqrt{q_{axis}}}{\sqrt{q_{wall} - q_{axis}}}
        \tan\bigg[ \dfrac{\sqrt{q_{axis}(q_{wall} - q_{axis})}}{\psi_{wall}}\psi_p \bigg]
    $$

    $$
    \dfrac{d\psi_p(\psi)}{d\psi} = ... = \dfrac{1}{q(\psi)} = \iota(\psi)
    $$

    $$
    \dfrac{d\psi(\psi_p)}{d\psi_p} =
        \dfrac{q_{axis}}{\cos^2 \bigg[
        \dfrac{\sqrt{q_{axis}(q_{wall} - q_{axis})}}{\psi_{wall}}\psi_p
        \bigg]}
        \overset{*}{=}
        q(\psi_p)
    $$

    $$
    q(\psi_p) = q_{axis} + q_{axis} \tan^2
        \bigg[ \dfrac{\sqrt{q_{axis}(q_{wall}-q_{axis})}}{\psi_{wall}} \psi_p \bigg]
    $$

    $^*$ Identity: $\dfrac{1}{\cos^2\theta} = 1 + \tan^2\theta$

    Note
    ----

    A ParabolicQfactor is defined with the help of the [`FluxWall`](types.md/#dexter.types.FluxWall)
    helper type, which changes the position where qwall is met.

    Parameters
    ----------
    qaxis
        The value of $q$ on the magnetic axis.
    qwall
        The value of $q$ on the wall.
    flux_wall
        Helper type to define a ParabolicQfactor with respect to one of the two fluxes’ values at the wall.

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
    >>> # Define q(ψ=ψwall=0.45) = qwall = 3.8
    >>> flux_wall: dex.types.FluxWall = ("Toroidal", 0.45)
    >>> qfactor = dex.ParabolicQfactor(qaxis=1.1, qwall=3.8, flux_wall=flux_wall)

    ```
    """


class NcQfactor(_PyNcQfactor, _FluxPlotter, _QfactorPlotter):
    r"""Numerical q-factor reconstructed from a NetCDF file.

    Related quantities are computed by interpolating over the data arrays.

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
    >>> current = dex.LarCurrent()

    ```
    """


class NcCurrent(_PyNcCurrent, _CurrentPlotter):
    r"""Numerical plasma current reconstructed from a NetCDF file.

    Related quantities are computed by interpolating over the data arrays.

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


# ================================================================================================


class LarBfield(_PyLarBfield):
    """Analytical Large Aspect Ratio magnetic field with B(ψ, θ) = 1 - sqrt(2ψ)cos(θ).

    Attributes
    ----------
    equilibrium_type
        The Equilibrium's type.

    Example
    -------
    ```python title="LarBfield creation"
    >>> bfield = dex.LarBfield()

    ```
    """


class NcBfield(_PyNcBfield):
    r"""Numerical magnetic field profile reconstructed from a netCDF file.

    Related quantities are computed by interpolating over the data arrays.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The type of 2D Interpolation.

    Attributes
    ----------
    path
        The path of the netCDF file.
    netcdf_version
        The netCDF convention version (SemVer).
    equilibrium_type
        The Equilibrium's type.
    interp_type
        The 2D Interpolation type.
    baxis
        The magnetic field strength on the magnetic axis $B_0$ in $[T]$.
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
    b_array
        The NetCDF $B$ data,

    Example
    -------

    ```python title="NcBfield creation"
    >>> bfield = dex.NcBfield(path, "Bicubic")

    ```
    """


# ================================================================================================


class CosHarmonic(_PyCosHarmonic, _HarmonicPlotter):
    r"""A simple analytical Harmonic of the form:

    $$
    \alpha_{\{m,n\}}(\psi/\psi_p, \theta, \zeta, t)
    $$

    where $\alpha$ and $\phi$ are constants.

    with '$\psi/\psi_p$' meaning it can be expressed as a function of either/both flux coordinates.

    Parameters
    ----------
    alpha
        The harmonic's constant amplitude $\alpha$ in Normalized Units.
    m
        The poloidal mode number $m$.
    n
        The poloidal mode number $n$.
    phase
        The harmonic's constant phase $\phi$ in $[rads]$.

    Attributes
    ----------
    equilibrium_type
        The Equilibrium's type.
    alpha
        The harmonic's constant amplitude $\alpha$ in Normalized Units.
    m
        The poloidal mode number $m$.
    n
        The poloidal mode number $n$.
    phase
        The harmonic's constant phase $\phi$ in $[rads]$.

    Example
    -------
    ```python title="CosHarmonic creation"
    >>> harmonic = dex.CosHarmonic(alpha=1e-3, m=3, n=2, phase=0)

    ```
    """


class NcHarmonic(_PyNcHarmonic, _HarmonicPlotter):
    r"""Single perturbation harmonic from a netCDF file.

    The harmonic has the form of:

    $$
    \alpha_{\{m, n\}}(\psi/\psi_p) \cos\big(m\theta - n\zeta + \phi(\psi/\psi_p) \big)
    $$

    with '$\psi/\psi_p$' meaning it can be expressed as a function of either/both flux coordinates.

    where $\alpha$ and $\phi$ can be expressed as functions of either or both $\psi, \psi_p$, and
    are calculated by interpolation over the numerical data.

    !!! info "Phase calculation configuration"

        $\phi$ calculation can be further configured with the `phase_method` optional parameter, which
        defaults to `Resonance`, meaning that a constant value equal to the value of $\phi$ at the
        resonance is used. If no valid value can be found, it falls back to `Zero`. See
        [`PhaseMethod`](types.md/#dexter.types.PhaseMethod)) for available configurations.

        !!! example

            ```python title="Phase calculation configuration"
            >>> harmonic = dex.NcHarmonic(path, "cubic", 3, 2, phase_method = "Average")
            >>> harmonic = dex.NcHarmonic(path, "cubic", 3, 2, phase_method = "Interpolation")
            >>> harmonic = dex.NcHarmonic(path, "cubic", 3, 2, phase_method = ("Custom", 3.1415))

            ...
            ```

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The type of 1D Interpolation.
    m
        The poloidal mode number $m$.
    n
        The poloidal mode number $n$.
    phase_method
        The method of calculation of the phase $\phi(\psi/\psi_p)$. Defaults to "Resonance".


    Attributes
    ----------
    path
        The path of the netCDF file.
    netcdf_version
        The netCDF convention version (SemVer).
    equilibrium_type
        The Equilibrium's type.
    interp_type
        The 1D Interpolation type.
    m
        The poloidal mode number $m$.
    n
        The poloidal mode number $n$.
    phase_method
        The method of calculation of the phase $\phi(\psi/\psi_p)$.
    phase_average
        The average value of the phase arrays, if `phase_method` is `Average`.
    psi_phase_resonance
        The toroidal flux’s value where the resonance is met, if `phase_method` is `Resonance` and
        the resonance is in bounds.
    psip_phase_resonance
        The poloidal flux’s value where the resonance is met, if `phase_method` is `Resonance` and
        the resonance is in bounds.
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
    alpha_array
        The NetCDF $\alpha$ data,
    alpha_array
        The NetCDF $\phi$ data,
    """


# ================================================================================================


class CosPerturbation(_PyCosPerturbation):
    r"""A sum of an arbitrary number of `CosHarmonics`.

    It has the general form of:

    $$
    \sum_{\{m,n\}} \bigg[ \alpha_{\{m,n\}}(\psi/\psi_p, \theta, \zeta, t) \bigg]
    $$

    with '$\psi/\psi_p$' meaning it can be expressed as a function of either/both flux coordinates.

    Parameters
    ----------
    harmonics
        List of the contained harmonics.

    Attributes
    ----------
    harmonics
        List of the contained harmonics.

    Example
    -------
    ```python title="CosPerturbation creation with specific harmonics"
    >>> perturbation = dex.CosPerturbation(
    ...     [
    ...         dex.CosHarmonic(1e-3, 1, 2, 0.0),
    ...         dex.CosHarmonic(1e-3, 1, 3, 0.0),
    ...         dex.CosHarmonic(1e-3, 1, 4, 0.0),
    ...         dex.CosHarmonic(1e-3, 1, 5, 0.0),
    ...     ]
    ... )

    ```

    ```python title="CosPerturbation creation with list iteration"
    >>> perturbation = dex.CosPerturbation(
    ...     [dex.CosHarmonic(1e-3, 1, n, 0.0) for n in range(1, 8)] # modes with m=1 and n=1-7
    ... )

    ```
    """

    harmonics: list[CosHarmonic]

    def __init__(self, harmonics: list[CosHarmonic]) -> None: ...


class NcPerturbation(_PyNcPerturbation):
    r"""A sum of an arbitrary number of `NcHarmonics`.

    It has the general form of:

    $$
    \sum_{\{m,n\}} \bigg[ \Phi_{\{m,n\}}(\psi/\psi_p, \theta, \zeta, t) \bigg] =
    \sum_{\{m,n\}} \bigg[ \alpha_{\{m,n\}}(\psi/\psi_p)
        \cos\big( m\theta - n\zeta + \phi(\psi/\psi_p) \big)\bigg]
    $$


    with '$\psi/\psi_p$' meaning it can be expressed as a function of either/both flux coordinates.

    Parameters
    ----------
    harmonics
        List of the contained harmonics.

    Attributes
    ----------
    harmonics
        List of the contained harmonics.

    Example
    -------
    ```python title="NcPerturbation creation with specific harmonics"
    >>> perturbation = dex.NcPerturbation(
    ...     [
    ...         dex.NcHarmonic(path, "cubic", 2, 1, phase_method="Interpolation"),
    ...         dex.NcHarmonic(path, "cubic", 2, 2, phase_method="Interpolation"),
    ...         dex.NcHarmonic(path, "cubic", 3, 2, phase_method="Interpolation"),
    ...     ]
    ... )

    ```
    """

    harmonics: list[NcHarmonic]

    def __init__(self, harmonics: list[NcHarmonic]) -> None: ...


# ================================================================================================


def perturbation(
    harmonics: list[CosHarmonic] | list[NcHarmonic],
) -> CosPerturbation | NcPerturbation:
    """Helper function to create a `Perturbation` object.

    The type of `Perturbation` is determined by the type of the passed harmonics.

    Parameters
    ----------
    harmonics
        The harmonics that comprise the perturbation. **Must be of the same type.**

    Returns
    -------
    The corresponding `Perturbation` type.
    """
    match harmonics:
        case []:
            return CosPerturbation([])
        case [*cos] if all([isinstance(harmonic, CosHarmonic) for harmonic in cos]):
            return CosPerturbation(harmonics)  # type: ignore
        case [*nc] if all([isinstance(harmonic, NcHarmonic) for harmonic in nc]):
            return NcPerturbation(harmonics)  # type: ignore
        case _:
            raise TypeError("All harmonics must be of the same type")
