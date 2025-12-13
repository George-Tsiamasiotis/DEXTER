"""This file mirrors all the definitions made in the dexter-python Rust API."""

from dexter.types import (
    NDArray1D,
    NDArray2D,
    Interp1DType,
    Interp2DType,
    NDArrayShape,
    PhaseMethod,
)

class NcGeometry:
    r"""Object describing the geometry of the device.

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
    psip_wall
        The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units.
    psi_wall
        The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units.
    r_wall
        The device's minor radius $r_{wall}$ in $[m]$.
    shape
        The shape of the 2 dimensinal data arrays, as in $(len(\psi_p), len(\theta_B))$.
    theta_data
        The $\theta_B$ data array.
    psip_data
        The $\psi_p$ data array.
    psi_data
        The $\psi$ data array.
    r_data
        The $r(\psi_p)$ data array.
    rlab_data
        The $R_{lab}(\psi_p, \theta_B)$ data array.
    zlab_data
        The $Z_{lab}(\psi_p, \theta_B)$ data array.
    jacobian_data
        The $J(\psi_p, \theta_B)$ data array.

    """

    path: str
    typ1d: Interp1DType
    typ2d: Interp2DType
    baxis: float
    raxis: float
    zaxis: float
    rgeo: float
    psip_wall: float
    psi_wall: float
    r_wall: float
    shape: NDArrayShape
    theta_data: NDArray1D
    psip_data: NDArray1D
    psi_data: NDArray1D
    r_data: NDArray1D
    rlab_data: NDArray2D
    zlab_data: NDArray2D
    jacobian_data: NDArray2D

    def __init__(self, path: str, typ1d: Interp1DType, typ2d: Interp2DType):
        """Constructs a `Geometry`.

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
        >>> r = geometry.r(psip)

        ```
        """

    def r(self, psip: float) -> float:
        r"""The $r(\psi_p)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def psip(self, r: float) -> float:
        r"""The $\psi_p(r)$ value in Normalized Units.

        Parameters
        ----------
        r
            The distance from the (geometrical) magnetic axis $r$ in $[m]$.
        """

    def psi(self, psip: float) -> float:
        r"""The $\psi(\psi_p)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def rlab(self, psip: float, theta: float) -> float:
        r"""The $R_{lab}(\psi_p, \theta_B)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta_B$ angle in $[rads]$.
        """

    def zlab(self, psip: float, theta: float) -> float:
        r"""The $Z_{lab}(\psi_p, \theta_B)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta_B$ angle in $[rads]$.
        """

    def jacobian(self, psip: float, theta: float) -> float:
        r"""The $J(\psi_p, \theta_B)$ value in $[m/T]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta_B$ angle in $[rads]$.
        """

class NcQfactor:
    r"""q-factor from a NetCDF file.

    Provides methods for calculating $q(\psi_p)$, $\psi(\psi_p)$ and $d\psi/d\psi_p$.

    Attributes
    ----------
    path
        The path to the NetCDF file.
    typ
        The 1D Interpolation type.
    psip_data
        The NetCDF $\psi_p$ data used to construct the $q(\psi_p)$ and $\psi(\psi_p)$ splines.
    q_data
        The NetCDF $q$ data used to construct the $q(\psi_p)$ spline.
    psi_data
        The NetCDF $\psi$ data used to construct the $\psi(\psi_p)$ spline.
    """

    path: str
    typ: Interp1DType
    psip_data: NDArray1D
    q_data: NDArray1D
    psi_data: NDArray1D

    def __init__(self, path: str, typ: Interp1DType) -> None:
        """Constructs a `Qfactor`.

        Parameters
        ----------
        path
            The path to the NetCDF file.
        typ
            The 1D Interpolation type.

        Example
        -------

        ```python title="Qfactor creation"
        >>> qfactor = dex.NcQfactor(path, "Steffen")
        >>>
        >>> # ψp->ψ interpolation
        >>> psip = 0.003
        >>> psi = qfactor.psi(psip)

        ```

        ```python title="dψ/dψp check:"
        >>> from math import isclose
        >>>
        >>> psip = 0.4*geometry.psip_wall
        >>> assert isclose(qfactor.q(psip), qfactor.dpsi_dpsip(psip), rel_tol=1e-4)

        ```
        """

    def q(self, psip: float) -> float:
        r"""The $q(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def psi(self, psip: float) -> float:
        r"""The $\psi(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def dpsi_dpsip(self, psip: float) -> float:
        r"""The $q(\psi_p)$ value, as calculated from $d\psi/d\psi_p$ by interpolation
        with the $\psi(\psi_p)$ spline.

        It's a good check that the values coincide with `Qfactor.q(psip)`.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def __len__(self) -> int:
        """Returns the number of ψp data points."""

class UnityQfactor:
    r"""q-factor of $q=1$."""

    def __init__(self) -> None:
        """Constructs a `Qfactor`.

        Example
        -------

        ```python title="UnityQfactor creation"
        >>> qfactor = dex.UnityQfactor()
        >>>
        >>> # ψp->ψ interpolation
        >>> psip = 0.003
        >>> qfactor.q(psip)
        1.0
        >>> qfactor.psi(psip)
        0.003

        ```

        """

    def q(self, psip: float) -> float:
        r"""Always returns 1."""

    def psi(self, psip: float) -> float:
        r"""Always returns `psip`."""

    def dpsi_dpsip(self, psip: float) -> float:
        r"""Always returns 1."""

    def __len__(self) -> int:
        """Returns the number of ψp data points."""

class NcCurrent:
    r"""Plasma reconstructed from a NetCDF file.

    Provides methods for calculating $g(\psi_p)$, $I(\psi_p)$, $dg/d\psi_p$ and
    $dI/d\psi_p$.

    Attributes
    ----------
    path
        The path to the NetCDF file.
    typ
        The 1D Interpolation type.
    psip_data
        The NetCDF $\psi_p$ data used to construct the $q(\psi_p)$ and $\psi(\psi_p)$ splines.
    g_data
        The NetCDF $g$ data used to construct the $g(\psi_p)$ spline.
    i_data
        The NetCDF $I$ data used to construct the $I(\psi_p)$ spline.
    """

    path: str
    typ: Interp1DType
    psip_data: NDArray1D
    g_data: NDArray1D
    i_data: NDArray1D

    def __init__(self, path: str, typ: Interp1DType) -> None:
        """Constructs a `Currents`.

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
        >>> g = current.g(psip)

        ```
        """

    def g(self, psip: float) -> float:
        r"""The $g(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def i(self, psip: float) -> float:
        r"""The $I(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def dg_dpsip(self, psip: float) -> float:
        r"""The $dg(\psi_p)/d\psi_p$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def di_dpsip(self, psip: float) -> float:
        r"""The $dI(\psi_p)/d\psi_p$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def __len__(self) -> int:
        """Returns the number of ψp data points."""

class LarCurrent:
    r"""Large Aspect Ratio plasma current approximation, with $g=1$ and $I=0$."""

    def __init__(self) -> None:
        """Constructs a `LarCurrent`.

        Example
        -------

        ```python
        >>> current = dex.LarCurrent()
        >>>
        >>> # ψp->g interpolation
        >>> psip = 0.015
        >>> current.g(psip)
        1.0
        >>> current.i(psip)
        0.0
        >>> current.dg_dpsip(psip)
        0.0
        >>> current.di_dpsip(psip)
        0.0

        ```
        """

    def g(self, psip: float) -> float:
        """Always returns 1."""

    def i(self, psip: float) -> float:
        """Always returns 0."""

    def dg_dpsip(self, psip: float) -> float:
        """Always returns 0."""

    def di_dpsip(self, psip: float) -> float:
        """Always returns 0."""

class NcBfield:
    r"""Magnetic field from a NetCDF file.

    Provides methods for calculating $B(\psi_p, \theta_B)$, and its derivatives.

    Attributes
    ----------
    path
        The path to the NetCDF file.
    typ
        The 2D Interpolation type.
    psip_data
        The NetCDF $\psi_p$ data used to construct the $B(\psi_p, \theta_B)$ spline.
    theta_data
        The NetCDF $\theta_B$ data used to construct the $B(\psi_p, \theta_B)$ spline.
    b_data
        The NetCDF $B$ data used to construct the $B(\psi_p, \theta_B)$ spline.
    db_dpsip_data
        The $dB/d\psi_p$ values evaluated at the $(\psi_p, \theta_B)$ data points,
        through interpolation.
    db_dtheta_data
        The $dB/d\theta_B$ values evaluated at the $(\psi_p, \theta_B)$ data points,
        through interpolation.
    """

    path: str
    typ: Interp2DType
    shape: NDArrayShape
    psip_data: NDArray1D
    theta_data: NDArray1D
    b_data: NDArray2D
    db_dpsip_data: NDArray2D
    db_dtheta_data: NDArray2D

    def __init__(self, path: str, typ: Interp2DType) -> None:
        """Constructs a `Bfield`.

        Parameters
        ----------
        path
            The path to the NetCDF file.
        typ
            The 2D Interpolation type.

        Example
        -------

        ```python title="NcBfield creation"
        >>> bfield = dex.NcBfield(path, "Bicubic")
        >>>
        >>> psip = 0.015
        >>> theta = 3.1415
        >>>
        >>> b = bfield.b(psip, theta)
        >>> db_dpsip = bfield.db_dpsip(psip, theta)

        ```
        """

    def b(self, psip: float, theta: float) -> float:
        r"""The $B(\psi_p, \theta_B)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        """

    def db_dpsip(self, psip: float, theta: float) -> float:
        r"""The $dB/d\psi_p$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        """

    def db_dtheta(self, psip: float, theta: float) -> float:
        r"""The $dB/d\theta_B$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        """

class NcHarmonic:
    r"""A single Harmonic from a NetCDF file.

    A Harmonic has the form:

    $$
    h_{m,n} = \alpha_{m,n}(\psi_p) \cdot \cos\big( m\theta - n\zeta + \phi_{m,n}(\psi_p) \big)
    $$

    where $m, n$ are integers.

    Attributes
    ----------
    path
        The path to the NetCDF file.
    typ
        The 1D Interpolation type.
    m
        The $\theta$ frequency mode number.
    n
        The $\zeta$ frequency mode number.
    phase_method
        The method for calculating the harmonic's phase.
    phase_average
        The average_value of the phase data array, if `phase_method="average"`.
    phase_resonance
        The value of the phase at the resonance $m/n$, if `phase_method="resonance"`.
    psip_data
        The NetCDF $\psi_p$ data used to construct the $\alpha(\psi_p)$ and $\phi(\psi_p)$ splines.
    a_data
        The NetCDF $\alpha$ data used to construct the $\alpha(\psi_p)$ spline.
    phase_data
        The NetCDF $\phi$ data used to construct the $\phi(\psi_p)$ spline.
    """

    path: str
    typ: Interp1DType
    m: int
    n: int
    phase_average: float | None
    phase_resonance: float | None
    phase_method: str
    psip_data: NDArray1D
    a_data: NDArray1D
    phase_data: NDArray1D

    def __init__(
        self,
        path: str,
        typ: Interp1DType,
        m: int,
        n: int,
        phase_method: PhaseMethod = "Zero",
    ) -> None:
        r"""Constructs a `Harmonic`.

        Parameters
        ----------
        path
            The path to the NetCDF file.
        typ
            The 1D Interpolation type.
        m
            The `θ` frequency mode number.
        n
            The `ζ` frequency mode number.

        Other Parameters
        ----------------
        phase_method
            Availiable methods for the calculation of an NcHarmonic's phase $\phi(\psi_p)$.

        Example
        -------

        ```python title="NcHarmonic creation"
        >>> harmonic1 = dex.NcHarmonic(path, "Steffen", m=2, n=1)
        >>> harmonic2 = dex.NcHarmonic(path, "Steffen", m=3, n=2, phase_method="Resonance")
        >>>
        >>> psip = 0.015
        >>> theta = 3.1415
        >>> zeta = 6.28
        >>>
        >>> h = harmonic1.h(psip, theta, zeta)
        >>> alpha = harmonic2.a(psip)

        ```
        """

    def h(self, psip: float, theta: float, zeta: float) -> float:
        r"""The $h(\psi_p, \theta_B, \zeta)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        zeta
            The angle $\zeta$ in $[rads]$.
        """

    def dh_dpsip(self, psip: float, theta: float, zeta: float) -> float:
        r"""The $dh/d\psi_p$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        zeta
            The angle $\zeta$ in $[rads]$.
        """

    def dh_dtheta(self, psip: float, theta: float, zeta: float) -> float:
        r"""The $dh/d\theta_B$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        zeta
            The angle $\zeta$ in $[rads]$.
        """

    def dh_dzeta(self, psip: float, theta: float, zeta: float) -> float:
        r"""The $dh/d\zeta$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        zeta
            The angle $\zeta$ in $[rads]$.
        """

    def dh_dt(self, psip: float, theta: float, zeta: float) -> float:
        r"""The $dh/dt$ value.

        This method always returs `0` at the moment.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        zeta
            The angle $\zeta$ in $[rads]$.
        """

    def a(self, psip: float) -> float:
        r"""The amplitude $\alpha(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def da_dpsip(self, psip: float) -> float:
        r"""The harmonic's amplitude $d\alpha/d\psi_p$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def phase(self, psip: float) -> float:
        r"""The harmonic's phase $\phi(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """

    def __len__(self) -> int:
        """Returns the number of ψp data points."""

class NcPerturbation:
    r"""A sum of different netCDF perturbation harmonics.

    A Perturbation has the form:

    $$
    \alpha = \alpha(\psi_p, \theta, \zeta) = \sum h_{m,n} =
    \sum\alpha_{m, n}(\psi_p)\cdot \cos\big(n\zeta - m\theta + \phi_{m,n}(\psi_p)\big)
    $$

    To avoid confusion with a Harmonic's amplitude $\alpha$, we denote the total sum of
    the harmonics as $p$.

    Attributes
    ----------
    harmonics
        The list of harmonics that appear in the perturbation.
    """

    harmonics: list[NcHarmonic]

    def __init__(self, harmonics: list[NcHarmonic]) -> None:
        """Constructs a `Perturbation`.

        Parameters
        ----------
        harmonics
            The list of harmonics that appear in the perturbation.

        Example
        -------

        ```python title="NcPerturbation creation with specific NcHarmonics"
        >>> perturbation = dex.NcPerturbation(
        ...     [
        ...         dex.NcHarmonic(path, "steffen", m=2, n=1),
        ...         dex.NcHarmonic(path, "steffen", m=3, n=1),
        ...         dex.NcHarmonic(path, "steffen", m=3, n=2),
        ...     ]
        ... )

        ```

        ```python title="NcPerturbation creation by list comprehension"
        >>> perturbation = dex.NcPerturbation(
        ...     [dex.NcHarmonic(path, "steffen", m=2, n=n) for n in range(1, 3)]
        ... )

        ```

        """
        # TODO: example

    def p(self, psip: float, theta: float, zeta: float) -> float:
        r"""The $p(\psi_p, \theta_B, \zeta)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        zeta
            The angle $\zeta$ in $[rads]$.

        """

    def dp_dpsip(self, psip: float, theta: float, zeta: float) -> float:
        r"""The $dp/d\psi_p$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        zeta
            The angle $\zeta$ in $[rads]$.
        """

    def dp_dtheta(self, psip: float, theta: float, zeta: float) -> float:
        r"""The $dp/d\theta_B$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        zeta
            The angle $\zeta$ in $[rads]$.
        """

    def dp_dzeta(self, psip: float, theta: float, zeta: float) -> float:
        r"""The $dp/d\zeta$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        zeta
            The angle $\zeta$ in $[rads]$.
        """

    def dp_dt(self, psip: float, theta: float, zeta: float) -> float:
        r"""The $dp/dt$ value.

        This method always returs `0` at the moment.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        zeta
            The angle $\zeta$ in $[rads]$.
        """

    def __getitem__(self, n: int) -> NcHarmonic:
        """Returns the n-th harmonic"""

    def __len__(self) -> int:
        """Returns the number of ψp data points."""

# ==========================================================================
