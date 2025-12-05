"""This file mirrors all the definitions made in the dexter-python Rust API.

Example
-------
>>> import dexter as dx

"""

from dexter.types import (
    CalculatedFrequency,
    NDArray1D,
    NDArray2D,
    Interp1DType,
    Interp2DType,
    NDArrayShape,
    OrbitType,
    IntegrationStatus,
    PoincareSection,
    SingePeriodIntersections,
)

class Geometry:
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
    zlab_data
        The $J_{lab}(\psi_p, \theta_B)$ data array.

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
        Creating a `Geometry`:

        ```python
        >>> geometry = dx.Geometry("./data.nc", "cubic", "bicubic")
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
        r"""The $J_{lab}(\psi_p, \theta_B)$ value in $[m/T]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta_B$ angle in $[rads]$.
        """

class Qfactor:
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
        Creating a `Qfactor`:

        ```python
        >>> qfactor = dx.Qfactor("./data.nc", "steffen")
        >>>
        >>> # ψp->ψ interpolation
        >>> psip = 0.003
        >>> psi = qfactor.psi(psip)

        ```
        dψ/dψp check:

        ```python
        >>> from math import isclose
        >>>
        >>> psip = geometry.psip_wall/2
        >>> assert isclose(qfactor.q(psip), qfactor.dpsi_dpsip(psip))

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

class Currents:
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
        Creating a `Currents`:

        ```python
        >>> currents = dx.Currents("./data.nc", "steffen")
        >>>
        >>> # ψp->g interpolation
        >>> psip = 0.015
        >>> g = currents.g(psip)

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

class Bfield:
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
        Creating a `Bfield`:

        ```python
        >>> bfield = dx.Bfield("./data.nc", "bicubic")
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

    def d2b_dpsip2(self, psip: float, theta: float) -> float:
        r"""The $d^2B/d\psi_p^2$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        """

    def d2b_dtheta2(self, psip: float, theta: float) -> float:
        r"""The $d^2B/d\theta_B^2$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        """

    def d2b_dpsip_dtheta(self, psip: float, theta: float) -> float:
        r"""The $d^2B/d\psi_p d\theta_B$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The angle $\theta$ in $[rads]$.
        """

class Harmonic:
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
    phase_average
        The average_value of the phase data array.
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
    phase_average: float
    psip_data: NDArray1D
    a_data: NDArray1D
    phase_data: NDArray1D

    def __init__(self, path: str, typ: Interp1DType, m: int, n: int) -> None:
        """Constructs a `Harmonic`.

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

        Example
        -------
        Creating a `Harmonic`:

        ```python
        >>> harmonic = dx.Harmonic("./data.nc", "steffen", m=1, n=7)
        >>>
        >>> psip = 0.015
        >>> theta = 3.1415
        >>> zeta = 6.28
        >>>
        >>> h = harmonic.h(psip, theta, zeta)
        >>> alpha = harmonic.a(psip)

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

class Perturbation:
    r"""A sum of different perturbation harmonics.

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

    harmonics: list[Harmonic]

    def __init__(self, harmonics: list[Harmonic]) -> None:
        """Constructs a `Perturbation`.

        Parameters
        ----------
        harmonics
            The list of harmonics that appear in the perturbation.
        """

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


        Example
        -------
        Creating a `Perturbation` by picking Harmonics:

        ```python
        >>> perturbation = dx.Perturbation(
        ...     [
        ...         dx.Harmonic("./data.nc", "steffen", m=1, n=7),
        ...         dx.Harmonic("./data.nc", "steffen", m=1, n=8),
        ...         dx.Harmonic("./data.nc", "steffen", m=1, n=9),
        ...     ]
        ... )

        ```

        Creating a `Perturbation` by list comprehension:

        ```python
        >>> perturbation = dx.Perturbation(
        ...     [dx.Harmonic("./data.nc", "steffen", m=0, n=n) for n in range(1, 5)]
        ... )

        ```
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

    def __getitem__(self, n: int) -> Harmonic:
        """Returns the n-th harmonic"""

    def __len__(self) -> int:
        """Returns the number of ψp data points."""

# ==========================================================================

class InitialConditions:
    r"""A set of initial conditions.

    The initial conditions are defined on the
    $(t, \theta, \psi_p, \rho, \zeta, \mu)$ space.

    Example
    -------
    Creating an `InitialConditions` set:

    ```python
    >>> initial = dx.InitialConditions(
    ...     time0=0,
    ...     theta0=0,
    ...     psip0=0.3 * geometry.psip_wall,
    ...     rho0=1e-5,
    ...     zeta0=0,
    ...     mu=1e-6,
    ... )

    ```
    """

    time0: float
    theta0: float
    psip0: float
    rho0: float
    zeta0: float
    mu: float

    def __init__(
        self,
        time0: float,
        theta0: float,
        psip0: float,
        rho0: float,
        zeta0: float,
        mu: float,
    ) -> None:
        r"""
        Parameters
        ----------
        time0
            The initial time.
        theta0
            The initial $\theta$ angle.
        psip0
            The initial poloidal magnetic flux $\psi_p$.
        rho0
            The initial parallel gyro radius $\rho$.
        zeta0
            The initial $\zeta$ angle.
        mu
            The magnetic moment $\mu$.
        """

class MappingParameters:
    """Defines all the necessary parameters of a Mapping.

    Example
    -------
    Creating an `MappingParameters`:

    ```python
    >>> params = dx.MappingParameters(
    ...     section="ConstTheta",
    ...     alpha=0, # θ=0 interection
    ...     intersections=2000,
    ... )

    ```
    """

    section: PoincareSection
    alpha: float
    intersections: int

    def __init__(
        self,
        section: PoincareSection,
        alpha: float,
        intersections: int,
    ) -> None:
        r"""
        Parameters
        ----------
        section
            The surface of section $\Sigma$, defined by an equation $x_i = \alpha$,
            where $x_i = \theta$ or $\zeta$.
        alpha
            The constant that defines the surface of section (modulo $2π$).
        intersections
            The number of interections to calculate.
        """

class Particle:
    r"""A Particle.

    By taking $\mu = 0$ and $\rho \rightarrow 0$, the particle traces magnetic field
    lines.

    Attributes
    ----------
    initial_conditions
        The initial conditions set.
    evolution
        The evolution time series of the particle's integration.
    status
        The particle's integration status.
    orbit_type
        The particle's orbit type, calculated form its $\theta$-span.
    frequencies
        The particle's calculated frequencies.
    initial_energy
        The particle's initial energy, calculated from its initial state, in Normalized Units.
    final_energy
        The particle's final energy, calculated from its final state, in Normalized Units.
    """

    initial_conditions: InitialConditions
    evolution: Evolution
    status: IntegrationStatus
    orbit_type: OrbitType
    frequencies: Frequencies
    initial_energy: float
    final_energy: float

    def __init__(self, initial_conditions: InitialConditions) -> None:
        """
        Parameters
        ----------
        initial_conditions: InitialConditions
            The initial conditions set.

        Example
        -------
        Creating a `Particle` from an `InitialConditions` set:

        ```python
        >>> initial = dx.InitialConditions(
        ...     time0=0,
        ...     theta0=0,
        ...     psip0=0.3 * geometry.psip_wall,
        ...     rho0=1e-5,
        ...     zeta0=0,
        ...     mu=1e-6,
        ... )
        >>>
        >>> particle = dx.Particle(initial)

        ```
        """

    def integrate(
        self,
        qfactor: Qfactor,
        bfield: Bfield,
        currents: Currents,
        perturbation: Perturbation,
        t_eval: tuple[float, float],
    ) -> None:
        """Integrates the particle, storing its evolution.

        Parameters
        ----------
        qfactor
            The equilibrium's qfactor.
        currents
            The equilibrium's plasma current.
        bfield
            The equilibrium's magnetic field.
        perturbation
            The equilibrium's perturbation.
        t_eval: tuple[float, float]
            The time span $(t_0, t_f)$ in which to integrate the particle, in Normalized Units.

        Example
        -------

        ```python
        >>> particle.integrate(
        ...     qfactor=qfactor,
        ...     currents=currents,
        ...     bfield=bfield,
        ...     perturbation=perturbation,
        ...     t_eval=(0, 500),
        ... )

        ```
        """

    def map(
        self,
        qfactor: Qfactor,
        currents: Currents,
        bfield: Bfield,
        perturbation: Perturbation,
        params: MappingParameters,
    ) -> None:
        """Integrates the particle, storing its intersections with the Poincare
        surface defined by `params`.

        Parameters
        ----------
        qfactor
            The equilibrium's qfactor.
        currents
            The equilibrium's plasma current.
        bfield
            The equilibrium's magnetic field.
        perturbation
            The equilibrium's perturbation.
        params
            The parameters of the Poincare mapping.

        Example
        -------

        ```python
        >>> params = dx.MappingParameters(
        ...     section="ConstTheta",
        ...     alpha=0, # θ=0 interection
        ...     intersections=100,
        ... )
        >>>
        >>> particle.map(
        ...     qfactor=qfactor,
        ...     currents=currents,
        ...     bfield=bfield,
        ...     perturbation=perturbation,
        ...     params=params,
        ... )

        ```
        """

    def calculate_frequencies(
        self,
        qfactor: Qfactor,
        currents: Currents,
        bfield: Bfield,
        perturbation: Perturbation,
    ) -> None:
        r"""Calculates $\omega_\theta$, $\omega_\zeta$ and $q_{kinetic}$.

        !!! note

            $\omega_\theta$ is calculated by integrating for a single period with respect
            to the $\theta-\psi_p$ variables.

            The orbit frequency $\omega_\zeta$ corresponds to the bounce/transit averaged
            rate of toroidal precession $\Delta\zeta / T_\omega$.

            Finally, $q_{kinetic} = \omega_\zeta / \omega_\theta$.

        !!! tip

            Avoid using values like $\pi, 0$ for $\theta$. The solver looks for
            intersections with the initial $\theta-\psi_p$ values, but $\pi$ and $0$
            tend to be local minima/maxima, and therefore no intersections can be detected.
            Use 'random' intermediate values instead.

        Parameters
        ----------
        qfactor
            The equilibrium's qfactor.
        currents
            The equilibrium's plasma current.
        bfield
            The equilibrium's magnetic field.
        perturbation: Perturbation
            The equilibrium's perturbation.

        Example
        -------

        ```python
        >>> initial = dx.InitialConditions(
        ...     time0=0,
        ...     theta0=1,
        ...     psip0=0.8 * geometry.psip_wall,
        ...     rho0=1e-5,
        ...     zeta0=0,
        ...     mu=1e-6,
        ... )
        >>>
        >>> particle = dx.Particle(initial)
        >>>
        >>> particle.calculate_frequencies(
        ...     qfactor=qfactor,
        ...     currents=currents,
        ...     bfield=bfield,
        ...     perturbation=dx.Perturbation([]),
        ... )
        >>> print(particle.frequencies)  # doctest: +SKIP
        Frequencies {
            omega_theta: "0.0000357",
            omega_zeta: "0.0000033",
            qkinetic: "0.0911628",
        }

        ```
        """

class Evolution:
    r"""Stores the time series of the particle's orbit.

    Not meant to be constructed. It is stored as a particle's attribute.

    Attributes
    ----------
    time
        The time evolution
    theta
        The $\theta$ time series.
    psip
        The $\psi_p$ time series.
    rho
        The $\rho_{||}$ time series.
    zeta
        The $\zeta$ time series.
    psi
        The $\psi$ time series.
    ptheta
        The $P_\theta$ time series.
    pzeta
        The $P_\zeta$ time series.
    energy
        The energy time series.
    energy_std
        The relative standard deviation ($\sigma/\mu$) of the energy time series.
    steps_taken
        The total number of steps taken to complete the integration.
    steps_stored
        The number of points stored. This can different from the `steps_taken` attribute,
        for example when `mapping` a particle, in which case only the intersections are stored.
    final_time
        The final time stored in the `time` array. Returns None if the particle has not
        been integrated.

    Example
    -------

    ```python
    >>> particle.evolution.theta # doctest: +SKIP
    array([  3.14      ,   3.14293393,   3.14607403, ..., 504.04624189,
        504.05828496, 504.07118786], shape=(41270,))
    ```
    """

    time: NDArray1D
    theta: NDArray1D
    psip: NDArray1D
    rho: NDArray1D
    zeta: NDArray1D
    psi: NDArray1D
    ptheta: NDArray1D
    pzeta: NDArray1D
    energy: NDArray1D
    energy_std: float
    steps_taken: int
    steps_stored: int
    final_time: float | None

class Frequencies:
    r"""Stores the Particle's calculated $\omega_\theta$, $\omega_\zeta$
    and $q_{kinetic}$.

    The values are calculated from `Particle.calculate_frequencies` since they
    must be calculated in a specific order. Otherwise they are absent.

    Attributes
    ----------
    omega_theta
        The $\omega_\theta$ calculated frequency.
    omega_zeta
        The $\omega_\zeta$ calculated frequency.
    qkinetic
        The calculated $q_{kinetic}$.
    theta_intersections
        The $\theta$ intersections found during the integration.
    psip_intersections
        The $\psi_p$ intersections found during the integration.

    Example
    -------
    For a simple passing particle:

    ```python
    >>> particle.frequencies  # doctest: +SKIP
    Frequencies: Frequencies {
        omega_theta: "0.0000360582",
        omega_zeta : "0.0000032916",
        qkinetic   : "0.0912849868",
        ψp-intersections: "(2, [284, 718])",
        θ-intersections: "(1, [718])",
    },
    ```
    """

    omega_theta: CalculatedFrequency
    omega_zeta: CalculatedFrequency
    qkinetic: CalculatedFrequency
    theta_intersections: SingePeriodIntersections
    psip_intersections: SingePeriodIntersections

# ==========================================================================

class HeapInitialConditions:
    """Sets up the initial conditions of the Particles.

    Example
    -------

    ```python
    >>> import numpy as np
    >>>
    >>> # set initial conditions on the `r` space
    >>> rs = np.linspace(0, geometry.r_wall, 40)
    >>> psips = np.asarray([geometry.psip(r) for r in rs])
    >>>
    >>> initials = dx.HeapInitialConditions(
    ...     psips=psips,
    ...     zetas=np.zeros(len(psips)),
    ...     thetas=np.zeros(len(psips)),
    ...     rhos=1e-5 * np.ones(len(psips)),
    ...     mus=np.zeros(len(psips)),
    ... )
    >>> print(len(initials))
    40

    ```
    """

    thetas: NDArray1D
    psips: NDArray1D
    rhos: NDArray1D
    zetas: NDArray1D
    mus: NDArray1D

    def __init__(
        self,
        thetas: NDArray1D,
        psips: NDArray1D,
        rhos: NDArray1D,
        zetas: NDArray1D,
        mus: NDArray1D,
    ) -> None:
        r"""
        Parameters
        ----------
        thetas
            The initial $\theta$ angles.
        psips
            The initial poloidal magnetic fluxes $\psi_p$.
        rhos
            The initial parallel gyro radii $\rho$.
        zetas
            The initial $\zeta$ angles.
        mus
            The magnetic moments $mu$.
        """

    def __len__(self) -> int:
        """Returns the number of InitialConditions sets."""

class Heap:
    """A collection of multiple Particles constructed from sets of Initial Conditions.

    The particles are stored in arbitrary order.

    Example
    -------

    ```python
    >>> initials = dx.HeapInitialConditions(
    ...     psips=psips,
    ...     zetas=np.zeros(len(psips)),
    ...     thetas=np.zeros(len(psips)),
    ...     rhos=1e-5 * np.ones(len(psips)),
    ...     mus=np.zeros(len(psips)),
    ... )
    >>>
    >>> heap = dx.Heap(initials)

    ```
    """

    zetas: NDArray2D
    psips: NDArray2D
    thetas: NDArray2D
    psis: NDArray2D
    routine: str

    def __init__(self, initials: HeapInitialConditions) -> None:
        """Constructs a Heap.

        Parameters
        ----------
        initials
            The Particles' initial conditions sets.
        """

    def poincare(
        self,
        qfactor: Qfactor,
        currents: Currents,
        bfield: Bfield,
        perturbation: Perturbation,
        params: MappingParameters,
    ) -> None:
        """Calculates the Poincare map.

        Parameters
        ----------
        qfactor
            The equilibrium's qfactor.
        currents
            The equilibrium's plasma current.
        bfield
            The equilibrium's magnetic field.
        perturbation
            The equilibrium's perturbation.
        params
            The parameters of the Poincare mapping.

        Example
        -------

        ```python
        >>> heap = dx.Heap(initials)
        >>> params = dx.MappingParameters("ConstTheta", np.pi, 100)
        >>> heap.poincare(qfactor, currents, bfield, perturbation, params)

        ```
        """

    def calculate_frequencies(
        self,
        qfactor: Qfactor,
        currents: Currents,
        bfield: Bfield,
        perturbation: Perturbation,
    ) -> None:
        r"""Calculates the $\omega_\theta$, $\omega_\zeta$ and $q_{kinetic}$ of the particles.

        !!! tip

            This method appears to be a bit sensitive to the solver's arithmetic error, so
            it's recommended to use the energy-adaptive-step method, and maybe tighten the
            tolerance a bit.

        Parameters
        ----------
        qfactor
            The equilibrium's qfactor.
        currents
            The equilibrium's plasma current.
        bfield
            The equilibrium's magnetic field.
        perturbation
            The equilibrium's perturbation.

        Example
        -------

        ```python
        >>> num = 5
        >>> initials = dx.HeapInitialConditions(
        ...     thetas=np.ones(num),
        ...     psips=0.8 * geometry.psip_wall * np.ones(num),
        ...     rhos=np.linspace(1e-5, 5e-5, num),
        ...     zetas=np.zeros(num),
        ...     mus=np.zeros(num),
        ... )
        >>>
        >>> heap = dx.Heap(initials)
        >>>
        >>> heap.calculate_frequencies(
        ...     qfactor,
        ...     currents,
        ...     bfield,
        ...     perturbation=dx.Perturbation([])
        ... )

        ```
        """

    def __getitem__(self, n: int) -> Particle:
        """Returns the n-th particle."""

    def __len__(self) -> int:
        """Returns the number of particles."""
