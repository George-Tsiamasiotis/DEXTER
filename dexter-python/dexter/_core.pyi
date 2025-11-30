"""This file mirrors all the definitions made in the dexter-python Rust API."""

from dexter.types import (
    CalculatedFrequency,
    NDArray1D,
    NDArray2D,
    Interp1DType,
    Interp2DType,
    NDArrayShape,
    PoincareSection,
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
        The device's major radius $R$ in $[m]$.
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

    """

    path: str
    typ1d: Interp1DType
    typ2d: Interp2DType
    baxis: float
    raxis: float
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

    def __init__(self, path: str, typ1d: Interp1DType, typ2d: Interp2DType):
        """Constructs a `Geometry`.

        Parameters
        ----------
        path:
            The path to the NetCDF file.
        typ1d
            The type of 1D Interpolation.
        typ2d
            The type of 2D Interpolation.
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

class Qfactor:
    r"""q-factor from a NetCDF file.

    Provides methods for calculating $q(\psi_p)$, $\psi(\psi_p)$ and $d\psi/d\psi_p$.

    Attributes
    ----------
    path
        The path to the NetCDF file.
    typ
        The 1D Interpolation type.
    psip_wall
        The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units.
    psi_wall
        The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units.
    psip_data
        The NetCDF $\psi_p$ data used to construct the $q(\psi_p)$ and $\psi(\psi_p)$ splines.
    q_data
        The NetCDF $q$ data used to construct the $q(\psi_p)$ spline.
    psi_data
        The NetCDF $\psi$ data used to construct the $\psi(\psi_p)$ spline.
    """

    path: str
    typ: Interp1DType
    psip_wall: float
    psi_wall: float
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
    psip_wall
        The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units.
    psi_wall
        The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units.
    psip_data
        The NetCDF $\psi_p$ data used to construct the $q(\psi_p)$ and $\psi(\psi_p)$ splines.
    g_data
        The NetCDF $g$ data used to construct the $g(\psi_p)$ spline.
    i_data
        The NetCDF $I$ data used to construct the $I(\psi_p)$ spline.
    """

    path: str
    typ: Interp1DType
    psip_wall: float
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
    psip_wall
        The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units.
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
    psip_wall: float
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
    h_{m,n} = \alpha_{m,n}(\psi_p) \cdot \cos\big( n\zeta - m\theta + \phi_{m,n}(\psi_p) \big)
    $$

    where $m, n$ are integers.

    Attributes
    ----------
    path
        The path to the NetCDF file.
    typ
        The 1D Interpolation type.
    psip_wall
        The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units.
    m
        The `θ` frequency mode number.
    n
        The `ζ` frequency mode number.
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
    psip_wall: float
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
    """A set of initial conditions."""

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
        """Creates a set of initial conditions.

        Parameters
        ----------
        time0: float
            The initial time.
        theta0: float
            The initial `θ` angle.
        psip0: float
            The initial poloidal magnetic flux `ψp`.
        rho0: float
            The initial parallel gyro radius `ρ`.
        zeta0: float
            The initial `ζ` angle.
        mu: float
            The magnetic moment `μ`.
        """

class MappingParameters:
    """Defines all the necessary parameters of a Mapping.

    Attributes
    ----------
    section: Literal["ConstTheta", "ConstZeta"]
        The surface of section Σ, defined by an equation xᵢ= α, where xᵢ= θ or ζ.
    alpha: float
        The constant that defines the surface of section (modulo 2π).
    intersections: int
        The number of interections to calculate.
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
        """Defines all the necessary parameters of a Poincare Map.

        Parameters
        ----------
        section: Literal["ConstTheta", "ConstZeta"]
            The surface of section Σ, defined by an equation xᵢ= α, where xᵢ= θ or ζ.
        alpha: float
            The constant that defines the surface of section (modulo 2π).
        intersections: int
            The number of interections to calculate.
        """

class Particle:
    """A particle.

    Attributes
    ----------
    initial_conditions: InitialConditions
        The initial conditions set.
    evolution: Evolution
        The evolution time series of the particle's integration.
    status: str
        The particle's integration status.
    frequencies: Frequencies
        The particle's calculated frequencies.
    """

    initial_conditions: InitialConditions
    evolution: Evolution
    status: str
    frequencies: Frequencies

    def __init__(self, initial_conditions: InitialConditions) -> None:
        """Creates a Particle from an `InitialConditions` set.

        Parameters
        ----------
        initial_conditions: InitialConditions
            The initial conditions set.
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
        qfactor: Qfactor
            The equilibrium's qfactor.
        currents: Currents
            The equilibrium's plasma current.
        bfield: Bfield
            The equilibrium's magnetic field.
        per: Perturbation
            The equilibrium's perturbation.
        t_eval: tuple[float, float]
            The time span (t0, tf) in which to integrate the particle [Normalized].
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
        surface defined by `MappingParameters`.

        Parameters
        ----------
        qfactor: Qfactor
            The equilibrium's qfactor.
        currents: Currents
            The equilibrium's plasma current.
        bfield: Bfield
            The equilibrium's magnetic field.
        per: Perturbation
            The equilibrium's perturbation.
        params: MappingParameters
            The parameters of the Poincare mapping.
        """

    def calculate_frequencies(
        self,
        qfactor: Qfactor,
        currents: Currents,
        bfield: Bfield,
        perturbation: Perturbation,
    ) -> None:
        """Integrates the particle for 1 period, calculating ωθ, ωζ and qkinetic.

        Parameters
        ----------
        qfactor: Qfactor
            The equilibrium's qfactor.
        currents: Currents
            The equilibrium's plasma current.
        bfield: Bfield
            The equilibrium's magnetic field.
        per: Perturbation
            The equilibrium's perturbation.
        """

class Evolution:
    """Time series of the particle's orbit.

    Not meant to be constructed. It is stored as a particle's attribute.
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

class Frequencies:
    """Stores the Particle's calculated ωθ, ωζ and qkinetic."""

    omega_theta: CalculatedFrequency
    omega_zeta: CalculatedFrequency
    qkinetic: CalculatedFrequency

# ==========================================================================

class HeapInitialConditions:
    """Sets up the initial conditions for a poincare map."""

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
        """Constructs a HeapInitialConditions.

        Parameters
        ----------
        thetas: NDArray1D
            The initial `θ` angles.
        psips: NDArray1D
            The initial poloidal magnetic fluxes `ψp`.
        rhos: NDArray1D
            The initial parallel gyro radii `ρ`.
        zetas: NDArray1D
            The initial `ζ` angles.
        mus: NDArray1D
            The magnetic moments `μ`.
        """

    def __len__(self) -> int:
        """Returns the number of InitialConditions sets."""

class Heap:
    """A collection of multiple Particles constructed from sets of Initial Conditions."""

    zetas: NDArray2D
    psips: NDArray2D
    thetas: NDArray2D
    psis: NDArray2D
    routine: str

    def __init__(self, initial: HeapInitialConditions) -> None:
        """Constructs a Heap.

        Parameters
        ----------
        initials: HeapInitialConditions
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
        qfactor: Qfactor
            The equilibrium's qfactor.
        currents: Currents
            The equilibrium's plasma current.
        bfield: Bfield
            The equilibrium's magnetic field.
        per: Perturbation
            The equilibrium's perturbation.
        params: MappingParameters
            The parameters of the Poincare mapping.
        """

    def __getitem__(self, n: int) -> Harmonic:
        """Returns the n-th particle."""

    def __len__(self) -> int:
        """Returns the number of particles."""
