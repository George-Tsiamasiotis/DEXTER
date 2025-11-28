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
    """Stores fluxes, angles and lab variables' data, and provides interpolation methods between them."""

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

    def __init__(self, path: str, typ1d: Interp1DType, typ2d: Interp2DType) -> None:
        """Constructs a Geometry.

        Parameters
        ----------
        path: str
            The path to the NetCDF file.
        typ1d: Interp1DType
            The type of 1D Interpolation.
        typ2d: Interp2DType
            The type of 2D Interpolation.
        """

    def r(self, psip: float) -> float:
        """The r(ψp) value in [m]."""

    def psip(self, r: float) -> float:
        """The ψp(r) value in Normalized Units."""

    def rlab(self, psip: float, theta: float) -> float:
        """The R(ψp, θ) value in [m]."""

    def zlab(self, psip: float, theta: float) -> float:
        """The Z(ψp, θ) value in [m]."""

class Qfactor:
    """q-factor reconstructed from a NetCDF file.

    Attributes
    ----------
    path: str
        The path to the NetCDF file.
    typ: Interp1DType
        The type of Interpolation.
    psip_wall: float
        The value of the poloidal flux ψp at the wall.
    psi_wall: float
        The value of the toroidal flux ψ at the wall.
    psip_data: NDArray1D
        The NetCDF ψp data used to construct the q(ψp) and ψ(ψp) splines.
    q_data: NDArray1D
        The NetCDF q data used to construct the q(ψp) spline.
    psi_data: NDArray1D
        The NetCDF ψp data used to construct the ψ(ψp) spline.
    q_data_derived: NDArray1D
        The q values, as calculated from dψ/dψp, at the ψp data.
    """

    path: str
    typ: Interp1DType
    psip_wall: float
    psi_wall: float
    psip_data: NDArray1D
    q_data: NDArray1D
    psi_data: NDArray1D
    r_data: NDArray1D

    def __init__(self, path: str, typ: Interp1DType) -> None:
        """q-factor reconstructed from a NetCDF file.

        Parameters
        ----------
        path: str
            The path to the NetCDF file.
        typ: Interp1DType
            The type of Interpolation.
        """

    def q(self, psip: float) -> float:
        """The q value evaluated at ψp"""

    def psi(self, psip: float) -> float:
        """The ψ value evaluated at ψp"""

    def dpsi_dpsip(self, psip: float) -> float:
        """The dψ/dψp value evaluated at ψp.

        This value should always equal `q(ψp)`.
        """

    def __len__(self) -> int:
        """Returns the number of ψp data points."""

class Currents:
    """Plasma current reconstructed from a NetCDF file.

    Attributes
    ----------
    path: str
        The path to the NetCDF file.
    typ: Interp1DType
        The type of Interpolation.
    psip_wall: float
        The value of the poloidal flux ψp at the wall.
    psip_data: NDArray1D
        The NetCDF ψp data used to construct the g(ψp) and I(ψp) splines.
    g_data: NDArray1D
        The NetCDF g data used to construct the g(ψp) spline.
    i_data: NDArray1D
        The NetCDF I data used to construct the I(ψp) spline.
    """

    path: str
    typ: Interp1DType
    psip_wall: float
    psip_data: NDArray1D
    g_data: NDArray1D
    i_data: NDArray1D

    def __init__(self, path: str, typ: Interp1DType) -> None:
        """Plasma current reconstructed from a NetCDF file.

        Parameters
        ----------
        path: str
            The path to the NetCDF file.
        typ: Interp1DType
            The type of Interpolation.
        """

    def g(self, psip: float) -> float:
        """The g value evaluated at ψp"""

    def i(self, psip: float) -> float:
        """The I value evaluated at ψp"""

    def dg_dpsip(self, psip: float) -> float:
        """The dg/dψp value evaluated at ψp"""

    def di_dpsip(self, psip: float) -> float:
        """The dI/dψp value evaluated at ψp"""

    def __len__(self) -> int:
        """Returns the number of ψp data points."""

class Bfield:
    """Magnetic field reconstructed from a NetCDF file.

    Attributes
    ----------
    path: str
        The path to the NetCDF file.
    typ: Interp2DType
        The type of Interpolation.
    psip_wall: float
        The value of the poloidal flux ψp at the wall.
    psip_data: NDArray1D
        The NetCDF ψp data used to construct the b(ψp, θ) spline.
    theta_data: NDArray1D
        The NetCDF θ data used to construct the b(ψp, θ) spline.
    b_data: NDArray2D
        The NetCDF b data used to construct the b(ψp, θ) spline.
    db_dpsip_data: NDArray2D
        The db/dψp values evaluated at the (ψp, θ) data, through interpolation.
    db_dtheta_data: NDArray2D
        The db/dψp values evaluated at the (ψp, θ) data, through interpolation.
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
        """Magnetic field reconstructed from a NetCDF file.

        Parameters
        ----------
        path: str
            The path to the NetCDF file.
        typ: Interp2DType
            The type of Interpolation.
        """

    def b(self, psip: float, theta: float) -> float:
        """The b value evaluated at (ψp, θ)"""

    def db_dpsip(self, psip: float, theta: float) -> float:
        """The db/dψp value evaluated at (ψp, θ)"""

    def db_dtheta(self, psip: float, theta: float) -> float:
        """The db/dθ value evaluated at (ψp, θ)"""

    def d2b_dpsip2(self, psip: float, theta: float) -> float:
        """The d2b/dψp2 value evaluated at (ψp, θ)"""

    def d2b_dtheta2(self, psip: float, theta: float) -> float:
        """The d2b/dθ2 value evaluated at (ψp, θ)"""

    def d2b_dpsip_dtheta(self, psip: float, theta: float) -> float:
        """The d2b/dψpdθ value evaluated at (ψp, θ)"""

class Harmonic:
    """A single perturbation harmonic.

    Attributes
    ----------
    path: str
        The path to the NetCDF file.
    typ: Interp1DType
        The type of Interpolation.
    psip_wall: float
        The value of the poloidal flux ψp at the wall.
    m: int
        The `θ` frequency number.
    n: int
        The `ζ` frequency number.
    phase_average: float
        The average_value of the phase data array.
    psip_data: float
        The NetCDF ψp data used to construct the a(ψp) spline.
    a_data: float
        The NetCDF a data used to construct the α(ψp) spline.
    phase_data: NDArray1D
        The NetCDF a data used to construct the φ(ψp) spline.
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
        """Creates a single perturbation harmonic.

        Parameters
        ----------
        path: str
            The path to the NetCDF file.
        typ: Interp1DType
            The type of Interpolation.
        m: int
            The `θ` frequency number.
        n: int
            The `ζ` frequency number.
        """

    def h(self, psip: float, theta: float, zeta: float) -> float:
        """The h value (value of the whole harmonic) evaluated at (ψp, θ, ζ)."""

    def dh_dpsip(self, psip: float, theta: float, zeta: float) -> float:
        """The dh/dψp value evaluated at (ψp, θ, ζ)."""

    def dh_dtheta(self, psip: float, theta: float, zeta: float) -> float:
        """The dh/dθ value evaluated at (ψp, θ, ζ)."""

    def dh_dzeta(self, psip: float, theta: float, zeta: float) -> float:
        """The dh/dζ value evaluated at (ψp, θ, ζ)."""

    def dh_dt(self, psip: float, theta: float, zeta: float) -> float:
        """The dh/dt value evaluated at (ψp, θ, ζ)."""

    def a(self, psip: float) -> float:
        """The `α(ψp) value (value of the harmonic's amplitude only)."""

    def da_dpsip(self, psip: float) -> float:
        """The `dα(ψp)/da_dpsip value (value of the harmonic's amplitude's derivative only)."""

    def phase(self, psip: float) -> float:
        """The `φ(ψp) value (value of the harmonic's phase only)."""

    def __len__(self) -> int:
        """Returns the number of ψp data points."""

class Perturbation:
    """A sum of different perturbation harmonics."""

    harmonics: list[Harmonic]

    def __init__(self, harmonics: list[Harmonic]) -> None:
        """Creates a Perturbation.

        Parameters
        ----------
        harmonics: list[Harmonics]
            The list of harmonics that appear in the perturbation.
        """

    def p(self, psip: float, theta: float, zeta: float) -> float:
        """The p value (value of the whole harmonic) evaluated at (ψp, θ, ζ)."""

    def dp_dpsip(self, psip: float, theta: float, zeta: float) -> float:
        """The dp/dψp value evaluated at (ψp, θ, ζ)."""

    def dp_dtheta(self, psip: float, theta: float, zeta: float) -> float:
        """The dp/dθ value evaluated at (ψp, θ, ζ)."""

    def dp_dzeta(self, psip: float, theta: float, zeta: float) -> float:
        """The dp/dζ value evaluated at (ψp, θ, ζ)."""

    def dp_dt(self, psip: float, theta: float, zeta: float) -> float:
        """The dp/dt value evaluated at (ψp, θ, ζ)."""

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
