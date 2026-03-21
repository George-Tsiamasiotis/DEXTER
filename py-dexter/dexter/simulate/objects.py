"""Final wrappers of the exported rust types.

Note
----

Python methods that wrap rust 'generic' methods, assemble the monomorphized method's name on call,
using the generic object's `_dyn` attribute.
"""

from typing import Callable, Optional

from dexter._core import (
    _PyInitialFlux,
    _PyInitialConditions,
    _PyIntersectParams,
    _PyParticle,
    _PyInitialFluxArray,
    _PyQueueInitialConditions,
    _PyQueue,
)

from ..equilibrium import Equilibrium
from ._plotters import _ParticlePlotter, _QueuePlotter

from dexter.types import (
    Array1,
    FluxCoordinate,
    Intersection,
    SteppingMethod,
    IntegrationStatus,
    Routine,
)

# ================================================================================================
# ================================================================================================


class InitialFlux:
    """Defines the flux coordinate to be used in the initial conditions of a `Particle`.

    Parameters
    ----------
    kind
        The kind of initial flux.
    value
        The flux' value

    Example
    -------
    ```python title="InitialFlux definition"
    >>> flux0=("Toroidal", 0.1), # ψ0 = 0.1
    >>> flux0=("Poloidal", 0.5), # ψp0 = 0.5

    ```
    """

    _rust: _PyInitialFlux

    def __init__(self, kind: FluxCoordinate, value: float) -> None:
        self._rust = _PyInitialFlux(kind=kind, value=value)

    @property
    def kind(self) -> FluxCoordinate:
        r"""The kind of initial flux ($\psi$ or $\psi_p$)"""
        return self._rust.kind

    @property
    def value(self) -> float:
        r"""The contained value, regardless of kind."""
        return self._rust.value

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


class InitialConditions:
    r"""Initial Conditions set for a Particle.

    The set can be in either `Boozer` or `Mixed` coordinates.

    Use [`InitialConditions.boozer`][dexter.InitialConditions.boozer]
    and [`InitialConditions.mixed`][dexter.InitialConditions.mixed] class methods to
    construct the set.
    """

    _rust: _PyInitialConditions

    def __init__(self) -> None:
        raise RuntimeError(
            "Use 'InitialConditions.boozer' or InitialConditions.mixed' "
            + "to construct an InitialConditions set."
        )

    @classmethod
    def boozer(
        cls,
        t0: float,
        flux0: InitialFlux,
        theta0: float,
        zeta0: float,
        rho0: float,
        mu0: float,
    ) -> InitialConditions:
        r"""Initial conditions for a Particle in Boozer coordinates.

        The initial conditions are defined on the
        $(t, \psi, \theta, \zeta, \rho, \mu)$ or
        $(t, \psi_p, \theta, \zeta, \rho, \mu)$
        space, depending on the value of `flux0`.

        Parameters
        ----------
        t0
            The initial time, in Normalized Units.
        flux0
            The initial $\psi / \psi_p$, in Normalized Units.
        theta0
            The initial $\theta$ angle, in rads.
        zeta0
            The initial $\zeta$ angle, in rads.
        rho0
            The initial $\rho_{||}$, in Normalized Units.
        mu0
            The initial magnetic moment $\mu$, in Normalized Units.

        Example
        -------
        ```python title="InitialConditions definition in Boozer coordinates"
        >>> initial_conditions = dex.InitialConditions.boozer(
        ...     t0=0,
        ...     flux0=dex.InitialFlux("Toroidal", 0.1), # ψ0 = 0.1
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-4,
        ...     mu0=7e-6,
        ... )

        ```
        """
        initial = InitialConditions.__new__(InitialConditions)
        initial._rust = _PyInitialConditions.boozer(
            t0=t0, flux0=flux0._rust, theta0=theta0, zeta0=zeta0, rho0=rho0, mu0=mu0
        )
        return initial

    @classmethod
    def mixed(
        cls,
        t0: float,
        flux0: InitialFlux,
        theta0: float,
        zeta0: float,
        pzeta0: float,
        mu0: float,
    ) -> InitialConditions:
        r"""Initial conditions for a Particle in Mixed coordinates.

        The initial conditions are defined on the
        $(t, \psi, \theta, \zeta, P\zeta, \mu)$ or
        $(t, \zeta, \psi_p, \theta, \zeta, P\zeta, \mu)$
        space, depending on the value of `flux0`.

        Parameters
        ----------
        t0
            The initial time, in Normalized Units.
        pzeta0
            The initial $P_\zeta$, in Normalized Units.
        flux0
            The initial $\psi / \psi_p$, in Normalized Units.
        theta0
            The initial $\theta$ angle, in rads.
        zeta0
            The initial $\zeta$ angle, in rads.
        mu0
            The initial magnetic moment $\mu$, in Normalized Units.

        Example
        -------
        ```python title="nitialConditions definition in Mixed Coordinates"
        >>> initial_conditions = dex.InitialConditions.mixed(
        ...     t0=0,
        ...     flux0=dex.InitialFlux("Toroidal", 0.1), # ψ0 = 0.1
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     pzeta0=-0.025,
        ...     mu0=7e-6,
        ... )

        ```
        """
        initial = InitialConditions.__new__(InitialConditions)
        initial._rust = _PyInitialConditions.mixed(
            t0=t0, flux0=flux0._rust, theta0=theta0, zeta0=zeta0, pzeta0=pzeta0, mu0=mu0
        )
        return initial

    @property
    def t0(self) -> float:
        """The initial time, in Normalized Units."""
        return self._rust.t0

    @property
    def flux0(self) -> InitialFlux:
        r"""The initial $\psi / \psi_p$, in Normalized Units."""
        initial_flux = InitialFlux.__new__(InitialFlux)
        initial_flux._rust = self._rust.flux0
        return initial_flux

    @property
    def theta0(self) -> float:
        r"""The initial $\theta$ angle, in Normalized Units."""
        return self._rust.theta0

    @property
    def zeta0(self) -> float:
        r"""The initial $\zeta$ angle, in Normalized Units."""
        return self._rust.zeta0

    @property
    def mu0(self) -> float:
        r"""The initial magnetic moment $\mu$, in Normalized Units."""
        return self._rust.mu0

    @property
    def rho0(self) -> Optional[float]:
        r"""The initial parallel radius $\rho_{||}$, in Normalized Units."""
        return self._rust.rho0

    @property
    def pzeta0(self) -> Optional[float]:
        r"""The initial canonical momentum $P_\zeta$, in Normalized Units."""
        return self._rust.pzeta0

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


class IntersectParams:
    r"""Parameters for the Particle's `intersect()` method.

    Parameters
    ----------
    intersection
        The surface of section Σ, defined by an equation $\chi_i = \alpha$, where $\chi_i = \theta$ or
        $\zeta$.
    angle
        The constant that defines the surface of section.
    turns
        The number of intersections to calculate.

    Example
    -------
    ```python title="IntersectParams definition"
    >>> intersect_params = dex.IntersectParams(
    ...     intersection = "ConstTheta",
    ...     angle = 3.1415,
    ...     turns = 100,
    ... )

    ```
    """

    _rust: _PyIntersectParams

    def __init__(
        self,
        intersection: Intersection,
        angle: float,
        turns: int,
    ):
        self._rust = _PyIntersectParams(
            intersection=intersection, angle=angle, turns=turns
        )

    @property
    def intersection(self) -> Intersection:
        """The intersection surface."""
        return self._rust.intersection

    @property
    def angle(self) -> float:
        """The intersection surface's angle."""
        return self._rust.angle

    @property
    def turns(self) -> int:
        """The number of intersections to calculate."""
        return self._rust.turns

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


class Particle(_ParticlePlotter):
    r"""A Particle.

    By taking $\mu = 0$ and $\rho \rightarrow 0$, the particle traces magnetic field
    lines.

    Parameters
    ----------
    initial_conditions
        The initial conditions set.

    Example
    -------
    ```python title="Particle creation from a Boozer Coordinates set"
    >>> initial_conditions = dex.InitialConditions.boozer(
    ...     t0=0,
    ...     flux0=dex.InitialFlux("Toroidal", 0.1),
    ...     theta0=3.14,
    ...     zeta0=0,
    ...     rho0=1e-4,
    ...     mu0=7e-6,
    ... )
    >>> particle = dex.Particle(initial_conditions)

    ```
    ```python title="Particle creation from a Mixed Coordinates set"
    >>> initial_conditions = dex.InitialConditions.mixed(
    ...     t0=0,
    ...     flux0=dex.InitialFlux("Toroidal", 0.1),
    ...     theta0=3.14,
    ...     zeta0=0,
    ...     pzeta0=-0.025,
    ...     mu0=7e-6,
    ... )
    >>> particle = dex.Particle(initial_conditions)

    ```
    """

    _initial_conditions: InitialConditions
    _rust: _PyParticle

    def __init__(self, initial_conditions: InitialConditions) -> None:
        self._initial_conditions = initial_conditions
        self._rust = _PyParticle(initial_conditions=initial_conditions._rust)

    @classmethod
    def _from_rust_pyparticle(cls, _rust: _PyParticle) -> Particle:
        """Creates a new empty Particle and manually sets its `_rust: _PyParticle` field to
        the passed `_rust: _PyParticle`, essentially creating a copy.

        This is necessary to make `Queue` iterable.
        """
        particle = Particle.__new__(Particle)
        particle._rust = _rust
        return particle

    @property
    def initial_conditions(self) -> InitialConditions:
        """The initial conditions set."""
        return self._initial_conditions

    @property
    def integration_status(self) -> IntegrationStatus:
        """The particle's integration status."""
        return self._rust.integration_status

    @property
    def steps_taken(self) -> int:
        """The total number of steps taken during the integration.

        This number is not necessarily the same as the number of steps stored."""
        return self._rust.steps_taken

    @property
    def steps_stored(self) -> int:
        """The number of steps stored in the the time series arrays."""
        return self._rust.steps_stored

    @property
    def initial_energy(self) -> float | None:
        """The particle's initial energy before the integration in Normalized Units."""
        return self._rust.initial_energy

    @property
    def final_energy(self) -> float | None:
        """The particle's final energy after the integration in Normalized Units."""
        return self._rust.final_energy

    @property
    def energy_var(self) -> float | None:
        """The variance of energy throughout the integration."""
        return self._rust.energy_var

    @property
    def t_array(self) -> Array1:
        """The times of the integration steps."""
        return self._rust.t_array

    @property
    def psi_array(self) -> Array1:
        r"""The $\psi$ time series."""
        return self._rust.psi_array

    @property
    def psip_array(self) -> Array1:
        r"""The $\psi_p$ time series."""
        return self._rust.psip_array

    @property
    def theta_array(self) -> Array1:
        r"""The $\theta$ time series."""
        return self._rust.theta_array

    @property
    def zeta_array(self) -> Array1:
        r"""The $\zeta$ time series."""
        return self._rust.zeta_array

    @property
    def rho_array(self) -> Array1:
        r"""The $\rho$ time series."""
        return self._rust.rho_array

    @property
    def mu_array(self) -> Array1:
        r"""The $\mu$ time series."""
        return self._rust.mu_array

    @property
    def ptheta_array(self) -> Array1:
        r"""The $P_\theta$ time series."""
        return self._rust.ptheta_array

    @property
    def pzeta_array(self) -> Array1:
        r"""The $P_\zeta$ time series."""
        return self._rust.pzeta_array

    @property
    def energy_array(self) -> Array1:
        r"""The energy time series."""
        return self._rust.energy_array

    def integrate(
        self,
        /,
        equilibrium: Equilibrium,
        teval: tuple[float, float],
        *,
        stepping_method: Optional[SteppingMethod] = "EnergyAdaptiveStep",
        max_steps: Optional[int] = 1_000_000,
        first_step: Optional[float] = 1e-1,
        safety_factor: Optional[float] = 0.9,
        energy_rel_tol: Optional[float] = 1e-12,
        energy_abs_tol: Optional[float] = 1e-14,
        error_rel_tol: Optional[float] = 1e-12,
        error_abs_tol: Optional[float] = 1e-14,
    ):
        r"""Integrates the particle for a certain time interval.

        The time interval is in Normalized Units (inverse gyro-frequency).

        Parameters
        ----------
        equilibrium
            The equilibrium in which to integrate the particle.
        teval
            The time span $(t_0, t_f)$ in which to integrate the particle, in Normalized Units.

        Other Parameters
        ----------------
        stepping_method
            The optimal step calculation method. Defaults to "EnergyAdaptiveStep".
        max_steps
            The maximum amount of steps a particle can make before terminating its integration. Defaults to
            1.000.000.
        first_step
            The initial time step for the RKF45 adaptive step method. The value is empirical. Defaults to
            1e-1.
        safety_factor
            The safety factor of the solver. Should be less than 1.0. Defaults to 0.9.
        energy_rel_tol
            The relative tolerance of the energy difference in every step. Defaults to 1e-12.
        energy_abs_tol
            The absolute tolerance of the energy difference in every step. Defaults to 1e-14.
        error_rel_tol
            The relative tolerance of the local truncation error in every step. Defaults to 1e-12.
        error_abs_tol
            The absolute tolerance of the local truncation error in every step. Defaults to 1e-14.

        Example
        -------
        ```python title="Particle integration"
        >>> equilibrium = dex.Equilibrium(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qwall=4.1, flux_wall=("Toroidal", 0.45)),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ...     perturbation=dex.Perturbation(
        ...         [
        ...             dex.CosHarmonic(alpha=1e-3, m=1, n=3, phase=0),
        ...             dex.CosHarmonic(alpha=1e-3, m=2, n=3, phase=0),
        ...         ]
        ...     )
        ... )
        >>> initial_conditions = dex.InitialConditions.boozer(
        ...     t0=0,
        ...     flux0=dex.InitialFlux("Toroidal", 0.1),
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-4,
        ...     mu0=7e-6,
        ... )
        >>> particle = dex.Particle(initial_conditions)
        >>> particle.integrate(
        ...     equilibrium=equilibrium,
        ...     teval=(0, 1e2),
        ...     first_step=1e-3,
        ...     energy_rel_tol=1e-11,
        ... )

        ```
        """
        prefix = "__integrate"
        q = equilibrium.qfactor._dyn
        c = equilibrium.current._dyn
        b = equilibrium.bfield._dyn
        p = equilibrium.perturbation._dyn
        method_name: Callable = getattr(self._rust, f"{prefix}_{q}_{c}_{b}_{p}")
        method_name(
            equilibrium.qfactor._rust,
            equilibrium.current._rust,
            equilibrium.bfield._rust,
            equilibrium.perturbation._rust,
            teval,
            stepping_method,
            max_steps,
            first_step,
            safety_factor,
            energy_rel_tol,
            energy_abs_tol,
            error_rel_tol,
            error_abs_tol,
        )

    def intersect(
        self,
        /,
        equilibrium: Equilibrium,
        intersect_params: IntersectParams,
        *,
        stepping_method: Optional[SteppingMethod] = "EnergyAdaptiveStep",
        max_steps: Optional[int] = 1_000_000,
        first_step: Optional[float] = 1e-1,
        safety_factor: Optional[float] = 0.9,
        energy_rel_tol: Optional[float] = 1e-12,
        energy_abs_tol: Optional[float] = 1e-14,
        error_rel_tol: Optional[float] = 1e-12,
        error_abs_tol: Optional[float] = 1e-14,
    ):
        r"""Integrates the particle, calculating its intersections with a constant $\theta$ or $\zeta$ surface.

        Using the method described by Hénon we can force the solver to step exactly on the intersection surface.

        The differences between two consecutive values of the corresponding angle variable are guaranteed to
        be $2\pi \pm \epsilon$, where $\epsilon$ a number smaller than the solver’s relative tolerance.

        Parameters
        ----------
        equilibrium
            The equilibrium in which to integrate the particle.
        intersect_params
            The parameters of the integration.

        Other Parameters
        ----------------
        stepping_method
            The optimal step calculation method. Defaults to "EnergyAdaptiveStep".
        max_steps
            The maximum amount of steps a particle can make before terminating its integration. Defaults to
            1.000.000.
        first_step
            The initial time step for the RKF45 adaptive step method. The value is empirical. Defaults to
            1e-1.
        safety_factor
            The safety factor of the solver. Should be less than 1.0. Defaults to 0.9.
        energy_rel_tol
            The relative tolerance of the energy difference in every step. Defaults to 1e-12.
        energy_abs_tol
            The absolute tolerance of the energy difference in every step. Defaults to 1e-14.
        error_rel_tol
            The relative tolerance of the local truncation error in every step. Defaults to 1e-12.
        error_abs_tol
            The absolute tolerance of the local truncation error in every step. Defaults to 1e-14.

        Example
        -------
        ```python title="Particle intersection integration"
        >>> equilibrium = dex.Equilibrium(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qwall=4.1, flux_wall=("Toroidal", 0.45)),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ...     perturbation=dex.Perturbation(
        ...         [
        ...             dex.CosHarmonic(alpha=1e-3, m=1, n=3, phase=0),
        ...             dex.CosHarmonic(alpha=1e-3, m=2, n=3, phase=0),
        ...         ]
        ...     )
        ... )
        >>> initial_conditions = dex.InitialConditions.boozer(
        ...     t0=0,
        ...     flux0=dex.InitialFlux("Toroidal", 0.1),
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-4,
        ...     mu0=7e-6,
        ... )
        >>> intersect_params = dex.IntersectParams(
        ...     intersection = "ConstTheta",
        ...     angle = 3.1415,
        ...     turns = 5,
        ... )
        >>>
        >>> particle = dex.Particle(initial_conditions)
        >>> particle.intersect(
        ...     equilibrium=equilibrium,
        ...     intersect_params=intersect_params,
        ... )

        ```
        """
        prefix = "__intersect"
        q = equilibrium.qfactor._dyn
        c = equilibrium.current._dyn
        b = equilibrium.bfield._dyn
        p = equilibrium.perturbation._dyn
        method_name: Callable = getattr(self._rust, f"{prefix}_{q}_{c}_{b}_{p}")
        method_name(
            equilibrium.qfactor._rust,
            equilibrium.current._rust,
            equilibrium.bfield._rust,
            equilibrium.perturbation._rust,
            intersect_params._rust,
            stepping_method,
            max_steps,
            first_step,
            safety_factor,
            energy_rel_tol,
            energy_abs_tol,
            error_rel_tol,
            error_abs_tol,
        )

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


class InitialFluxArray:
    """A 1D array of initial "Toroidal" or "Poloidal" fluxes.

    Useful when creating a [`QueueInitialConditions`][dexter.QueueInitialConditions].

    Example
    -------
    ```python title="InitialFluxArray definition"
    >>> psi0s = InitialFluxArray("Toroidal", np.linspace(0, 0.5, 10))
    >>> psip0s = InitialFluxArray("Poloidal", np.linspace(0, 0.8, 20))

    ```
    """

    _rust: _PyInitialFluxArray

    def __init__(self, kind: FluxCoordinate, values: Array1):
        self._rust = _PyInitialFluxArray(kind=kind, values=values)

    @property
    def kind(self) -> FluxCoordinate:
        """The kind of the contained fluxes ("Toroidal" or "Poloidal")."""
        return self._rust.kind

    @property
    def values(self) -> Array1:
        """The flux values."""
        return self._rust.values

    def __str__(self) -> str:
        return f"kind: {self.kind}\nvalues: {self.values}"

    def __repr__(self) -> str:
        return self.__str__()


class QueueInitialConditions:
    r"""Sets of initial conditions for initializing a Queue.

    Use the [`InitialFluxArray`][dexter.InitialFluxArray] helper object to initialize the $\psi$/$\psi_p$ variables.

    Example
    -------
    ```python title="QueueInitialConditions definition"
    >>> num = 10
    >>> psi0s = InitialFluxArray("Toroidal", np.linspace(0, 0.5, num))
    >>>
    >>> initial_conditions = QueueInitialConditions(
    ...     t0=np.zeros(num),
    ...     flux0=psi0s,
    ...     theta0=np.zeros(num),
    ...     zeta0=np.zeros(num),
    ...     rho0=np.logspace(1e-6, 1e-5, num),
    ...     mu0=np.full(num, 1e-6),
    ... )

    ```
    """

    _rust: _PyQueueInitialConditions

    def __init__(
        self,
        t0: Array1,
        flux0: InitialFluxArray,
        theta0: Array1,
        zeta0: Array1,
        rho0: Array1,
        mu0: Array1,
    ) -> None:
        self._rust = _PyQueueInitialConditions(
            t0=t0,
            flux0=flux0._rust,
            theta0=theta0,
            zeta0=zeta0,
            rho0=rho0,
            mu0=mu0,
        )

    @property
    def t_array(self) -> Array1:
        """The initial times array."""
        return self._rust.t_array

    @property
    def theta_array(self) -> Array1:
        r"""The initial $\theta$ array."""
        return self._rust.theta_array

    @property
    def zeta_array(self) -> Array1:
        r"""The initial $\zeta$ array."""
        return self._rust.zeta_array

    @property
    def rho_array(self) -> Array1:
        r"""The initial $\rho$ array."""
        return self._rust.rho_array

    @property
    def mu_array(self) -> Array1:
        r"""The initial $\mu$ array."""
        return self._rust.mu_array

    def __len__(self) -> int:
        return self._rust.__len__()

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


class Queue(_QueuePlotter):
    """
    A collection of multiple [`Particles`][dexter.Particle] , constructed from a
    [`QueueInitialConditions`][dexter.QueueInitialConditions].

    Offers the ability to batch [`Particle::integrate`][dexter.Particle.integrate] or
    [`Particle::intersect`][dexter.Particle.intersect] all the contained particles,
    using multiple threads.

    Parameters
    ----------
    initial_conditions
        The initial conditions sets.

    Example
    -------
    ```python title="Queue integration"
    >>> # Initial Conditions setup
    >>> num = 10
    >>> psi0s = InitialFluxArray("Toroidal", np.linspace(0, 0.5, num))
    >>>
    >>> initial_conditions = QueueInitialConditions(
    ...     t0=np.zeros(num),
    ...     flux0=psi0s,
    ...     theta0=np.zeros(num),
    ...     zeta0=np.zeros(num),
    ...     rho0=np.logspace(1e-6, 1e-5, num),
    ...     mu0=np.full(num, 1e-6),
    ... )
    >>>
    >>> # Queue setup
    >>> queue = dex.Queue(initial_conditions)

    ```
    """

    _initial_conditions: QueueInitialConditions
    _intersect_params: IntersectParams  # for plot use
    _rust: _PyQueue

    def __init__(self, initial_conditions: QueueInitialConditions) -> None:
        self._initial_conditions = initial_conditions
        self._rust = _PyQueue(initial_conditions=initial_conditions._rust)

    @property
    def initial_conditions(self) -> QueueInitialConditions:
        """The Queue's initial conditions."""
        return self._initial_conditions

    @property
    def particle_count(self) -> int:
        """The number of contained Particles."""
        return self._rust.particle_count

    @property
    def particles(self) -> list[Particle]:
        """The contained Particles."""
        return [self[index] for index in range(self.particle_count)]

    @property
    def routine(self) -> Routine:
        """The routine run by the Queue."""
        return self._rust.routine

    def integrate(
        self,
        /,
        equilibrium: Equilibrium,
        teval: tuple[float, float],
        *,
        stepping_method: Optional[SteppingMethod] = "EnergyAdaptiveStep",
        max_steps: Optional[int] = 1_000_000,
        first_step: Optional[float] = 1e-1,
        safety_factor: Optional[float] = 0.9,
        energy_rel_tol: Optional[float] = 1e-12,
        energy_abs_tol: Optional[float] = 1e-14,
        error_rel_tol: Optional[float] = 1e-12,
        error_abs_tol: Optional[float] = 1e-14,
    ):
        r"""Integrates all the contained particles for a specific time interval.

        The time interval is in Normalized Units (inverse gyro-frequency).

        Parameters
        ----------
        equilibrium
            The equilibrium in which to integrate the particle.
        teval
            The time span $(t_0, t_f)$ in which to integrate the particles, in Normalized Units.

        Other Parameters
        ----------------
        stepping_method
            The optimal step calculation method. Defaults to "EnergyAdaptiveStep".
        max_steps
            The maximum amount of steps a particle can make before terminating its integration. Defaults to
            1.000.000.
        first_step
            The initial time step for the RKF45 adaptive step method. The value is empirical. Defaults to
            1e-1.
        safety_factor
            The safety factor of the solver. Should be less than 1.0. Defaults to 0.9.
        energy_rel_tol
            The relative tolerance of the energy difference in every step. Defaults to 1e-12.
        energy_abs_tol
            The absolute tolerance of the energy difference in every step. Defaults to 1e-14.
        error_rel_tol
            The relative tolerance of the local truncation error in every step. Defaults to 1e-12.
        error_abs_tol
            The absolute tolerance of the local truncation error in every step. Defaults to 1e-14.

        Example
        -------
        ```python title="Queue integration"
        >>> # Equilibrium setup
        >>> equilibrium = dex.Equilibrium(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qwall=4.1, flux_wall=("Toroidal", 0.45)),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ...     perturbation = dex.Perturbation(
        ...         [
        ...             dex.CosHarmonic(alpha=1e-3, m=1, n=3, phase=0),
        ...             dex.CosHarmonic(alpha=1e-3, m=2, n=3, phase=0),
        ...         ]
        ...     )
        ... )
        >>>
        >>> # Initial Conditions setup
        >>> num = 10
        >>> psi0s = InitialFluxArray("Toroidal", np.linspace(0, 0.5, num))
        >>>
        >>> initial_conditions = QueueInitialConditions(
        ...     t0=np.zeros(num),
        ...     flux0=psi0s,
        ...     theta0=np.zeros(num),
        ...     zeta0=np.zeros(num),
        ...     rho0=np.logspace(1e-6, 1e-5, num),
        ...     mu0=np.full(num, 1e-6),
        ... )
        >>>
        >>> # Queue setup
        >>> queue = dex.Queue(initial_conditions)
        >>>
        >>> # Run
        >>> queue.integrate(
        ...     equilibrium=    equilibrium,
        ...     teval=(0, 1e2),
        ...     first_step=1e-3,
        ...     energy_rel_tol=1e-11,
        ... )

        ```
        """
        prefix = "__integrate"
        q = equilibrium.qfactor._dyn
        c = equilibrium.current._dyn
        b = equilibrium.bfield._dyn
        p = equilibrium.perturbation._dyn
        method_name: Callable = getattr(self._rust, f"{prefix}_{q}_{c}_{b}_{p}")
        method_name(
            equilibrium.qfactor._rust,
            equilibrium.current._rust,
            equilibrium.bfield._rust,
            equilibrium.perturbation._rust,
            teval,
            stepping_method,
            max_steps,
            first_step,
            safety_factor,
            energy_rel_tol,
            energy_abs_tol,
            error_rel_tol,
            error_abs_tol,
        )

    def intersect(
        self,
        /,
        equilibrium: Equilibrium,
        intersect_params: IntersectParams,
        *,
        stepping_method: Optional[SteppingMethod] = "EnergyAdaptiveStep",
        max_steps: Optional[int] = 1_000_000,
        first_step: Optional[float] = 1e-1,
        safety_factor: Optional[float] = 0.9,
        energy_rel_tol: Optional[float] = 1e-12,
        energy_abs_tol: Optional[float] = 1e-14,
        error_rel_tol: Optional[float] = 1e-12,
        error_abs_tol: Optional[float] = 1e-14,
    ):
        r"""Integrates all the contained particles, calculating their intersections with a constant $\theta$ or $\zeta$ surface.

        Otherwise known as a Poincare map.

        The intersection surface, angle, and number of turns are configured with the helper class
        [`IntersectParams`][dexter.IntersectParams].

        Parameters
        ----------
        equilibrium
            The equilibrium in which to integrate the particle.
        intersect_params
            The parameters of the integration.

        Other Parameters
        ----------------
        stepping_method
            The optimal step calculation method. Defaults to "EnergyAdaptiveStep".
        max_steps
            The maximum amount of steps a particle can make before terminating its integration. Defaults to
            1.000.000.
        first_step
            The initial time step for the RKF45 adaptive step method. The value is empirical. Defaults to
            1e-1.
        safety_factor
            The safety factor of the solver. Should be less than 1.0. Defaults to 0.9.
        energy_rel_tol
            The relative tolerance of the energy difference in every step. Defaults to 1e-12.
        energy_abs_tol
            The absolute tolerance of the energy difference in every step. Defaults to 1e-14.
        error_rel_tol
            The relative tolerance of the local truncation error in every step. Defaults to 1e-12.
        error_abs_tol
            The absolute tolerance of the local truncation error in every step. Defaults to 1e-14.

        Example
        -------
        ```python title="Queue intersection integration"
        >>> # Equilibrium setup
        >>> equilibrium = dex.Equilibrium(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qwall=4.1, flux_wall=("Toroidal", 0.45)),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ...     perturbation = dex.Perturbation(
        ...         [
        ...             dex.CosHarmonic(alpha=1e-3, m=1, n=3, phase=0),
        ...             dex.CosHarmonic(alpha=1e-3, m=2, n=3, phase=0),
        ...         ]
        ...     )
        ... )
        >>>
        >>> # Initial Conditions setup
        >>> num = 10
        >>> psi0s = InitialFluxArray("Toroidal", np.linspace(0.01, 0.4, num))
        >>>
        >>> initial_conditions = QueueInitialConditions(
        ...     t0=np.zeros(num),
        ...     flux0=psi0s,
        ...     theta0=np.zeros(num),
        ...     zeta0=np.zeros(num),
        ...     rho0=np.linspace(1e-6, 2e-6, num),
        ...     mu0=np.full(num, 7e-7),
        ... )
        >>>
        >>> # IntersectParams setup
        >>> intersect_params = dex.IntersectParams(
        ...     intersection = "ConstTheta",
        ...     angle = 0.0,
        ...     turns = 5,
        ... )
        >>>
        >>> # Queue setup
        >>> queue = dex.Queue(initial_conditions)
        >>>
        >>> # Run
        >>> queue.intersect(
        ...     equilibrium=equilibrium,
        ...     intersect_params=intersect_params,
        ...     energy_rel_tol=1e-11,
        ... )

        ```
        """
        self._intersect_params = intersect_params
        prefix = "__intersect"
        q = equilibrium.qfactor._dyn
        c = equilibrium.current._dyn
        b = equilibrium.bfield._dyn
        p = equilibrium.perturbation._dyn
        method_name: Callable = getattr(self._rust, f"{prefix}_{q}_{c}_{b}_{p}")
        method_name(
            equilibrium.qfactor._rust,
            equilibrium.current._rust,
            equilibrium.bfield._rust,
            equilibrium.perturbation._rust,
            intersect_params._rust,
            stepping_method,
            max_steps,
            first_step,
            safety_factor,
            energy_rel_tol,
            energy_abs_tol,
            error_rel_tol,
            error_abs_tol,
        )

    def __getitem__(self, index: int):
        return Particle._from_rust_pyparticle(self._rust[index])

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()
