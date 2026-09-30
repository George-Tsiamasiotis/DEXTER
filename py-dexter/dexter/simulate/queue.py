"""Definition of the `Queue` object."""

from collections.abc import Collection, Iterable, Sequence

from dexter.types import Array1, EnergyPzetaPosition, Intersection, OrbitType, Routine

from dexter.machine.machine import Machine
from dexter.simulate.params import SolverParams, IntersectParams
from dexter.simulate.initial import QueueInitialConditions
from dexter.simulate.particle import Particle

from dexter._utils import _ReprStrImpl
from dexter._core import _PySolverParams, _PyQueue, _PyIntersectParams


class Queue(_ReprStrImpl):
    r"""A collection of multiple [`Particles`][dexter.Particle], constructed from a
    [`QueueInitialConditions`][dexter.QueueInitialConditions].

    Example
    -------

    ``` py title="Queue Instantiation"
    >>> num = 100
    >>> psi0s = np.linspace(0, 0.15, num)
    >>> initial_conditions = dex.QueueInitialConditions.Boozer(
    ...     t0=np.zeros(num),
    ...     flux0=dex.MagneticFluxArray.Toroidal(psi0s),
    ...     theta0=np.zeros(num),
    ...     zeta0=np.zeros(num),
    ...     rho0=np.full(num, 1e-5),
    ...     mu0=np.full(num, 1e-6),
    ... )
    >>>
    >>> # Queue setup and integration
    >>> queue = dex.Queue(initial_conditions)

    ```

    Attributes
    ----------
    initial_conditions
        The Queue's initial conditions
    routines
        The routines Queue has run
    energy_array
        A 1D array with the particle's energies.
    energy_var_array
        A 1D array with the particle's energy variances.
    omega_theta_array
        A 1D array with the particle's $\omega_\theta$.
    omega_zeta_array
        A 1D array with the particle's $\omega_\zeta$.
    qkinetic_array
        A 1D array with the particle's $q_{kin}$.

    Note
    ----

    In arrays containing particle quantities, the particles are visited in order of instantiation.
    """

    _r: _PyQueue

    def __init__(self, initial_conditions: QueueInitialConditions) -> None:
        self._r = _PyQueue(initial_conditions._r)

    @classmethod
    def FromParticles(cls, particles: Collection[Particle]) -> Queue:
        """Creates a `Queue` from a collection of `Particles`, copying all their attributes.

        The particles do not have to be initialized, but they must be defined on the same coordinate
        set (boozer/mixed).
        """
        _particles = [particle._r for particle in particles]

        obj = Queue.__new__(Queue)
        obj._r = _PyQueue.from_particles(_particles)
        return obj

    def integrate(
        self,
        machine: Machine,
        teval: tuple[float, float],
        solver_params: SolverParams = SolverParams(),
    ):
        r"""Integrates the particles for a specific time interval.

        The time interval is in Normalized Units (inverse gyro-frequency).

        Parameters
        ----------
        machine
            The machine in which to integrate the particles.
        teval
            The time span $(t_0, t_f)$ in which to integrate the particles, in Normalized Units.
        solver_params
            The parameters passed to the solver.

        Example
        -------
        ```python title="Queue integration"
        >>> # Machine setup
        >>> LCFS = dex.MagneticFlux.Toroidal(0.45)
        >>> machine = dex.Machine(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.1, lcfs=LCFS),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ...     perturbation=dex.Perturbation(
        ...         [
        ...             dex.FluteMode(epsilon=1e-3, lcfs=LCFS, m=1, n=3, phase=0),
        ...             dex.FluteMode(epsilon=1e-3, lcfs=LCFS, m=2, n=3, phase=0),
        ...         ]
        ...     )
        ... )
        >>>
        >>> # Initial conditions setup
        >>> num = 10
        >>> psi0s = dex.MagneticFluxArray.Toroidal(np.linspace(0, machine.psi_last.value, num))
        >>> initial_conditions = dex.QueueInitialConditions.Boozer(
        ...     t0=np.zeros(num),
        ...     flux0=psi0s,
        ...     theta0=np.zeros(num),
        ...     zeta0=np.zeros(num),
        ...     rho0=np.full(num, 1e-5),
        ...     mu0=np.full(num, 1e-6),
        ... )
        >>>
        >>> # Queue setup and integration
        >>> queue = dex.Queue(initial_conditions)
        >>>
        >>> solver_params = dex.SolverParams(
        ...     method=dex.SteppingMethod.EnergyAdaptiveStep(rel_tol=1e-7, abs_tol=1e-9),
        ...     max_steps=100_000,
        ... )
        >>> queue.integrate(machine, teval=(0, 1e2), solver_params=solver_params)

        ```
        """
        self._r.integrate(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            perturbation=machine.perturbation._r,
            teval=teval,
            solver_params=solver_params._r,
        )

    def intersect(
        self,
        machine: Machine,
        intersect_params: IntersectParams,
        solver_params: SolverParams = SolverParams(),
    ):
        r"""Integrates the particles, calculating their intersections with a constant $\theta$ or $\zeta$ surface.

        Using the method described by Hénon we can force the solver to step exactly on the intersection surface.

        The differences between two consecutive values of the corresponding angle variable are guaranteed to
        be $2\pi \pm \epsilon$, where $\epsilon$ a number smaller than the solver’s relative tolerance.

        Parameters
        ----------
        machine
            The machine in which to integrate the particles.
        intersect_params
            The intersection parameters.
        solver_params
            The parameters passed to the solver.

        Example
        -------
        ```python title="Queue intersection integration"
        >>> # Machine setup
        >>> LCFS = dex.MagneticFlux.Toroidal(0.45)
        >>> machine = dex.Machine(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.1, lcfs=LCFS),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ...     perturbation=dex.Perturbation(
        ...         [
        ...             dex.FluteMode(epsilon=1e-3, lcfs=LCFS, m=1, n=3, phase=0),
        ...             dex.FluteMode(epsilon=1e-3, lcfs=LCFS, m=2, n=3, phase=0),
        ...         ]
        ...     )
        ... )
        >>>
        >>> # Initial conditions and Intersection Parameters setup
        >>> num = 10
        >>> psi0s = dex.MagneticFluxArray.Toroidal(np.linspace(0, machine.psi_last.value, num))
        >>> initial_conditions = dex.QueueInitialConditions.Boozer(
        ...     t0=np.zeros(num),
        ...     flux0=psi0s,
        ...     theta0=np.zeros(num),
        ...     zeta0=np.zeros(num),
        ...     rho0=np.full(num, 1e-6),
        ...     mu0=np.full(num, 1e-4),
        ... )
        >>>
        >>> # Queue setup and integration
        >>> queue = dex.Queue(initial_conditions)
        >>>
        >>> intersect_params = dex.IntersectParams("ConstZeta", angle=3.1415, turns=5)
        >>> queue.intersect(machine, intersect_params)

        ```
        """
        self._r.intersect(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            perturbation=machine.perturbation._r,
            intersect_params=intersect_params._r,
            solver_params=solver_params._r,
        )

    def close(
        self,
        machine: Machine,
        periods: int = 1,
        solver_params: SolverParams = SolverParams(),
    ):
        r"""Integrates the particles for a certain amount of $\theta-\psi$ periods.

        Parameters
        ----------
        machine
            The machine in which to integrate the particles.
        periods
            The amount of periods to integrate.
        solver_params
            The parameters passed to the solver.

        Example
        -------
        ```python title="Queue orbit closing"
        >>> # Machine setup
        >>> LCFS = dex.MagneticFlux.Toroidal(0.45)
        >>> machine = dex.Machine(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.1, lcfs=LCFS),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ... )
        >>>
        >>> # Initial conditions and Intersection Parameters setup
        >>> num = 10
        >>> psi0s = dex.MagneticFluxArray.Toroidal(np.linspace(0, machine.psi_last.value, num))
        >>> initial_conditions = dex.QueueInitialConditions.Boozer(
        ...     t0=np.zeros(num),
        ...     flux0=psi0s,
        ...     theta0=np.zeros(num),
        ...     zeta0=np.zeros(num),
        ...     rho0=np.full(num, 1e-5),
        ...     mu0=np.full(num, 1e-6),
        ... )
        >>>
        >>> # Queue setup and integration
        >>> queue = dex.Queue(initial_conditions)
        >>> queue.close(machine)

        ```
        """
        self._r.close(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            perturbation=machine.perturbation._r,
            periods=periods,
            solver_params=solver_params._r,
        )

    def classify(
        self,
        machine: Machine,
    ):
        r"""Classifies the particles’ orbit using their position on the $(E, P_\zeta, \mu=const)$
        plane without integrating.

        This routine only works for LAR-like equilibria.

        Parameters
        ----------
        machine
            The machine in which to classify the particles.

        Example
        -------
        ```python title="Queue classification"
        >>> # Machine setup
        >>> LCFS = dex.MagneticFlux.Toroidal(0.45)
        >>> machine = dex.Machine(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.1, lcfs=LCFS),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ... )
        >>>
        >>> # Initial conditions and Intersection Parameters setup
        >>> num = 10
        >>> psi0s = dex.MagneticFluxArray.Toroidal(np.linspace(0, machine.psi_last.value, num))
        >>> initial_conditions = dex.QueueInitialConditions.Boozer(
        ...     t0=np.zeros(num),
        ...     flux0=psi0s,
        ...     theta0=np.zeros(num),
        ...     zeta0=np.zeros(num),
        ...     rho0=np.full(num, 1e-5),
        ...     mu0=np.full(num, 1e-6),
        ... )
        >>>
        >>> # Queue setup and integration
        >>> queue = dex.Queue(initial_conditions)
        >>> queue.classify(machine)

        ```
        """
        # Whether or not to use `classify_common_mu` is handled on the rust side.
        self._r.classify(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
        )

    def particles(self) -> list[Particle]:
        r"""Returns a list with all the contained objects.

        Note that calling this method **copies** all particles.
        """
        _particles = self._r.particles()
        return [Particle._wrap(_particle) for _particle in _particles]

    def retain_pzeta(self, start: float, end: float) -> None:
        r"""Iterates through `self`'s particles, keeping only the ones with an initial $P_\zeta$
        within the given span, discarding the rest.

        Parameters
        ----------
        start
            The $P_\zeta$ starting value.
        end
            The $P_\zeta$ ending value.
        """
        self._r.retain_pzeta(start, end)

    def retain_energy(self, start: float, end: float) -> None:
        r"""Iterates through `self`'s particles, keeping only the ones with an initial energy
        within the given span, discarding the rest.

        Parameters
        ----------
        start
            The energy starting value.
        end
            The energy ending value.
        """
        self._r.retain_energy(start, end)

    def retain_energy_pzeta_positions(self, positions: Iterable[EnergyPzetaPosition]):
        r"""Iterates through `self`’s particles, keeping only the ones with an
        [`EnergyPzetaPosition`][dexter.EnergyPzetaPosition] that matches an element of `positions`.

        Parameters
        ----------
        positions
            The $E-P_\zeta$ positions.
        """
        self._r.retain_energy_pzeta_positions(list(positions))

    def retain_orbit_types(self, orbit_types: Iterable[OrbitType]):
        r"""Iterates through `self`’s particles, keeping only the ones with an
        [`OrbitType`][dexter.OrbitType] that matches an element of `orbit_types`.

        Parameters
        ----------
        orbit_types
            The orbit types.
        """
        self._r.retain_orbit_types(list(orbit_types))

    @property
    def initial_conditions(self) -> QueueInitialConditions:
        return QueueInitialConditions._wrap(self._r.initial_conditions)

    @property
    def routines(self) -> list[Routine]:
        return self._r.routines

    @property
    def energy_array(self) -> Array1:
        return self._r.get_array("energy_array")

    @property
    def energy_rsd_array(self) -> Array1:
        return self._r.get_array("energy_rsd_array")

    @property
    def omega_theta_array(self) -> Array1:
        return self._r.get_array("omega_theta_array")

    @property
    def omega_zeta_array(self) -> Array1:
        return self._r.get_array("omega_zeta_array")

    @property
    def qkinetic_array(self) -> Array1:
        return self._r.get_array("qkinetic_array")

    def __getitem__(self, index: int) -> Particle:
        r"""Returns a **copy** of the `index`-th particle."""
        if -len(self) < index < 0:
            index = len(self) + index
        elif index <= -len(self) or index >= len(self):
            raise IndexError(
                f"Particle index '{index}' is out of range for Queue with length '{len(self)}'"
            )

        _particle = self._r[index]
        assert _particle is not None
        return Particle._wrap(_particle)

    def __len__(self) -> int:
        """Returns the particle count."""
        return self._r.__len__()
