"""Definition of the `Queue` object."""

from dexter.types import SteppingMethod, Intersection

from dexter.machine.machine import Machine
from dexter.simulate.initial import QueueInitialConditions

from dexter._utils import _ReprStrImpl
from dexter._core import _PySolverParams, _PyQueue, _PyIntersectParams


class Queue(_ReprStrImpl):
    r"""A collection of multiple [`Particles`][dexter.Particle], constructed from a
    [`QueueInitialConditions`][dexter.QueueInitialConditions].
    """

    _r: _PyQueue

    def __init__(self, initial_conditions: QueueInitialConditions) -> None:
        self._r = _PyQueue(initial_conditions._r)

    def integrate(
        self,
        machine: Machine,
        teval: tuple[float, float],
        *,
        stepping_method: SteppingMethod | None = "EnergyAdaptiveStep",
        max_steps: int | None = 1_000_000,
        first_step: float | None = 1e-1,
        safety_factor: float | None = 0.9,
        energy_rel_tol: float | None = 1e-12,
        energy_abs_tol: float | None = 1e-14,
        error_rel_tol: float | None = 1e-12,
        error_abs_tol: float | None = 1e-14,
    ):
        r"""Integrates the particles for a specific time interval.

        The time interval is in Normalized Units (inverse gyro-frequency).

        Parameters
        ----------
        machine
            The machine in which to integrate the particles.
        teval
            The time span $(t_0, t_f)$ in which to integrate the particles, in Normalized Units.

        Other Parameters
        ----------------
        stepping_method
            The optimal step calculation method.
        max_steps
            The maximum amount of steps each particle can make before terminating its integration.
        first_step
            The initial time step for the RKF45 adaptive step method.
        safety_factor
            The safety factor of the solver.
        energy_rel_tol
            The relative tolerance of the energy difference in every step.
        energy_abs_tol
            The absolute tolerance of the energy difference in every step.
        error_rel_tol
            The relative tolerance of the local truncation error in every step.
        error_abs_tol
            The absolute tolerance of the local truncation error in every step.

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
        >>> queue.integrate(
        ...     machine=machine,
        ...     teval=(0, 1e2),
        ...     energy_rel_tol=1e-11,
        ...     energy_abs_tol=1e-13,
        ... )

        ```
        """
        solver_params = _PySolverParams(
            stepping_method=stepping_method,
            max_steps=max_steps,
            first_step=first_step,
            safety_factor=safety_factor,
            energy_rel_tol=energy_rel_tol,
            energy_abs_tol=energy_abs_tol,
            error_rel_tol=error_rel_tol,
            error_abs_tol=error_abs_tol,
        )
        self._r.integrate(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            perturbation=machine.perturbation._r,
            teval=teval,
            solver_params=solver_params,
        )

    def intersect(
        self,
        machine: Machine,
        intersection: Intersection,
        angle: float,
        turns: int,
        *,
        stepping_method: SteppingMethod | None = "EnergyAdaptiveStep",
        max_steps: int | None = 1_000_000,
        first_step: float | None = 1e-1,
        safety_factor: float | None = 0.9,
        energy_rel_tol: float | None = 1e-12,
        energy_abs_tol: float | None = 1e-14,
        error_rel_tol: float | None = 1e-12,
        error_abs_tol: float | None = 1e-14,
    ):
        r"""Integrates the particles, calculating their intersections with a constant $\theta$ or $\zeta$ surface.

        Using the method described by Hénon we can force the solver to step exactly on the intersection surface.

        The differences between two consecutive values of the corresponding angle variable are guaranteed to
        be $2\pi \pm \epsilon$, where $\epsilon$ a number smaller than the solver’s relative tolerance.

        Parameters
        ----------
        machine
            The machine in which to integrate the particles.
        intersection
            The surface of section Σ, defined by an equation $\chi_i = \alpha$, where $\chi_i = \theta$ or
            $\zeta$.
        angle
            The constant that defines the surface of section.
        turns
            The number of intersections to calculate.

        Other Parameters
        ----------------
        stepping_method
            The optimal step calculation method.
        max_steps
            The maximum amount of steps each particle can make before terminating its integration.
        first_step
            The initial time step for the RKF45 adaptive step method.
        safety_factor
            The safety factor of the solver.
        energy_rel_tol
            The relative tolerance of the energy difference in every step.
        energy_abs_tol
            The absolute tolerance of the energy difference in every step.
        error_rel_tol
            The relative tolerance of the local truncation error in every step.
        error_abs_tol
            The absolute tolerance of the local truncation error in every step.

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
        >>> queue.intersect(
        ...     machine=machine,
        ...     intersection="ConstZeta",
        ...     angle=3.1415,
        ...     turns=5,
        ... )

        ```
        """
        solver_params = _PySolverParams(
            stepping_method=stepping_method,
            max_steps=max_steps,
            first_step=first_step,
            safety_factor=safety_factor,
            energy_rel_tol=energy_rel_tol,
            energy_abs_tol=energy_abs_tol,
            error_rel_tol=error_rel_tol,
            error_abs_tol=error_abs_tol,
        )
        intersect_params = _PyIntersectParams(
            intersection=intersection,
            angle=angle,
            turns=turns,
        )

        self._r.intersect(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            perturbation=machine.perturbation._r,
            intersect_params=intersect_params,
            solver_params=solver_params,
        )

    def close(
        self,
        machine: Machine,
        periods: int | None = 1,
        *,
        stepping_method: SteppingMethod | None = "EnergyAdaptiveStep",
        max_steps: int | None = 1_000_000,
        first_step: float | None = 1e-1,
        safety_factor: float | None = 0.9,
        energy_rel_tol: float | None = 1e-12,
        energy_abs_tol: float | None = 1e-14,
        error_rel_tol: float | None = 1e-12,
        error_abs_tol: float | None = 1e-14,
    ):
        r"""Integrates the particles for a certain amount of $\theta-\psi$ periods.

        Parameters
        ----------
        machine
            The machine in which to integrate the particles.
        periods
            The amount of periods to integrate.

        Other Parameters
        ----------------
        stepping_method
            The optimal step calculation method.
        max_steps
            The maximum amount of steps each particle can make before terminating its integration.
        first_step
            The initial time step for the RKF45 adaptive step method.
        safety_factor
            The safety factor of the solver.
        energy_rel_tol
            The relative tolerance of the energy difference in every step.
        energy_abs_tol
            The absolute tolerance of the energy difference in every step.
        error_rel_tol
            The relative tolerance of the local truncation error in every step.
        error_abs_tol
            The absolute tolerance of the local truncation error in every step.

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
        solver_params = _PySolverParams(
            stepping_method=stepping_method,
            max_steps=max_steps,
            first_step=first_step,
            safety_factor=safety_factor,
            energy_rel_tol=energy_rel_tol,
            energy_abs_tol=energy_abs_tol,
            error_rel_tol=error_rel_tol,
            error_abs_tol=error_abs_tol,
        )

        self._r.close(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            perturbation=machine.perturbation._r,
            periods=periods if periods is not None else 1,
            solver_params=solver_params,
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
