"""Definition of the `Particle` object.

Note
----

Many of the Particle's attributes are stored as properties as they must be pulled from the `_r`
object to be up-to-date.
"""

from collections.abc import Sequence

from dexter.machine.machine import Machine
from dexter.simulate.initial import InitialConditions
from dexter.simulate.params import SolverParams, IntersectParams
from dexter.simulate import plot

from dexter.types import (
    Array1,
    EnergyPzetaPosition,
    IntegrationStatus,
    Intersection,
    OrbitType,
)

from dexter._core import _PyParticle, _PySolverParams, _PyIntersectParams
from dexter._utils import _ReprStrImpl, _RustTypeWrapper


class Particle(_ReprStrImpl, _RustTypeWrapper):
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
    >>> initial_conditions = dex.InitialConditions.Boozer(
    ...     t0=0,
    ...     flux0=dex.MagneticFlux.Toroidal(0.1),
    ...     theta0=3.14,
    ...     zeta0=0,
    ...     rho0=1e-4,
    ...     mu0=7e-6,
    ... )
    >>> particle = dex.Particle(initial_conditions)

    ```
    ```python title="Particle creation from a Mixed Coordinates set"
    >>> initial_conditions = dex.InitialConditions.Mixed(
    ...     t0=0,
    ...     flux0=dex.MagneticFlux.Poloidal(0.1),
    ...     theta0=3.14,
    ...     zeta0=0,
    ...     pzeta0=-0.025,
    ...     mu0=7e-6,
    ... )
    >>> particle = dex.Particle(initial_conditions)

    ```

    Attributes
    ----------
    initial_conditions
        The Particle's [`InitialConditions`][dexter.InitialConditions] set.
    integration_status
        The Particle's [`IntegrationStatus`][dexter.IntegrationStatus].
    steps_taken
        The total number of steps taken during the integration.
        This number is not necessarily the same as the number of steps stored.
    steps_stored
        The total number of steps stored in the time series arrays.
    duration
        The duration of the integration routine.
    initial_energy
        The Particle's initial energy in Normalized Units.
    final_energy
        The Particle's final energy in Normalized Units.
    energy_var
        The variance of the particles `energy_array`.
    energy_pzeta_position
        The Particle's [`EnergyPzetaPosition`][dexter.EnergyPzetaPosition].
    orbit_type
        The Particle's [`OrbitType`][dexter.OrbitType].
    omega_theta
        The Particle's $\omega_\theta$ frequency in Normalized Units.
    omega_zeta
        The Particle's $\omega_\zeta$ frequency in Normalized Units.
    qkinetic
        The Particle's $q_{kinetic}$.
    flux_cache_hits
        The magnetic flux' Accelerator cache hits.
    flux_cache_misses
        The magnetic flux' Accelerator cache misses.
    theta_cache_hits
        The $\theta$ coordinate's Accelerator cache misses.
    theta_cache_misses
        The $\theta$ coordinate's Accelerator cache misses.
    mode_cache_hits
        The modes' Accelerator cache hits.
    mode_cache_misses
        The modes' Accelerator cache misses.
    t_array
        The $t$ array.
    psi_array
        The $\psi$ array.
    psip_array
        The $\psi_p$ array.
    theta_array
        The $\theta$ array.
    zeta_array
        The $\zeta$ array.
    rho_array
        The $\rho_{||}$ array.
    mu_array
        The $\mu$ array.
    ptheta_array
        The $P_\theta$ array.
    pzeta_array
        The $P_\zeta$ array.
    energy_array
        The energy array.
    """

    _r: _PyParticle

    def __init__(self, initial_conditions: InitialConditions) -> None:
        self._r = _PyParticle(initial_conditions._r)
        pass

    def integrate(
        self,
        machine: Machine,
        teval: tuple[float, float],
        solver_params: SolverParams = SolverParams(),
    ):
        r"""Integrates the particle for a specific time interval.

        The time interval is in Normalized Units (inverse gyro-frequency).

        Parameters
        ----------
        machine
            The machine in which to integrate the particle.
        teval
            The time span $(t_0, t_f)$ in which to integrate the particle, in Normalized Units.
        solver_params
            The parameters passed to the solver.

        Example
        -------
        ```python title="Particle integration"
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
        >>> initial_conditions = dex.InitialConditions.Boozer(
        ...     t0=0,
        ...     flux0=dex.MagneticFlux.Toroidal(0.1),
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-4,
        ...     mu0=7e-6,
        ... )
        >>>
        >>> # Particle setup and integration
        >>> particle = dex.Particle(initial_conditions)
        >>>
        >>> solver_params = dex.SolverParams(
        ...     method=dex.SteppingMethod.EnergyAdaptiveStep(rel_tol=1e-7, abs_tol=1e-9),
        ...     max_steps=100_000,
        ... )
        >>> particle.integrate(machine, teval=(0, 1e2), solver_params=solver_params)

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
        r"""Integrates the particle, calculating its intersections with a constant $\theta$ or $\zeta$ surface.

        Using the method described by Hénon we can force the solver to step exactly on the intersection surface.

        The differences between two consecutive values of the corresponding angle variable are guaranteed to
        be $2\pi \pm \epsilon$ or $0 \pm \epsilon$, where $\epsilon$ a number smaller than the solver’s
        relative tolerance.

        Parameters
        ----------
        machine
            The machine in which to integrate the particle.
        intersect_params
            The intersection parameters.
        solver_params
            The parameters passed to the solver.

        Example
        -------
        ```python title="Particle intersection integration"
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
        >>> initial_conditions = dex.InitialConditions.Boozer(
        ...     t0=0,
        ...     flux0=dex.MagneticFlux.Toroidal(0.3),
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-3,
        ...     mu0=7e-5,
        ... )
        >>>
        >>> # Particle setup and intersection
        >>> particle = dex.Particle(initial_conditions)
        >>>
        >>> intersect_params = dex.IntersectParams("ConstZeta", angle=3.1415, turns=5)
        >>> particle.intersect(machine, intersect_params)

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
        r"""Integrates the particle for a certain amount of $\theta-\psi$ periods.

        Parameters
        ----------
        machine
            The machine in which to integrate the particle.
        periods
            The amount of periods to integrate.
        solver_params
            The parameters passed to the solver.

        Example
        -------
        ```python title="Particle orbit closing"
        >>> # Machine setup
        >>> LCFS = dex.MagneticFlux.Toroidal(0.45)
        >>> machine = dex.Machine(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.1, lcfs=LCFS),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ... )
        >>>
        >>> # Initial conditions setup
        >>> initial_conditions = dex.InitialConditions.Boozer(
        ...     t0=0,
        ...     flux0=dex.MagneticFlux.Toroidal(0.1),
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-5,
        ...     mu0=7e-6,
        ... )
        >>>
        >>> # Particle setup and integration
        >>> particle = dex.Particle(initial_conditions)
        >>> particle.close(machine)

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
        r"""Classifies the particle’s orbit using its position on the $(E, P_\zeta, \mu=const)$
        plane without integrating.

        This routine only works for LAR-like equilibria.

        Parameters
        ----------
        machine
            The machine in which to classify the particle.

        Example
        -------
        ```python title="Particle classification"
        >>> # Machine setup
        >>> LCFS = dex.MagneticFlux.Toroidal(0.45)
        >>> machine = dex.Machine(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.1, lcfs=LCFS),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ... )
        >>>
        >>> # Initial conditions setup
        >>> initial_conditions = dex.InitialConditions.Boozer(
        ...     t0=0,
        ...     flux0=dex.MagneticFlux.Toroidal(0.1),
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-4,
        ...     mu0=7e-6,
        ... )
        >>>
        >>> # Particle setup and classification
        >>> particle = dex.Particle(initial_conditions)
        >>> particle.classify(machine)

        ```
        """
        self._r.classify(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
        )

    @property
    def initial_conditions(self) -> InitialConditions:
        return InitialConditions._wrap(self._r.initial_conditions)

    @property
    def integration_status(self) -> IntegrationStatus:
        return self._r.integration_status

    @property
    def steps_taken(self) -> int:
        return self._r.steps_taken

    @property
    def steps_stored(self) -> int:
        return self._r.steps_stored

    @property
    def duration(self) -> str:
        return self._r.duration

    @property
    def initial_energy(self) -> float:
        if self._r.initial_energy is None:
            raise AttributeError("`initial_energy` has not been calculated")
        return self._r.initial_energy

    @property
    def final_energy(self) -> float:
        if self._r.final_energy is None:
            raise AttributeError("`final_energy` has not been calculated")
        return self._r.final_energy

    @property
    def energy_var(self) -> float:
        if self._r.energy_var is None:
            raise AttributeError("`energy_var` has not been calculated")
        return self._r.energy_var

    @property
    def energy_pzeta_position(self) -> EnergyPzetaPosition:
        if self._r.energy_pzeta_position is None:
            raise AttributeError("`energy_pzeta_position` has not been calculated")
        return self._r.energy_pzeta_position

    @property
    def orbit_type(self) -> OrbitType:
        if self._r.orbit_type is None:
            raise AttributeError("`orbit_type` has not been calculated")
        return self._r.orbit_type

    @property
    def omega_theta(self) -> float:
        if self._r.omega_theta is None:
            raise AttributeError("`omega_theta` has not been calculated")
        return self._r.omega_theta

    @property
    def omega_zeta(self) -> float:
        if self._r.omega_zeta is None:
            raise AttributeError("`omega_zeta` has not been calculated")
        return self._r.omega_zeta

    @property
    def qkinetic(self) -> float:
        if self._r.qkinetic is None:
            raise AttributeError("`qkinetic` has not been calculated")
        return self._r.qkinetic

    def print_caches(self):
        self._r.print_caches()

    def discard_arrays(self):
        self._r.discard_arrays()

    @property
    def flux_cache_hits(self) -> int:
        return self._r.flux_cache_hits

    @property
    def flux_cache_misses(self) -> int:
        return self._r.flux_cache_misses

    @property
    def theta_cache_hits(self) -> int:
        return self._r.theta_cache_hits

    @property
    def theta_cache_misses(self) -> int:
        return self._r.theta_cache_misses

    @property
    def mode_cache_hits(self) -> int:
        return self._r.mode_cache_hits

    @property
    def mode_cache_misses(self) -> int:
        return self._r.mode_cache_misses

    @property
    def t_array(self) -> Array1:
        return self._r.get_array("t_array")

    @property
    def psi_array(self) -> Array1:
        return self._r.get_array("psi_array")

    @property
    def psip_array(self) -> Array1:
        return self._r.get_array("psip_array")

    @property
    def theta_array(self) -> Array1:
        return self._r.get_array("theta_array")

    @property
    def zeta_array(self) -> Array1:
        return self._r.get_array("zeta_array")

    @property
    def rho_array(self) -> Array1:
        return self._r.get_array("rho_array")

    @property
    def mu_array(self) -> Array1:
        return self._r.get_array("mu_array")

    @property
    def ptheta_array(self) -> Array1:
        return self._r.get_array("ptheta_array")

    @property
    def pzeta_array(self) -> Array1:
        return self._r.get_array("pzeta_array")

    @property
    def energy_array(self) -> Array1:
        return self._r.get_array("energy_array")

    def plot_evolution(self, machine: Machine, **kwargs):
        """Wrapper around [`dexter.plot_evolution`][dexter.plot_evolution].


        Parameters
        ----------
        machine
            The machine in which the particle was integrated. It is used to convert specific
            quantities to SI units.
        **kwargs
            Extra arguments passed to [`dexter.plot_evolution`][dexter.plot_evolution].
        """
        plot.plot_evolution(machine, self, **kwargs)

    def plot_poloidal_drift(self, machine: Machine, **kwargs):
        """Wrapper around [`dexter.plot_poloidal_drift`][dexter.plot_poloidal_drift].


        Parameters
        ----------
        machine
            The machine in which the particle was integrated.
        **kwargs
            Extra arguments passed to [`dexter.plot_poloidal_drift`][dexter.plot_poloidal_drift].
        """
        plot.plot_poloidal_drift(machine, self, **kwargs)
