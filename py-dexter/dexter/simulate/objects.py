"""Final wrappers of the exported rust types.

Note
----

Python methods that wrap rust 'generic' methods, assemble the monomorphized method's name on call,
using the generic object's `_dyn` attribute.
"""

from typing import Callable, Optional, Any
from dexter._core import _PyInitialConditions, _PyIntersectParams, _PyParticle

from ..equilibrium import Qfactor, Current, Bfield, Perturbation
from ._plotters import _ParticlePlotter

from dexter.types import (
    Array1,
    InitialFlux,
    Intersection,
    SteppingMethod,
    IntegrationStatus,
)

# ================================================================================================
# ================================================================================================


class InitialConditions:
    r"""Initial conditions for a Particle.

    The initial conditions are defined on the
    $(t, \theta, \psi, \rho, \zeta, \mu)$ or
    $(t, \theta, \psi_p, \rho, \zeta, \mu)$
    space, depending on the value of `flux0`.

    Example
    -------
    ```python title="InitialConditions definition"
    >>> initial_conditions = dex.InitialConditions(
    ...     t0=0,
    ...     flux0=("Toroidal", 0.1), # ψ0 = 0.1
    ...     theta0=3.14,
    ...     zeta0=0,
    ...     rho0=1e-4,
    ...     mu0=7e-6,
    ... )

    ```
    """

    _rust: _PyInitialConditions

    def __init__(
        self,
        t0: float,
        flux0: InitialFlux,
        theta0: float,
        zeta0: float,
        rho0: float,
        mu0: float,
    ):
        self._rust = _PyInitialConditions(
            t0=t0, flux0=flux0, theta0=theta0, zeta0=zeta0, rho0=rho0, mu0=mu0
        )

    @property
    def t0(self) -> float:
        """The initial time, in Normalized Units."""
        return self._rust.t0

    @property
    def flux0(self) -> InitialFlux:
        r"""The initial $\psi / \psi_p$, in Normalized Units."""
        return self._rust.flux0

    @property
    def theta0(self) -> float:
        r"""The initial $\theta$ angle, in Normalized Units."""
        return self._rust.theta0

    @property
    def zeta0(self) -> float:
        r"""The initial $\zeta$ angle, in Normalized Units."""
        return self._rust.zeta0

    @property
    def rho0(self) -> float:
        r"""The initial $\rho_{||}$, in Normalized Units."""
        return self._rust.rho0

    @property
    def mu0(self) -> float:
        r"""The initial magnetic moment $\mu$, in Normalized Units."""
        return self._rust.mu0

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
    ```python title="Particle creation"
    >>> initial_conditions = dex.InitialConditions(
    ...     t0=0,
    ...     flux0=("Toroidal", 0.1),
    ...     theta0=3.14,
    ...     zeta0=0,
    ...     rho0=1e-4,
    ...     mu0=7e-6,
    ... )
    >>> particle = dex.Particle(initial_conditions)

    ```
    """

    _rust: _PyParticle

    def __init__(self, initial_conditions: InitialConditions) -> None:
        self._initial_conditions = initial_conditions
        self._rust = _PyParticle(initial_conditions=initial_conditions._rust)

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
        qfactor: Qfactor,
        current: Current,
        bfield: Bfield,
        perturbation: Perturbation,
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
        qfactor
            The equilibrium's qfactor.
        current
            The equilibrium's plasma current.
        bfield
            The equilibrium's magnetic field.
        perturbation
            The equilibrium's perturbation.
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
        >>> qfactor = dex.ParabolicQfactor(qaxis=1.1, qwall=4.1, flux_wall=("Toroidal", 0.45))
        >>> current = dex.LarCurrent()
        >>> bfield = dex.LarBfield()
        >>> perturbation = dex.Perturbation(
        ...     [
        ...         dex.CosHarmonic(alpha=1e-3, m=1, n=3, phase=0),
        ...         dex.CosHarmonic(alpha=1e-3, m=2, n=3, phase=0),
        ...     ]
        ... )
        >>> initial_conditions = dex.InitialConditions(
        ...     t0=0,
        ...     flux0=("Toroidal", 0.1),
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-4,
        ...     mu0=7e-6,
        ... )
        >>> particle = dex.Particle(initial_conditions)
        >>> particle.integrate(
        ...     qfactor=qfactor,
        ...     current=current,
        ...     bfield=bfield,
        ...     perturbation=perturbation,
        ...     teval=(0, 1e2),
        ...     first_step=1e-3,
        ...     energy_rel_tol=1e-11,
        ... )

        ```
        """
        prefix = "__integrate"
        q = qfactor._dyn
        c = current._dyn
        b = bfield._dyn
        p = perturbation._dyn
        method_name: Callable = getattr(self._rust, f"{prefix}_{q}_{c}_{b}_{p}")
        method_name(
            qfactor._rust,
            current._rust,
            bfield._rust,
            perturbation._rust,
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
        qfactor: Qfactor,
        current: Current,
        bfield: Bfield,
        perturbation: Perturbation,
        intersect_params: Any,
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
        qfactor
            The equilibrium's qfactor.
        current
            The equilibrium's plasma current.
        bfield
            The equilibrium's magnetic field.
        perturbation
            The equilibrium's perturbation.
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
        >>> qfactor = dex.ParabolicQfactor(qaxis=1.1, qwall=4.1, flux_wall=("Toroidal", 0.45))
        >>> current = dex.LarCurrent()
        >>> bfield = dex.LarBfield()
        >>> perturbation = dex.Perturbation(
        ...     [
        ...         dex.CosHarmonic(alpha=1e-3, m=1, n=3, phase=0),
        ...         dex.CosHarmonic(alpha=1e-3, m=2, n=3, phase=0),
        ...     ]
        ... )
        >>> initial_conditions = dex.InitialConditions(
        ...     t0=0,
        ...     flux0=("Toroidal", 0.1),
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
        ...     qfactor=qfactor,
        ...     current=current,
        ...     bfield=bfield,
        ...     perturbation=perturbation,
        ...     intersect_params=intersect_params,
        ... )

        ```
        """
        prefix = "__intersect"
        q = qfactor._dyn
        c = current._dyn
        b = bfield._dyn
        p = perturbation._dyn
        method_name: Callable = getattr(self._rust, f"{prefix}_{q}_{c}_{b}_{p}")
        method_name(
            qfactor._rust,
            current._rust,
            bfield._rust,
            perturbation._rust,
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
