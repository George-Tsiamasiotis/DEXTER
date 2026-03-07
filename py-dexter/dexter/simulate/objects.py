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

from dexter.types import Array1, InitialFlux, SteppingMethod, IntegrationStatus

# ================================================================================================
# ================================================================================================


class InitialConditions(_PyInitialConditions):
    r"""Initial conditions for a Particle.

    The initial conditions are defined on the
    $(t, \theta, \psi, \rho, \zeta, \mu)$ or
    $(t, \theta, \psi_p, \rho, \zeta, \mu)$
    space, depending on the value of `flux0`.

    Parameters
    ----------
    t0
        The initial time, in Normalized Units.
    flux0
        The initial $\psi / \psi_p$, in Normalized Units.
    theta0
        The initial $\theta$ angle, in Normalized Units.
    zeta0
        The initial $\zeta$ angle, in Normalized Units.
    rho0
        The initial $\rho_{||}$, in Normalized Units.
    mu0
        The initial magnetic moment $\mu$, in Normalized Units.

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

    t0: float
    flux0: InitialFlux
    theta0: float
    zeta0: float
    rho0: float
    mu0: float


class IntersectParams(_PyIntersectParams):
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

    intersection: str
    angle: float
    turns: int


class Particle(_PyParticle, _ParticlePlotter):
    r"""A Particle.

    By taking $\mu = 0$ and $\rho \rightarrow 0$, the particle traces magnetic field
    lines.

    Parameters
    ----------
    initial_conditions
        The initial conditions set.

    Attributes
    ----------
    initial_conditions
        The initial conditions set.
    status
        The particle's integration status.
    steps_taken
        The total number of steps taken during the integration. This number is not necessarily
        the same as the number of steps stored.
    steps_stored
        The number of steps stored in the the time series arrays.
    initial_energy
        The particle's initial energy before the integration in Normalized Units.
    final_energy
        The particle's final energy after the integration in Normalized Units.
    energy_var
        The variance of energy throughout the integration.
    t_array
        The times of the integration steps.
    psi_array
        The $ψ$ time series.
    psip_array
        The $ψp$ time series.
    theta_array
        The $\theta$ time series.
    zeta_array
        The $\zeta$ time series.
    rho_array
        The $\rho$ time series.
    mu_array
        The $\mu$ time series.
    ptheta_array
        The $P_\theta$ time series.
    pzeta_array
        The $P_\zeta$ time series.
    energy_array
        The energy time series.

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

    initial_conditions: InitialConditions
    integration_status: IntegrationStatus
    steps_taken: int
    steps_stored: int
    initial_energy: float | None
    final_energy: float | None
    energy_var: float | None
    t_array: Array1
    psi_array: Array1
    psip_array: Array1
    theta_array: Array1
    zeta_array: Array1
    rho_array: Array1
    mu_array: Array1
    ptheta_array: Array1
    pzeta_array: Array1
    energy_array: Array1

    def __init__(self, initial_conditions: InitialConditions): ...

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
        >>> perturbation = dex.CosPerturbation(
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
        method_name: Callable = getattr(self, f"{prefix}_{q}_{c}_{b}_{p}")
        method_name(
            qfactor,
            current,
            bfield,
            perturbation,
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
        >>> perturbation = dex.CosPerturbation(
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
        ...     stepping_method="ErrorAdaptiveStep",
        ... )

        ```
        """
        prefix = "__intersect"
        q = qfactor._dyn
        c = current._dyn
        b = bfield._dyn
        p = perturbation._dyn
        method_name: Callable = getattr(self, f"{prefix}_{q}_{c}_{b}_{p}")
        method_name(
            qfactor,
            current,
            bfield,
            perturbation,
            intersect_params,
            stepping_method,
            max_steps,
            first_step,
            safety_factor,
            energy_rel_tol,
            energy_abs_tol,
            error_rel_tol,
            error_abs_tol,
        )
