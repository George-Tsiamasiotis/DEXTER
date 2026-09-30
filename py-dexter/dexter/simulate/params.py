"""Helper types to pass parameters to the solver.

Classes
-------
SteppingMethod
    Sets the solver's stepping method.
SolverParams
    Container for the solver parameters.
IntersectParams
    Parameters for the calculation of Poincare maps.
"""

from dexter._core import _PySteppingMethod, _PySolverParams, _PyIntersectParams
from dexter._utils import _ReprStrImpl

from dexter.types import Directionality, Intersection


class SteppingMethod(_ReprStrImpl):
    r"""The method used to calculate the next optimal step.

    See `classmethods` for the different stepping methods.
    """

    _r: _PySteppingMethod

    def __init__(self) -> None:
        raise RuntimeError("Cannot instantiate class")

    @classmethod
    def EnergyAdaptiveStep(
        cls, rel_tol: float = 1e-8, abs_tol: float = 1e-10
    ) -> SteppingMethod:
        r"""Forces the step size to be small enough so that the Energy difference from step
        to step is under a certain threshold.

        Parameters
        ----------
        rel_tol
            The relative tolerance to compare the relative energy error with.
        abs_tol
            The absolute error tolerance. Prevents the relative energy error from becoming too
            small, causing the particle to get stuck.

        Example
        -------

        ``` py
        >>> method = dex.SteppingMethod.EnergyAdaptiveStep(1e-11, 1e-14)

        ```
        """
        obj = SteppingMethod.__new__(SteppingMethod)
        obj._r = _PySteppingMethod.energy_adaptive_step(
            rel_tol=rel_tol, abs_tol=abs_tol
        )
        return obj

    @classmethod
    def ErrorAdaptiveStep(cls, rel_tol: float, abs_tol: float) -> SteppingMethod:
        r"""Classic RK error estimation : Adjust the step size to minimize the local truncation error.

        Note that the errors are not normalized, therefore the tolerances must be set according
        to the scale of the system’s time derivatives. A good starting point is
        `#!python rel_tol=1e-17` and `#!python abs_tol=1e-19`.

        Parameters
        ----------
        rel_tol
            The relative tolerance to compare the relative error with.
        abs_tol
            The absolute error tolerance. Prevents the relative error from becoming too small,
            causing the particle to get stuck.

        Example
        -------

        ``` py
        >>> method = dex.SteppingMethod.ErrorAdaptiveStep(1e-17, 1e-19)

        ```
        """
        obj = SteppingMethod.__new__(SteppingMethod)
        obj._r = _PySteppingMethod.error_adaptive_step(rel_tol=rel_tol, abs_tol=abs_tol)
        return obj

    @classmethod
    def FixedStep(cls, step: float) -> SteppingMethod:
        r"""Fixed time step.

        Parameters
        ----------
        step
            The time step in Normalized units.

        Example
        -------

        ``` py
        >>> method = dex.SteppingMethod.FixedStep(0.01)

        ```
        """
        obj = SteppingMethod.__new__(SteppingMethod)
        obj._r = _PySteppingMethod.fixed_step(step)
        return obj


class SolverParams(_ReprStrImpl):
    r"""Container for the solver parameters.

    Parameters
    ----------
    method:
        The optimal step calculation method.
    max_steps
        The maximum amount of steps a particle can make before terminating its integration.
    first_step
        The initial time step for the RKF45 adaptive step method in Normalized Units. The value is
        empirical.
    safety_factor
        The safety factor of the solver. Should be less than $1.0$.

    Example
    -------
    ``` py
    >>> solver_params = dex.SolverParams(
    ...     method=dex.SteppingMethod.EnergyAdaptiveStep(rel_tol=1e-7, abs_tol=1e-9),
    ...     max_steps=100_000,
    ...     first_step=1e-3,
    ... )

    ```

    """

    _r: _PySolverParams

    def __init__(
        self,
        method: SteppingMethod = SteppingMethod.EnergyAdaptiveStep(),
        max_steps: int = 1_000_000,
        first_step: float = 1e-1,
        safety_factor: float = 0.9,
    ) -> None:
        self._r = _PySolverParams(
            method=method._r,
            max_steps=max_steps,
            first_step=first_step,
            safety_factor=safety_factor,
        )


class IntersectParams(_ReprStrImpl):
    r"""Defines all necessary parameters for the [`Particle.intersect`][dexter.Particle.intersect] routine.

    Parameters
    ----------
    intersection
        The surface of section $\Sigma$, defined by an equation $x_i=\alpha$,
        where $x_i = \theta$ or $\zeta$.
    angle
        The constant that defines the surface of section.
    turns
        The number of intersections to calculate.
    directionality
        The method with which Poincare intersections are recorded.

    Example
    -------

    ``` py
    >>> intersect_params = dex.IntersectParams("ConstTheta", 0, 800)
    >>> intersect_params = dex.IntersectParams("ConstZeta", np.pi, 1000)

    ```

    Attributes
    ----------
    intersection
        The surface of section $\Sigma$, defined by an equation $x_i=\alpha$,
        where $x_i = \theta$ or $\zeta$.
    angle
        The constant that defines the surface of section.
    turns
        The number of intersections to calculate.
    directionality
        The method with which Poincare intersections are recorded.

    """

    _r: _PyIntersectParams

    def __init__(
        self,
        intersection: Intersection,
        angle: float,
        turns: int,
        directionality: Directionality = "Initial",
    ) -> None:
        self._r = _PyIntersectParams(
            intersection=intersection,
            angle=angle,
            turns=turns,
            directionality=directionality,
        )

    @property
    def intersection(self) -> Intersection:
        return self._r.intersection

    @property
    def angle(self) -> float:
        return self._r.angle

    @property
    def turns(self) -> int:
        return self._r.turns

    @property
    def directionality(self) -> Directionality:
        return self._r.directionality
