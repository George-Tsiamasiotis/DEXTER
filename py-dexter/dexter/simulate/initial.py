"""Defines types associated with Particle and Queue initialization."""

from dexter.machine.utils import MagneticFlux
from dexter.types import CoordinateSet

from dexter._utils import _ReprStrImpl, _RustTypeWrapper
from dexter._core import _PyInitialConditions


class InitialConditions(_ReprStrImpl, _RustTypeWrapper):
    r"""Initial conditions set for a Particle.

    This type is instantiated through the [`Boozer`][dexter.InitialConditions.Boozer] and
    [`Mixed`][dexter.InitialConditions.Mixed] class methods.

    Attributes
    ----------
    t0
        The initial time $t_0$.
    flux0
        The initial magnetic flux $\psi_0$ or $\psi_{p,0}$.
    theta0
        The initial $\theta$ angle.
    theta0
        The initial $\zeta$ angle.
    rho0
        The initial $\rho_{||,0}$. If the set was initialized from [`Mixed`][dexter.InitialConditions.Mixed]
        and no particle routines have run, then `rho0` is `None`.
    pzeta0
        The initial $P_{\zeta,0}$. If the set was initialized from [`Boozer`][dexter.InitialConditions.Boozer]
        and no particle routines have run, then `pzeta0` is `None`.
    coordinate_set
        The kind of initial conditions set.
    """

    _r: _PyInitialConditions
    t0: float
    flux0: MagneticFlux
    theta0: float
    zeta0: float
    mu0: float
    coordinate_set: CoordinateSet

    def __init__(self) -> None:
        raise RuntimeError("Cannot instantiate class")

    @classmethod
    def Boozer(
        cls,
        t0: float,
        flux0: MagneticFlux,
        theta0: float,
        zeta0: float,
        rho0: float,
        mu0: float,
    ) -> InitialConditions:
        r"""Creates initial conditions for a Particle in Boozer coordinates.

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
        ```python title="InitialConditions definition in Boozer-Toroidal coordinates"
        >>> initial_conditions = dex.InitialConditions.Boozer(
        ...     t0=0,
        ...     flux0=dex.MagneticFlux.Toroidal(0.01),  # ψ0 = 0.01
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-4,
        ...     mu0=7e-6,
        ... )

        ```
        """
        obj = InitialConditions.__new__(InitialConditions)
        obj._r = _PyInitialConditions.boozer(t0, flux0._r, theta0, zeta0, rho0, mu0)
        obj.t0 = obj._r.t0
        obj.flux0 = MagneticFlux._wrap(obj._r.flux0)
        obj.theta0 = obj._r.theta0
        obj.zeta0 = obj._r.zeta0
        obj.mu0 = obj._r.mu0
        obj.coordinate_set = obj._r.coordinate_set

        return obj

    @classmethod
    def Mixed(
        cls,
        t0: float,
        flux0: MagneticFlux,
        theta0: float,
        zeta0: float,
        pzeta0: float,
        mu0: float,
    ) -> InitialConditions:
        r"""Creates initial conditions for a Particle in Mixed coordinates.

        The initial conditions are defined on the
        $(t, \psi, \theta, \zeta, P_\zeta, \mu)$ or
        $(t, \psi_p, \theta, \zeta, P_\zeta, \mu)$
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
        pzeta0
            The initial $P_\zeta$, in Normalized Units.
        mu0
            The initial magnetic moment $\mu$, in Normalized Units.

        Example
        -------
        ```python title="InitialConditions definition in Mixed-Poloidal coordinates"
        >>> flux0=dex.MagneticFlux.Poloidal(0.02)  # ψ0 = 0.02
        >>> initial_conditions = dex.InitialConditions.Mixed(
        ...     t0=0,
        ...     flux0=flux0,
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     pzeta0=-0.4 * flux0.value,
        ...     mu0=7e-6,
        ... )

        ```
        """
        obj = InitialConditions.__new__(InitialConditions)
        obj._r = _PyInitialConditions.mixed(t0, flux0._r, theta0, zeta0, pzeta0, mu0)
        obj.t0 = obj._r.t0
        obj.flux0 = MagneticFlux._wrap(obj._r.flux0)
        obj.theta0 = obj._r.theta0
        obj.zeta0 = obj._r.zeta0
        obj.mu0 = obj._r.mu0
        obj.coordinate_set = obj._r.coordinate_set

        return obj

    @property
    def rho0(self) -> float:
        r"""The initial parallel radius $\rho_{||}$, in Normalized Units."""
        if self._r.rho0 is None:
            raise AttributeError("'rho0' has not been defined")
        else:
            return self._r.rho0

    @property
    def pzeta0(self) -> float:
        r"""The initial canonical momentum $P_\zeta$, in Normalized Units."""
        if self._r.pzeta0 is None:
            raise AttributeError("'pzeta0' has not been defined")
        else:
            return self._r.pzeta0
