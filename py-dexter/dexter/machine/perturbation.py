"""Defines the Perturbation wrapper over `_PyPerturbation`"""

from collections.abc import Collection
from typing import TypeAlias

from .modes import FluteMode, NcFluteMode
from .base import _flux_eval_wrap4d
from dexter._utils import _ReprStrImpl
from dexter._core import _PyPerturbation
from dexter.types import Array, ArrayLike

ModeObject: TypeAlias = FluteMode | NcFluteMode


class Perturbation(_ReprStrImpl):
    """A container type for [`ModeObjects`](modes), representing the total perturbation of the system.

    A `Perturbation` can consist of an arbitrary amount of modes, not necessarily of the same type.

    Parameters
    ----------
    modes
        Collection containing the comprising modes. If `None`, the perturbation is 0.

    Example
    -------
    ``` py title="Perturbation consisting of analytical flute modes"
    >>> LCFS = dex.MagneticFlux.Toroidal(0.05)
    >>> per = dex.Perturbation(
    ...     [
    ...         dex.FluteMode(1e-4, LCFS, 3, 2, 0),
    ...         dex.FluteMode(2e-4, LCFS, 4, 3, 0),
    ...         dex.FluteMode(3e-4, LCFS, 5, 4, 0),
    ...     ]
    ... )

    ```

    Example
    -------
    ``` py title="Perturbation consisting of flute modes from a netCDF file"
    >>> LCFS = dex.MagneticFlux.Toroidal(0.05)
    >>> per = dex.Perturbation(
    ...     [
    ...         dex.NcFluteMode(path, "Cubic", 2, 1),
    ...         dex.NcFluteMode(path, "Cubic", 3, 2),
    ...     ]
    ... )

    ```

    Example
    -------
    ``` py title="Perturbation consisting of mixed types of flute modes"
    >>> LCFS = dex.MagneticFlux.Toroidal(0.05)
    >>> per = dex.Perturbation(
    ...     [
    ...         dex.FluteMode(1e-4, LCFS, 3, 2, 0),
    ...         dex.FluteMode(2e-4, LCFS, 5, 2, 0),
    ...         dex.NcFluteMode(path, "Cubic", 3, 2),
    ...     ]
    ... )

    ```

    """

    _r: _PyPerturbation

    def __init__(self, modes: Collection[ModeObject] | None = None) -> None:
        if modes is not None:
            _modes = [mode._r for mode in modes]
        else:
            _modes = []
        self._r = _PyPerturbation(_modes)
        self._eval_ms = _flux_eval_wrap4d(self._r.eval_ms)
        self._eval_dms_dflux = _flux_eval_wrap4d(self._r.eval_dms_dflux)
        self._eval_dms_dtheta = _flux_eval_wrap4d(self._r.eval_dms_dtheta)
        self._eval_dms_dzeta = _flux_eval_wrap4d(self._r.eval_dms_dzeta)
        self._eval_dms_dt = _flux_eval_wrap4d(self._r.eval_dms_dt)
        self._eval_phis = _flux_eval_wrap4d(self._r.eval_phis)
        self._eval_dphis_dflux = _flux_eval_wrap4d(self._r.eval_dphis_dflux)
        self._eval_dphis_dtheta = _flux_eval_wrap4d(self._r.eval_dphis_dtheta)
        self._eval_dphis_dzeta = _flux_eval_wrap4d(self._r.eval_dphis_dzeta)
        self._eval_dphis_dt = _flux_eval_wrap4d(self._r.eval_dphis_dt)

    def eval_ms(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the perturbation's value as a function of $(\psi, \theta, \zeta, t)$."""
        return self._eval_ms(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_dms_dflux(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the perturbation's derivative with respect to $\psi$."""
        return self._eval_dms_dflux(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_dms_dtheta(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the perturbation's derivative with respect to $\theta$."""
        return self._eval_dms_dtheta(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_dms_dzeta(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the perturbation's derivative with respect to $\zeta$."""
        return self._eval_dms_dzeta(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_dms_dt(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the perturbation's derivative with respect to $t$."""
        return self._eval_dms_dt(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_phis(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the value of the Perturbation’s induced electric potential."""
        return self._eval_phis(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_dphis_dflux(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the derivative of the Perturbation’s induced electric potential with respect to $\psi$."""
        return self._eval_dphis_dflux(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_dphis_dtheta(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the derivative of the Perturbation’s induced electric potential with respect to $\theta$."""
        return self._eval_dphis_dtheta(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_dphis_dzeta(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the derivative of the Perturbation’s induced electric potential with respect to $\zeta$."""
        return self._eval_dphis_dzeta(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_dphis_dt(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the derivative of the Perturbation’s induced electric potential with respect to $t$."""
        return self._eval_dphis_dt(theta, zeta, t, psit=psit, psip=psip)[()]

    def __len__(self) -> int:
        """Returns the total number of modes."""
        return self._r.__len__()
