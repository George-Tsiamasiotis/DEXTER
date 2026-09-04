"""Defines the Perturbation wrapper over `_PyPerturbation`"""

import numpy as np
from collections.abc import Collection
from typing import TypeAlias

from .modes import FluteMode, NcFluteMode
from .base import _flux_eval_wrap4d
from dexter._utils import _ReprStrImpl
from dexter._core import _PyPerturbation
from dexter.types import ArrayLike, Array

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
        self._eval_p = _flux_eval_wrap4d(self._r.eval_p)
        self._eval_deriv_flux = _flux_eval_wrap4d(self._r.eval_deriv_flux)
        self._eval_deriv_theta = _flux_eval_wrap4d(self._r.eval_deriv_theta)
        self._eval_deriv_zeta = _flux_eval_wrap4d(self._r.eval_deriv_zeta)
        self._eval_deriv_t = _flux_eval_wrap4d(self._r.eval_deriv_t)

    def eval_p(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psi: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the perturbation's value $p(\psi/\psi_p, \theta, \zeta, t)$, in Normalized Units."""
        return self._eval_p(theta, zeta, t, psi=psi, psip=psip)[()]

    def eval_deriv_flux(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psi: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi/\psi_p, \theta, \zeta, t)/d(\psi/\psi_p)$, in Normalized Units."""
        return self._eval_deriv_flux(theta, zeta, t, psi=psi, psip=psip)[()]

    def eval_deriv_theta(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psi: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi/\psi_p, \theta, \zeta, t)/d\theta$, in Normalized Units."""
        return self._eval_deriv_theta(theta, zeta, t, psi=psi, psip=psip)[()]

    def eval_deriv_zeta(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psi: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi/\psi_p, \theta, \zeta, t)/d\zeta$, in Normalized Units."""
        return self._eval_deriv_zeta(theta, zeta, t, psi=psi, psip=psip)[()]

    def eval_deriv_t(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psi: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi/\psi_p, \theta, \zeta, t)/dt$, in Normalized Units."""
        return self._eval_deriv_t(theta, zeta, t, psi=psi, psip=psip)[()]

    def __len__(self) -> int:
        """Returns the total number of modes."""
        return self._r.__len__()
