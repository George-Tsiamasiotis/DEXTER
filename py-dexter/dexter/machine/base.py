"""Machine objects' base classes.

Machine objects define evaluations over machine quantities, provide information about the
[`state`][dexter.types.FluxCoordinateState] of each magnetic flux coordinate, as well as useful
scalar quantities and data arrays.

Each parent class corresponds to an evaluation Trait on the Rust API.

Classes
-------
MachineObject
    Common attributes in all machine objects.
Qfactor
    q-factor related quantities and evaluation methods.
Current
    Plasma current related evaluation methods.
Bfield
    Magnetic field related evaluation methods.
Efield
    Electric field related evaluation methods.
Geometry
    Device geometry related evaluation methods.
Mode
    Single perturbation mode related evaluation methods.

"""

import numpy as np
from numpy import nan as NAN
from typing import TypeAlias, Any
from collections.abc import Callable

from dexter._core import (
    _PyQfactor,
    _PyCurrent,
    _PyBfield,
    _PyEfield,
    _PyGeometry,
    _PyMode,
)
from dexter._utils import _ReprStrImpl
from dexter.machine.flux import MagneticFlux
from dexter.types import ArrayLike, Array, Array1, FluxCoordinateState, MachineType

# Evaluation method signatures as defined in `_core.pyi`
_FluxEval1dMethod: TypeAlias = Callable[[float, float], float]
_FluxEval2dMethod: TypeAlias = Callable[[float, float, float], float]
_FluxEval4dMethod: TypeAlias = Callable[[float, float, float, float, float], float]


def _flux_eval_wrap1d(method: _FluxEval1dMethod) -> np.vectorize:
    """Wraps and vectorizes methods with signature `func(psit, psip)`."""

    def new_func(psit: float | None, psip: float | None) -> float:
        if psit is not None and psip is None:
            return method(psit, NAN)
        elif psip is not None and psit is None:
            return method(NAN, psip)
        else:
            raise TypeError("One of `psit` or `psip` must be passed")

    return np.vectorize(new_func)


def _flux_eval_wrap2d(method: _FluxEval2dMethod) -> np.vectorize:
    """Wraps and vectorizes methods with signature `func(theta, psit, psip)`."""

    def new_func(theta: float, psit: float | None, psip: float | None) -> float:
        if psit is not None and psip is None:
            return method(theta, psit, NAN)
        elif psip is not None and psit is None:
            return method(theta, NAN, psip)
        else:
            raise TypeError("One of `psit` or `psip` must be passed")

    return np.vectorize(new_func)


def _flux_eval_wrap4d(method: _FluxEval4dMethod) -> np.vectorize:
    """Wraps and vectorizes methods with signature `func(theta, zeta, t, psit, psip)`."""

    def new_func(
        theta: float,
        zeta: float,
        t: float,
        psit: float | None,
        psip: float | None,
    ) -> float:
        if psit is not None and psip is None:
            return method(theta, zeta, t, psit, NAN)
        elif psip is not None and psit is None:
            return method(theta, zeta, t, NAN, psip)
        else:
            raise TypeError("One of `psit` or `psip` must be passed")

    return np.vectorize(new_func)


class MachineObject(_ReprStrImpl):
    r"""Common attributes in all machine objects.

    Attributes
    ----------
    machine_type
        The type of the machine.
    psit_state
        The state of the toroidal flux coordinate $\psi_t$.
    psip_state
        The state of the toroidal flux coordinate $\psi_p$.

    """

    _r: Any
    machine_type: MachineType
    psit_state: FluxCoordinateState
    psip_state: FluxCoordinateState

    def __init__(self) -> None:
        self.machine_type = self._r.machine_type
        self.psit_state = self._r.psit_state
        self.psip_state = self._r.psip_state


class Geometry(_ReprStrImpl):
    r"""Geometry related evaluation methods.

    Attributes
    ----------
    baxis
        The magnetic field strength on the axis $B_0$ in $[T]$.
    raxis
        The horizontal position of the magnetic axis $R_0$ in $[m]$.
    zaxis
        The vertical position of the magnetic axis in $[m]$.
    rgeo
        The horizontal position of the geometric axis (device major radius) in $[m]$.
    rlast
        The $r$ coordinate's value at the last closed flux surface in $[m]$.
    rlab_last
        The last $R$ values that correspond to the device's last closed flux surface, in $[m]$.
    zlab_last
        The last $Z$ values that correspond to the device's last closed flux surface, in $[m]$.
    """

    _r: _PyGeometry
    baxis: float
    raxis: float
    zaxis: float
    rgeo: float
    rlast: float
    # `psit_last` and `psip_last` must be defined on the children classes

    def __init__(self) -> None:
        self.baxis = self._r.baxis
        self.raxis = self._r.raxis
        self.zaxis = self._r.zaxis
        self.rgeo = self._r.rgeo
        self.rlast = self._r.rlast
        self._eval_r = _flux_eval_wrap1d(self._r.eval_r)
        self._eval_psit_of_r = np.vectorize(self._r.eval_psit_of_r)
        self._eval_psip_of_r = np.vectorize(self._r.eval_psip_of_r)
        self._eval_rlab = _flux_eval_wrap2d(self._r.eval_rlab)
        self._eval_zlab = _flux_eval_wrap2d(self._r.eval_zlab)
        self._eval_jacobian = _flux_eval_wrap2d(self._r.eval_jacobian)

    @property
    def rlab_last(self) -> Array1:
        return self._r.rlab_last

    @property
    def zlab_last(self) -> Array1:
        return self._r.zlab_last

    def eval_r(
        self,
        psi: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $r(\psi)$, where $r$ in $[m]$."""
        return self._eval_r(psi, psip)[()]

    def eval_psit_of_r(self, r: ArrayLike) -> Array:
        r"""Calculates $\psi_t(r)$, where $r$ in $[m]$."""
        return self._eval_psit_of_r(r)[()]

    def eval_psip_of_r(self, r: ArrayLike) -> Array:
        r"""Calculates $\psi_p(r)$, where $r$ in $[m]$."""
        return self._eval_psip_of_r(r)[()]

    def eval_rlab(
        self,
        theta: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $R(\psi, \theta)$, where $R$ in $[m]$."""
        return self._eval_rlab(theta, psit, psip)[()]

    def eval_zlab(
        self,
        theta: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $Z(\psi, \theta)$, where $Z$ in $[m]$."""
        return self._eval_zlab(theta, psit, psip)[()]

    def eval_jacobian(
        self,
        theta: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array | None:
        r"""Calculates the Jacobian $J(\psi, \theta)$."""
        return self._eval_jacobian(theta, psit, psip)[()]


class Qfactor(_ReprStrImpl):
    r"""q-factor related quantities and evaluation methods.

    Attributes
    ----------
    psit_last
        The value of the last closed toroidal flux $\psi_{t,LCFS}$.
    psip_last
        The value of the last closed poloidal flux $\psi_{p,LCFS}$.
    qlast
        The q-factor's value at the last closed flux surface, $q_{LCFS}$.
    qaxis
        The q-factor's value at the magnetic axis, $q_{axis}$.
    """

    _r: _PyQfactor
    qlast: float
    qaxis: float
    psit_last: MagneticFlux
    psip_last: MagneticFlux

    def __init__(self) -> None:
        self.qlast = self._r.qlast
        self.qaxis = self._r.qaxis
        self.psit_last = MagneticFlux._wrap(self._r.psit_last)
        self.psip_last = MagneticFlux._wrap(self._r.psip_last)
        self._eval_q = _flux_eval_wrap1d(self._r.eval_q)
        self._eval_q_deriv = _flux_eval_wrap1d(self._r.eval_q_deriv)
        self._eval_other = _flux_eval_wrap1d(self._r.eval_other)
        self._eval_psit_of_q = np.vectorize(self._r.eval_psit_of_q)
        self._eval_psip_of_q = np.vectorize(self._r.eval_psip_of_q)
        self._eval_deriv_of_other = _flux_eval_wrap1d(self._r.eval_deriv_of_other)
        self._eval_deriv_wrt_other = _flux_eval_wrap1d(self._r.eval_deriv_wrt_other)
        self._eval_iota = _flux_eval_wrap1d(self._r.eval_iota)
        self._eval_iota_deriv = _flux_eval_wrap1d(self._r.eval_iota_deriv)

    def eval_q(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $q(\psi)$."""
        return self._eval_q(psit, psip)[()]

    def eval_q_deriv(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $dq(\psi)/d\psi$."""
        return self._eval_q_deriv(psit, psip)[()]

    def eval_other(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Converts a [`MagneticFlux`][dexter.MagneticFlux] to the other variant."""
        return self._eval_other(psit, psip)[()]

    def eval_psit_of_q(self, q: ArrayLike) -> Array:
        r"""Calculates $\psi_t(q)$."""
        return self._eval_psit_of_q(q)[()]

    def eval_psip_of_q(self, q: ArrayLike) -> Array:
        r"""Calculates $\psi_p(q)$."""
        return self._eval_psip_of_q(q)[()]

    def eval_deriv_of_other(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the derivative of the other magnetic flux with respect to the passed flux.

        + If `psit` is passed, then $d\psi_p/d\psi_t$ is calculated.
        + If `psip` is passed, then $d\psi_t/d\psi_p$ is calculated.

        In contrast to `Qfactor.eval_deriv_wrt_other()`, this method only requires one of the
        fluxes to be in a “good” state (the one corresponding to the passed flux argument).

        This method is useful for ensuring that $d\psi_t/d\psi_p = q$ and $d\psi_p/d\psi_t = \iota$.
        The corresponding methods `Qfactor.eval_q` and `Qfactor.eval_iota` should be used in
        calculations as they are faster and more accurate.
        """
        return self._eval_deriv_of_other(psit, psip)[()]

    def eval_deriv_wrt_other(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the derivative of the magnetic flux with respect to the other.

        + If `psit` is passed, then $d\psi_t/d\psi_p$ is calculated.
        + If `psip` is passed, then $d\psi_p/d\psi_t$ is calculated.

        This method requires both fluxes to be in a “good” state. If this is not true,
        `Qfactor.eval_deriv_of_other` should be used.

        This method is useful for ensuring that $d\psi_t/d\psi_p = q$ and $d\psi_p/d\psi_t = \iota$.
        The corresponding methods `Qfactor.eval_q` and `Qfactor.eval_iota` should be used in
        calculations as they are faster and more accurate.
        """
        return self._eval_deriv_wrt_other(psit, psip)[()]

    def eval_iota(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $\iota(\psi) = 1/q$."""
        return self._eval_iota(psit, psip)[()]

    def eval_iota_deriv(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $d\iota(\psi)/d\psi = -q'/q^2$."""
        return self._eval_iota_deriv(psit, psip)[()]


class Current(_ReprStrImpl):
    """Plasma current related evaluation methods."""

    _r: _PyCurrent

    def __init__(self) -> None:
        self._eval_g = _flux_eval_wrap1d(self._r.eval_g)
        self._eval_i = _flux_eval_wrap1d(self._r.eval_i)
        self._eval_g_deriv = _flux_eval_wrap1d(self._r.eval_g_deriv)
        self._eval_i_deriv = _flux_eval_wrap1d(self._r.eval_i_deriv)
        self._eval_g_deriv2 = _flux_eval_wrap1d(self._r.eval_g_deriv2)
        self._eval_i_deriv2 = _flux_eval_wrap1d(self._r.eval_i_deriv2)

    def eval_g(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $g(\psi)$."""
        return self._eval_g(psit, psip)[()]

    def eval_i(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $I(\psi)$."""
        return self._eval_i(psit, psip)[()]

    def eval_g_deriv(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $dg(\psi)/d\psi$."""
        return self._eval_g_deriv(psit, psip)[()]

    def eval_i_deriv(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $dI(\psi)/d\psi$."""
        return self._eval_i_deriv(psit, psip)[()]

    def eval_g_deriv2(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $d^2g(\psi)/d\psi^2$."""
        return self._eval_g_deriv2(psit, psip)[()]

    def eval_i_deriv2(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $d^2I(\psi)/d\psi^2$."""
        return self._eval_i_deriv2(psit, psip)[()]


class Bfield(_ReprStrImpl):
    """Magnetic field related evaluation methods."""

    _r: _PyBfield

    def __init__(self) -> None:
        self._eval_b = _flux_eval_wrap2d(self._r.eval_b)
        self._eval_deriv_flux = _flux_eval_wrap2d(self._r.eval_deriv_flux)
        self._eval_deriv_theta = _flux_eval_wrap2d(self._r.eval_deriv_theta)
        self._eval_deriv_flux2 = _flux_eval_wrap2d(self._r.eval_deriv_flux2)
        self._eval_deriv_theta2 = _flux_eval_wrap2d(self._r.eval_deriv_theta2)
        self._eval_deriv_mixed = _flux_eval_wrap2d(self._r.eval_deriv_mixed)

    def eval_b(
        self,
        theta: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $B(\psi, \theta)$."""
        return self._eval_b(theta, psit, psip)[()]

    def eval_deriv_flux(
        self,
        theta: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $dB(\psi, \theta)/d\psi$."""
        return self._eval_deriv_flux(theta, psit, psip)[()]

    def eval_deriv_theta(
        self,
        theta: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $dB(\psi, \theta)/d\theta$."""
        return self._eval_deriv_theta(theta, psit, psip)[()]

    def eval_deriv_flux2(
        self,
        theta: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $d^2B(\psi, \theta)/d\psi^2$."""
        return self._eval_deriv_flux2(theta, psit, psip)[()]

    def eval_deriv_theta2(
        self,
        theta: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $d^2B(\psi, \theta)/d\theta^2$."""
        return self._eval_deriv_theta2(theta, psit, psip)[()]

    def eval_deriv_mixed(
        self,
        theta: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $d^2B(\psi, \theta)/d\psi d\theta$."""
        return self._eval_deriv_mixed(theta, psit, psip)[()]


class Efield(_ReprStrImpl):
    """Electric field related evaluation methods."""

    _r: _PyEfield

    def __init__(self) -> None:
        self._eval_e = _flux_eval_wrap1d(self._r.eval_e)
        self._eval_phi = _flux_eval_wrap1d(self._r.eval_phi)
        self._eval_phi_deriv = _flux_eval_wrap1d(self._r.eval_phi_deriv)
        self._eval_phi_deriv2 = _flux_eval_wrap1d(self._r.eval_phi_deriv2)

    def eval_e(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $E(\psi)$."""
        return self._eval_e(psit, psip)[()]

    def eval_phi(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $\Phi(\psi)$."""
        return self._eval_phi(psit, psip)[()]

    def eval_phi_deriv(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $d\Phi(\psi)/d\psi$.

        Since the field is radial, this method is equivalent to $-E(\psi)$.
        """
        return self._eval_phi_deriv(psit, psip)[()]

    def eval_phi_deriv2(
        self,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates $d^2\Phi(\psi)/d\psi^2$."""
        return self._eval_phi_deriv2(psit, psip)[()]


class Mode(_ReprStrImpl):
    r"""Single perturbation mode related evaluation methods.

    Attributes
    ----------
    m
        The poloidal number $m$.
    n
        The toroidal number $n$.
    """

    _r: _PyMode
    m: int
    n: int

    def __init__(self) -> None:
        self.m = self._r.m
        self.n = self._r.n

        self._eval_amplitude = _flux_eval_wrap4d(self._r.eval_amplitude)
        self._eval_phase = _flux_eval_wrap4d(self._r.eval_phase)
        self._eval_m = _flux_eval_wrap4d(self._r.eval_m)
        self._eval_deriv_flux = _flux_eval_wrap4d(self._r.eval_deriv_flux)
        self._eval_deriv_theta = _flux_eval_wrap4d(self._r.eval_deriv_theta)
        self._eval_deriv_zeta = _flux_eval_wrap4d(self._r.eval_deriv_zeta)
        self._eval_deriv_t = _flux_eval_wrap4d(self._r.eval_deriv_t)

    def eval_amplitude(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates **the amplitude** $\alpha_{m,n}(\psi, \theta, \zeta, t)$, in Normalized Units."""
        return self._eval_amplitude(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_phase(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates **the phase** $\phi(\psi, \theta, \zeta, t)$."""
        return self._eval_phase(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_m(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the mode's value $m(\psi, \theta, \zeta, t)$, in Normalized Units."""
        return self._eval_m(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_deriv_flux(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi, \theta, \zeta, t)/d\psi$, in Normalized Units."""
        return self._eval_deriv_flux(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_deriv_theta(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi, \theta, \zeta, t)/d\theta$, in Normalized Units."""
        return self._eval_deriv_theta(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_deriv_zeta(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi, \theta, \zeta, t)/d\zeta$, in Normalized Units."""
        return self._eval_deriv_zeta(theta, zeta, t, psit=psit, psip=psip)[()]

    def eval_deriv_t(
        self,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
        psit: ArrayLike | None = None,
        psip: ArrayLike | None = None,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi, \theta, \zeta, t)/dt$, in Normalized Units."""
        return self._eval_deriv_t(theta, zeta, t, psit=psit, psip=psip)[()]
