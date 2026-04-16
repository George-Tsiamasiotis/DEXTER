"""Evaluation traits implementation and documentation.

To vectorize the provided methods, each Trait's constructor must be called from the child object
*after* the `_rust` object has been created.
"""

import numpy as np
from numpy.typing import ArrayLike, NDArray

from dexter._core import _PyLarGeometry, _PyNcGeometry
from dexter._core import _PyUnityQfactor, _PyParabolicQfactor, _PyNcQfactor
from dexter._core import _PyLarCurrent, _PyNcCurrent
from dexter._core import _PyLarBfield, _PyNcBfield
from dexter._core import _PyCosHarmonic, _PyNcHarmonic
from dexter.types import FluxState


class _FluxCommuteTrait:
    """Documents the methods provided by the 'FluxCommute' trait."""

    _rust: _PyNcGeometry | _PyUnityQfactor | _PyParabolicQfactor | _PyNcQfactor
    """`FluxCommute` implementors"""

    def __init__(self) -> None:
        self._psi_of_psip = np.vectorize(self._rust.psi_of_psip)
        self._psip_of_psi = np.vectorize(self._rust.psip_of_psi)

    def psip_of_psi(self, psi: ArrayLike) -> NDArray:
        r"""The $\psi_p(\psi)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._psip_of_psi(psi)[()]

    def psi_of_psip(self, psip: ArrayLike) -> NDArray:
        r"""The $\psi(\psi_p)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._psi_of_psip(psip)[()]


class _GeometryTrait:
    """Documents the methods provided by the 'Geometries' trait."""

    _rust: _PyLarGeometry | _PyNcGeometry
    """`Geometry` implementors"""

    def __init__(self) -> None:
        self._r_of_psi = np.vectorize(self._rust.r_of_psi)
        self._r_of_psip = np.vectorize(self._rust.r_of_psip)
        self._psi_of_r = np.vectorize(self._rust.psi_of_r)
        self._psip_of_r = np.vectorize(self._rust.psip_of_r)
        self._rlab_of_psi = np.vectorize(self._rust.rlab_of_psi)
        self._rlab_of_psip = np.vectorize(self._rust.rlab_of_psip)
        self._zlab_of_psi = np.vectorize(self._rust.zlab_of_psi)
        self._zlab_of_psip = np.vectorize(self._rust.zlab_of_psip)
        self._jacobian_of_psi = np.vectorize(self._rust.jacobian_of_psi)
        self._jacobian_of_psip = np.vectorize(self._rust.jacobian_of_psip)

    @property
    def psi_state(self) -> FluxState:
        """The state of the toroidal flux coordinate."""
        return self._rust.psi_state

    @property
    def psip_state(self) -> FluxState:
        """The state of the poloidal flux coordinate."""
        return self._rust.psip_state

    def r_of_psi(self, psi: ArrayLike) -> NDArray:
        r"""The $r(\psi)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._r_of_psi(psi)[()]

    def r_of_psip(self, psip: ArrayLike) -> NDArray:
        r"""The $r(\psi_p)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._r_of_psip(psip)[()]

    def psi_of_r(self, r: ArrayLike) -> NDArray:
        r"""The $\psi(r)$ value in Normalized Units.

        Parameters
        ----------
        r
            The radial distance $r$ in $[m]$.
        """
        return self._psi_of_r(r)[()]

    def psip_of_r(self, r: ArrayLike) -> NDArray:
        r"""The $\psi_p(r)$ value in Normalized Units.

        Parameters
        ----------
        r
            The radial distance $r$ in $[m]$.
        """
        return self._psip_of_r(r)[()]

    def rlab_of_psi(self, psi: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $R_{lab}(\psi, \theta)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rlab_of_psi(psi, theta)[()]

    def rlab_of_psip(self, psip: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $R_{lab}(\psi_p, \theta)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rlab_of_psip(psip, theta)[()]

    def zlab_of_psi(self, psi: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $Z_{lab}(\psi, \theta)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._zlab_of_psi(psi, theta)[()]

    def zlab_of_psip(self, psip: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $Z_{lab}(\psi_p, \theta)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._zlab_of_psip(psip, theta)[()]

    def jacobian_of_psi(self, psi: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $J(\psi, \theta)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._jacobian_of_psi(psi, theta)[()]

    def jacobian_of_psip(self, psip: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $J(\psi_p, \theta)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._jacobian_of_psip(psip, theta)[()]

    # Implemented as properties
    # def rlab_last(self) -> Array1:
    # def zlab_last(self) -> Array1:


class _QfactorTrait:
    """Documents the methods provided by the 'Qfactor' trait."""

    _rust: _PyUnityQfactor | _PyParabolicQfactor | _PyNcQfactor
    """`Qfactor` implementors"""

    def __init__(self) -> None:
        self._q_of_psi = np.vectorize(self._rust.q_of_psi)
        self._q_of_psip = np.vectorize(self._rust.q_of_psip)
        self._iota_of_psi = np.vectorize(self._rust.iota_of_psi)
        self._iota_of_psip = np.vectorize(self._rust.iota_of_psip)
        self._dpsip_dpsi = np.vectorize(self._rust.dpsip_dpsi)
        self._dpsi_dpsip = np.vectorize(self._rust.dpsi_dpsip)
        self._psi_of_q = np.vectorize(self._rust.psi_of_q)
        self._psip_of_q = np.vectorize(self._rust.psip_of_q)

    @property
    def psi_state(self) -> FluxState:
        """The state of the toroidal flux coordinate."""
        return self._rust.psi_state

    @property
    def psip_state(self) -> FluxState:
        """The state of the poloidal flux coordinate."""
        return self._rust.psip_state

    def q_of_psi(self, psi: ArrayLike) -> NDArray:
        r"""The $q(\psi)$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._q_of_psi(psi)[()]

    def q_of_psip(self, psip: ArrayLike) -> NDArray:
        r"""The $q(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._q_of_psip(psip)[()]

    def dpsip_dpsi(self, psi: ArrayLike) -> NDArray:
        r"""The derivative $d\psi_p(\psi)/d\psi$ value in Normalized Units.

        It's a good check that the values coincide with `qfactor.iota_of_psi(psi)`.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._dpsip_dpsi(psi)[()]

    def dpsi_dpsip(self, psip: ArrayLike) -> NDArray:
        r"""The derivative $d\psi(\psi_p)/d\psi_p$ value in Normalized Units.

        It's a good check that the values coincide with `qfactor.q_of_psip(psip)`.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._dpsi_dpsip(psip)[()]

    def iota_of_psi(self, psi: ArrayLike) -> NDArray:
        r"""The $\iota(\psi) = \dfrac{1}{q(\psi)}$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._iota_of_psi(psi)[()]

    def iota_of_psip(self, psip: ArrayLike) -> NDArray:
        r"""The $\iota(\psi_p) = \dfrac{1}{q(\psi_p)}$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._iota_of_psip(psip)[()]

    def psi_of_q(self, q: ArrayLike) -> NDArray:
        r"""The toroidal flux $\psi$ in Normalized Units.

        Parameters
        ----------
        q
            The q-factor value $q$.
        """
        return self._psi_of_q(q)[()]

    def psip_of_q(self, q: ArrayLike) -> NDArray:
        r"""The poloidal flux $\psi_p$ in Normalized Units.

        Parameters
        ----------
        q
            The q-factor value $q$.
        """
        return self._psip_of_q(q)[()]


class _CurrentTrait:
    """Documents the methods provided by the 'Current' trait."""

    _rust: _PyLarCurrent | _PyNcCurrent
    """`Current` implementors"""

    def __init__(self) -> None:
        self._g_of_psi = np.vectorize(self._rust.g_of_psi)
        self._g_of_psip = np.vectorize(self._rust.g_of_psip)
        self._i_of_psi = np.vectorize(self._rust.i_of_psi)
        self._i_of_psip = np.vectorize(self._rust.i_of_psip)
        self._dg_dpsi = np.vectorize(self._rust.dg_dpsi)
        self._dg_dpsip = np.vectorize(self._rust.dg_dpsip)
        self._di_dpsi = np.vectorize(self._rust.di_dpsi)
        self._di_dpsip = np.vectorize(self._rust.di_dpsip)

    @property
    def psi_state(self) -> FluxState:
        """The state of the toroidal flux coordinate."""
        return self._rust.psi_state

    @property
    def psip_state(self) -> FluxState:
        """The state of the poloidal flux coordinate."""
        return self._rust.psip_state

    def g_of_psi(self, psi: ArrayLike) -> NDArray:
        r"""The $g(\psi)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._g_of_psi(psi)[()]

    def g_of_psip(self, psip: ArrayLike) -> NDArray:
        r"""The $g(\psi_p)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._g_of_psip(psip)[()]

    def i_of_psi(self, psi: ArrayLike) -> NDArray:
        r"""The $I(\psi)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._i_of_psi(psi)[()]

    def i_of_psip(self, psip: ArrayLike) -> NDArray:
        r"""The $I(\psi_p)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._i_of_psip(psip)[()]

    def dg_dpsi(self, psi: ArrayLike) -> NDArray:
        r"""The $dg(\psi)/d\psi$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._dg_dpsi(psi)[()]

    def dg_dpsip(self, psip: ArrayLike) -> NDArray:
        r"""The $dg(\psi_p)/d\psi_p$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi$ in Normalized Units.
        """
        return self._dg_dpsip(psip)[()]

    def di_dpsi(self, psi: ArrayLike) -> NDArray:
        r"""The $dg(\psi)/d\psi$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._di_dpsi(psi)[()]

    def di_dpsip(self, psip: ArrayLike) -> NDArray:
        r"""The $dg(\psi_p)/d\psi_p$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi$ in Normalized Units.
        """
        return self._di_dpsip(psip)[()]


class _BfieldTrait:
    """Documents the methods provided by the 'Bfield' trait."""

    _rust: _PyLarBfield | _PyNcBfield
    """`Bfield` implementors"""

    def __init__(self) -> None:
        self._b_of_psi = np.vectorize(self._rust.b_of_psi)
        self._b_of_psip = np.vectorize(self._rust.b_of_psip)
        self._db_dpsi = np.vectorize(self._rust.db_dpsi)
        self._db_dpsip = np.vectorize(self._rust.db_dpsip)
        self._db_of_psi_dtheta = np.vectorize(self._rust.db_of_psi_dtheta)
        self._db_of_psip_dtheta = np.vectorize(self._rust.db_of_psip_dtheta)

    @property
    def psi_state(self) -> FluxState:
        """The state of the toroidal flux coordinate."""
        return self._rust.psi_state

    @property
    def psip_state(self) -> FluxState:
        """The state of the poloidal flux coordinate."""
        return self._rust.psip_state

    def b_of_psi(self, psi: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $B(\psi, \theta)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._b_of_psi(psi, theta)[()]

    def b_of_psip(self, psip: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $B(\psi_p, \theta)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._b_of_psip(psip, theta)[()]

    def db_dpsi(self, psi: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $dB(\psi, \theta)/d\psi$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._db_dpsi(psi, theta)[()]

    def db_dpsip(self, psip: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $dB(\psi_p, \theta)/d\psi_p$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._db_dpsip(psip, theta)[()]

    def db_of_psi_dtheta(self, psi: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $dB(\psi, \theta)/d\theta$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._db_of_psi_dtheta(psi, theta)[()]

    def db_of_psip_dtheta(self, psip: ArrayLike, theta: ArrayLike) -> NDArray:
        r"""The $dB(\psi_p, \theta)/d\theta$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._db_of_psip_dtheta(psip, theta)[()]


class _HarmonicTrait:
    """Documents the methods provided by the 'Harmonic' trait."""

    _rust: _PyCosHarmonic | _PyNcHarmonic
    """`Harmonic` implementors"""

    def __init__(self) -> None:
        self._alpha_of_psi = np.vectorize(self._rust.alpha_of_psi)
        self._alpha_of_psip = np.vectorize(self._rust.alpha_of_psip)
        self._phase_of_psi = np.vectorize(self._rust.phase_of_psi)
        self._phase_of_psip = np.vectorize(self._rust.phase_of_psip)

        self._h_of_psi = np.vectorize(self._rust.h_of_psi)
        self._h_of_psip = np.vectorize(self._rust.h_of_psip)
        self._dh_dpsi = np.vectorize(self._rust.dh_dpsi)
        self._dh_dpsip = np.vectorize(self._rust.dh_dpsip)
        self._dh_of_psi_dtheta = np.vectorize(self._rust.dh_of_psi_dtheta)
        self._dh_of_psip_dtheta = np.vectorize(self._rust.dh_of_psip_dtheta)
        self._dh_of_psi_dzeta = np.vectorize(self._rust.dh_of_psi_dzeta)
        self._dh_of_psip_dzeta = np.vectorize(self._rust.dh_of_psip_dzeta)
        self._dh_of_psi_dt = np.vectorize(self._rust.dh_of_psi_dt)
        self._dh_of_psip_dt = np.vectorize(self._rust.dh_of_psip_dt)

    @property
    def psi_state(self) -> FluxState:
        """The state of the toroidal flux coordinate."""
        return self._rust.psi_state

    @property
    def psip_state(self) -> FluxState:
        """The state of the poloidal flux coordinate."""
        return self._rust.psip_state

    def alpha_of_psi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's amplitude $\alpha(\psi, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._alpha_of_psi(psi, theta, zeta, t)[()]

    def alpha_of_psip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's amplitude $\alpha(\psi_p, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._alpha_of_psip(psip, theta, zeta, t)[()]

    def phase_of_psi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's phase $\phi(\psi, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._phase_of_psi(psi, theta, zeta, t)[()]

    def phase_of_psip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's phase $\phi(\psi_p, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._phase_of_psip(psip, theta, zeta, t)[()]

    def h_of_psi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The full harmonic's value $h(\psi, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._h_of_psi(psi, theta, zeta, t)[()]

    def h_of_psip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The full harmonic's value $h(\psi_p, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._h_of_psip(psip, theta, zeta, t)[()]

    def dh_dpsi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's derivative with respect to $\psi$, $\partial h(\psi, \theta, \zeta, t)/\partial\psi$
        in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dh_dpsi(psi, theta, zeta, t)[()]

    def dh_dpsip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's derivative with respect to $\psi_p$, $\partial h(\psi_p, \theta, \zeta, t)/\partial \psi_p$
        in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dh_dpsip(psip, theta, zeta, t)[()]

    def dh_of_psi_dtheta(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's derivative with respect to $\theta$, $\partial h(\psi, \theta, \zeta, t)/\partial \theta$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dh_of_psi_dtheta(psi, theta, zeta, t)[()]

    def dh_of_psip_dtheta(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's derivative with respect to $\theta$, $\partial h(\psi_p, \theta, \zeta, t)/\partial \theta$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dh_of_psip_dtheta(psip, theta, zeta, t)[()]

    def dh_of_psi_dzeta(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's derivative with respect to $\zeta$, $\partial h(\psi, \theta, \zeta, t)/\partial \zeta$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dh_of_psi_dzeta(psi, theta, zeta, t)[()]

    def dh_of_psip_dzeta(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's derivative with respect to $\zeta$, $\partial h(\psi_p, \theta, \zeta, t)/\partial \zeta$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dh_of_psip_dzeta(psip, theta, zeta, t)[()]

    def dh_of_psi_dt(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's derivative with respect to the time $t$, $\partial h(\psi, \theta, \zeta, t)/\partial t$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dh_of_psi_dt(psi, theta, zeta, t)[()]

    def dh_of_psip_dt(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The harmonic's derivative with respect to the time $t$, $\partial h(\psi_p, \theta, \zeta, t)/\partial t$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dh_of_psip_dt(psip, theta, zeta, t)[()]
