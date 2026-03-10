"""Evaluation traits implementation and documentation"""

from dexter._core import _PyLarGeometry, _PyNcGeometry
from dexter._core import _PyUnityQfactor, _PyParabolicQfactor, _PyNcQfactor
from dexter._core import _PyLarCurrent, _PyNcCurrent
from dexter._core import _PyLarBfield, _PyNcBfield
from dexter._core import _PyCosHarmonic, _PyNcHarmonic


class _FluxCommuteTrait:
    """Documents the methods provided by the 'FluxCommute' trait."""

    _rust: _PyNcGeometry | _PyUnityQfactor | _PyParabolicQfactor | _PyNcQfactor
    """`FluxCommute` implementors"""

    def psip_of_psi(self, psi: float) -> float:
        r"""The $\psi_p(\psi)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._rust.psip_of_psi(psi)

    def psi_of_psip(self, psip: float) -> float:
        r"""The $\psi(\psi_p)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._rust.psi_of_psip(psip)


class _GeometryTrait:
    """Documents the methods provided by the 'Geometries' trait."""

    _rust: _PyLarGeometry | _PyNcGeometry
    """`Geometry` implementors"""

    def r_of_psi(self, psi: float) -> float:
        r"""The $r(\psi)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._rust.r_of_psi(psi)

    def r_of_psip(self, psip: float) -> float:
        r"""The $r(\psi_p)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._rust.r_of_psip(psip)

    def psi_of_r(self, r: float) -> float:
        r"""The $\psi(r)$ value in Normalized Units.

        Parameters
        ----------
        r
            The radial distance $r$ in $[m]$.
        """
        return self._rust.psi_of_r(r)

    def psip_of_r(self, r: float) -> float:
        r"""The $\psi_p(r)$ value in Normalized Units.

        Parameters
        ----------
        r
            The radial distance $r$ in $[m]$.
        """
        return self._rust.psip_of_r(r)

    def rlab_of_psi(self, psi: float, theta: float) -> float:
        r"""The $R_{lab}(\psi, \theta)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.rlab_of_psi(psi, theta)

    def rlab_of_psip(self, psip: float, theta: float) -> float:
        r"""The $R_{lab}(\psi_p, \theta)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.rlab_of_psip(psip, theta)

    def zlab_of_psi(self, psi: float, theta: float) -> float:
        r"""The $Z_{lab}(\psi, \theta)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.zlab_of_psi(psi, theta)

    def zlab_of_psip(self, psip: float, theta: float) -> float:
        r"""The $Z_{lab}(\psi_p, \theta)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.zlab_of_psip(psip, theta)

    def jacobian_of_psi(self, psi: float, theta: float) -> float:
        r"""The $J(\psi, \theta)$ value in $[m]$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.jacobian_of_psi(psi, theta)

    def jacobian_of_psip(self, psip: float, theta: float) -> float:
        r"""The $J(\psi_p, \theta)$ value in $[m]$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.jacobian_of_psip(psip, theta)

    # Implemented as properties
    # def rlab_wall(self) -> Array1:
    # def zlab_wall(self) -> Array1:


class _QfactorTrait:
    """Documents the methods provided by the 'Qfactor' trait."""

    _rust: _PyUnityQfactor | _PyParabolicQfactor | _PyNcQfactor
    """`Qfactor` implementors"""

    def q_of_psi(self, psi: float) -> float:
        r"""The $q(\psi)$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._rust.q_of_psi(psi)

    def q_of_psip(self, psip: float) -> float:
        r"""The $q(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._rust.q_of_psip(psip)

    def dpsip_dpsi(self, psi: float) -> float:
        r"""The derivative $d\psi_p(\psi)/d\psi$ value.

        It's a good check that the values coincide with `qfactor.iota_of_psi(psi)`.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._rust.dpsip_dpsi(psi)

    def dpsi_dpsip(self, psip: float) -> float:
        r"""The derivative $d\psi(\psi_p)/d\psi_p$ value.

        It's a good check that the values coincide with `qfactor.q_of_psip(psip)`.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._rust.dpsi_dpsip(psip)

    def iota_of_psi(self, psi: float) -> float:
        r"""The $\iota(\psi) = \dfrac{1}{q(\psi)}$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._rust.iota_of_psi(psi)

    def iota_of_psip(self, psip: float) -> float:
        r"""The $\iota(\psi_p) = \dfrac{1}{q(\psi_p)}$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._rust.iota_of_psip(psip)


class _CurrentTrait:
    """Documents the methods provided by the 'Current' trait."""

    _rust: _PyLarCurrent | _PyNcCurrent
    """`Current` implementors"""

    def g_of_psi(self, psi: float) -> float:
        r"""The $g(\psi)$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._rust.g_of_psi(psi)

    def g_of_psip(self, psip: float) -> float:
        r"""The $g(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._rust.g_of_psip(psip)

    def i_of_psi(self, psi: float) -> float:
        r"""The $I(\psi)$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._rust.i_of_psi(psi)

    def i_of_psip(self, psip: float) -> float:
        r"""The $I(\psi_p)$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        """
        return self._rust.i_of_psip(psip)

    def dg_dpsi(self, psi: float) -> float:
        r"""The $dg(\psi)/d\psi$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._rust.dg_dpsi(psi)

    def dg_dpsip(self, psip: float) -> float:
        r"""The $dg(\psi_p)/d\psi_p$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi$ in Normalized Units.
        """
        return self._rust.dg_dpsip(psip)

    def di_dpsi(self, psi: float) -> float:
        r"""The $dg(\psi)/d\psi$ value.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        """
        return self._rust.di_dpsi(psi)

    def di_dpsip(self, psip: float) -> float:
        r"""The $dg(\psi_p)/d\psi_p$ value.

        Parameters
        ----------
        psip
            The poloidal flux $\psi$ in Normalized Units.
        """
        return self._rust.di_dpsip(psip)


class _BfieldTrait:
    """Documents the methods provided by the 'Bfield' trait."""

    _rust: _PyLarBfield | _PyNcBfield
    """`Bfield` implementors"""

    def b_of_psi(self, psi: float, theta: float) -> float:
        r"""The $B(\psi, \theta)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.b_of_psi(psi, theta)

    def b_of_psip(self, psip: float, theta: float) -> float:
        r"""The $B(\psi_p, \theta)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.b_of_psip(psip, theta)

    def db_dpsi(self, psi: float, theta: float) -> float:
        r"""The $dB(\psi, \theta)/d\psi$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.db_dpsi(psi, theta)

    def db_dpsip(self, psip: float, theta: float) -> float:
        r"""The $dB(\psi_p, \theta)/d\psi_p$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.db_dpsip(psip, theta)

    def db_of_psi_dtheta(self, psi: float, theta: float) -> float:
        r"""The $dB(\psi, \theta)/d\theta$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.db_of_psi_dtheta(psi, theta)

    def db_of_psip_dtheta(self, psip: float, theta: float) -> float:
        r"""The $dB(\psi_p, \theta)/d\theta$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        """
        return self._rust.db_of_psip_dtheta(psip, theta)


class _HarmonicTrait:
    """Documents the methods provided by the 'Harmonic' trait."""

    _rust: _PyCosHarmonic | _PyNcHarmonic
    """`Harmonic` implementors"""

    def alpha_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float:
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
        return self._rust.alpha_of_psi(psi, theta, zeta, t)

    def alpha_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float:
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
        return self._rust.alpha_of_psip(psip, theta, zeta, t)

    def phase_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float:
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
        return self._rust.phase_of_psi(psi, theta, zeta, t)

    def phase_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float:
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
        return self._rust.phase_of_psip(psip, theta, zeta, t)

    def h_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float:
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
        return self._rust.h_of_psi(psi, theta, zeta, t)

    def h_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float:
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
        return self._rust.h_of_psip(psip, theta, zeta, t)

    def dh_dpsi(self, psi: float, theta: float, zeta: float, t: float) -> float:
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
        return self._rust.dh_dpsi(psi, theta, zeta, t)

    def dh_dpsip(self, psip: float, theta: float, zeta: float, t: float) -> float:
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
        return self._rust.dh_dpsip(psip, theta, zeta, t)

    def dh_of_psi_dtheta(
        self, psi: float, theta: float, zeta: float, t: float
    ) -> float:
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
        return self._rust.dh_of_psi_dtheta(psi, theta, zeta, t)

    def dh_of_psip_dtheta(
        self, psip: float, theta: float, zeta: float, t: float
    ) -> float:
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
        return self._rust.dh_of_psip_dtheta(psip, theta, zeta, t)

    def dh_of_psi_dzeta(self, psi: float, theta: float, zeta: float, t: float) -> float:
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
        return self._rust.dh_of_psi_dzeta(psi, theta, zeta, t)

    def dh_of_psip_dzeta(
        self, psip: float, theta: float, zeta: float, t: float
    ) -> float:
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
        return self._rust.dh_of_psip_dzeta(psip, theta, zeta, t)

    def dh_of_psi_dt(self, psi: float, theta: float, zeta: float, t: float) -> float:
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
        return self._rust.dh_of_psi_dt(psi, theta, zeta, t)

    def dh_of_psip_dt(self, psip: float, theta: float, zeta: float, t: float) -> float:
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
        return self._rust.dh_of_psip_dt(psip, theta, zeta, t)
