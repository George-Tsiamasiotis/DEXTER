r"""Defines the `GCHamiltonian` unit type."""

from warnings import warn

from dexter.machine.machine import Machine
from dexter.machine.flux import MagneticFlux, MagneticFluxArray
from dexter.types import Array1, Array2
from dexter.comspace.curve import StationaryCurve

from dexter._core import _PyGCHamiltonian
from dexter._utils import _ReprStrImpl, _into_pyarray1f64


class GCHamiltonian(_ReprStrImpl):
    r"""Guiding center Hamiltonian.

    The guiding center Hamiltonian is defined as:

    $$
    \mathcal{H}(\theta, \psi/\psi_p; P_\zeta, \mu)
    = \dfrac{(P_\zeta+\psi_p)^2B^2}{2g^2} + \mu B.
    $$

    This object has no internal state, it only provides evaluation methods of the Hamiltonian and
    its derived quantities.
    """

    _r: type[_PyGCHamiltonian]

    def __init__(self) -> None:
        self._r = _PyGCHamiltonian

    @classmethod
    def eval(
        cls,
        machine: Machine,
        theta: float,
        flux: MagneticFlux,
        pzeta: float | MagneticFlux,
        mu: float,
    ) -> Array1:
        r"""Evaluates the Hamiltonian $\mathcal{H}(\theta, \psi/\psi_p; P_\zeta, \mu)$.

        Parameters
        ----------
        machine
            The machine on which the Hamiltonian is defined.
        theta
            The $\theta$ value.
        flux
            The magnetic flux value $\psi$ or $\psi_p$.
        pzeta
            The $P_\zeta$ constant of motion.
        mu
            The $\mu$ constant of motion.

        Example
        -------

        ``` py
        >>> geometry = dex.LarGeometry(baxis=3, raxis=1.75, rlast=0.5)
        >>> LCFS = geometry.psi_last
        >>> machine = dex.Machine(
        ...     geometry=geometry,
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.9, lcfs=LCFS),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ... )
        >>>
        >>> psi = dex.MagneticFlux.Toroidal(0.01)
        >>> pzeta = -0.1 * machine.psip_last
        >>> h = dex.GCHamiltonian.eval(machine, 3.14, psi, pzeta, 1e-4)

        ```

        """
        if isinstance(pzeta, MagneticFlux):
            if pzeta.kind == "Toroidal":
                warn("'MagneticFlux.Toroidal' value encountered for 'pzeta0'")
            pzeta = pzeta.value

        return _PyGCHamiltonian.eval(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            theta=theta,
            flux=flux._r,
            pzeta=pzeta,
            mu=mu,
        )

    @classmethod
    def eval_array(
        cls,
        machine: Machine,
        theta_array: Array1,
        flux_array: MagneticFluxArray,
        pzeta: float | MagneticFlux,
        mu: float,
    ) -> Array1:
        r"""Evaluates the Hamiltonian $\mathcal{H}(\theta, \psi/\psi_p; P_\zeta, \mu)$ over a 1D array.

        !!! note

            The produced array has the form:

            $$| \mathcal{H}(\theta_0, \psi_0),\
                \mathcal{H}(\theta_1, \psi_1),\
                \dotsc,\
                \mathcal{H}(\theta_N, \psi_N)\ |$$

        This method uses the [`#!python flux_array.kind`][dexter.MagneticFluxKind] flux coordinate
        to perform all evaluations. Since only $\psi_p$ appears explicitly in the Hamiltonian, if
        `#!python flux_array.kind == "Toroidal"` then an extra evaluation is needed to calculate $\psi_p$.

        Parameters
        ----------
        machine
            The machine on which the Hamiltonian is defined.
        theta_array
            The $\theta$ values.
        flux_array
            The $\psi/\psi_p$ values.
        pzeta
            The $P_\zeta$ constant of motion.
        mu
            The $\mu$ constant of motion.

        Example
        -------

        ``` py
        >>> geometry = dex.LarGeometry(baxis=3, raxis=1.75, rlast=0.5)
        >>> LCFS = geometry.psi_last
        >>> machine = dex.Machine.FromNetcdf(path, "Akima", "Bicubic")
        >>>
        >>> pzeta, mu = -0.1 * machine.psip_last, 1e-4
        >>> theta_array = np.linspace(0, 6.28, 10)
        >>>
        >>> # H(ψ, θ)
        >>> psi_array = dex.MagneticFluxArray.Toroidal(np.linspace(0, machine.psi_last.value, 10))
        >>> h = dex.GCHamiltonian.eval_array(machine, theta_array, psi_array, pzeta, mu)
        >>>
        >>> # H(ψp, θ)
        >>> psip_array = dex.MagneticFluxArray.Poloidal(np.linspace(0, machine.psip_last.value, 10))
        >>> h = dex.GCHamiltonian.eval_array(machine, theta_array, psip_array, pzeta, mu)

        ```

        """
        if isinstance(pzeta, MagneticFlux):
            if pzeta.kind == "Toroidal":
                warn("'MagneticFlux.Toroidal' value encountered for 'pzeta0'")
            pzeta = pzeta.value

        return _PyGCHamiltonian.eval_array(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            theta_array=_into_pyarray1f64(theta_array),
            flux_array=flux_array._r,
            pzeta=pzeta,
            mu=mu,
        )

    @classmethod
    def eval_grid(
        cls,
        machine: Machine,
        theta_array: Array1,
        flux_array: MagneticFluxArray,
        pzeta: float | MagneticFlux,
        mu: float,
    ) -> Array2:
        r"""Evaluates the Hamiltonian $\mathcal{H}(\theta, \psi/\psi_p; P_\zeta, \mu)$ over a 2D grid.

        !!! note

            The $\theta$ and $\psi/\psi_p$ arrays do not need to have the same length. The produced
            grid has the form:

            $$|\ \mathcal{H}(\theta_0, \psi_0),
               \ \mathcal{H}(\theta_1, \psi_0),
               \ \dotsc,
               \ \mathcal{H}(\theta_M, \psi_0)\ | \\\
              |\ \mathcal{H}(\theta_0, \psi_1),
               \ \mathcal{H}(\theta_1, \psi_1),
               \ \dotsc,
               \ \mathcal{H}(\theta_M, \psi_1)\ | \\\
              ... \\\
              |\mathcal{H}(\theta_0, \psi_N),
               \mathcal{H}(\theta_1, \psi_N),
               \dotsc,
               \mathcal{H}(\theta_M, \psi_N) | $$

        This method uses the [`#!python flux_array.kind`][dexter.MagneticFluxKind] flux coordinate
        to perform all evaluations. Since only $\psi_p$ appears explicitly in the Hamiltonian, if
        `#!python flux_array.kind == "Toroidal"` then an extra evaluation is needed to calculate $\psi_p$.

        Parameters
        ----------
        machine
            The machine on which the Hamiltonian is defined.
        theta_array
            The $\theta$ values.
        flux_array
            The $\psi/\psi_p$ values.
        pzeta
            The $P_\zeta$ constant of motion.
        mu
            The $\mu$ constant of motion.

        Example
        -------

        ``` py
        >>> geometry = dex.LarGeometry(baxis=3, raxis=1.75, rlast=0.5)
        >>> LCFS = geometry.psi_last
        >>> machine = dex.Machine.FromNetcdf(path, "Akima", "Bicubic")
        >>>
        >>> pzeta, mu = -0.1 * machine.psip_last, 1e-4
        >>> theta_array = np.linspace(0, 6.28, 4)
        >>>
        >>> # H(ψ, θ)
        >>> psi_array = dex.MagneticFluxArray.Toroidal(np.linspace(0, machine.psi_last.value, 5))
        >>> h = dex.GCHamiltonian.eval_grid(machine, theta_array, psi_array, pzeta, mu)
        >>>
        >>> # H(ψp, θ)
        >>> psip_array = dex.MagneticFluxArray.Poloidal(np.linspace(0, machine.psip_last.value, 5))
        >>> h = dex.GCHamiltonian.eval_grid(machine, theta_array, psip_array, pzeta, mu)

        ```

        """
        if isinstance(pzeta, MagneticFlux):
            if pzeta.kind == "Toroidal":
                warn("'MagneticFlux.Toroidal' value encountered for 'pzeta0'")
            pzeta = pzeta.value

        return _PyGCHamiltonian.eval_grid(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            theta_array=_into_pyarray1f64(theta_array),
            flux_array=flux_array._r,
            pzeta=pzeta,
            mu=mu,
        )

    @classmethod
    def eval_theta_dot_array(
        cls,
        machine: Machine,
        theta_array: Array1,
        flux_array: MagneticFluxArray,
        pzeta: float | MagneticFlux,
        mu: float,
    ) -> Array1:
        r"""Evaluates $\dot\theta = \partial\mathcal{H}(\theta, \psi/\psi_p; P_\zeta, \mu)/\partial P_\theta$
        over a 1D array

        This method uses the [`#!python flux_array.kind`][dexter.MagneticFluxKind] flux coordinate
        to perform all evaluations. Since only $\psi_p$ appears explicitly in the equation (inside the
        $\rho_{||}$ term), if `#!python flux_array.kind == "Toroidal"` then an extra evaluation is
        needed to calculate $\psi_p$.

        Parameters
        ----------
        machine
            The machine on which the Hamiltonian is defined.
        theta_array
            The $\theta$ values.
        flux_array
            The $\psi/\psi_p$ values.
        pzeta
            The $P_\zeta$ constant of motion.
        mu
            The $\mu$ constant of motion.

        Example
        -------

        ``` py
        >>> geometry = dex.LarGeometry(baxis=3, raxis=1.75, rlast=0.5)
        >>> LCFS = geometry.psi_last
        >>> machine = dex.Machine.FromNetcdf(path, "Akima", "Bicubic")
        >>>
        >>> pzeta, mu = -0.1 * machine.psip_last, 1e-4
        >>> theta_array = np.linspace(0, 6.28, 10)
        >>>
        >>> # dot(θ(ψ, θ))
        >>> psi_array = dex.MagneticFluxArray.Toroidal(np.linspace(0, machine.psi_last.value, 10))
        >>> theta_dot = dex.GCHamiltonian.eval_theta_dot_array(
        ...     machine,
        ...     theta_array,
        ...     psi_array,
        ...     pzeta,
        ...     mu
        ... )
        >>>
        >>> # dot(θ(ψp, θ))
        >>> psip_array = dex.MagneticFluxArray.Poloidal(np.linspace(0, machine.psip_last.value, 10))
        >>> theta_dot = dex.GCHamiltonian.eval_theta_dot_array(
        ...     machine,
        ...     theta_array,
        ...     psip_array,
        ...     pzeta,
        ...     mu
        ... )

        ```

        """
        if isinstance(pzeta, MagneticFlux):
            if pzeta.kind == "Toroidal":
                warn("'MagneticFlux.Toroidal' value encountered for 'pzeta0'")
            pzeta = pzeta.value

        return _PyGCHamiltonian.eval_theta_dot_array(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            theta_array=_into_pyarray1f64(theta_array),
            flux_array=flux_array._r,
            pzeta=pzeta,
            mu=mu,
        )

    @classmethod
    def build_stationary_curve(cls, machine: Machine) -> StationaryCurve:
        r"""Builds the [`StationaryCurve`][dexter.StationaryCurve] of `machine`.

        This routines tries to construct the curve w.r.t $\psi$ first. If $\psi$ is not a good
        coordinate, then $\psi_p$ is used.
        """
        return StationaryCurve._wrap(
            _PyGCHamiltonian.build_stationary_curve(
                qfactor=machine.qfactor._r,
                current=machine.current._r,
                bfield=machine.bfield._r,
            )
        )
