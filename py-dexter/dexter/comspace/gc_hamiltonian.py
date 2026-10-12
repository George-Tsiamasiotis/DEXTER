r"""Defines the `GCHamiltonian` type and some helper free functions."""

from warnings import warn

from dexter.machine.flux import MagneticFlux
from dexter.types import Array1, Array2

from dexter._core import _PyGCHamiltonian, _py_create_poloidal_grid
from dexter._utils import _ReprStrImpl

from dexter.machine.machine import Machine


def create_poloidal_grid(
    theta_array: Array1,
    flux_array: Array1,
) -> tuple[Array2, Array2]:
    r"""Creates a $(\theta, \psi/\psi_p)$ meshgrid from the two 1D arrays.

    This method is equivalent to `#!python np.meshgrid(theta_array, flux_array, indexing="xy")`.

    Parameters
    ----------
    theta_array
        The 1D array representing the $\theta$ coordinate of the grid.
    flux_array
        The 1D array representing the $\psi$ or $\psi_p$ coordinate of the grid.
    """
    return _py_create_poloidal_grid(theta_array=theta_array, flux_array=flux_array)


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
    def eval_psi_array(
        cls,
        machine: Machine,
        theta_array: Array1,
        psi_array: Array1,
        pzeta: float | MagneticFlux,
        mu: float,
    ) -> Array1:
        r"""Evaluates the Hamiltonian $\mathcal{H}(\theta, \psi; P_\zeta, \mu)$ over a 1D array.

        Use this method when you *know* that `psi` is a [`Good`][dexter.FluxCoordinateState] coordinate.

        Parameters
        ----------
        machine
            The machine on which the Hamiltonian is defined.
        theta_array
            The $\theta$ values.
        psi_array
            The $\psi$ values.
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
        >>> theta_array = np.linspace(0, 6.28, 10)
        >>> psi_array = np.linspace(0, machine.psi_last.value, 10)
        >>> pzeta = -0.1 * machine.psip_last
        >>> h = dex.GCHamiltonian.eval_psi_array(machine, theta_array, psi_array, pzeta, 1e-4)

        ```

        """
        if isinstance(pzeta, MagneticFlux):
            if pzeta.kind == "Toroidal":
                warn("'MagneticFlux.Toroidal' value encountered for 'pzeta0'")
            pzeta = pzeta.value

        return _PyGCHamiltonian.eval_psi_array(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            theta_array=theta_array,
            psi_array=psi_array,
            pzeta=pzeta,
            mu=mu,
        )

    @classmethod
    def eval_psip_array(
        cls,
        machine: Machine,
        theta_array: Array1,
        psip_array: Array1,
        pzeta: float | MagneticFlux,
        mu: float,
    ) -> Array1:
        r"""Evaluates the Hamiltonian $\mathcal{H}(\theta, \psi_p; P_\zeta, \mu)$ over a 1D array.

        Use this method when you *know* that `psip` is a [`Good`][dexter.FluxCoordinateState] coordinate.

        Parameters
        ----------
        machine
            The machine on which the Hamiltonian is defined.
        theta_array
            The $\theta$ values.
        psip_array
            The $\psi_p$ values.
        pzeta
            The $P_\zeta$ constant of motion.
        mu
            The $\mu$ constant of motion.

        Example
        -------

        ``` py
        >>> machine = dex.Machine.FromNetcdf(path, "Akima", "Bicubic")
        >>>
        >>> theta_array = np.linspace(0, 6.28, 10)
        >>> psip_array = np.linspace(0, machine.psip_last.value, 10)
        >>> pzeta = -0.1 * machine.psip_last
        >>> h = dex.GCHamiltonian.eval_psip_array(machine, theta_array, psip_array, pzeta, 1e-4)

        ```

        """
        if isinstance(pzeta, MagneticFlux):
            if pzeta.kind == "Toroidal":
                warn("'MagneticFlux.Toroidal' value encountered for 'pzeta0'")
            pzeta = pzeta.value

        return _PyGCHamiltonian.eval_psip_array(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            theta_array=theta_array,
            psip_array=psip_array,
            pzeta=pzeta,
            mu=mu,
        )

    @classmethod
    def eval_psi_grid(
        cls,
        machine: Machine,
        theta_grid: Array2,
        psi_grid: Array2,
        pzeta: float | MagneticFlux,
        mu: float,
    ) -> Array2:
        r"""Evaluates the Hamiltonian $\mathcal{H}(\theta, \psi; P_\zeta, \mu)$ over a 2D grid.

        Use this method when you *know* that `psi` is a [`Good`][dexter.FluxCoordinateState] coordinate.

        Parameters
        ----------
        machine
            The machine on which the Hamiltonian is defined.
        theta_grid
            The $\theta$ values.
        psi_grid
            The $\psi$ values.
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
        >>> theta_array = np.linspace(0, 6.28, 10)
        >>> psi_array = np.linspace(0, machine.psi_last.value, 10)
        >>> theta_grid, psi_grid = dex.create_poloidal_grid(theta_array, psi_array)
        >>> pzeta = -0.1 * machine.psip_last
        >>> h = dex.GCHamiltonian.eval_psi_grid(machine, theta_grid, psi_grid, pzeta, 1e-4)

        ```

        """
        if isinstance(pzeta, MagneticFlux):
            if pzeta.kind == "Toroidal":
                warn("'MagneticFlux.Toroidal' value encountered for 'pzeta0'")
            pzeta = pzeta.value

        return _PyGCHamiltonian.eval_psi_grid(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            theta_grid=theta_grid,
            psi_grid=psi_grid,
            pzeta=pzeta,
            mu=mu,
        )

    @classmethod
    def eval_psip_grid(
        cls,
        machine: Machine,
        theta_grid: Array2,
        psip_grid: Array2,
        pzeta: float | MagneticFlux,
        mu: float,
    ) -> Array2:
        r"""Evaluates the Hamiltonian $\mathcal{H}(\theta, \psi_p; P_\zeta, \mu)$ over a 2D grid.

        Use this method when you *know* that `psip` is a [`Good`][dexter.FluxCoordinateState] coordinate.

        Parameters
        ----------
        machine
            The machine on which the Hamiltonian is defined.
        theta_grid
            The $\theta$ values.
        psip_grid
            The $\psi_p$ values.
        pzeta
            The $P_\zeta$ constant of motion.
        mu
            The $\mu$ constant of motion.

        Example
        -------

        ``` py
        >>> machine = dex.Machine.FromNetcdf(path, "Akima", "Bicubic")
        >>>
        >>> theta_array = np.linspace(0, 6.28, 10)
        >>> psip_array = np.linspace(0, machine.psip_last.value, 10)
        >>> theta_grid, psip_grid = dex.create_poloidal_grid(theta_array, psip_array)
        >>> pzeta = -0.1 * machine.psip_last
        >>> h = dex.GCHamiltonian.eval_psip_grid(machine, theta_grid, psip_grid, pzeta, 1e-4)

        ```

        """
        if isinstance(pzeta, MagneticFlux):
            if pzeta.kind == "Toroidal":
                warn("'MagneticFlux.Toroidal' value encountered for 'pzeta0'")
            pzeta = pzeta.value

        return _PyGCHamiltonian.eval_psip_grid(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            theta_grid=theta_grid,
            psip_grid=psip_grid,
            pzeta=pzeta,
            mu=mu,
        )
