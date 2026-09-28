r"""Calculations of the Energy in a 2D grids.

Functions
---------
energy_of_psi_grid
    Calculates the energy on a 2D meshgrid of the $\theta$ and $\psi$ arrays, in Normalized Units.
energy_of_psip_grid
    Calculates the energy on a 2D meshgrid of the $\theta$ and $\psi_p$ arrays, in Normalized Units.
create_poloidal_grid
    Creates a $(\theta, \psi/\psi_p)$ meshgrid from the two 1D arrays.
"""

from dexter._core import (
    _py_create_poloidal_grid,
    _py_energy_of_psi_grid,
    _py_energy_of_psip_grid,
)

from dexter.machine.machine import Machine
from dexter.types import Array1, Array2


def create_poloidal_grid(
    theta_array: Array1,
    flux_array: Array1,
) -> tuple[Array2, Array2]:
    r"""Creates a $(\theta, \psi/\psi_p)$ meshgrid from the two 1D arrays.

    This method is equivalent to `#!python np.meshgrid(theta_array, flux_array)`.

    Parameters
    ----------
    theta_array
        The 1D array representing the $\theta$ coordinate of the grid.
    flux_array
        The 1D array representing the $\psi$ or $\psi_p$ coordinate of the grid.
    """
    return _py_create_poloidal_grid(theta_array=theta_array, flux_array=flux_array)


def energy_of_psi_grid(
    machine: Machine,
    pzeta: float,
    mu: float,
    theta_array: Array2,
    psi_array: Array2,
) -> Array2:
    r"""Calculates the energy on a 2D meshgrid of the $\theta$ and $\psi$ arrays, in Normalized Units.

    Use [`create_poloidal_grid`][dexter.create_poloidal_grid] to construct the meshgrid arrays.

    Parameters
    ----------
    machine
        The machine on which to evaluate the Hamiltonian.
    pzeta
        The $P_\zeta$ constant of motion.
    mu
        The $\mu$ constant of motion.
    theta_array
        The 2D $\theta$ array
    psi_array
        The 2D $\psi$ array
    """
    return _py_energy_of_psi_grid(
        qfactor=machine.qfactor._r,
        current=machine.current._r,
        bfield=machine.bfield._r,
        perturbation=machine.perturbation._r,
        pzeta=pzeta,
        mu=mu,
        theta_array=theta_array,
        psi_array=psi_array,
    )
    pass


def energy_of_psip_grid(
    machine: Machine,
    pzeta: float,
    mu: float,
    theta_array: Array2,
    psip_array: Array2,
) -> Array2:
    r"""Calculates the energy on a 2D meshgrid of the $\theta$ and $\psi_p$ arrays, in Normalized Units.

    Use [`create_poloidal_grid`][dexter.create_poloidal_grid] to construct the meshgrid arrays.

    Parameters
    ----------
    machine
        The machine on which to evaluate the Hamiltonian.
    pzeta
        The $P_\zeta$ constant of motion.
    mu
        The $\mu$ constant of motion.
    theta_array
        The 2D $\theta$ array
    psip_array
        The 2D $\psi_p$ array
    """
    return _py_energy_of_psip_grid(
        qfactor=machine.qfactor._r,
        current=machine.current._r,
        bfield=machine.bfield._r,
        perturbation=machine.perturbation._r,
        pzeta=pzeta,
        mu=mu,
        theta_array=theta_array,
        psip_array=psip_array,
    )
    pass
