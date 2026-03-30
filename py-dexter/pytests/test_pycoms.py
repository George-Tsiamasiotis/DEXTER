import numpy as np
from dexter import COMs, Equilibrium


def test_instantiation_and_fields():
    coms = COMs(energy=1e-5, pzeta=0.01, mu=1e-5)
    assert coms.energy == 1e-5
    assert coms.pzeta == 0.01
    assert coms.mu == 1e-5


def test_energy_of_psi_grid(toroidal_nc_equilibrium: Equilibrium):
    coms = COMs(pzeta=0.01, mu=1e-5)
    psi_array = np.array([0.01, 0.02])
    theta_array = np.linspace(-np.pi, np.pi, 10)
    energy_grid = coms.energy_of_psi_grid(
        toroidal_nc_equilibrium,
        psi_array,
        theta_array,
    )
    assert energy_grid.shape == (2, 10)
    assert not np.any(np.isnan(energy_grid))


def test_energy_of_psip_grid(poloidal_nc_equilibrium: Equilibrium):
    coms = COMs(pzeta=0.01, mu=1e-5)
    psip_array = np.array([0.01, 0.02])
    theta_array = np.linspace(-np.pi, np.pi, 10)
    energy_grid = coms.energy_of_psip_grid(
        poloidal_nc_equilibrium,
        psip_array,
        theta_array,
    )
    assert energy_grid.shape == (2, 10)
    assert not np.any(np.isnan(energy_grid))
