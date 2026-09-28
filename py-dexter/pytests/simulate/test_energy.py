import numpy as np
import dexter as dex


def test_energy_of_psi_grid(lar_machine: dex.Machine):
    theta_array = np.linspace(0, 2 * np.pi, 30)
    psi_array = np.linspace(0, lar_machine.psi_last.value, 20)
    theta_grid, psi_grid = dex.create_poloidal_grid(theta_array, psi_array)
    energy_grid = dex.energy_of_psi_grid(lar_machine, -0.03, 1e-4, theta_grid, psi_grid)
    assert energy_grid.shape == (20, 30)
    assert np.all(np.isfinite(energy_grid))


def test_energy_of_psip_grid(nc_machine: dex.Machine):
    theta_array = np.linspace(0, 2 * np.pi, 30)
    psip_array = np.linspace(0, nc_machine.psip_last.value, 20)
    theta_grid, psip_grid = dex.create_poloidal_grid(theta_array, psip_array)
    energy_grid = dex.energy_of_psip_grid(
        nc_machine, -0.03, 1e-4, theta_grid, psip_grid
    )
    assert energy_grid.shape == (20, 30)
    assert np.all(np.isfinite(energy_grid))
