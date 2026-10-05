from math import isclose, isfinite

import dexter as dex
import numpy as np


def test_gc_hamiltonian_eval(lar_machine: dex.Machine):
    h = dex.GCHamiltonian

    h1 = h.eval(lar_machine, np.pi, dex.MagneticFlux.Toroidal(0.01), -0.01, 1e-4)
    h2 = h.eval(lar_machine, np.pi, dex.MagneticFlux.Poloidal(0.01), -0.01, 1e-4)
    assert isfinite(h1)
    assert isfinite(h2)
    assert not isclose(h1, h2)


def test_gc_hamiltonian_eval_array(nc_machine: dex.Machine):
    h = dex.GCHamiltonian()

    theta = np.linspace(0, np.pi, 10)
    flux = np.linspace(0, 0.01, 10)
    h1 = h.eval_psi_array(nc_machine, theta, flux, -0.01, 1e-4)
    h2 = h.eval_psip_array(nc_machine, theta, flux, -0.01, 1e-4)
    assert np.all(np.isfinite(h1))
    assert np.all(np.isfinite(h2))
    assert not np.all(np.isclose(h1, h2))


def test_gc_hamiltonian_eval_grid(nc_machine: dex.Machine):
    h = dex.GCHamiltonian()

    theta = np.linspace(0, np.pi, 4)
    flux = np.linspace(0, 0.01, 4)
    theta, flux = dex.create_poloidal_grid(theta, flux)
    h1 = h.eval_psi_grid(nc_machine, theta, flux, -0.01, 1e-4)
    h2 = h.eval_psip_grid(nc_machine, theta, flux, -0.01, 1e-4)
    assert np.all(np.isfinite(h1))
    assert np.all(np.isfinite(h2))
    assert not np.all(np.isclose(h1, h2))
