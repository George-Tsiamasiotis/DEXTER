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

    theta_array = np.linspace(0, np.pi, 10)
    psit_array = dex.MagneticFluxArray.Toroidal(np.linspace(0, 0.01, 10))

    psip_array = dex.MagneticFluxArray.Poloidal(
        [nc_machine.qfactor.eval_other(psit.value) for psit in psit_array]
    )

    h_psit = h.eval_array(nc_machine, theta_array, psit_array, -0.01, 1e-4)
    h_psip = h.eval_array(nc_machine, theta_array, psip_array, -0.01, 1e-4)
    assert np.all(np.isfinite(h_psit))
    assert np.all(np.isfinite(h_psip))
    assert np.all(np.isclose(h_psit, h_psip))


def test_gc_hamiltonian_eval_grid(nc_machine: dex.Machine):
    h = dex.GCHamiltonian()

    theta_array = np.linspace(0, np.pi, 4)
    psit_array = dex.MagneticFluxArray.Toroidal(np.linspace(0, 0.01, 5))

    psip_array = dex.MagneticFluxArray.Poloidal(
        [nc_machine.qfactor.eval_other(psit.value) for psit in psit_array]
    )

    h_psit = h.eval_grid(nc_machine, theta_array, psit_array, -0.01, 1e-4)
    h_psip = h.eval_grid(nc_machine, theta_array, psip_array, -0.01, 1e-4)
    assert np.all(np.isfinite(h_psit))
    assert np.all(np.isfinite(h_psip))
    assert np.all(np.isclose(h_psit, h_psip))


def test_gc_hamiltonian_eval_theta_dot_array(nc_machine: dex.Machine):
    h = dex.GCHamiltonian()

    theta_array = np.linspace(0, np.pi, 10)
    psit_array = dex.MagneticFluxArray.Toroidal(np.linspace(0, 0.01, 10))

    psip_array = dex.MagneticFluxArray.Poloidal(
        [nc_machine.qfactor.eval_other(psit.value) for psit in psit_array]
    )

    h_psit = h.eval_theta_dot_array(nc_machine, theta_array, psit_array, -0.01, 1e-4)
    h_psip = h.eval_theta_dot_array(nc_machine, theta_array, psip_array, -0.01, 1e-4)
    assert np.all(np.isfinite(h_psit))
    assert np.all(np.isfinite(h_psip))
    assert np.all(np.isclose(h_psit, h_psip))
