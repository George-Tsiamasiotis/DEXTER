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
    psi_array = dex.MagneticFluxArray.Toroidal(np.linspace(0, 0.01, 10))

    psip_array = dex.MagneticFluxArray.Poloidal(
        [nc_machine.qfactor.eval_other(psi.value) for psi in psi_array]
    )

    h_psi = h.eval_array(nc_machine, theta_array, psi_array, -0.01, 1e-4)
    h_psip = h.eval_array(nc_machine, theta_array, psip_array, -0.01, 1e-4)
    assert np.all(np.isfinite(h_psi))
    assert np.all(np.isfinite(h_psip))
    assert np.all(np.isclose(h_psi, h_psip))


def test_gc_hamiltonian_eval_grid(nc_machine: dex.Machine):
    h = dex.GCHamiltonian()

    theta_array = np.linspace(0, np.pi, 4)
    psi_array = dex.MagneticFluxArray.Toroidal(np.linspace(0, 0.01, 5))

    psip_array = dex.MagneticFluxArray.Poloidal(
        [nc_machine.qfactor.eval_other(psi.value) for psi in psi_array]
    )

    h_psi = h.eval_grid(nc_machine, theta_array, psi_array, -0.01, 1e-4)
    h_psip = h.eval_grid(nc_machine, theta_array, psip_array, -0.01, 1e-4)
    assert np.all(np.isfinite(h_psi))
    assert np.all(np.isfinite(h_psip))
    assert np.all(np.isclose(h_psi, h_psip))


def test_gc_hamiltonian_eval_theta_dot_array(nc_machine: dex.Machine):
    h = dex.GCHamiltonian()

    theta_array = np.linspace(0, np.pi, 10)
    psi_array = dex.MagneticFluxArray.Toroidal(np.linspace(0, 0.01, 10))

    psip_array = dex.MagneticFluxArray.Poloidal(
        [nc_machine.qfactor.eval_other(psi.value) for psi in psi_array]
    )

    h_psi = h.eval_theta_dot_array(nc_machine, theta_array, psi_array, -0.01, 1e-4)
    h_psip = h.eval_theta_dot_array(nc_machine, theta_array, psip_array, -0.01, 1e-4)
    assert np.all(np.isfinite(h_psi))
    assert np.all(np.isfinite(h_psip))
    assert np.all(np.isclose(h_psi, h_psip))


def test_parametrized_segment(lar_machine: dex.Machine):
    segment = dex.ParametrizedSegment(
        dex.MagneticFluxArray.Toroidal([0.001, 0.002]),
        theta=np.asarray([0.3, 0.4]),
    )
    assert np.all(segment.t_array == np.array([0.0, 1.0]))
    assert np.all(segment.flux_array.values == np.array([0.001, 0.002]))
    assert np.all(segment.theta_array == np.array([0.3, 0.4]))
    assert segment.theta_dot_array is None

    segment.calculate_theta_dot(lar_machine, -0.0123, 1e-5)

    assert segment.theta_dot_array is not None
    assert np.all(np.isfinite(segment.theta_dot_array))
