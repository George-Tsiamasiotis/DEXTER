import dexter as dex
import pytest


from dexter.utils import get_max_threads
from math import isclose


def test_threads():
    dex.set_num_threads(get_max_threads())


def test_magnetic_flux():
    with pytest.raises(RuntimeError):
        dex.MagneticFlux()

    flux = dex.MagneticFlux.Toroidal(0.1)
    flux.__repr__()
    flux.__str__()
    assert flux.value == 0.1
    assert flux.kind == "Toroidal"
    flux = dex.MagneticFlux.Poloidal(0.2)
    assert flux.value == 0.2
    assert flux.kind == "Poloidal"

    assert dex.MagneticFlux.Toroidal(0.1) == dex.MagneticFlux.Toroidal(0.1)
    assert dex.MagneticFlux.Toroidal(0.1) != dex.MagneticFlux.Toroidal(0.2)
    assert dex.MagneticFlux.Toroidal(0.1) != dex.MagneticFlux.Poloidal(0.1)


def test_magnetic_flux_add():
    flux = dex.MagneticFlux.Toroidal(1)
    flux2 = 2 + flux + 2
    assert flux.kind == "Toroidal"
    assert isclose(flux2.value, 5)

    flux = dex.MagneticFlux.Poloidal(1)
    flux2 = 2 + flux + 2
    assert flux.kind == "Poloidal"
    assert isclose(flux2.value, 5)


def test_magnetic_flux_sub():
    flux = dex.MagneticFlux.Toroidal(2)
    flux2 = 6 - flux - 1
    assert flux.kind == "Toroidal"
    assert isclose(flux2.value, 3)

    flux = dex.MagneticFlux.Poloidal(2)
    flux2 = 6 - flux - 1
    assert flux.kind == "Poloidal"
    assert isclose(flux2.value, 3)


def test_magnetic_flux_mul():
    flux = dex.MagneticFlux.Toroidal(3)
    flux2 = 2 * flux * 2
    assert flux.kind == "Toroidal"
    assert isclose(flux2.value, 12)

    flux = dex.MagneticFlux.Poloidal(3)
    flux2 = 2 * flux * 2
    assert flux.kind == "Poloidal"
    assert isclose(flux2.value, 12)


def test_magnetic_flux_div():
    flux = dex.MagneticFlux.Toroidal(5)
    flux2 = flux / 2
    assert flux.kind == "Toroidal"
    assert isclose(flux2.value, 2.5)

    flux = dex.MagneticFlux.Poloidal(5)
    flux2 = flux / 2
    assert flux.kind == "Poloidal"
    assert isclose(flux2.value, 2.5)


def test_flux_eval_wrappers(nc_machine_perturbed: dex.Machine):
    qfactor = nc_machine_perturbed.qfactor
    bfield = nc_machine_perturbed.bfield
    perturbation = nc_machine_perturbed.perturbation

    flux = 0.0002
    theta = 1.2
    zeta = 2.4
    t = 10

    # `_flux_eval_wrap1d`
    qfactor.eval_q(psi=flux)
    qfactor.eval_q(psip=flux)
    with pytest.raises(TypeError):
        qfactor.eval_q(psi=flux, psip=flux)
    with pytest.raises(TypeError):
        qfactor.eval_q()

    # `_flux_eval_wrap2d`
    bfield.eval_b(psi=flux, theta=theta)
    bfield.eval_b(psip=flux, theta=theta)
    with pytest.raises(TypeError):
        bfield.eval_b(psi=flux, psip=flux, theta=theta)
    with pytest.raises(TypeError):
        bfield.eval_b(theta=theta)

    # `_flux_eval_wrap4d`
    perturbation.eval_p(psi=flux, theta=theta, zeta=zeta, t=t)
    perturbation.eval_p(psip=flux, theta=theta, zeta=zeta, t=t)
    with pytest.raises(TypeError):
        perturbation.eval_p(psi=flux, psip=flux, theta=theta, zeta=zeta, t=t)
    with pytest.raises(TypeError):
        perturbation.eval_p(theta=theta, zeta=zeta, t=t)
