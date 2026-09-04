import dexter as dex
import pytest


def test_lcfs():
    with pytest.raises(RuntimeError):
        dex.MagneticFlux()

    lcfs = dex.MagneticFlux.Toroidal(0.1)
    lcfs.__repr__()
    lcfs.__str__()
    assert lcfs.value == 0.1
    assert lcfs.kind == "Toroidal"
    lcfs = dex.MagneticFlux.Poloidal(0.2)
    assert lcfs.value == 0.2
    assert lcfs.kind == "Poloidal"

    assert dex.MagneticFlux.Toroidal(0.1) == dex.MagneticFlux.Toroidal(0.1)
    assert dex.MagneticFlux.Toroidal(0.1) != dex.MagneticFlux.Toroidal(0.2)
    assert dex.MagneticFlux.Toroidal(0.1) != dex.MagneticFlux.Poloidal(0.1)


def test_flux_eval_wrappers(nc_machine_perturbed: dex.Machine):
    qfactor = nc_machine_perturbed.qfactor
    bfield = nc_machine_perturbed.bfield
    geometry = nc_machine_perturbed.geometry
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
