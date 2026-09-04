import numpy as np
import dexter as dex

from math import isfinite

LCFS = dex.MagneticFlux.Toroidal(0.05)


def test_empty_perturbation():
    per = dex.Perturbation()
    assert len(per) == 0
    per.__repr__()
    per.__str__()
    _test_evals(per)


def test_flute_perturbation():
    mode1 = dex.FluteMode(1e-4, LCFS, 3, 2, 0)
    mode2 = dex.FluteMode(2e-4, LCFS, 4, 3, 0)
    mode3 = dex.FluteMode(3e-4, LCFS, 5, 4, 0)
    per = dex.Perturbation([mode1, mode2, mode3])
    assert len(per) == 3
    per.__repr__()
    per.__str__()
    _test_evals(per)


def test_nc_perturbation(nc_flute_mode: dex.NcFluteMode):
    per = dex.Perturbation([nc_flute_mode, nc_flute_mode])
    assert len(per) == 2
    per.__repr__()
    per.__str__()
    _test_evals(per)


def _test_evals(per: dex.Perturbation):
    flux = 0.02
    fluxes = np.linspace(0.01, 0.04, 10)
    theta = zeta = t = 1.2
    thetas = zetas = ts = np.linspace(0.01, np.pi, 10)

    try:

        per.eval_p(psi=flux, theta=theta, zeta=zeta, t=t)
        per.eval_p(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_p(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_p(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_deriv_flux(psi=flux, theta=theta, zeta=zeta, t=t)
        per.eval_deriv_flux(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_deriv_flux(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_deriv_flux(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_deriv_theta(psi=flux, theta=theta, zeta=zeta, t=t)
        per.eval_deriv_theta(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_deriv_theta(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_deriv_theta(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_deriv_zeta(psi=flux, theta=theta, zeta=zeta, t=t)
        per.eval_deriv_zeta(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_deriv_zeta(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_deriv_zeta(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_deriv_t(psi=flux, theta=theta, zeta=zeta, t=t)
        per.eval_deriv_t(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_deriv_t(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_deriv_t(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )
