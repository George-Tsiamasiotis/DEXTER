import numpy as np
import dexter as dex

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

        per.eval_ms(psit=flux, theta=theta, zeta=zeta, t=t)
        per.eval_ms(psit=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_ms(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_ms(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_dms_dflux(psit=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dms_dflux(psit=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_dms_dflux(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dms_dflux(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_dms_dtheta(psit=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dms_dtheta(psit=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_dms_dtheta(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dms_dtheta(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_dms_dzeta(psit=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dms_dzeta(psit=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_dms_dzeta(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dms_dzeta(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_dms_dt(psit=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dms_dt(psit=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_dms_dt(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dms_dt(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_phis(psit=flux, theta=theta, zeta=zeta, t=t)
        per.eval_phis(psit=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_phis(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_phis(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_dphis_dflux(psit=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dphis_dflux(psit=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_dphis_dflux(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dphis_dflux(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_dphis_dtheta(psit=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dphis_dtheta(psit=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_dphis_dtheta(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dphis_dtheta(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_dphis_dzeta(psit=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dphis_dzeta(psit=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_dphis_dzeta(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dphis_dzeta(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        per.eval_dphis_dt(psit=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dphis_dt(psit=fluxes, theta=thetas, zeta=zetas, t=ts)
        per.eval_dphis_dt(psip=flux, theta=theta, zeta=zeta, t=t)
        per.eval_dphis_dt(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )
