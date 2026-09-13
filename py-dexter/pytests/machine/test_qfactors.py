import numpy as np
import pytest
import dexter as dex
from math import isclose, isfinite

from semver import Version

LCFS = dex.MagneticFlux.Toroidal(0.05)


def test_unity():
    qfactor = dex.UnityQfactor(LCFS)
    assert qfactor.machine_type == "Analytical"
    assert qfactor.psi_state == "Good"
    assert qfactor.psip_state == "Good"
    assert qfactor.psi_last == LCFS
    assert qfactor.psip_last == dex.MagneticFlux.Poloidal(LCFS.value)
    assert qfactor.qaxis == 1
    assert qfactor.qlast == 1
    qfactor.__repr__()
    qfactor.__str__()
    _test_qfactor_base(qfactor)
    with pytest.raises(Exception):
        qfactor.eval_psi_of_q(1)
    with pytest.raises(Exception):
        qfactor.eval_psip_of_q(1)


def test_parabolic():
    qfactor = dex.ParabolicQfactor(1.1, 3.9, LCFS)
    _test_qfactor_base(qfactor)
    assert qfactor.psi_last == LCFS
    assert isclose(
        qfactor.psip_last.value, qfactor.eval_other(psi=qfactor.psi_last.value)
    )
    assert qfactor.qaxis == 1.1
    assert qfactor.qlast == 3.9


def test_nc(nc_qfactor: dex.NcQfactor):
    _test_qfactor_base(nc_qfactor)
    assert nc_qfactor.interp_type == "Cubic"
    assert isinstance(nc_qfactor.path, str)
    assert isinstance(nc_qfactor.netcdf_version, Version)
    assert isinstance(nc_qfactor.psi_array, np.ndarray)
    assert isinstance(nc_qfactor.psip_array, np.ndarray)
    assert isinstance(nc_qfactor.q_array, np.ndarray)


def _test_qfactor_base(qfactor: dex.QfactorObject):

    qfactor.__repr__()
    qfactor.__str__()

    assert qfactor.machine_type in ["Numerical", "Analytical"]
    assert qfactor.psi_state in ["Good", "Bad"]
    assert qfactor.psip_state in ["Good", "Bad"]

    assert isfinite(qfactor.psi_last.value)
    assert isfinite(qfactor.psip_last.value)
    assert isfinite(qfactor.qlast)
    assert isfinite(qfactor.qaxis)

    flux = 0.02
    fluxes = np.linspace(0, 0.04, 10)
    q = 1.2
    qs = np.linspace(1.1, 3, 10)

    try:

        qfactor.eval_q(psi=flux)
        qfactor.eval_q(psi=fluxes)
        qfactor.eval_q(psip=flux)
        qfactor.eval_q(psip=fluxes)

        qfactor.eval_other(psi=flux)
        qfactor.eval_other(psi=fluxes)
        qfactor.eval_other(psip=flux)
        qfactor.eval_other(psip=fluxes)

        qfactor.eval_psi_of_q(q=q)
        qfactor.eval_psi_of_q(q=qs)
        qfactor.eval_psip_of_q(q=q)
        qfactor.eval_psip_of_q(q=qs)

        qfactor.eval_deriv_of_other(psi=flux)
        qfactor.eval_deriv_of_other(psi=fluxes)
        qfactor.eval_deriv_of_other(psip=flux)
        qfactor.eval_deriv_of_other(psip=fluxes)

        qfactor.eval_deriv_wrt_other(psi=flux)
        qfactor.eval_deriv_wrt_other(psi=fluxes)
        qfactor.eval_deriv_wrt_other(psip=flux)
        qfactor.eval_deriv_wrt_other(psip=fluxes)

        qfactor.eval_iota(psi=flux)
        qfactor.eval_iota(psi=fluxes)
        qfactor.eval_iota(psip=flux)
        qfactor.eval_iota(psip=fluxes)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )
