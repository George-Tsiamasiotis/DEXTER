import pytest
from math import isfinite
from dexter import CosHarmonic, NcHarmonic
from dexter import CosPerturbation, NcPerturbation, Perturbation, perturbation


def test_helper_constructor(nc_harmonic: NcHarmonic):
    emptyp = perturbation([])
    cosp = perturbation([CosHarmonic(0, 0, 0, 0)])
    ncp = perturbation([nc_harmonic])

    assert isinstance(emptyp, CosPerturbation)
    assert isinstance(cosp, CosPerturbation)
    assert isinstance(ncp, NcPerturbation)

    with pytest.raises(TypeError):
        perturbation([CosHarmonic(0, 0, 0, 0), nc_harmonic])  # type: ignore


def test_empty_perturbation():
    p = NcPerturbation([])
    assert p.harmonics == []
    assert len(p) == 0
    with pytest.raises(BaseException):
        p[0]

    args = (0.01, 1, 2, 3)
    assert p.p_of_psi(*args) == 0
    assert p.p_of_psip(*args) == 0
    assert p.dp_dpsi(*args) == 0
    assert p.dp_dpsip(*args) == 0
    assert p.dp_of_psi_dtheta(*args) == 0
    assert p.dp_of_psip_dzeta(*args) == 0
    assert p.dp_of_psi_dtheta(*args) == 0
    assert p.dp_of_psip_dzeta(*args) == 0
    assert p.dp_of_psi_dt(*args) == 0
    assert p.dp_of_psip_dt(*args) == 0


def test_cos_perturbation():
    p = CosPerturbation(
        [
            CosHarmonic(1e-3, 1, 2, 0.0),
            CosHarmonic(1e-3, 1, 3, 0.0),
            CosHarmonic(1e-3, 1, 4, 0.0),
            CosHarmonic(1e-3, 1, 5, 0.0),
        ]
    )
    assert len(p) == 4
    _test_perturbation_evals(p)


def test_nc_perturbation(nc_perturbation: NcPerturbation):
    assert len(nc_perturbation) == 3
    _test_perturbation_evals(nc_perturbation)


def _test_perturbation_evals(p: Perturbation):
    args = (0.01, 1, 2, 3)
    assert isfinite(p.p_of_psi(*args))
    assert isfinite(p.p_of_psip(*args))
    assert isfinite(p.dp_dpsi(*args))
    assert isfinite(p.dp_dpsip(*args))
    assert isfinite(p.dp_of_psi_dtheta(*args))
    assert isfinite(p.dp_of_psip_dzeta(*args))
    assert isfinite(p.dp_of_psi_dtheta(*args))
    assert isfinite(p.dp_of_psip_dzeta(*args))
    assert isfinite(p.dp_of_psi_dt(*args))
    assert isfinite(p.dp_of_psip_dt(*args))
