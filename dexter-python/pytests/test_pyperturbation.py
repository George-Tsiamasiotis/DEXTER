import pytest
from dexter import NcHarmonic, NcPerturbation


def test_nc_perturbation_eval(nc_perturbation: NcPerturbation):
    psip = 0.015
    theta = 1
    zeta = 2
    assert isinstance(nc_perturbation.p(psip, theta, zeta), float)
    assert isinstance(nc_perturbation.dp_dpsip(psip, theta, zeta), float)
    assert isinstance(nc_perturbation.dp_dtheta(psip, theta, zeta), float)
    assert isinstance(nc_perturbation.dp_dzeta(psip, theta, zeta), float)
    assert isinstance(nc_perturbation.dp_dt(psip, theta, zeta), float)
    assert len(nc_perturbation) == 2
    assert isinstance(nc_perturbation.__str__(), str)
    assert isinstance(nc_perturbation.__repr__(), str)


def test_nc_perturbation_getitem(nc_perturbation: NcPerturbation):
    assert isinstance(nc_perturbation[0], NcHarmonic)
    assert isinstance(nc_perturbation[1], NcHarmonic)
    assert len(nc_perturbation) == 2
    with pytest.raises(match="NcHarmonic index out of bounds"):
        nc_perturbation[20]
