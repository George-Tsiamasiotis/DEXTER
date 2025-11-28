import pytest
from dexter import Harmonic, Perturbation


def test_eval(perturbation: Perturbation):
    psip = 0.015
    theta = 1
    zeta = 2
    assert isinstance(perturbation.p(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dpsip(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dtheta(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dzeta(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dt(psip, theta, zeta), float)


def test_getitem(perturbation: Perturbation):
    assert isinstance(perturbation[0], Harmonic)
    assert isinstance(perturbation[1], Harmonic)
    assert len(perturbation) == 2
    with pytest.raises(match="Harmonic index out of bounds"):
        perturbation[20]


def test_magic(perturbation: Perturbation):
    assert len(perturbation) == 2
    assert isinstance(perturbation.__str__(), str)
    assert isinstance(perturbation.__repr__(), str)
