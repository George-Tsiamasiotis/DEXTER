import pytest
from dexter import InitialConditions


def test_fields(initial_conditions: InitialConditions):
    assert isinstance(initial_conditions.time0, float)
    assert isinstance(initial_conditions.theta0, float)
    assert isinstance(initial_conditions.psip0, float)
    assert isinstance(initial_conditions.rho0, float)
    assert isinstance(initial_conditions.zeta0, float)
    assert isinstance(initial_conditions.mu, float)


def test_immutability(initial_conditions: InitialConditions):
    with pytest.raises(AttributeError):
        initial_conditions.time0 *= 2
    with pytest.raises(AttributeError):
        initial_conditions.theta0 *= 2
    with pytest.raises(AttributeError):
        initial_conditions.psip0 *= 2
    with pytest.raises(AttributeError):
        initial_conditions.rho0 *= 2
    with pytest.raises(AttributeError):
        initial_conditions.zeta0 *= 2
    with pytest.raises(AttributeError):
        initial_conditions.mu *= 2


def test_repr(initial_conditions: InitialConditions):
    str(initial_conditions)
