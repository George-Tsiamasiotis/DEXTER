import pytest
import numpy as np
import dexter as dex
from math import isclose


def test_initial_conditions_boozer():
    with pytest.raises(RuntimeError):
        dex.InitialConditions()

    flux0 = dex.MagneticFlux.Toroidal(0.01)
    init = dex.InitialConditions.Boozer(0, flux0, 1, 2, 3, 4)
    assert init.coordinate_set == "BoozerToroidal"
    assert isinstance(init.flux0, dex.MagneticFlux)
    assert init.flux0.kind == "Toroidal"
    assert isclose(init.t0, 0)
    assert isclose(init.flux0.value, 0.01)
    assert isclose(init.theta0, 1)
    assert isclose(init.zeta0, 2)
    assert isclose(init.rho0, 3)
    assert isclose(init.mu0, 4)
    with pytest.raises(AttributeError):
        init.pzeta0
    init.__repr__()
    init.__str__()


def test_initial_conditions_mixed():
    flux0 = dex.MagneticFlux.Poloidal(0.01)
    init = dex.InitialConditions.Mixed(0, flux0, 1, 2, 3, 4)
    assert init.coordinate_set == "MixedPoloidal"
    assert isinstance(init.flux0, dex.MagneticFlux)
    assert init.flux0.kind == "Poloidal"
    assert isclose(init.t0, 0)
    assert isclose(init.flux0.value, 0.01)
    assert isclose(init.theta0, 1)
    assert isclose(init.zeta0, 2)
    assert isclose(init.pzeta0, 3)
    assert isclose(init.mu0, 4)
    with pytest.raises(AttributeError):
        init.rho0
    init.__repr__()
    init.__str__()


def test_magnetic_flux_array():
    with pytest.raises(RuntimeError):
        dex.MagneticFluxArray()

    dex.MagneticFluxArray.Toroidal(np.linspace(0, 1, 10))
    dex.MagneticFluxArray.Poloidal(np.linspace(0, 1, 10))
    dex.MagneticFluxArray.Toroidal(1)
    dex.MagneticFluxArray.Toroidal((1, 2))
    with pytest.raises(TypeError):
        dex.MagneticFluxArray.Toroidal(np.zeros((2, 2)))
    with pytest.raises(ValueError):
        dex.MagneticFluxArray.Toroidal(np.nan)


def test_queue_initial_conditions():
    with pytest.raises(RuntimeError):
        dex.QueueInitialConditions()

    num = 10
    psi0s = dex.MagneticFluxArray.Toroidal(np.linspace(0, 0.05, num))
    initial = dex.QueueInitialConditions.Boozer(
        t0=np.zeros(num),
        flux0=psi0s,
        theta0=np.zeros(num),
        zeta0=np.zeros(num),
        rho0=np.full(num, 1e-5),
        mu0=np.full(num, 1e-6),
    )
    psip0s = dex.MagneticFluxArray.Poloidal(np.linspace(0, 0.05, num))
    initial = dex.QueueInitialConditions.Mixed(
        t0=np.zeros(num),
        flux0=psip0s,
        theta0=np.zeros(num),
        zeta0=np.zeros(num),
        pzeta0=np.full(num, -0.02),
        mu0=np.full(num, 1e-6),
    )
