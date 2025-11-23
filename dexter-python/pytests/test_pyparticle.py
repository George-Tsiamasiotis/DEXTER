import numpy as np
from dexter import Evolution, InitialConditions, Particle


def test_initialization(initial_conditions: InitialConditions):
    Particle(initial_conditions)


def test_initial_condtions(initial_conditions: InitialConditions):
    assert isinstance(initial_conditions.time0, float)
    assert isinstance(initial_conditions.theta0, float)
    assert isinstance(initial_conditions.psip0, float)
    assert isinstance(initial_conditions.rho0, float)
    assert isinstance(initial_conditions.zeta0, float)
    assert isinstance(initial_conditions.mu, float)
    assert isinstance(
        Particle(initial_conditions).initial_conditions, InitialConditions
    )


def test_initial_evolution(initial_conditions: InitialConditions):
    particle = Particle(initial_conditions)
    evolution = particle.evolution
    assert isinstance(evolution, Evolution)
    assert evolution.steps_taken == 0
    assert evolution.steps_stored == 1  # initial state
    assert np.isnan(evolution.energy_std)
    evolution_arrays = [
        "time",
        "theta",
        "psip",
        "rho",
        "zeta",
        "psi",
        "ptheta",
        "pzeta",
        "energy",
    ]
    for array_name in evolution_arrays:
        assert hasattr(evolution, array_name)
        array = getattr(evolution, array_name)
        assert isinstance(array, np.ndarray)
        assert array.shape == (1,)


def test_initial_status(initial_conditions: InitialConditions):
    particle = Particle(initial_conditions)
    assert particle.status == "Initialized"


def test_initial_frequencies(initial_conditions: InitialConditions):
    particle = Particle(initial_conditions)
    assert particle.frequencies.omega_theta is None
    assert particle.frequencies.omega_zeta is None
    assert particle.frequencies.qkinetic is None
