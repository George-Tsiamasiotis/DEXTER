import pytest
import numpy as np
from dexter import (
    Bfield,
    Currents,
    Evolution,
    InitialConditions,
    Particle,
    MappingParameters,
    Perturbation,
    Qfactor,
)


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


# ==========================================================================


def test_mapping_parameters():
    params = MappingParameters("ConstTheta", alpha=3.1415, intersections=10)
    assert params.section == "ConstTheta"
    assert params.alpha == 3.1415
    assert params.intersections == 10
    params = MappingParameters("ConstZeta", alpha=1, intersections=100)
    assert params.section == "ConstZeta"
    assert params.alpha == 1
    assert params.intersections == 100


# ==========================================================================


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


# ==========================================================================


def test_initial_status(initial_conditions: InitialConditions):
    particle = Particle(initial_conditions)
    assert particle.status == "Initialized"


# ==========================================================================


def test_initial_frequencies(initial_conditions: InitialConditions):
    particle = Particle(initial_conditions)
    assert particle.frequencies.omega_theta is None
    assert particle.frequencies.omega_zeta is None
    assert particle.frequencies.qkinetic is None


# ==========================================================================
# ==========================================================================


def test_particle_integrate(
    initial_conditions: InitialConditions,
    qfactor: Qfactor,
    currents: Currents,
    bfield: Bfield,
    perturbation: Perturbation,
):
    particle = Particle(initial_conditions)

    assert particle.status == "Initialized"
    particle.integrate(
        qfactor=qfactor,
        currents=currents,
        bfield=bfield,
        perturbation=perturbation,
        t_eval=(0, 100),
    )
    assert particle.status == "Integrated"
    assert particle.evolution.steps_taken > 2
    assert particle.evolution.steps_stored > 2

    str(particle)  # Inlcudes all nested fields


def test_particle_map(
    initial_conditions: InitialConditions,
    qfactor: Qfactor,
    currents: Currents,
    bfield: Bfield,
    perturbation: Perturbation,
):
    params = MappingParameters(section="ConstTheta", alpha=3.14, intersections=5)

    particle = Particle(initial_conditions)

    assert particle.status == "Initialized"
    particle.map(
        qfactor=qfactor,
        currents=currents,
        bfield=bfield,
        perturbation=perturbation,
        params=params,
    )
    assert particle.status == "Mapped"
    assert particle.evolution.steps_taken > 6
    assert particle.evolution.steps_stored == 6  # + initial point


def test_particle_calculate_frequencies(
    initial_conditions: InitialConditions,
    qfactor: Qfactor,
    currents: Currents,
    bfield: Bfield,
):
    particle = Particle(initial_conditions)

    assert particle.status == "Initialized"
    assert particle.frequencies.omega_theta is None
    assert particle.frequencies.omega_zeta is None
    assert particle.frequencies.qkinetic is None

    particle.calculate_frequencies(
        qfactor=qfactor,
        currents=currents,
        bfield=bfield,
        perturbation=Perturbation([]),  # Unperturbed system for frequencies
    )
    assert particle.status == "SinglePeriodIntegrated"
    assert particle.evolution.steps_taken > 2
    assert particle.evolution.steps_stored > 2

    assert particle.frequencies.omega_theta is not None
    assert particle.frequencies.omega_zeta is not None
    assert particle.frequencies.qkinetic is not None
