import pytest
import numpy as np

from dexter import (
    LastClosedFluxSurface,
    InitialConditions,
    InitialFluxArray,
    QueueInitialConditions,
    Queue,
    Particle,
    Equilibrium,
    ParabolicQfactor,
    LarCurrent,
    LarBfield,
    Perturbation,
    CosHarmonic,
    IntersectParams,
)


def test_initial_flux_array1():
    i = InitialFluxArray("Toroidal", np.linspace(0, 0.5, 10))
    assert i.kind == "Toroidal"
    assert isinstance(i.values, np.ndarray)
    i = InitialFluxArray("Poloidal", np.linspace(0, 0.5, 10))
    assert i.kind == "Poloidal"
    assert isinstance(i.values, np.ndarray)
    assert isinstance(i.__str__(), str)
    assert isinstance(i.__repr__(), str)


def test_queue_boozer_initial_conditions():
    particle_count = 10
    psi0s = InitialFluxArray("Toroidal", np.linspace(0, 0.5, particle_count))
    init = QueueInitialConditions.boozer(
        t0=np.zeros(particle_count),
        flux0=psi0s,
        theta0=np.zeros(particle_count),
        zeta0=np.zeros(particle_count),
        rho0=np.full(particle_count, 1e-5),
        mu0=np.full(particle_count, 1e-6),
    )
    for i in init:
        assert isinstance(i, InitialConditions)
        assert i.coordinate_set == "BoozerToroidal"
    assert len(init) == particle_count
    assert init.t_array.shape == (particle_count,)
    assert init.theta_array.shape == (particle_count,)
    assert init.zeta_array.shape == (particle_count,)
    assert init.mu_array.shape == (particle_count,)
    assert init.rho_array is not None and init.rho_array.shape == (particle_count,)
    with pytest.raises(AttributeError):
        init.pzeta_array
    assert isinstance(init.__str__(), str)
    assert isinstance(init.__repr__(), str)


def test_queue_mixed_initial_conditions():
    particle_count = 10
    psi0s = InitialFluxArray("Poloidal", np.linspace(0, 0.5, particle_count))
    init = QueueInitialConditions.mixed(
        t0=np.zeros(particle_count),
        flux0=psi0s,
        theta0=np.zeros(particle_count),
        zeta0=np.zeros(particle_count),
        pzeta0=np.full(particle_count, -0.025),
        mu0=np.full(particle_count, 1e-6),
    )
    for i in init:
        assert isinstance(i, InitialConditions)
        assert i.coordinate_set == "MixedPoloidal"
    assert len(init) == particle_count
    assert init.t_array.shape == (particle_count,)
    assert init.theta_array.shape == (particle_count,)
    assert init.zeta_array.shape == (particle_count,)
    assert init.mu_array.shape == (particle_count,)
    with pytest.raises(AttributeError):
        init.rho_array
    assert init.pzeta_array is not None and init.pzeta_array.shape == (particle_count,)
    assert isinstance(init.__str__(), str)
    assert isinstance(init.__repr__(), str)


def test_queue_boozer_instantiation():
    particle_count = 10
    psi0s = InitialFluxArray("Toroidal", np.linspace(0, 0.5, particle_count))
    initial_conditions = QueueInitialConditions.boozer(
        t0=np.zeros(particle_count),
        flux0=psi0s,
        theta0=np.zeros(particle_count),
        zeta0=np.zeros(particle_count),
        rho0=np.full(particle_count, 1e-5),
        mu0=np.full(particle_count, 1e-6),
    )
    queue = Queue(initial_conditions)
    particles = queue.particles
    assert queue.particle_count == len(particles) == particle_count
    assert isinstance(queue.initial_conditions, QueueInitialConditions)
    assert isinstance(queue[0], Particle)


def test_queue_mixed_instantiation():
    particle_count = 10
    psi0s = InitialFluxArray("Toroidal", np.linspace(0, 0.5, particle_count))
    initial_conditions = QueueInitialConditions.mixed(
        t0=np.zeros(particle_count),
        flux0=psi0s,
        theta0=np.zeros(particle_count),
        zeta0=np.zeros(particle_count),
        pzeta0=np.full(particle_count, -0.025),
        mu0=np.full(particle_count, 1e-6),
    )
    queue = Queue(initial_conditions)
    particles = queue.particles
    assert queue.particle_count == len(particles) == particle_count
    assert isinstance(queue.initial_conditions, QueueInitialConditions)
    assert isinstance(queue[0], Particle)


class TestQueue:
    @classmethod
    def setup_class(cls) -> None:
        cls.lcfs = LastClosedFluxSurface("Toroidal", 0.45)
        cls.equilibrium = Equilibrium(
            qfactor=ParabolicQfactor(1.1, 3.9, cls.lcfs),
            current=LarCurrent(),
            bfield=LarBfield(),
            perturbation=Perturbation(
                [
                    CosHarmonic(1e-3, cls.lcfs, 1, 3, 0),
                    CosHarmonic(1e-3, cls.lcfs, 2, 3, 0),
                ]
            ),
        )
        cls.unperturbed_equilibrium = Equilibrium(
            qfactor=ParabolicQfactor(1.1, 3.9, cls.lcfs),
            current=LarCurrent(),
            bfield=LarBfield(),
        )
        particle_count = 10
        psi0s = InitialFluxArray("Toroidal", np.linspace(0.01, 0.4, particle_count))
        cls.initial_conditions = QueueInitialConditions.boozer(
            t0=np.zeros(particle_count),
            flux0=psi0s,
            theta0=np.zeros(particle_count),
            zeta0=np.zeros(particle_count),
            rho0=np.full(particle_count, 1e-6),
            mu0=np.full(particle_count, 1e-7),
        )

    def test_queue_integration(self):
        queue = Queue(self.initial_conditions)
        queue.integrate(
            self.equilibrium,
            (0, 1e4),
        )
        for particle in queue.particles:
            assert particle.steps_taken > 10
            assert particle.integration_status == "Integrated"
        assert np.all(np.isfinite(queue.energy_array))
        assert np.all(np.isfinite(queue.durations))

    def test_queue_intersection(self):
        queue = Queue(self.initial_conditions)
        intersect_params = IntersectParams("ConstTheta", 0.0, 5)
        queue.intersect(
            self.equilibrium,
            intersect_params,
        )
        for particle in queue.particles:
            assert isinstance(particle, Particle)
            assert particle.steps_taken > 10
            assert particle.integration_status == "Intersected"

    def test_queue_close(self):
        queue = Queue(self.initial_conditions)
        queue.close(self.unperturbed_equilibrium)
        for particle in queue.particles:
            assert isinstance(particle, Particle)
            assert particle.steps_taken > 10
            assert particle.integration_status == "ClosedPeriods(1)"
        assert np.all(np.isfinite(queue.omega_theta_array))
        assert np.all(np.isfinite(queue.omega_zeta_array))
        assert np.all(np.isfinite(queue.qkinetic_array))
        assert queue.steps_taken_array.ndim == 1
        assert queue.steps_stored_array.ndim == 1
