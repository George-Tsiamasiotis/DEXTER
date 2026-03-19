import numpy as np

from dexter import (
    InitialFluxArray1,
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
    i = InitialFluxArray1("Toroidal", np.linspace(0, 0.5, 10))
    assert i.kind == "Toroidal"
    assert isinstance(i.values, np.ndarray)
    i = InitialFluxArray1("Poloidal", np.linspace(0, 0.5, 10))
    assert i.kind == "Poloidal"
    assert isinstance(i.values, np.ndarray)
    assert isinstance(i.__str__(), str)
    assert isinstance(i.__repr__(), str)


def test_queue_initial_conditions():
    particle_count = 10
    psi0s = InitialFluxArray1("Toroidal", np.linspace(0, 0.5, particle_count))
    q = QueueInitialConditions(
        t0=np.zeros(particle_count),
        flux0=psi0s,
        theta0=np.zeros(particle_count),
        zeta0=np.zeros(particle_count),
        rho0=np.full(particle_count, 1e-5),
        mu0=np.full(particle_count, 1e-6),
    )
    assert len(q) == particle_count
    assert q.t_array.shape == (particle_count,)
    assert q.theta_array.shape == (particle_count,)
    assert q.zeta_array.shape == (particle_count,)
    assert q.rho_array.shape == (particle_count,)
    assert q.mu_array.shape == (particle_count,)
    assert isinstance(q.__str__(), str)
    assert isinstance(q.__repr__(), str)


def test_queue_intstantiation():
    particle_count = 10
    psi0s = InitialFluxArray1("Toroidal", np.linspace(0, 0.5, particle_count))
    initial_conditions = QueueInitialConditions(
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


class TestQueue:
    @classmethod
    def setup_class(cls) -> None:
        cls.equilibrium = Equilibrium(
            qfactor=ParabolicQfactor(1.1, 3.9, ("Toroidal", 0.45)),
            current=LarCurrent(),
            bfield=LarBfield(),
            perturbation=Perturbation(
                [
                    CosHarmonic(1e-3, 1, 3, 0),
                    CosHarmonic(1e-3, 2, 3, 0),
                ]
            ),
        )
        particle_count = 10
        psi0s = InitialFluxArray1("Toroidal", np.linspace(0.01, 0.4, particle_count))
        cls.initial_conditions = QueueInitialConditions(
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

    def test_queue_intersection(self):
        queue = Queue(self.initial_conditions)
        intersect_params = IntersectParams("ConstTheta", 0.0, 5)
        queue.intersect(
            self.equilibrium,
            intersect_params,
        )
        for particle in queue.particles:
            assert particle.steps_taken > 10
            assert particle.integration_status == "Intersected"
