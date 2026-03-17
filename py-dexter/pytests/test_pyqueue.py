import numpy as np

from dexter import InitialFluxArray1, QueueInitialConditions


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
