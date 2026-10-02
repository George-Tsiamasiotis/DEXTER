import pytest
import numpy as np
import dexter as dex


def test_queue_from_particles(lar_machine: dex.Machine):
    psi0 = dex.MagneticFlux.Toroidal(0.01)
    i = dex.InitialConditions.Boozer(0, psi0, 0, 0, 1e-4, 1e-6)
    p1 = dex.Particle(i)
    p2 = dex.Particle(i)
    p1.close(lar_machine)
    queue = dex.Queue.FromParticles([p1, p2])
    particles = queue.particles()
    assert len(particles) == 2
    assert isinstance(particles[0], dex.Particle)
    assert isinstance(particles[1], dex.Particle)
    assert "Closed" in particles[0].integration_status
    assert particles[1].integration_status == "Initialized"


def test_queue_routines_analytical(lar_machine: dex.Machine):
    num = 4
    psi0s = dex.MagneticFluxArray.Toroidal(
        np.linspace(1e-6, lar_machine.psi_last.value, num)
    )
    initial = dex.QueueInitialConditions.Boozer(
        t0=np.zeros(num),
        flux0=psi0s,
        theta0=np.zeros(num),
        zeta0=np.zeros(num),
        rho0=np.full(num, 1e-5),
        mu0=np.full(num, 1e-6),
    )
    queue = dex.Queue(initial)

    queue.integrate(lar_machine, (0, 100))
    intersect_params = dex.IntersectParams("ConstZeta", 0, 4)
    queue.intersect(lar_machine, intersect_params)
    queue.close(lar_machine)
    queue.classify(lar_machine)


def test_queue_routines_numerical(nc_machine: dex.Machine):
    num = 4
    psi0s = dex.MagneticFluxArray.Toroidal(
        np.linspace(1e-6, nc_machine.psi_last.value, num)
    )
    initial = dex.QueueInitialConditions.Boozer(
        t0=np.zeros(num),
        flux0=psi0s,
        theta0=np.zeros(num),
        zeta0=np.zeros(num),
        rho0=np.full(num, 1e-3),
        mu0=np.full(num, 1e-5),
    )
    queue = dex.Queue(initial)

    queue.integrate(nc_machine, (0, 100))
    intersect_params = dex.IntersectParams("ConstZeta", 0, 4)
    queue.intersect(nc_machine, intersect_params)
    queue.close(nc_machine)
    queue.classify(nc_machine)


def test_queue_slicing(lar_machine: dex.Machine):
    num = 4
    psi0s = dex.MagneticFluxArray.Toroidal(
        np.linspace(1e-6, lar_machine.psi_last.value, num)
    )
    initial = dex.QueueInitialConditions.Boozer(
        t0=np.zeros(num),
        flux0=psi0s,
        theta0=np.zeros(num),
        zeta0=np.zeros(num),
        rho0=np.full(num, 1e-5),
        mu0=np.full(num, 1e-6),
    )
    queue = dex.Queue(initial)

    queue.retain_pzeta(-0.1, 0)
    queue.retain_energy(-0.1, 0)
    queue.retain_energy_pzeta_positions(["Alpha", "Iota"])
    queue.retain_orbit_types(["TrappedConfined", "Potato"])


def test_queue_getters(lar_machine: dex.Machine):
    num = 4
    psi0s = dex.MagneticFluxArray.Toroidal(
        np.linspace(1e-6, lar_machine.psi_last.value, num)
    )
    initial = dex.QueueInitialConditions.Mixed(
        t0=np.zeros(num),
        flux0=psi0s,
        theta0=np.zeros(num),
        zeta0=np.zeros(num),
        pzeta0=np.full(num, 1e-5),
        mu0=np.full(num, 1e-6),
    )
    queue = dex.Queue(initial)
    queue.close(lar_machine)

    assert isinstance(queue.initial_conditions, dex.QueueInitialConditions)
    assert queue.routines == ["Close"]
    assert len(queue) == 4

    assert isinstance(queue[0], dex.Particle)
    assert isinstance(queue[-1], dex.Particle)
    with pytest.raises(IndexError):
        queue[10]
    with pytest.raises(IndexError):
        queue[-5]
    for p in queue:
        p.steps_stored

    assert queue.energy_array is not None and queue.energy_array.shape == (num,)
    assert queue.energy_rsd_array is not None and queue.energy_rsd_array.shape == (num,)
    assert queue.omega_theta_array is not None and queue.omega_theta_array.shape == (
        num,
    )
    assert queue.omega_zeta_array is not None and queue.omega_zeta_array.shape == (num,)
    assert queue.qkinetic_array is not None and queue.qkinetic_array.shape == (num,)

    assert np.all(np.isclose(queue._r._initial_pzetas(), np.full(num, 1e-5)))
    assert isinstance(queue._r._initial_energies(), np.ndarray)
    assert isinstance(queue._r._initial_energies()[0], float)
    assert isinstance(queue._r._orbit_types(), list)
    assert isinstance(queue._r._orbit_types()[0], str)
    assert queue._r._orbit_types()[0] in dex.OrbitType.__args__
